#############################################################################################################
# pipeline_template
#
# This program is free software that contains third party software subject to various licenses, 
# namely, the GNU General Public License version 3 (or later), the GNU Affero General Public License 
# version 3 (or later), and the LaTeX Project Public License v.1.3(c). A list of the software contained 
# in this program, including the applicable licenses, can be accessed here: 
#
# You can redistribute and/or modify this program, including its components, only under the terms of 
# the applicable license(s).  
#
# This program is distributed in the hope that it will be useful, but "as is," WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#
# Program:  init-qdnaseq-cnv-analysis.r
# Version:  0.1
# Author:   Sami R. Cherikh
# Purpose:  Quantitative DNA sequencing for chromosomal aberrations. The genome is divided into non-overlapping
#           fixed-sized bins, number of sequence reads in each counted, adjusted with a simultaneous two-dimensional
#           loess correction for sequence mappability and GC content, and filtered to remove spurious regions in the genome.
# Input:    <COMMAND LINE INPUTS>
# Output:
#############################################################################################################

libs = c('QDNAseq', 'QDNAseq.hg19');
invisible(suppressPackageStartupMessages(lapply(libs, require, character.only=T)))


## get path to config.xlsx and source dir command line arguments
args = commandArgs(trailingOnly = TRUE)
if(length(args)>0) {
	predir = args[1]
	sample = args[2]
}else{
    stop("ERROR: No predir and/or sample name supplied.")
}

#setwd(paste0(predir,'/cnv'))

## obtain bin annotations
bins <- getBinAnnotations(binSize=15)

## Read in sample bam files
readCounts <- binReadCounts(bins, bamfiles=paste0(predir,'/bam/',sample,'.bam'))

## Plot raw copy number profile (read counts across the genome) and highlight bins to remove with default filtering
pdf(file = paste0(predir,'/cnv/',sample,'_cn_profile_pre.pdf'))
plot(readCounts, logTransform=FALSE, ylim=c(-50, 200))
highlightFilters(readCounts, logTransform=FALSE, residual=TRUE, blacklist=TRUE)
garb=dev.off()

## Apply filters
readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE)

## plot median read counts as a function of GC content and mappability
pdf(file = paste0(predir,'/cnv/',sample,'_med_read_counts.pdf'))
isobarPlot(readCountsFiltered)
garb=dev.off()

## Estimate the correction for GC content and mappability
readCountsFiltered <- estimateCorrection(readCountsFiltered)

## plot for the relationship between the observed standard deviation in the data and its read depth
pdf(file = paste0(predir,'/cnv/',sample,'_std_read_depth.pdf'))
noisePlot(readCountsFiltered)
garb=dev.off()

## apply the correction for GC content and mappability then normalize, smooth outliers
copyNumbers <- correctBins(readCountsFiltered)
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)

## Copy number profile after correcting for GC content and mappability
pdf(file = paste0(predir,'/cnv/',sample,'_cn_profile_post.pdf'))
plot(copyNumbersSmooth)
garb=dev.off()

## export as BED file
invisible(exportBins(copyNumbersSmooth, file=paste0(predir,'/cnv/',sample,'_cnv.bed'), format="bed"))


