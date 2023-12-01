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

libs = c('QDNAseq', 'GenomicRanges', 'rtracklayer', 'remotes');
invisible(suppressPackageStartupMessages(lapply(libs, require, character.only=T)))
## qdnaseq uses hg19 default with package
## installing hg38 qdnaseq bin annot from package asntech/github - package is not in conda currently
if (!require("QDNAseq.hg38")){
	remotes::install_github("asntech/QDNAseq.hg38@main", upgrade="never")
}
library(QDNAseq.hg38)


## cmd line args
args = commandArgs(trailingOnly = TRUE)
if(length(args)>0) {
	predir = args[1]
	sample = args[2]
	annot.gtf.file = args[3]
}else{
    stop("ERROR: No predir, sample name, and/or annot gtf supplied.")
}



###########
#
# Call CNV
#
###########
## get bin annots
bins <- getBinAnnotations(binSize=15, genome="hg38")

## Read in sample bam files
readCounts <- binReadCounts(bins, bamfiles=paste0(predir,'/bam/',sample,'.bam'))

## Plot raw copy number profile (read counts across the genome) and highlight bins to remove with default filtering
png(file = paste0(predir,'/cnv/',sample,'_cn_profile_pre.png'))
plot(readCounts, logTransform=FALSE)
highlightFilters(readCounts, logTransform=FALSE, residual=TRUE, blacklist=TRUE)
garb=dev.off()

## Apply filters
readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE, chromosomes = c("Y", "MT"))

## plot median read counts as a function of GC content and mappability
png(file = paste0(predir,'/cnv/',sample,'_med_read_counts.png'))
isobarPlot(readCountsFiltered)
garb=dev.off()

## Estimate the correction for GC content and mappability
readCountsFiltered <- estimateCorrection(readCountsFiltered)

## plot for the relationship between the observed standard deviation in the data and its read depth
png(file = paste0(predir,'/cnv/',sample,'_std_read_depth.png'))
noisePlot(readCountsFiltered)
garb=dev.off()

## apply the correction for GC content and mappability then normalize, smooth outliers
copyNumbers <- correctBins(readCountsFiltered)
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)

## Copy number profile after correcting for GC content and mappability
png(file = paste0(predir,'/cnv/',sample,'_cn_profile_post.png'))
plot(copyNumbersSmooth)
garb=dev.off()

## export as BED file
invisible(exportBins(copyNumbersSmooth, file=paste0(predir,'/cnv/',sample,'_cnv.bed'), format="bed"))
invisible(exportBins(copyNumbersSmooth, file=paste0(predir,'/cnv/',sample,'_cnv.igv'), format="igv"))

## Segment
copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)

## Copy number profile of segments
png(file = paste0(predir,'/cnv/',sample,'_cn_profile_segment.png'))
plot(copyNumbersSegmented)
garb=dev.off()

## Call copy numbers
copyNumbersCalled <- callBins(copyNumbersSegmented)

## Copy number profile of segments
png(file = paste0(predir,'/cnv/',sample,'_cn_profile_calls.png'))
plot(copyNumbersCalled)
garb=dev.off()

## Export segmented and called copy numbers as IGV and BED
invisible(exportBins(copyNumbersSegmented, file=paste0(predir,'/cnv/',sample,"_cnv_segmented.bed"), format="bed", type="segments"))
invisible(exportBins(copyNumbersSegmented, file=paste0(predir,'/cnv/',sample,"_cnv_segmented.igv"), format="igv", type="segments"))
invisible(exportBins(copyNumbersCalled, file=paste0(predir,'/cnv/',sample,"_cnv_calls.bed"), format="bed", logTransform=FALSE, type="calls"))
invisible(exportBins(copyNumbersCalled, file=paste0(predir,'/cnv/',sample,"_cnv_calls.igv"), format="igv", logTransform=FALSE, type="calls"))



############
#
# Annot CNV
#
############
## sample cnv segmented calls
cnv.file = paste0(predir,'/cnv/',sample,'_cnv_segmented.bed')
cnv = read.csv(cnv.file, header=F, sep='')
cnv = cnv[-1,]
colnames(cnv) = c('chr', 'start', 'end', 'feature', 'score', 'strand')

## reference annot gtf
annot = rtracklayer::readGFF(annot.gtf.file)
annot = annot[annot$type=='gene', c('seqid','start', 'end', 'gene_id', 'description', 'gene_biotype', 'strand')]
colnames(annot) = c('chr', 'start', 'end', 'gene.id', 'description', 'gene.type', 'strand')
annot$chr = sub('chr','',annot$chr)

## create genomic ranges objects
cnv.ranges = GRanges(seqnames=cnv$chr, ranges=IRanges(start=as.numeric(cnv$start), end=as.numeric(cnv$end)))
annot.ranges = GRanges(seqnames=annot$chr, ranges=IRanges(start=as.numeric(annot$start), end=as.numeric(annot$end)))

## find overlaps in ranges
cnv.annot = findOverlaps(cnv.ranges, annot.ranges)

## add annot info to cnv data based on overlap hits
cnv$gene.id=NA
cnv$description=NA
cnv$gene.type=NA
cnv[queryHits(cnv.annot),'gene.id'] = annot[subjectHits(cnv.annot), 'gene.id']
cnv[queryHits(cnv.annot),'description'] = annot[subjectHits(cnv.annot), 'description']
cnv[queryHits(cnv.annot),'gene.type'] = annot[subjectHits(cnv.annot), 'gene.type']

## write annotated cnv
write.table(cnv, paste0(predir,'/cnv/',sample,'_cnv_annot.csv'), sep=',', quote=F, row.names=F, col.names=T)

