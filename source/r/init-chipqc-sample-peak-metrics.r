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
# Program:  init-chipqc-sample-peak-metrics.r
# Version:  0.1
# Author:   Sami R. Cherikh
# Purpose:  Constructing a ChIPQsample object and retrieving peak metrics per sample
# Input:    <COMMAND LINE INPUTS>
# Output:
#############################################################################################################

libs = c('ChIPQC', 'TxDb.Hsapiens.UCSC.hg38.knownGene','data.table', 'GenomicAlignments');
#libs = c('ChIPQC', 'data.table');
invisible(suppressPackageStartupMessages(lapply(libs, require, character.only=T)))


# #############################
# ##    For each sample      ##
# #############################
## get path to config.xlsx and source dir command line arguments
args = commandArgs(trailingOnly = TRUE)
if(length(args)>0) {
	predir = args[1]
	sample = args[2]
	bam = args[3]
	peak = args[4]
}else{
    stop("ERROR: No predir and/or needed sample info supplied.")
}

## quantify for one sample
#data(blacklist_hg19)
#chipqc.obj = ChIPQCsample(paste0(predir,'/',bam), paste0(predir,'/',peak), annotation="hg19", blacklist = blacklist.hg19)#, blacklist = blacklist.hg38)
chipqc.obj = ChIPQCsample(bam, peak, annotation="hg38")
#chipqc.obj = ChIPQCsample(paste0(predir,'/bam/SAM01.bam'), peaks=paste0(predir,'/peak/SAM01_peaks.narrowPeak'), annotation="hg19", blacklist = blacklist.hg19, chromosomes='chr22') ## test

## create report
#ChIPQCreport(chipqc.obj, reportName=paste0(sample,'_chipqc_report'), reportFolder=paste0(predir,'/chipqc/',sample))

## convert sample QC metrics to data frame
res.sample = data.frame(as.list(QCmetrics(chipqc.obj)))
colnames(res.sample) = names(QCmetrics(chipqc.obj))
res.sample$SampleID = sample
res.sample = res.sample[,c(ncol(res.sample),1:(ncol(res.sample)-1))]

## read in genomic alignments from bam
align.obj = readGAlignments(bam)

## convert GAlignments object to data.table for fast aggregation 
aln = data.table(strand=as.factor(BiocGenerics::as.vector(strand(align.obj))),
                 seqnames=as.factor(BiocGenerics::as.vector(seqnames(align.obj))),
                 pos=ifelse(strand(align.obj) == "+", start(align.obj), end(align.obj))
                )

## aggregate reads by position and count them
readsPerPosition = aln[,list(count=.N), by=list(strand, seqnames, pos)]$count

## Number of genomic locations to which EXACTLY one unique mapping read maps
n1 = sum(readsPerPosition == 1)
## Number of genomic locations to which AT LEAST one unique mapping read maps, i.e. the number of non-redundant, unique mapping reads
nd = length(readsPerPosition)
## PBC = positions with exactly 1 read / positions with at least 1 read
pbc = n1 / nd

## add PBC values to the result dataframe
res.sample[c("N1", "Nd", "PBC")] = list(n1, nd, pbc)

## write out sample results
write.table(res.sample, paste0(predir,'/chipqc/',sample,'_chipqc.csv'), sep=',', quote=F, row.names=F, col.names=T)
#system(paste0('gzip ',predir,'/chipqc/',sample,'_chipqc.csv'))


