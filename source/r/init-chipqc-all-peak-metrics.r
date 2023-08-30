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
# Program:  init-chipqc-all-peak-metrics.r
# Version:  0.1
# Author:   Sami R. Cherikh
# Purpose:  Constructing a ChIPQsample object and retrieving peak metrics for all samples in mta
# Input:    <COMMAND LINE INPUTS>
# Output:
#############################################################################################################

libs = c('ChIPQC', 'TxDb.Hsapiens.UCSC.hg38.knownGene','data.table');
#libs = c('ChIPQC', 'data.table');
invisible(suppressPackageStartupMessages(lapply(libs, require, character.only=T)))


#############################
##   For samples in mta    ##
#############################
## get path to config.xlsx and source dir command line arguments
args = commandArgs(trailingOnly = TRUE)
if(length(args)>0) {
	predir = args[1]
	dtadir = args[2]
	mta.file = args[3]
	make.report = args[4]
}else{
    #stop("ERROR: No predir and/or sample metadata file supplied.")
    stop("ERROR: No predir and path to valid sample meta data file supplied.")
}
#predir = '~/atac_test_output'
#dtadir = paste0(predir,'/analysis/data')
#mta.file = paste0(predir,'/sample_metadata.csv')
#make.report = 'make-report'

## read in sample metadata
mta = read.csv(mta.file, h=T, stringsAsFactors=F)
## subset metadata to columns needed for ChIPQC
mta.chipqc = mta[,c('samid', 'spct', 'factor', 'replicate','trt')]
colnames(mta.chipqc) = c('SampleID', 'Tissue', 'Factor', 'Replicate','Condition')

mta.chipqc$bamReads=NA
mta.chipqc$Peaks=NA
## get bam and peak file locations for samples
for(s in 1:length(mta.chipqc$SampleID)){
    samid = mta.chipqc$SampleID[s]

    ## get bam file
    bams = list.files(paste0(predir,'/bam'))
    bam = bams[bams == paste0(samid,'.bam')]

    ## get peak file
    peaks = list.files(paste0(predir,'/peak'))
    peak = peaks[peaks %like% paste0(samid,'_peaks\\..+r+.+Peak$')]

    ## add to mta
    mta.chipqc[mta.chipqc$SampleID==samid, 'bamReads'] = paste0(predir,'/bam/',bam)
    mta.chipqc[mta.chipqc$SampleID==samid, 'Peaks'] = paste0(predir,'/peak/',peak)

}

## quantify for each sample listed in metadata ## TODO error maybe due being able to test with only one sample currently
#data(blacklist_hg19)
#chipqc.obj = ChIPQC(mta.chipqc, annotation="hg19", blacklist = blacklist.hg19)
chipqc.obj = ChIPQC(mta.chipqc, annotation="hg38", chromosomes=NULL)

## if report
if(make.report=='make-report'){
    ChIPQCreport(chipqc.obj, reportName=paste0('all_samples_chipqc_report'), reportFolder=paste0(dtadir,'/../report'), facetBy=c("Condition"))

    png(paste0(predir,'/analysis/report/Regi.png'))
    plotRegi(chipqc.obj, facetBy=c("Condition"))
    dev.off()
}

## convert sample QC metrics to data frame
res.all = as.data.frame(QCmetrics(chipqc.obj))
#colnames(res.all) = names(QCmetrics(chipqc.obj))
#res.sample$SampleID = sample
#res.sample = res.sample[,c(ncol(res.sample),1:(ncol(res.sample)-1))]

## write out sample results
write.table(res.all, paste0(dtadir,'/all_chipqc.csv'), sep=',', quote=F, row.names=T, col.names=T)
system(paste0('gzip -f ',dtadir,'/all_chipqc.csv'))
