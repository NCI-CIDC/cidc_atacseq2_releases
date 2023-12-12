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
    peak.mode = args[4]
	make.report = args[5]
}else{
    stop("ERROR: No dirs, path to valid sample meta data file, and/or peak mode supplied.")
}


## read in sample metadata
mta = read.csv(mta.file, h=T, stringsAsFactors=F)
## subset metadata to columns needed for ChIPQC
mta.chipqc = mta[,c('samid', 'spct', 'factor', 'replicate','trt')]
colnames(mta.chipqc) = c('SampleID', 'Tissue', 'Factor', 'Replicate','Condition')

## add bam and peak files to mta
mta.chipqc$bamReads = paste0(predir,'/bam/',mta.chipqc$SampleID,'_sampled.bam')
mta.chipqc$Peaks = paste0(predir,'/peak/',mta.chipqc$SampleID,'_peaks.',peak.mode,'Peak')


## quantify for each sample listed in metadata
chipqc.obj = ChIPQC(mta.chipqc, annotation="hg38", chromosomes=NULL)

## if report
if(make.report=='make-report'){
    ChIPQCreport(chipqc.obj, reportName=paste0('all_samples_chipqc_report'), reportFolder=paste0(dtadir,'/../report'), facetBy=c("Condition"))

    png(paste0(predir,'/analysis/report/Regi.png'))
    plotRegi(chipqc.obj, facetBy=c("Condition"))
    garb=dev.off()
}

## convert sample QC metrics to data frame
res.all = as.data.frame(QCmetrics(chipqc.obj))

## write out sample results
write.table(res.all, paste0(dtadir,'/all_chipqc.csv'), sep=',', quote=F, row.names=T, col.names=T)
