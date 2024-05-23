# QDNAseq uses hg19 by default with the package. This script installs QDNAseq.hg38 in the conda environment for future use.

args <- commandArgs(trailingOnly = TRUE)
output_file = args[1]

library(QDNAseq)
library(GenomicRanges)
library(rtracklayer)
library(remotes)

remotes::install_github("asntech/QDNAseq.hg38@main", upgrade="never")

library(QDNAseq.hg38)

args <- commandArgs(trailingOnly = TRUE)
output_file = args[1]

file.create(output_file)
