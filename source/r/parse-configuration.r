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
# Program:  parse-configuration.r
# Version:  0.1
# Author:   Travis L. Jensen and Sami R. Cherikh
# Purpose:  Run sanity checks on configuration, download and parse data sets,
#				parse workflow and analysis configurations, create work space directories.
# Input:    <COMMAND LINE INPUTS>
# Output:
#############################################################################################################

#libs = c('openxlsx','stringr','biomaRt','yaml');
libs = c('openxlsx','stringr','yaml');
invisible(suppressPackageStartupMessages(lapply(libs, require, character.only=T)))


## get path to config.xlsx and source dir command line arguments
args = commandArgs(trailingOnly = TRUE)
if(length(args)>0) {
	infile       = args[1]
	source.dir   = args[2]
}


#####################
#
# Import data
#
#####################

## import metadata
metadata = read.xlsx(infile,sheet=1)

## remove unused fields
metadata = metadata[c(3:nrow(metadata)),]

## if no samples provided
if (nrow(metadata) == 0) stop('No samples were found\n')

## import workflow configs
workflow.config = read.xlsx(infile,sheet=2)[,c('Name','Value')] # key value fields only
workflow.config[is.na(workflow.config)] = '' # replace missing values with empty string

## import analysis configs
analysis.config = read.xlsx(infile,sheet=3)[,c('Name','Value')] # key value fields only

## define any needed variables - Example: ensembl version
ensembl.version = workflow.config$Value[workflow.config$Name=='ensembl_version']


#################
#
# Sanity Checks
#
#################
cat('Performing sanity checks on configuration options\n')

## SKIP this since using conda
# ## ensure all programs exist - TODO auto populate program list from config
# progs = c(workflow.config$Value[match(c('openssl_prog','samtools_prog','fastqc_prog','bwa_prog','cloud_prog','cutadapt_prog','trimmomatic_prog','rseqc_dir'),workflow.config$Name)])
# for (i in 1:length(progs)) {
# 	if(!file.exists(progs[i])) {
# 		cat(paste('Program does not exist at this location!',progs[i],'\n'))
# 		q(status=25,save='no')
# 	}
# }

## ensure ensembl version is valid ## TODO turning off for now
# ensembl.test = tryCatch(is(useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=ensembl.version),'Mart'), error = function(e) FALSE)
# if(ensembl.test==F) {
# 	cat(paste('The specified Ensembl version is invalid:',ensembl.version,'\n'))
# 	q(status=25,save='no')
# }


############################
#
# Parse Configurations
#
############################
cat('Parsing workflow and analysis configuration\n')

results.dir = analysis.config$Value[analysis.config$Name=='report_dir']
pre.dir = workflow.config$Value[workflow.config$Name=='pre_dir']
sub.dirs = workflow.config$Value[workflow.config$Name=='sub_dirs']

## create pre-processing dir and subdirs
if(!dir.exists(pre.dir)) {
	system(paste('mkdir',pre.dir))
}

if(!dir.exists(paste(pre.dir,'genome',sep='/'))) {
	system(paste('mkdir',paste(pre.dir,'genome',sep='/')))
}

if(!dir.exists(paste(pre.dir,'annot',sep='/'))) {
	system(paste('mkdir',paste(pre.dir,'annot',sep='/')))
}


## create analysis results dir and subdirs
if(!dir.exists(results.dir)) {
	system(paste('mkdir',results.dir))
}

## create analysis data dir
dta.dir = paste(results.dir,'data',sep='/')
if(!dir.exists(dta.dir)) {
	system(paste('mkdir',dta.dir))
}

## create analysis report dir
if(!dir.exists(paste(results.dir,'report',sep='/'))) {
	system(paste('mkdir',paste(results.dir,'report',sep='/')))
}

## add any variables from workflow config to analysis config
analysis.config = rbind(analysis.config,c('ensembl_version',ensembl.version))

## write analysis configuration
write.table(analysis.config,paste(dta.dir,'analysis_config.csv',sep='/'),sep=',',quote=T,row.names=F)


## write CSV formatted metadata file
cat('Creating sample metadata\n')
write.table(metadata,paste(pre.dir,'sample_metadata.csv',sep='/'),na='',sep=',',quote=F,row.names=F)


##################################
#
# Download/parse additional files
#
#################################
cat('Downloading / parsing annotation files\n')

## download version of genome from ensembl FTP
genome.file = paste(pre.dir,'/genome/Homo_sapiens.ensembl.version',ensembl.version,'.genome.fa',sep='')
if(!file.exists(genome.file)) {
	cat('Downloading Genome annotations from Ensembl FTP site\n')
	system(paste(source.dir,'/shell/download-ensembl-genome.sh ',ensembl.version,' ',dirname(genome.file),sep=''))
}

## download version of GTF from ensembl FTP
gtf.file = paste(pre.dir,'/annot/Homo_sapiens.ensembl.version',ensembl.version,'.chr.gtf',sep='')
if(!file.exists(gtf.file)) {
	cat('Downloading GTF annotations from Ensembl FTP site\n')
	system(paste(source.dir,'/shell/download-ensembl-gtf.sh ',ensembl.version,' ',dirname(gtf.file),sep=''))
}



#############################
#
# Convert config to YAML
#
#############################

workflow.config.yaml = list()
for (x in 1:nrow(workflow.config)) {
	workflow.config.yaml[[workflow.config$Name[x]]] = workflow.config$Value[x]
} 

## add sub dirs
workflow.config.yaml[['subdirs']]  = sub.dirs

## add sample IDs, paths to fastqs
workflow.config.yaml[['samid']]  = metadata$samid	
workflow.config.yaml[['fastq1']] = metadata$fastq_file_1

if (all(is.na(metadata$fastq_file_2))) {
	workflow.config.yaml[['fastq2']] = ''
	cat('Working in single-end mode\n')
} else {
	workflow.config.yaml[['fastq2']] = metadata$fastq_file_2
	cat('Working in paired-end mode\n')
}

## check for and add adapter sequences
if (all(is.na(metadata$threep_adapter_seq))) {
	cat("Skipping 3' adapter trimming... some or all adapters omitted.\n")
} 
if (all(is.na(metadata$fivep_adapter_seq))) {
	cat("Skipping 5' adapter trimming... some or all adapters omitted.\n")
} 

metadata$threep_adapter_seq[is.na(metadata$threep_adapter_seq)] = 'NA'
metadata$fivep_adapter_seq[is.na(metadata$fivep_adapter_seq)]   = 'NA'
workflow.config.yaml[['tp_adapter_seq']]  = metadata$threep_adapter_seq
workflow.config.yaml[['fp_adapter_seq']]  = metadata$fivep_adapter_seq


## check to make sure an appropriate quality trimming value was specified
if (is.na(workflow.config.yaml[['quality_trim']])) {
	stop("No quality trimming value specified")
} else if (!as.numeric(workflow.config.yaml[['quality_trim']]) %in% 0:40) {
	stop("Quality trimming value specified is outside allowed range: ", workflow.config.yaml[['quality_trim']])
}


## write out preprocessing config in YAML
sink(paste(pre.dir,'preprocess_config.yaml',sep='/'))
cat(as.yaml(workflow.config.yaml))
sink()


## write dirs to file
## workflow directory, analysis directory, workflow configuration, and sample_metadata csv
tmp.dta = c(pre.dir, analysis.config[analysis.config$Name=='report_dir','Value'], paste(pre.dir,'preprocess_config.yaml',sep='/'), paste(pre.dir,'sample_metadata.csv',sep='/'))
write.table(tmp.dta,paste(source.dir,'/../dir.csv',sep=''),sep=',',quote=F,row.names=F,col.names=F)