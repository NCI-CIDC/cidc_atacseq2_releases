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
#
# Program:  Snakefile
# Version:  0.1
# Author:   Sami R. Cherikh, Travis L. Jensen
# Purpose:  Main snakefile for workflow template.
# Input:    Sample fastq files stored in a cloud location (google cloud, aws)
# Output:   'sample_metadata.csv','rseqc/bam_qc_parsed.tab', 'rseqc/bam_gc_parsed.tab'
#############################################################################################################

## Import modules
import shutil
import logging as _logging
import psutil
import os
from pathlib import Path
from box import Box
import yaml

wildcard_constraints:
    sample="[^_]+"


##############################
#       CONFIGURATION        #
##############################
# Specify YAML config file location
configfile: "config/config.yaml"
# Set workflow working (output) dir
workdir: config['predir']

# Directories
# source dir
SOURCEDIR  = config["srcdir"]
# working output dir
PREDIR     = config["predir"]
# analysis data results and reporting within working dir
DATADIR    = config["datadir"]
REPDIR     = config["repdir"]


def create_path_accessor(prefix: Path = Path(PREDIR)) -> Box:
    """Create a Box to provide '.' access to hierarchy of paths"""
    data = yaml.load(Path(config["file_layout"]).open(), Loader=yaml.SafeLoader)
    paths = {}
    for directory in data.keys():
        paths[directory] = {}
        for file_alias, file_name in data[directory].items():
            p = str(prefix / directory / file_name)
            paths[directory][file_alias] = str(p)
    return Box(paths, frozen_box=True)

## create file accessor
paths = create_path_accessor()


## Use the source dir to import helper modules
sys.path.append(SOURCEDIR+'/python')
import trimadapters
import getfile  
import putfile  
import utils 
import time



##############################
#       SET GLOBAL VARS      #
##############################
# Workflow info
## number of cores dedicated to run
NCORES  = int(config["ncores"])
## Logging config
LOG_FILE  = config["log_file"]
## init sub folders
SUBDIRS  = config["subdirs"]


# Reference genome ftp
GENOME_FA_FTP = config["genome_fa_ftp"]
GENOME_GTF_FTP = config["genome_gtf_ftp"]
GENOME_BWA_FTP = config["genome_gtf_ftp"]
GENOME_BLACKLIST_URL = config["genome_blacklist_url"]
GENOME_DHS_GC = config["genome_dhs_gc"]


# Sample info
## List of samples to process
SAMID        = utils.toList(config["samid"])
## List of input files
FASTQ_1      = utils.toList(config["fastq1"])
FASTQ_2      = utils.toList(config["fastq2"])
## Adapter sequences
FP_ADAPTERS   = [x.strip() for x in utils.toList(config["fp_adapter_seq"])]
TP_ADAPTERS   = [x.strip() for x in utils.toList(config["tp_adapter_seq"])]
## Set single or paired end
if (FASTQ_2 != ['']):
    ENDS  = ['1','2']
else:
    ENDS  = ['1']
## Determine whether adapters should be trimmed or not
TRIM_FP = sum([x == 'NA'  for x in FP_ADAPTERS]) == 0
TRIM_TP = sum([x == 'NA'  for x in TP_ADAPTERS]) == 0
TRIM_ADAPTERS_OUTPUT = '.fastq.gz' if (TRIM_FP or TRIM_TP) else '.skipped'


# Preprocessing options
## Configure quality trimming level
QUAL_CUTOFF = config["quality_trim"]


# Peak calling options
## Macs peak caller mode
PEAK_MODE = str(config["peak_mode"])
## Macs peak caller ext size for broad calls
EXTSIZE = str(config["ext_size"])
## Macs peak caller arbitrary shift in bp here
SHIFT = str(config["shift"])
## Macs min qvalue for sig peak
PEAK_FDR = str(config["peak_fdr"])
## Create peak mode str for macs
if PEAK_MODE == 'broad':
    PEAK_MODE_STR = '--broad --broad-cutoff %s --nomodel --shift %s --extsize %s' % (PEAK_FDR, SHIFT, EXTSIZE)
else:
    PEAK_MODE_STR = '--nomodel --shift %s --extsize %s' % (SHIFT, EXTSIZE)


# Genome track options
TRACK_REGION = str(config["track_region"])


# Download input / upload output info
DECRYPT_PASS = config["decrypt_pass"]
ARCHIVE = config["archive_bucket"]
DOARCHIVE = len(ARCHIVE) != 0
HASH = config["hash"]


## Program names/locations referenced in some python scripts
CLOUD        = config["cloud_prog"]
OPENSSL      = config["openssl_prog"]
FASTQDUMP    = config["fastqdump_prog"]
CUTADAPT     = config["cutadapt_prog"]
TRIMMOMATIC  = config["trimmomatic_prog"]


## Set up logging to log file
_logging.basicConfig(level=_logging.INFO,
                    format='[cidc-atac] %(asctime)s - %(levelname)s - %(message)s',
                    datefmt='%m/%d/%Y %I:%M:%S %p',
                    handlers=[_logging.FileHandler(LOG_FILE)])



################################
#     DEFINE TARGET OUTPUT     #
################################
OUTPUT = [expand(paths.rseqc.bamqc_txt, sample=SAMID),
          expand(paths.rseqc.bamgc_txt, sample=SAMID),
          expand(paths.fastqc.targz, sample=SAMID),
          expand(paths.cnv.csv, sample=SAMID),
          expand(paths.peak.bw, sample=SAMID),
          expand(paths.chipqc.csv, sample=SAMID),
          expand(paths.peak.annot_tab, sample=SAMID),
          expand(paths.ptw.gobp, sample=SAMID),
          expand(paths.ptw.kegg, sample=SAMID),
          expand(paths.track.png, sample=SAMID)]

if PEAK_MODE == "narrow":
    OUTPUT.append(expand(paths.motif.narrow_peak, sample=SAMID))
    OUTPUT.append(expand(paths.motif.summit, sample=SAMID))
    OUTPUT.append(expand(paths.targets.narrowPeak_1k, sample=SAMID))
    OUTPUT.append(expand(paths.targets.narrowPeak_10k, sample=SAMID))
    OUTPUT.append(expand(paths.targets.narrowPeak_100k, sample=SAMID))
    OUTPUT.append(expand(paths.peak.dhs_narrowPeak, sample=SAMID))
    OUTPUT.append(expand(paths.peak.csv_narrowPeak, sample=SAMID))
else:
    OUTPUT.append(expand(paths.motif.broad_peak, sample=SAMID))
    OUTPUT.append(expand(paths.targets.broadPeak_1k, sample=SAMID))
    OUTPUT.append(expand(paths.targets.broadPeak_10k, sample=SAMID))
    OUTPUT.append(expand(paths.targets.broadPeak_100k, sample=SAMID))
    OUTPUT.append(expand(paths.peak.dhs_broadPeak, sample=SAMID))
    OUTPUT.append(expand(paths.peak.csv_broadPeak, sample=SAMID))



#########################################
#    Define any onstart or onsuccess    #
#########################################
onsuccess:
    ## Merge sample rseqc results into single result files
    merged_results = utils.mergeRSEQC(SOURCEDIR)

    ## Copy some results to analysis data dir
    [shutil.copy2(x, DATADIR) for x in merged_results]

    ## call ChIPQC on samples in mta and output results to analysis data dir - set last arg to one of the options [make-report, no-report]
    shell('Rscript --vanilla '+SOURCEDIR+'/r/init-chipqc-all-peak-metrics.r '+PREDIR+' '+DATADIR+' sample_metadata.csv make-report')

    ## add chipqc output to results to enc and archive if needed
    merged_results.append(['analysis/report/'+f for f in os.listdir('analysis/report')])
    merged_results.append('analysis/data/all_chipqc.csv.gz')

    ## knit rmarkdown html report
    shell('Rscript --vanilla '+SOURCEDIR+'/r/run-report.r '+SOURCEDIR+'/r/cidc_atac-report-slidy.Rmd '+PREDIR+' '+DATADIR+'/../report')
    merged_results.append('analysis/report/cidc_atac-report-slidy.html')

    ## Upload main results if needed
    if DOARCHIVE:
        [putfile.upload(file=x, destination=ARCHIVE, prog=CLOUD) for x in merged_results]

    shell("echo 'Pipeline complete!'")



################################
#   PASS OUTPUT TO all RULE    #
################################
rule all:
    input:
        OUTPUT
    


################################
#        PIPELINE RULES        #
################################
include: "./rules/initialization.smk"
include: "./rules/ingest.smk"
include: "./rules/mapping.smk"
include: "./rules/analysis.smk"
include: "./rules/format_peaks.smk"
include: "./rules/track.smk"
include: "./rules/motif.smk"
include: "./rules/targets.smk"
