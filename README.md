# cidc_atac

*ATAC-seq data processing pipeline* 

*This framework is currently
set up to perform the following across sample files: read pre-processing (adapter trimming and quality filtering), read QC, reference
genome alignment, obtaining alignment metrics, CNV analysis, peak calling with MACS using settings optimal for detection of chromatin-accessibility signatures, peak QC and quantification with ChIPQC. 
Currently a basic report is output containing some alignment and peak metrics using ChIPQC's builtin reporting feature.*

## Software
* See SOFTWARE.xlsx to see more info on required base software
* Software dependencies for rule execution is handled by conda and snakemake. See conda environment files in env folder for software specifications

### Configure config excel file (config.xlsx)
* This is where the sample metadata is set (most importantly sample IDs and google cloud location of input fastq files), as well as workflow and analysis directories and options
* Please reference the test sra config file config/config-sra-test.xlsx as a basic example that can be used with the template framework. This lists what is needed at minimum for execution: info for test samples, the workflow dirs (pre_dir; sub_dirs) and options, analysis data dir (report_dir)


### Working:
* Refinement of ATAC-seq specific steps, etc.


## Please reach out to Sami Cherikh (cherikhsr@nih.gov) with any questions or to report any issues concerning this framework.
