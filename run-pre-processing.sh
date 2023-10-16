#!/bin/bash
#############################################################################################################
# pipeline_template
# 
# This program is free software that contains third party software subject to various licenses, 
# namely, the GNU General Public License version 3 (or later), the GNU Affero General Public License 
# version 3 (or later), and the LaTeX Project Public License v.1.3(c). A list of the software contained 
# in this program, including the applicable licenses, can be accessed here: 
#
# 
# You can redistribute and/or modify this program, including its components, only under the terms of 
# the applicable license(s).  
#
# This program is distributed in the hope that it will be useful, but "as is," WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#
# Program:  run-pre-processing.sh
# Version:  v0.1
# Author:   Sami R. Cherikh, Travis L. Jensen
# Purpose:  script to run pre-processing only steps
# Input:    N/A
# Output:   N/A
#############################################################################################################

## Locate Source Directory from root
SRCDIR="$(cd `dirname $0` && pwd)/source";

## add config script for fastq-dump to throw the cache in /tmp.
## otherwise .sra cache files are stored on the home directory
#mkdir $HOME/.ncbi -p;
#echo '/repository/user/main/public/root = "/tmp"' > $HOME/.ncbi/user-settings.mkfg

## Pass command line arguments to python program that calls snakemake
#python3 $SRCDIR/python/init-workflow --run-preprocessing $@ --config <path-to-config-file> --threads <vCPU-cores-dedicated-to-pipeline>
python3 $SRCDIR/python/init-workflow --run-preprocessing $@ --config $SRCDIR/../config/config-sra-test-subset.xlsx --threads 16 --continue-run

## Example command using the test config file with args that are useful when testing or debugging
#python3 $SRCDIR/python/init-workflow --run-preprocessing $@ --config $SRCDIR/../config/config-sra-test-subset.xlsx --threads 16 --save-int-local-files --continue-run --dryrun


