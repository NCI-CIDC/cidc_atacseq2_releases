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
# Program:  parse-rseqc-bam-gc-results.pl 
# Version:  v0.1
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:	Parse RSeqC read_GC.py results
# Input:	list of absolute file paths
# Output: 	N/A
#############################################################################################################

use strict;
use warnings;
use File::Basename;

my $isFirst=1;

my $header = "sample_id\tmin\tq1\tmedian\tmean\tq3\tmax";

while(<>) {
	chomp;
	my $file = $_;
	
	open(DATA, "<$file") or die "Couldn't open file $file, $!";
	
	if($isFirst) {
		print $header."\n";
		$isFirst=0;
	} 
	
	while(<DATA>){
		print;
	}
	close DATA;	
}