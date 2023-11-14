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
# Program:  init-peak-annot-ptw-enr.r
# Version:  0.1
# Author:   Sami R. Cherikh
# Purpose:  Annot peaks with chipseeker R package and perform func enr ptw analyis with clusterProfiler R package
# Input:    <COMMAND LINE INPUTS>
# Output:
#############################################################################################################

libs = c('ChIPseeker','clusterProfiler','org.Hs.eg.db');
invisible(suppressPackageStartupMessages(lapply(libs, require, character.only=T)))


#############################
##      cmd line args      ##
#############################
## get path to config.xlsx and source dir command line arguments
args = commandArgs(trailingOnly = TRUE)
if(length(args)>0) {
	predir = args[1]
	sample = args[2]
	annot.gtf.file = args[3]
	peak = args[4]
}else{
    stop("ERROR: No predir and/or needed sample info supplied.")
}


## read in genome annot gtf file
txdb = GenomicFeatures::makeTxDbFromGFF(annot.gtf.file)


##############
#
# Annot
#
#############
## annotate peaks with chipseeker
annot.obj = annotatePeak(peak, TxDb=txdb, annoDb="org.Hs.eg.db", tssRegion=c(-1000, 1000))
annot.stat = data.frame(annot.obj@annoStat)
annot.tab = data.frame(annot.obj@anno)

## write to file
write.table(annot.stat, paste0(predir,'/peak/',sample,'_peaks_annot_stat.csv'), sep=',', quote=F, row.names=F, col.names=T)
write.table(annot.tab, paste0(predir,'/peak/',sample,'_peaks_annot.csv'), sep=',', quote=F, row.names=F, col.names=T)

## barplot of genomic feature representation
pdf(file = paste0(predir,'/peak/',sample,'_peaks_annot_bar.pdf'))
plotAnnoBar(annot.obj, xlab=sample)
garb=dev.off()

## barplot of distribution of TF-binding loci relative to TSS
pdf(file = paste0(predir,'/peak/',sample,'_peaks_annot_bar_tss.pdf'))
plotDistToTSS(annot.obj, title="Distribution of transcription factor-binding loci \n relative to TSS", xlab=sample)
garb=dev.off()


##############
#
# Func enr
#
#############
## go bp ptw enr
ego = enrichGO(gene = annot.tab$ENTREZID,
                    keyType = "ENTREZID",
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05,
                    readable = TRUE)
## write to file
ego.summary = data.frame(ego)
write.table(ego.summary, paste0(predir,'/ptw/',sample,'_gobp_ptw.csv'), sep=',', quote=F, row.names=F, col.names=T)
## dot plot
pdf(file = paste0(predir,'/ptw/',sample,'_gobp_ptw_dot.pdf'))
dotplot(ego, showCategory=20, font.size=10, title = paste0("GO BP Enrichment Analysis - ",sample))
garb=dev.off()


## kegg ptw enr
ekegg = enrichKEGG(gene = annot.tab$ENTREZID,
                 organism = 'hsa',
                 pvalueCutoff = 0.05)
## write to file
ekegg.summary = data.frame(ekegg)
write.table(ekegg.summary, paste0(predir,'/ptw/',sample,'_kegg_ptw.csv'), sep=',', quote=F, row.names=F, col.names=T)
## dot plot
pdf(file = paste0(predir,'/ptw/',sample,'_kegg_ptw_dot.pdf'))
dotplot(ekegg, showCategory=20,  font.size=10, title = paste0("KEGG Enrichment Analysis - ",sample))
garb=dev.off()