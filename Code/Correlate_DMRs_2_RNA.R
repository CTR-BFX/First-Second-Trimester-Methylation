#!/usr/local/bin/Rscript

#------------------------------------------------------------------------------
# EPIC Analysis        
# Data from GenomeScan EPIC Array 
#
# Link to publication
# TO ADD ONCE AVAILABLE
#
# Script available from:
# https://github.com/CTR-BFX/First-Second-Trimester-Methylation
#
#
# Analysis Performed by Russell S. Hamilton
# CTR Bioinformatics Facility
# Centre for Trophoblast Reseach, University of Cambridge, UK
# Copyright Russell S. Hamilton (rsh46@cam.ac.uk)
#
#------------------------------------------------------------------------------
# License Information
# This program is free software: you can redistribute it and/or modify  
# it under the terms of the GNU General Public License as published by  
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#------------------------------------------------------------------------------



library("DESeq2")
library("biomaRt")
library("ggplot2")
library("ggrepel")
library("ComplexHeatmap")
library("circlize")
library("dplyr")

baseDir <- "~/Documents/CTR-Data/CTR_EPIC/GitHub_First_Second"
setwd(baseDir)

Project      <- "CTR_EPIC.First_Second"

message("+-------------------------------------------------------------------------------")
message("+ Use ensEMBL Annotations")
message("+-------------------------------------------------------------------------------")

ensembl    <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description','chromosome_name', 'entrezgene_id', 'gene_biotype'), mart = ensembl)  
head(ensEMBL2id)
nrow(ensEMBL2id)


DMRs <- read.table("~/Documents/CTR-Data/CTR_EPIC/test.DMRs.srt.ann.bed", sep="\t", header=F, stringsAsFactors=T)
colnames(DMRs) <- c("chr", "start", "end", "name", "value", "p.val", "strand", "ensembl_gene_id", "dist2gene")
head(DMRs)
nrow(DMRs)

RNAs <- read.csv("~/Documents/CTR-Data/CTR_EPIC/CTR_gjb2_0001_STAR_shr_BlockSex_deseq2_results.edited.csv", header=T, stringsAsFactors=F)
RNAs$padj           <- as.numeric(RNAs$padj)
RNAs$log2FoldChange <- as.numeric(RNAs$log2FoldChange)
head(RNAs)

RNA.counts            <- read.csv("~/Documents/CTR-Data/CTR_EPIC/CTR_gjb2_0001_STAR_DESeq2_shrinkage___bulk_vs_scRNA_seq_BlockSex_deseq2_counts_norm.csv", header=T, stringsAsFactors=F)
RNA.counts$ensembl_gene_id <- RNA.counts$X

RNA.counts[,c(2:15)] <- log2(RNA.counts[,c(2:15)]+1)

RNA.counts$FirstMean  <- rowMeans(RNA.counts[,c(2:9)])
RNA.counts$SecondMean <- rowMeans(RNA.counts[,c(10:15)])

head(RNA.counts)

RNA.counts.meancentered           <- RNA.counts[,c(2:15)] - rowMeans(RNA.counts[,c(2:15)])   
RNA.counts.meancentered$FT        <- rowMeans(RNA.counts.meancentered[,c(2:9)])
RNA.counts.meancentered$ST        <- rowMeans(RNA.counts.meancentered[,c(10:15)])
rownames(RNA.counts.meancentered) <- RNA.counts$ensembl_gene_id
RNA.counts.meancentered$ensembl_gene_id <- rownames(RNA.counts.meancentered)

head(RNA.counts.meancentered)

#RNA.counts$meanCentered_ST <- log2(abs(RNA.counts.meancentered$ST)+1)
#RNA.counts$meanCentered_FT <- log2(abs(RNA.counts.meancentered$FT)+1)


#RNA.counts$meanCentered_FT_z <- RNA.counts$meanCentered_FT / sd( RNA.counts.meancentered[,c(2:9)] )

RNA.counts.meancentered <- RNA.counts.meancentered %>% mutate(FirstStd  = apply(.[(2:9)],  1,sd))
RNA.counts.meancentered <- RNA.counts.meancentered %>% mutate(SecondStd = apply(.[(10:15)],1,sd))

RNA.counts.meancentered$FT_Z <- RNA.counts.meancentered$FT / RNA.counts.meancentered$FirstStd
RNA.counts.meancentered$ST_Z <- RNA.counts.meancentered$ST / RNA.counts.meancentered$SecondStd


head(RNA.counts.meancentered)


#RNA.counts$meanCentered_ST_z <- RNA.counts$meanCentered_ST / sd( RNA.counts.meancentered[,c(10:15)] )



#if(RNA.counts.meancentered$FT < 0){ RNA.counts$meanCentered_FT = -1*(RNA.counts$meanCentered_FT)}
#if(RNA.counts.meancentered$ST < 0){ RNA.counts$meanCentered_ST = -1*(RNA.counts$meanCentered_ST)}


#RNA.counts$meanCentered_FT[is.nan(RNA.counts$meanCentered_FT)] <- 0
#RNA.counts$meanCentered_ST[is.nan(RNA.counts$meanCentered_ST)] <- 0

#head(RNA.counts)



DMRs.RNAs <- merge(DMRs, RNAs, by="ensembl_gene_id")
DMRs.RNAs <- DMRs.RNAs[!duplicated(DMRs.RNAs[,c("name")]), ]
DMRs.RNAs <- merge(DMRs.RNAs, RNA.counts.meancentered, by="ensembl_gene_id")

# Remove pseudo-genes HLA-L and HLA-H
DMRs.RNAs <- subset(DMRs.RNAs, external_gene_name!="HLA-L")
DMRs.RNAs <- subset(DMRs.RNAs, external_gene_name!="HLA-H")


nrow( subset(DMRs.RNAs, abs(log2FoldChange > 1)  & abs(value) > 0.25 & padj <= 0.05))

meth_thresh <- 0.2

pdf(paste0(baseDir, "/Figures/", Project, ".Correlation_DMR_RNA.pdf"),width=4.5,height=4.5, onefile=FALSE)
par(bg=NA)
ggplot(DMRs.RNAs, aes(x=value, y=log2FoldChange, label=external_gene_name)) +
  geom_vline(xintercept = 0,  colour="black", linetype='solid', alpha=.5) +
  geom_hline(yintercept = 0,  colour="black", linetype='solid', alpha=.5) +
  geom_vline(xintercept = meth_thresh,  colour="red", linetype='dashed') +
  geom_vline(xintercept = -meth_thresh, colour="red", linetype='dashed') +
  geom_hline(yintercept = 1,  colour="red", linetype='dashed') +
  geom_hline(yintercept = -1, colour="red", linetype='dashed') +
  
  geom_point(data=subset(DMRs.RNAs, abs(log2FoldChange) < 1 | abs(value) < meth_thresh | padj > 0.05), 
             colour='grey', size=1.5) +
  geom_point(data=subset(DMRs.RNAs, log2FoldChange > 1  & value < -meth_thresh & padj <= 0.05 & dist2gene < 1000), 
             colour='purple', size=1.5) +
  geom_point(data=subset(DMRs.RNAs, log2FoldChange > 1  & value > meth_thresh  & padj <= 0.05 & dist2gene < 1000), 
             colour='blue', size=1.5) +
  geom_point(data=subset(DMRs.RNAs, log2FoldChange < -1 & value > meth_thresh  & padj <= 0.05 & dist2gene < 1000), 
             colour='red', size=1.5) +
  geom_point(data=subset(DMRs.RNAs, log2FoldChange < -1 & value < -meth_thresh & padj <= 0.05 & dist2gene < 1000), 
             colour='darkgreen', size=1.5) +
  
  geom_text_repel(data=subset(DMRs.RNAs, abs(log2FoldChange) > 1 & abs(value) > meth_thresh & padj <= 0.05 & dist2gene < 1000), 
                  segment.alpha=0.6, size=2.75, force=10 ) +
  
  xlab(bquote("Methylation Difference (DMRs, "*beta~"value)")  )+
  ylab(bquote("Expression Difference (log"[2]~"Fold Change)") ) +
  scale_x_continuous(breaks=seq(-2,1,0.2)) +
  scale_y_continuous(breaks=seq(-3,4,1)) +
  #theme_bw() +
  theme_cowplot(12) +
  theme(axis.text=element_text(size=10), 
        axis.text.x = element_text(angle = -45, vjust = 0.25, hjust=0.3),
        axis.title=element_text(size=12,face="bold"),
        aspect.ratio=1)
dev.off()


head(DMRs.RNAs)


functionMakeHeatmaps <- function(MATRIX,R1,R2,R3,R4) {

  DMR.mat                 <- MATRIX[,c(6), drop=F] 
  rownames(DMR.mat)       <- paste0(MATRIX$external_gene_name ) #, "_", MATRIX$name )
  colnames(DMR.mat)       <- c("Methylation")
  
  RNA.mat                 <- MATRIX[,c("FT_Z", "ST_Z"), drop=F] 
  rownames(RNA.mat)       <- paste0(MATRIX$external_gene_name ) #, "_", MATRIX$name )
  #RNA.mat[,2] <- log2(RNA.mat[,2]+1)
  #RNA.mat[,3] <- log2(RNA.mat[,3]+1)
  
  print( min(RNA.mat[,1]))
  print( max(RNA.mat[,2]))
  
  colnames(RNA.mat)       <- c("First", "Second")
  
  print(head(RNA.mat))
  
  col_DMR = colorRamp2(c(R1, 0, R2), c("darkgreen", "grey95", "purple"))
  col_RNA = colorRamp2(c(R3, 0, R4), c("darkgreen", "grey95", "purple"))
 # col_RNA = colorRamp2( c(-20,-10,0,10,20), c("darkgreen", "green3", "grey95", "slateblue", "darkorchid4"), space = "RGB") 
  
  hm.height <- nrow(RNA.mat)*.25
  
  ht1 = Heatmap(DMR.mat, name="DMR", col =col_DMR, cluster_columns=F, show_row_names=F, width=unit(1,"cm"))
  ht2 = Heatmap(RNA.mat, name="RNA", col =col_RNA, cluster_columns=F, show_row_names=T, width=unit(1,"cm"),
                row_names_gp = gpar(fontsize = 6))
  ht.both <- ht2 + ht1
  
  
 # draw(ht.both, row_title = "Three heatmaps, row title", row_title_gp = gpar(col = "red"),
#       column_title = "Three heatmaps, column title", column_title_gp = gpar(fontsize = 16))
  return(ht.both)
}

extra <- 1.0
rowsize <- 0.15





DMRs.RNAs.purple <- subset(DMRs.RNAs, log2FoldChange > 1 & value < -0.25 & padj <= 0.05 & dist2gene < 1000)
hm.purple        <- functionMakeHeatmaps(DMRs.RNAs.purple,-1,1,-1.5,1.5)
pdf(paste(Project, "_Heatmap_purple.pdf", sep=""),width=4,height=(extra+(nrow(DMRs.RNAs.purple)*rowsize)), onefile=FALSE)
par(bg=NA)
hm.purple
dev.off()

DMRs.RNAs.green <- subset(DMRs.RNAs, log2FoldChange < -1 & value < -0.25 & padj <= 0.05 & dist2gene < 1000)
hm.green        <- functionMakeHeatmaps(DMRs.RNAs.green,-1,1,-1.5,1.5)
pdf(paste(Project, "_Heatmap_green.pdf", sep=""),width=4,height=(extra+(nrow(DMRs.RNAs.green)*rowsize)), onefile=FALSE)
par(bg=NA)
hm.green
dev.off()

DMRs.RNAs.red <- subset(DMRs.RNAs, log2FoldChange < -1 & value > 0.25 & padj <= 0.05 & dist2gene < 1000)
hm.red        <- functionMakeHeatmaps(DMRs.RNAs.red,-1,1,-1.5,1.5)
pdf(paste(Project, "_Heatmap_red.pdf", sep=""),width=4,height=(extra+(nrow(DMRs.RNAs.red)*rowsize)), onefile=FALSE)
par(bg=NA)
hm.red
dev.off()

DMRs.RNAs.blue <- subset(DMRs.RNAs, log2FoldChange > 1 & value > 0.25 & padj <= 0.05 & dist2gene < 1000)
hm.blue        <- functionMakeHeatmaps(DMRs.RNAs.blue,-1,1,-1.5,1.5)
pdf(paste(Project, "_Heatmap_blue.pdf", sep=""),width=4,height=(extra+(nrow(DMRs.RNAs.blue)*rowsize)), onefile=FALSE)
par(bg=NA)
hm.blue
dev.off()





enrichR.list.purple <- ensEMBL2id[ensEMBL2id$ensembl_gene_id %in% DMRs.RNAs.purple$ensembl_gene_id,]
enrichR.list.purple[,c("ensembl_gene_id", "external_gene_name", "entrezgene_id")]

enrichR.list.green <- ensEMBL2id[ensEMBL2id$ensembl_gene_id %in% DMRs.RNAs.green$ensembl_gene_id,]
enrichR.list.green[,c("ensembl_gene_id", "external_gene_name", "entrezgene_id")]

enrichR.list.red <- ensEMBL2id[ensEMBL2id$ensembl_gene_id %in% DMRs.RNAs.red$ensembl_gene_id,]
enrichR.list.red[,c("ensembl_gene_id", "external_gene_name", "entrezgene_id")]

enrichR.list.blue <- ensEMBL2id[ensEMBL2id$ensembl_gene_id %in% DMRs.RNAs.blue$ensembl_gene_id,]
enrichR.list.blue[,c("ensembl_gene_id", "external_gene_name", "entrezgene_id")]



cat( paste( enrichR.list.purple[,c("ensembl_gene_id", "external_gene_name", "entrezgene_id")]$external_gene_name, collapse='\n' ) )

xx <- DMRs.RNAs.purple.summary[, c("external_gene_name", "value")]
rownames(xx) <- xx$external_gene_name
print( xx[,c("value"), drop=F] )


write.table(xx[,c("value"), drop=F], paste0(baseDir, "/Figures/", Project, ".Correlation_DMR_RNA.purple.txt"), 
            append = FALSE, sep = ",", row.names=T, col.names=F, quote=F)

#
#
#

DMRs.RNAs.purple.summary           <-  DMRs.RNAs.purple[, c("ensembl_gene_id", "external_gene_name", 
                                                      "value", "p.val", "log2FoldChange", "padj")]
DMRs.RNAs.purple.summary$quadrant  <- "purple"
colnames(DMRs.RNAs.purple.summary) <- c("ensembl_gene_id", "external_gene_name", "meth_diff", 
                                        "meth.p.val", "log2FoldChange", "padj", "quadrant")


DMRs.RNAs.green.summary            <-  DMRs.RNAs.green[, c("ensembl_gene_id", "external_gene_name", 
                                                           "value", "p.val", "log2FoldChange", "padj")]
DMRs.RNAs.green.summary$quadrant   <- "green"
colnames(DMRs.RNAs.green.summary)  <- c("ensembl_gene_id", "external_gene_name", "meth_diff", 
                                        "meth.p.val", "log2FoldChange", "padj", "quadrant")


DMRs.RNAs.red.summary              <-  DMRs.RNAs.red[, c("ensembl_gene_id", "external_gene_name", 
                                                         "value", "p.val", "log2FoldChange", "padj")]
DMRs.RNAs.red.summary$quadrant     <- "red"
colnames(DMRs.RNAs.red.summary)    <- c("ensembl_gene_id", "external_gene_name", "meth_diff", 
                                        "meth.p.val", "log2FoldChange", "padj", "quadrant")


DMRs.RNAs.blue.summary             <-  DMRs.RNAs.blue[, c("ensembl_gene_id", "external_gene_name", 
                                                          "value", "p.val", "log2FoldChange", "padj")]
DMRs.RNAs.blue.summary$quadrant    <- "blue"
colnames(DMRs.RNAs.blue.summary)   <- c("ensembl_gene_id", "external_gene_name", "meth_diff", 
                                        "meth.p.val", "log2FoldChange", "padj", "quadrant")



DMRs.RNAs.summary <- rbind( DMRs.RNAs.purple.summary, DMRs.RNAs.green.summary, DMRs.RNAs.red.summary, DMRs.RNAs.blue.summary  )

head(DMRs.RNAs.summary)
write.csv(DMRs.RNAs.summary, file=paste0(Project, "_Methylation_RNA_Correlations", ".csv"))

