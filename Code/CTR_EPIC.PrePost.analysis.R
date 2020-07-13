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

#BiocManager::install("ChAMP")


library("ggplot2")
library('gplots')
library("ggforce")
library('minfi')
library('lumi')
library('methylumi')
library('limma')
#library('VennDiagram');
require('FDb.InfiniumMethylation.hg19')
require('IlluminaHumanMethylationEPICmanifest')
library('ChAMP')
library("data.table")
library("cowplot")
library("reshape")
library("ggrepel")
library("ggdendro")
library("useful")
library("MASS")
library("viridis")
library("scales")
library("biomaRt")
require("plyr")
library("ggalt")
library("ggridges")
library("dplyr")
library("tidyverse")

message("+-------------------------------------------------------------------------------")
message("+ Set the colours for first and second trimesters")
message("+-------------------------------------------------------------------------------")

col_1st <- "firebrick2"
col_2nd <- "steelblue3"


message("+-------------------------------------------------------------------------------")
message("+ Load in the EPIC probe manifest")
message("+-------------------------------------------------------------------------------")

data(probe.features.epic)
str(probe.features)
probe.features$probe_id <- rownames(probe.features)
head(probe.features)


probe.locations.epic <- probe.features[, c("CHR","MAPINFO", "Strand")]
probe.locations.epic$probe_id <- rownames(probe.locations.epic)
head(probe.locations.epic)


message("+-------------------------------------------------------------------------------")
message("+ Set the directories and baseDir for the IDAT files")
message("+-------------------------------------------------------------------------------")


setwd("/Users/rhamilto/Documents/CTR-Data/CTR_EPIC/")
baseDir      <- getwd()
baseDirIDATs <- sprintf("%s/IDATs",baseDir)

Project      <- "CTR_EPIC.First_Second"
GitHubDir    <- "GitHub_First_Second"


message("+-------------------------------------------------------------------------------")
message("+ Load the IDATS into CHAMP for processing")
message("+-------------------------------------------------------------------------------")


myLoad <- champ.load(directory = baseDirIDATs, methValue='B', filterBeads=TRUE, filterXY=TRUE, 
                     detPcut = 0.01, arraytype="EPIC")
head(myLoad$pd)

sampleTable.champ     <- as.data.frame(myLoad$pd)

#
# Sex determined from paired RNA-Seq experiment
#
sampleTable.champ$Sex <- sampleTable.champ$Sample_Name
sampleTable.champ$Sex <- gsub("first_63_oxBS",  "M", sampleTable.champ$Sex)
sampleTable.champ$Sex <- gsub("first_64_oxBS",  "F", sampleTable.champ$Sex)
sampleTable.champ$Sex <- gsub("first_65_oxBS",  "F", sampleTable.champ$Sex)
sampleTable.champ$Sex <- gsub("first_66_oxBS",  "M", sampleTable.champ$Sex)
sampleTable.champ$Sex <- gsub("first_67_oxBS",  "M", sampleTable.champ$Sex)
sampleTable.champ$Sex <- gsub("first_69_oxBS",  "F", sampleTable.champ$Sex)
sampleTable.champ$Sex <- gsub("first_70_oxBS",  "F", sampleTable.champ$Sex)

sampleTable.champ$Sex <- gsub("second_71_oxBS", "M", sampleTable.champ$Sex)
sampleTable.champ$Sex <- gsub("second_72_oxBS", "M", sampleTable.champ$Sex)
sampleTable.champ$Sex <- gsub("second_73_oxBS", "F", sampleTable.champ$Sex)
sampleTable.champ$Sex <- gsub("second_74_oxBS", "M", sampleTable.champ$Sex)
sampleTable.champ$Sex <- gsub("second_75_oxBS", "M", sampleTable.champ$Sex)
sampleTable.champ$Sex <- gsub("second_76_oxBS", "M", sampleTable.champ$Sex)
sampleTable.champ$Sex <- gsub("second_77_oxBS", "X", sampleTable.champ$Sex)



message("+-------------------------------------------------------------------------------")
message("+ Run CHAMP QC and Normalisation")
message("+-------------------------------------------------------------------------------")


myQC   <- champ.QC()

myNorm <- champ.norm()


message("+-------------------------------------------------------------------------------")
message("+ Plot PCA and dendrogram")
message("+-------------------------------------------------------------------------------")

customPCADendro <- function(ProjectTitle, myNorm, TOPNUM, sampleTable.champ) {
  
  myNorm.df <- as.data.frame(myNorm)
  rv        <- rowVars(as.matrix(myNorm.df))
  select    <- order(rv, decreasing = TRUE)[seq_len(min(TOPNUM, length(rv)))]
  pca       <- prcomp(t(myNorm.df[select, ]))
  pc1var    <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var    <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab    <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab    <- paste0("PC2 (",as.character(pc2var),"%)")
  scores    <- data.frame(sampleName=sampleTable.champ$Sample_Name, pca$x, 
                          Sample_Group=sampleTable.champ$Sample_Group, Sex=sampleTable.champ$Sex)

  scores$sampleName   <- gsub("_oxBS", "", scores$sampleName)
  scores$Sample_Group <- gsub("_oxBS", "", scores$Sample_Group)
  scores$Label        <- scores$Sample_Group
  scores$Label        <- gsub("first",  "First Trimester", scores$Label)
  scores$Label        <- gsub("second", "Second Trimester", scores$Label)
  
  plt.pca <- ggplot(scores, aes(x = PC1, y = PC2, fill=Sample_Group, colour=Sample_Group, shape=Sex) ) +
             geom_mark_ellipse(aes(fill = NULL, group=Sample_Group, color=Sample_Group, label=Label), 
                               alpha=0.1, label.fontsize = 10, label.buffer = unit(1, 'mm')) +
          #   geom_text_repel(aes(label=sampleName), show.legend = FALSE, size=2, colour="black") +
             geom_point(size = 4, alpha=0.75 ) + 
             xlab(pc1lab) + 
             ylab(pc2lab) +
             scale_colour_manual(name="", values = c("first"="firebrick2", "second"="steelblue3"), 
                                          labels=c("first"="First Trimester", "second"="Second Trimester"), guide=F) +
             scale_fill_manual(  name="", values = c("first"="firebrick2", "second"="steelblue3"), 
                                          labels=c("first"="First Trimester", "second"="Second Trimester"), guide=F) +
             expand_limits(x = c(-5, 5), y = c(-4,4)) +
             scale_x_continuous(breaks = seq(-5,5, 2.5)) +
             scale_y_continuous(breaks = seq(-4,4, 2)) +
             #theme_bw() +
             theme_cowplot(12) +
             theme(axis.text=element_text(size=10), 
                   legend.title=element_text(size=10),
                   legend.text=element_text(size=8),
                   axis.title=element_text(size=10),
                   aspect.ratio=1, legend.position="right") 
  
  loadings        <- as.data.frame(pca$rotation)
  topX            <- 20
  pca.1           <-  loadings[ order(loadings$PC1,decreasing=TRUE), ]
  pca.1.topX      <-  pca.1[c(1:topX),]
  pca.2           <-  loadings[ order(loadings$PC2,decreasing=TRUE), ]
  pca.2.topX      <-  pca.2[c(1:topX),]

  dd.row          <- as.dendrogram(hclust(dist(t( myNorm.df[select, ]))))
  ddata           <- dendro_data(dd.row)
  labs            <- label(ddata)
  labs$group      <- labs$label
  labs$group      <- gsub("first.*",  "FT", labs$group)
  labs$group      <- gsub("second.*", "ST", labs$group)

  newLab          <- label(ddata)
  newLab$label    <- gsub("_oxBS", "", newLab$label)
  newLab$label    <- gsub("first", "FT", newLab$label)
  newLab$label    <- gsub("second", "ST", newLab$label)
  
  plt.dendro <- ggplot(segment(ddata)) +
                geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
                geom_text(data=newLab,aes(label=label, x=x, y=-0.25, angle=-90, colour=labs$group), hjust=0) +
                scale_colour_manual(name="", guide=FALSE,
                                    values = c("FT"="firebrick2", "ST"="steelblue3")) +
                ylim(-1.5, max(ddata$segments$y)) + ylab("distance") + xlab("") +
                ggtitle(paste0("Dendrogram Top ", TOPNUM, " MV")) +
                theme(axis.text.x = element_blank(),  axis.ticks.x = element_blank()) 
  
  myNorm.df <- ""
  
  return(list(plt.pca, plt.dendro))
}

plt.QC <- customPCADendro(Project, myNorm, 500, sampleTable.champ)

pdf(paste0(baseDir, "/", GitHubDir, "/Figures/", Project, ".PCA_oxBS.pdf"),width=4.5,height=4, onefile=FALSE)
par(bg=NA)
plt.QC[[1]]
dev.off()

pdf(paste0(baseDir, "/", GitHubDir, "/Figures/",Project, ".Dendro_oxBS.pdf"),width=5,height=5, onefile=FALSE)
par(bg=NA)
plt.QC[[2]]
dev.off()






message("+-------------------------------------------------------------------------------")
message("+ Calculate DMPs using true mC (oxBS)")
message("+-------------------------------------------------------------------------------")

first_oxBS_vs_second_oxBS_DMP     <- champ.DMP(beta=myNorm, pheno=myLoad$pd$Sample_Group, 
                                                compare.group=c("first_oxBS", "second_oxBS"), arraytype="EPIC")
first_oxBS_vs_second_oxBS_DMP.flt <- first_oxBS_vs_second_oxBS_DMP[[1]][first_oxBS_vs_second_oxBS_DMP[[1]]$deltaBeta>0,]

head(first_oxBS_vs_second_oxBS_DMP.flt)
nrow(first_oxBS_vs_second_oxBS_DMP.flt)

write.csv(first_oxBS_vs_second_oxBS_DMP.flt, 
          file=paste0(baseDir, "/", GitHubDir, "/Data/", Project, "_first_vs_second_DMPs_using_oxBS", ".csv"))

message("+-------------------------------------------------------------------------------")
message("+ Calculate DMRs using true mC (oxBS)")
message("+-------------------------------------------------------------------------------")

# bumpHunter
first_oxBS_vs_second_oxBS_DMR <- champ.DMR(beta=myNorm, pheno=myLoad$pd$Sample_Group, 
                                           compare.group=c("first_oxBS", "second_oxBS"), arraytype="EPIC")

head(first_oxBS_vs_second_oxBS_DMR$BumphunterDMR)
nrow(first_oxBS_vs_second_oxBS_DMR$BumphunterDMR)

write.csv(first_oxBS_vs_second_oxBS_DMR$BumphunterDMR, 
          file=paste0(baseDir, "/", GitHubDir, "/Data/", Project, "_first_vs_second_DMRs_using_oxBS", ".csv"))


# DMRcate
first_oxBS_vs_second_oxBS_DMRcate <- champ.DMR(beta=myNorm, pheno=myLoad$pd$Sample_Group, method="DMRcate",
                                               compare.group=c("first_oxBS", "second_oxBS"), arraytype="EPIC")

head(first_oxBS_vs_second_oxBS_DMRcate$DMRcateDMR)
nrow(first_oxBS_vs_second_oxBS_DMRcate$DMRcateDMR)

message("+-------------------------------------------------------------------------------")
message("+ Calculate DMBs using true mC (oxBS)")
message("+-------------------------------------------------------------------------------")


first_oxBS_vs_second_oxBS_DMB <- champ.Block(beta=myNorm,pheno=myLoad$pd$Sample_Group, arraytype="EPIC")

head(first_oxBS_vs_second_oxBS_DMB$Block)
nrow(first_oxBS_vs_second_oxBS_DMB$Block)

write.csv(first_oxBS_vs_second_oxBS_DMB$Block, 
          file=paste0(baseDir, "/", GitHubDir, "/Data/", Project, "_first_vs_second_DMBs_using_oxBS", ".csv"))






message("+-------------------------------------------------------------------------------")
message("+ Annotate the genomic positions")
message("+-------------------------------------------------------------------------------")

threshold <- 0.8

data(probe.features.epic)
str(probe.features)
probe.features$probe_id <- rownames(probe.features)
head(probe.features)
nrow(probe.features)

myNorm.probes.df                <- as.data.frame(myNorm)
myNorm.probes.df$probe_id       <- rownames(myNorm.probes.df)
myNorm.probes.ann               <- merge(myNorm.probes.df, probe.features, by="probe_id")
myNorm.probes.ann$Strand        <- gsub("F", "1", myNorm.probes.ann$Strand)
myNorm.probes.ann$Strand        <- gsub("R", "-1", myNorm.probes.ann$Strand)

myNorm.probes.ann$promoter      <- myNorm.probes.ann$feature
myNorm.probes.ann$promoter      <- gsub("TSS1500", "promoter", myNorm.probes.ann$promoter)
myNorm.probes.ann$promoter      <- gsub("TSS200",  "promoter", myNorm.probes.ann$promoter)

myNorm.probes.ann$first_oxBS_meanBeta  <- rowMeans( myNorm.probes.ann[,grep("first_.*_oxBS",colnames(myNorm.probes.ann))] )
myNorm.probes.ann$second_oxBS_meanBeta <- rowMeans( myNorm.probes.ann[,grep("second_.*_oxBS",colnames(myNorm.probes.ann))] )

head(myNorm.probes.ann)
nrow(myNorm.probes.ann)





message("+-------------------------------------------------------------------------------")
message("+ Per Sample Features")
message("+-------------------------------------------------------------------------------")


functionGetPerSampleFeature <- function(myNorm, Feature, SubFeature, sampleGroup) {
  
  myNorm           <- melt(myNorm[, c("probe_id", Feature, colnames(myNorm)[grep("meanBeta",colnames(myNorm))] )])
  colnames(myNorm) <- c("probe_id", "Feature","variable","value")
  myNorm           <- subset(myNorm, Feature==SubFeature)
  print( head(myNorm))
  return(myNorm) 
}

myNorm.probes.ann.summary           <- myNorm.probes.ann
rownames(myNorm.probes.ann.summary) <- myNorm.probes.ann.summary$probe_id
myNorm.probes.ann.summary           <- myNorm.probes.ann.summary[, c("CHR","MAPINFO", "Strand", 
                                                                     "first_oxBS_meanBeta", "second_oxBS_meanBeta")]
head(myNorm.probes.ann.summary)
nrow(myNorm.probes.ann.summary)
write.table(myNorm.probes.ann.summary, col.names=T, row.names=T, quote=F, sep=",",
            file=paste0(baseDir, "/", GitHubDir, "/Data/", Project, ".myNorm.probes.ann.summary", ".csv"))



message("+-------------------------------------------------------------------------------")
message("+ Plot mean beta value correlation plot 1st vs 2nd")
message("+-------------------------------------------------------------------------------")


meth_corr <- myNorm.probes.ann[, c("probe_id", "cgi", "gene", "feature", "feat.cgi", "promoter", colnames(myNorm.probes.ann)[grep("meanBeta",colnames(myNorm.probes.ann))] )]

threshold <- 0.2
plt.meth.corr <- ggplot(data=meth_corr, aes(x=first_oxBS_meanBeta, y=second_oxBS_meanBeta, 
                                            colour=abs(first_oxBS_meanBeta-second_oxBS_meanBeta), label=gene)) +
                 geom_point(data=subset(meth_corr, abs(first_oxBS_meanBeta-second_oxBS_meanBeta) >= 0 ), 
                            size=0.25, alpha=0.75) +
                 scale_colour_gradient2(name="Methylation\nDifference", 
                                        low="grey", mid="grey", high="darkblue", midpoint=0.02) +
                 geom_abline(intercept = 2*threshold, slope = 1, linetype="dashed") +
                 geom_abline(intercept = threshold, slope = 1, linetype="dashed") +
                 geom_abline(intercept = 0, slope = 1, linetype="dashed") +
                 geom_abline(intercept = -(threshold), slope = 1, linetype="dashed") +
                 geom_abline(intercept = -(2*threshold), slope = 1, linetype="dashed") +
            #     geom_label_repel(data=subset(meth_corr, 
            #                                  (first_oxBS_meanBeta-second_oxBS_meanBeta) <= -2*threshold), 
            #                      force=10, size=2.5, colour="black", nudge_x=-0.05, nudge_y=0.05 ) +
           #      geom_label_repel(data=subset(meth_corr, 
           #                                   (first_oxBS_meanBeta-second_oxBS_meanBeta) >=  2*threshold), 
           #                       force=10, size=2.5, colour="black", nudge_x=0.05,  nudge_y=-0.05) +
                 scale_x_continuous(name=bquote("First Trimester ("*beta~"value)"),  
                                    limits=c(0,1.01), expand=c(0,0), breaks=seq(0,1,0.2)) +
                 scale_y_continuous(name=bquote("Second Trimester ("*beta~"value)"), 
                                    limits=c(0,1.01), expand=c(0,0), breaks=seq(0,1,0.2)) +
  xlab(bquote("Methylation Difference (DMRs, "*beta~"value)")  )+
                 coord_fixed() +
                 #theme_bw()
                 theme_cowplot(12) +
                 theme(axis.text=element_text(size=10), 
                       #axis.text.x = element_text(angle = -45, vjust = 0.25, hjust=0.3),
                       axis.title=element_text(size=12,face="bold"),
                       aspect.ratio=1)


pdf(paste0(baseDir, "/", GitHubDir, "/Figures/", Project, ".MethylationCorrelation.pdf"),width=5,height=4.0, onefile=FALSE)
par(bg=NA)
plt.meth.corr
dev.off()



message("+-------------------------------------------------------------------------------")
message("+ Plot methyaltion for specific genomic features")
message("+-------------------------------------------------------------------------------")

meth_corr.feature <- myNorm.probes.ann[, c("probe_id", "feature", 
                                           colnames(myNorm.probes.ann)[grep("meanBeta",colnames(myNorm.probes.ann))] )]
head(meth_corr.feature)
nrow(meth_corr.feature)

meth_corr.feature.m <- melt(meth_corr.feature, id=c("probe_id","feature"))

plt.features <- ggplot(data=meth_corr.feature.m, aes(x=feature, y=value, fill=variable)) +
                geom_boxplot(aes(fill=variable), outlier.alpha = 0.0) +
                scale_fill_manual(name="", values = c("second_oxBS_meanBeta"=col_2nd, "first_oxBS_meanBeta"=col_1st), 
                                  limits=c("first_oxBS_meanBeta","second_oxBS_meanBeta"), 
                                  labels=c("First Trimester","Second Trimester")) +
                xlab("") +
                scale_y_continuous(name=bquote("Methylation ("*beta~"value)"), 
                                   limits=c(0,1.01), expand=c(0,0), breaks=seq(0,1,0.2)) +
                theme_bw() +
                theme(text=element_text(size=12,  family="sans"), axis.text.x = element_text(angle = 45, hjust = 1),
                axis.text.y = element_text(size=12),
                legend.position="top")

pdf(paste0(baseDir, "/", GitHubDir, "/Figures/", Project, ".GenomicFeature.BoxPlot.pdf"),width=6,height=5, onefile=FALSE)
par(bg=NA)
plt.features
dev.off()


message("+-------------------------------------------------------------------------------")
message("+ Plot methyaltion for selected  genomic features")
message("+-------------------------------------------------------------------------------")


cgi.island.first    <- functionGetPerSampleFeature(myNorm.probes.ann, "cgi", "island", "first_meanBeta")  
cgi.island.second   <- functionGetPerSampleFeature(myNorm.probes.ann, "cgi", "island", "second_meanBeta")  

feature.body.first  <- functionGetPerSampleFeature(myNorm.probes.ann, "feature", "Body", "first_meanBeta")
feature.body.second <- functionGetPerSampleFeature(myNorm.probes.ann, "feature", "Body", "second_meanBeta")

promoter.first      <- functionGetPerSampleFeature(myNorm.probes.ann, "promoter", "promoter", "first_meanBeta")  
promoter.second     <- functionGetPerSampleFeature(myNorm.probes.ann, "promoter", "promoter", "second_meanBeta")  


feature.summary.tbl <- rbind(cgi.island.first,   cgi.island.second,   
                             promoter.first,     promoter.second,    
                             feature.body.first, feature.body.second )
head(feature.summary.tbl)



pdf(paste0(baseDir, "/", GitHubDir, "/Figures/", Project, ".GenomicFeatureSelected.BoxPlot.pdf"), width=6,height=5)
par(bg=NA)
ggplot(data=feature.summary.tbl, aes(x=Feature, y=value)) +
  geom_boxplot(aes(fill=variable), outlier.alpha = 0.0) +
  scale_fill_manual(name="", values = c("second_oxBS_meanBeta"=col_2nd, "first_oxBS_meanBeta"=col_1st), 
                    limits=c("first_oxBS_meanBeta","second_oxBS_meanBeta"), 
                    labels=c("First Trimester","Second Trimester")) +
  scale_x_discrete(name="", labels=c("CpG Islands", "Promoters", "Gene Bodies")) +
  scale_y_continuous(name=bquote("Methylation ("*beta~"value)"), 
                     limits=c(0,1.01), expand=c(0,0), breaks=seq(0,1,0.2)) +
  theme_bw() +
  theme(text=element_text(size=12,  family="sans"), axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size=12),
        legend.position="top")
dev.off()




message("+-------------------------------------------------------------------------------")
message("+ Calculate some p-values")
message("+-------------------------------------------------------------------------------")

head(myNorm.probes.ann)

table.for.pearsons.EPIC <- subset(myNorm.probes.ann, cgi=="island")
table.for.pearsons.EPIC <- table.for.pearsons.EPIC[,c("probe_id","first_oxBS_meanBeta", "second_oxBS_meanBeta" )]
nrow(table.for.pearsons.EPIC)
island.first_second.cor  <- cor.test(table.for.pearsons.EPIC[,c("first_oxBS_meanBeta")],table.for.pearsons.EPIC[,c("second_oxBS_meanBeta")])
island.first_second.cor

table.for.pearsons.EPIC <- subset(myNorm.probes.ann, feature=="Body")
table.for.pearsons.EPIC <- table.for.pearsons.EPIC[,c("probe_id","first_oxBS_meanBeta", "second_oxBS_meanBeta" )]
nrow(table.for.pearsons.EPIC)
body.first_second.cor  <- cor.test(table.for.pearsons.EPIC[,c("first_oxBS_meanBeta")],table.for.pearsons.EPIC[,c("second_oxBS_meanBeta")])
body.first_second.cor

table.for.pearsons.EPIC <- subset(myNorm.probes.ann, promoter=="promoter")
table.for.pearsons.EPIC <- table.for.pearsons.EPIC[,c("probe_id","first_oxBS_meanBeta", "second_oxBS_meanBeta" )]
nrow(table.for.pearsons.EPIC)
promoter.first_second.cor  <- cor.test(table.for.pearsons.EPIC[,c("first_oxBS_meanBeta")],table.for.pearsons.EPIC[,c("second_oxBS_meanBeta")])
promoter.first_second.cor


message("+-------------------------------------------------------------------------------")
message("+ Plot Specific Genomic Feature methylation for  1st and 2nd")
message("+-------------------------------------------------------------------------------")


feature.summary.tbl2 <- feature.summary.tbl
head(feature.summary.tbl2)
feature.summary.tbl2$Feature <- gsub("island", "CpG Islands", feature.summary.tbl2$Feature)
feature.summary.tbl2$Feature <- gsub("Body", "Gene Bodies", feature.summary.tbl2$Feature)
feature.summary.tbl2$Feature <- gsub("promoter", "Promoters", feature.summary.tbl2$Feature)

feature.summary.tbl2$Feature_ordered = factor(feature.summary.tbl2$Feature, levels=c("Promoters","Gene Bodies","CpG Islands" ))
head(feature.summary.tbl2)


pdf(paste0(baseDir, "/", GitHubDir, "/Figures/", Project, ".GenomicFeatureSelected.BoxPlot.pdf"), width=10,height=5)
par(bg=NA)
ggplot(feature.summary.tbl2, aes(x = value, y = variable, fill=variable )) + 
  geom_vline(xintercept=0.2, linetype="dashed") +
  geom_vline(xintercept=0.8, linetype="dashed") +
  scale_fill_manual(name="", values = c("first_oxBS_meanBeta"=col_1st, "second_oxBS_meanBeta"=col_2nd), 
                    limits=c("first_oxBS_meanBeta","second_oxBS_meanBeta"), labels=c("First","Second")) +
  scale_x_continuous(breaks=seq(0,1,0.2)) +
  geom_density_ridges(scale=0.75, alpha=0.75) +
  facet_grid(. ~Feature_ordered) +
  xlab( bquote("Methylation ("*beta~"value)")) +
  theme(text=element_text(size=12,  family="sans"),
        strip.background = element_rect(colour=NULL, fill="white"),
        axis.text.x = element_text(size=10),
        legend.position="none")
dev.off()





message("+-------------------------------------------------------------------------------")
message("+ Plot Specific Genomic Feature methylation for  1st and 2nd with Densities")
message("+-------------------------------------------------------------------------------")


"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer( data = data, mapping = mapping, stat = stat, geom = GeomFlatViolin,
         position = position, show.legend = show.legend, inherit.aes = inherit.aes,
         params = list( trim = trim, scale = scale, ... )  )
}

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomFlatViolin <- ggproto("GeomFlatViolin", Geom,
          
            setup_data = function(data, params) {
                         data$width <- data$width %||%
                         params$width %||% (resolution(data$x, FALSE) * 0.9)
                         data %>% group_by(group) %>% mutate(ymin = min(y), ymax = max(y), xmin = x, xmax = x + width / 2)},
  
            draw_group = function(data, panel_scales, coord) {
                         data    <- transform(data, xminv = x, xmaxv = x + violinwidth * (xmax - x))
                         newdata <- rbind(plyr::arrange(transform(data, x = xminv), y), plyr::arrange(transform(data, x = xmaxv), -y))
                         newdata <- rbind(newdata, newdata[1,])
                         ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord)) },
  
            draw_key = draw_key_polygon, default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5, alpha = NA, linetype = "solid"), required_aes = c("x", "y")
  )
head(feature.summary.tbl2)

pdf(paste0(baseDir, "/", GitHubDir, "/Figures/", 
           Project, ".GenomicFeatureSelected.DensityBoxPlotFreeScale.pdf"), width=10,height=5)
par(bg=NA)
ggplot(feature.summary.tbl2, aes(y = value, x = variable, fill=variable)) + 
  geom_hline(yintercept=0.2, linetype="dashed", colour="darkgrey") +
  geom_hline(yintercept=0.8, linetype="dashed", colour="darkgrey") +
  geom_flat_violin(position = position_nudge(x = 0.1, y = 0), alpha=0.95, scale="area") +
  geom_boxplot( width=0.1, outlier.shape = NA, alpha=0.95) + 
  coord_flip() +
  scale_fill_manual(name="", values=c("first_oxBS_meanBeta"=col_1st, "second_oxBS_meanBeta"=col_2nd), 
                    limits=c("first_oxBS_meanBeta","second_oxBS_meanBeta") ) +
  scale_x_discrete(name="", labels=c("first_oxBS_meanBeta"="First Trimester",
                                     "second_oxBS_meanBeta"="Second Trimester")) +
  scale_y_continuous(breaks=seq(0,1,0.2)) +
  facet_grid(. ~Feature_ordered) +
  ylab( bquote("Methylation ("*beta~"value)")) +
  theme(text=element_text(size=12,  family="sans"),
        strip.background = element_rect(colour=NULL, fill="white"),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10), 
        legend.position="none")
dev.off()  
  

message("+-------------------------------------------------------------------------------")
message("+ Plot Specific Genomic Feature methylation for  1st and 2nd with Densities overlaid")
message("+-------------------------------------------------------------------------------")


head(feature.summary.tbl2)
#after_stat()

pdf(paste0(baseDir, "/", GitHubDir, "/Figures/", 
           Project, ".GenomicFeatureSelected.DensityOverlay.pdf"), width=3,height=4.5)
par(bg=NA)
ggplot(feature.summary.tbl2, aes(x =  (value), group=variable, colour=variable, fill=NULL)) + 
  geom_vline(xintercept=0.2, linetype="dashed", colour="darkgrey") +
  geom_vline(xintercept=0.8, linetype="dashed", colour="darkgrey") +
 # geom_density(outline.type="upper", show.legend=F, size=1.5, alpha=0.25 ) +
  stat_density(geom="line",position="identity", size=0.6, alpha=0.75) +
  scale_colour_manual(name="", 
                      values=c("first_oxBS_meanBeta"=col_1st, "second_oxBS_meanBeta"=col_2nd), 
                      limits=c("first_oxBS_meanBeta","second_oxBS_meanBeta"),
                      labels=c("first_oxBS_meanBeta"="FT","second_oxBS_meanBeta"="ST")) +
  scale_x_continuous(name=bquote("Methylation ("*beta~"value)"), breaks=seq(0,1,0.2)) +
  facet_wrap(~Feature_ordered, nrow=3, scales="free") +
  ylab( "" ) +
  guides(linetype = guide_legend(override.aes = list(size = 8))) +
#  guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
#  guides(colour=guide_legend(ncol=1,nrow=1,byrow=TRUE)) +
  theme_cowplot(12) +
  theme(text=element_text(size=12,  family="sans"),
        legend.text=element_text(size=8),
        strip.text.x=element_text(size=10, face="bold"), 
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10), 
        strip.background = element_rect(colour=NULL, fill="white"),
      #  legend.key.height=unit(2, "cm"),
    #  legend.box="vertical",
        legend.position="right")
dev.off()


message("+-------------------------------------------------------------------------------")
message("+ END OF SCRIPT")
message("+-------------------------------------------------------------------------------")