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


library("ggplot2")
library('gplots')
library("ggforce")
library('minfi')
library('lumi')
library('methylumi')
library('limma')
library('VennDiagram');
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


col_1st <- "firebrick2"
col_2nd <- "steelblue3"



#ensembl    <-  useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
#ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description'), mart = ensembl)

#grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
#listAttributes(grch37)
# http://zwdzwd.github.io/InfiniumAnnotation
#test <- read.table(gzfile("EPIC.hg38.manifest.gencode.v22.tsv.gz"), sep="\t", header=TRUE)
#head(test, 100)

#library("IlluminaHumanMethylationEPICmanifest")
#library("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

#
# Load in the EPIC probe manifest
#
data(probe.features.epic)
str(probe.features)
probe.features$probe_id <- rownames(probe.features)
head(probe.features)


probe.locations.epic <- probe.features[, c("CHR","MAPINFO", "Strand")]
probe.locations.epic$probe_id <- rownames(probe.locations.epic)
head(probe.locations.epic)


#anno <- IlluminaHumanMethylationEPICanno.ilm10b2.hg19@data
#names(anno)
#head(anno)
#names(IlluminaHumanMethylationEPICanno.ilm10b2.hg19@defaults)




#data(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
#ILMN.GR <- minfi::getLocations(IlluminaHumanMethylationEPICanno.ilm10b2.hg19::IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
#ILMN.GR <- ILMN.GR[substring(names(ILMN.GR), 1, 2) != "ch"]
#ILMN.GR$Index <- names(ILMN.GR)
#head(ILMN.GR, 25)

#promoters(genes(ILMN.GR), upstream = 1500, downstream = 500)


#------------------------------------------------------------------------------
# Set the directories and baseDir for the 450K experiment
setwd("/Users/rhamilto/Documents/CTR-Data/CTR_EPIC/")
baseDir <- getwd()
baseDirIDATs <- sprintf("%s/IDATs",baseDir)

#Project <- "CTR_EPIC.organoid"
Project <- "CTR_EPIC.pre-post"




#------------------------------------------------------------------------------

# Specifically for the Gender plot
#myLoad <- champ.load(directory = baseDirIDATs, methValue='B', filterBeads=FALSE, filterXY=FALSE, detPcut = 0.01, arraytype="EPIC")

myLoad <- champ.load(directory = baseDirIDATs, methValue='B', filterBeads=TRUE, filterXY=TRUE, detPcut = 0.01, arraytype="EPIC")
head(myLoad$pd)

#myLoad$pd$Sample_Name <- gsub("pre",  "first", myLoad$pd$Sample_Name)
#myLoad$pd$Sample_Name <- gsub("post", "second", myLoad$pd$Sample_Name)
#print(myLoad$pd)


sampleTable.champ <- as.data.frame(myLoad$pd)



sampleTable.champ$Sex <- sampleTable.champ$Sample_Name
sampleTable.champ$Sex <- gsub("first_63_oxBS", "M", sampleTable.champ$Sex)
sampleTable.champ$Sex <- gsub("first_64_oxBS", "F", sampleTable.champ$Sex)
sampleTable.champ$Sex <- gsub("first_65_oxBS", "F", sampleTable.champ$Sex)
sampleTable.champ$Sex <- gsub("first_66_oxBS", "M", sampleTable.champ$Sex)
sampleTable.champ$Sex <- gsub("first_67_oxBS", "M", sampleTable.champ$Sex)

sampleTable.champ$Sex <- gsub("first_69_oxBS", "F", sampleTable.champ$Sex)
sampleTable.champ$Sex <- gsub("first_70_oxBS", "F", sampleTable.champ$Sex)

sampleTable.champ$Sex <- gsub("second_71_oxBS", "M", sampleTable.champ$Sex)
sampleTable.champ$Sex <- gsub("second_72_oxBS", "M", sampleTable.champ$Sex)
sampleTable.champ$Sex <- gsub("second_73_oxBS", "F", sampleTable.champ$Sex)
sampleTable.champ$Sex <- gsub("second_74_oxBS", "M", sampleTable.champ$Sex)
sampleTable.champ$Sex <- gsub("second_75_oxBS", "M", sampleTable.champ$Sex)
sampleTable.champ$Sex <- gsub("second_76_oxBS", "M", sampleTable.champ$Sex)
sampleTable.champ$Sex <- gsub("second_77_oxBS", "X", sampleTable.champ$Sex)





#myLoad$pd$GeneScan_Name <- myLoad$pd$Sample_Name
#myLoad$pd$Sample_Name   <- myLoad$pd$Label


myQC <- champ.QC()

myNorm <- champ.norm()
## PCA Analysis
## PCAPlot(t(myNorm$beta),myLoad$pd$Sample_Group,ProjectName,multifigure=T)


#myLoad <- ""

customPCADendro <- function(ProjectTitle, myNorm, TOPNUM, sampleTable.champ) {
  
  myNorm.df <- as.data.frame(myNorm)
  
  rv        <- rowVars(as.matrix(myNorm.df))
  select    <- order(rv, decreasing = TRUE)[seq_len(min(TOPNUM, length(rv)))]
  pca       <- prcomp(t(myNorm.df[select, ]))
  pc1var    <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var    <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab    <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab    <- paste0("PC2 (",as.character(pc2var),"%)")
  scores    <- data.frame(sampleName=sampleTable.champ$Sample_Name, pca$x, Sample_Group=sampleTable.champ$Sample_Group, Sex=sampleTable.champ$Sex)

  scores$sampleName   <- gsub("_oxBS", "", scores$sampleName)
  scores$Sample_Group <- gsub("_oxBS", "", scores$Sample_Group)
  
  scores$Label <- scores$Sample_Group
  scores$Label <- gsub("first",  "First Trimester", scores$Label)
  scores$Label <- gsub("second", "Second Trimester", scores$Label)
  
  
  plt.pca <- ggplot(scores, aes(x = PC1, y = PC2, fill=Sample_Group, colour=Sample_Group, shape=Sex) ) +
             geom_mark_ellipse(aes(fill = NULL, group=Sample_Group, color=Sample_Group, label=Label), 
                               alpha=0.1, label.fontsize = 10, label.buffer = unit(1, 'mm')) +
          #   geom_text_repel(aes(label=sampleName), show.legend = FALSE, size=2, colour="black") +
             geom_point(size = 4, alpha=0.75 ) + 
             xlab(pc1lab) + ylab(pc2lab) + #coord_fixed() +
             scale_colour_manual(name="", values = c("first"="firebrick2", "second"="steelblue3"), labels=c("first"="First Trimester", "second"="Second Trimester")) +
             scale_fill_manual(  name="", values = c("first"="firebrick2", "second"="steelblue3"), labels=c("first"="First Trimester", "second"="Second Trimester")) +
             expand_limits(x = c(-5, 5), y = c(-4,4)) +
             scale_x_continuous(breaks = seq(-5,5, 2.5)) +
             scale_y_continuous(breaks = seq(-4,4, 2)) +
    
             theme_bw() +
             theme(text = element_text(size=12), aspect.ratio=1, legend.position="none") 

  png(paste0(Project, ".QC.PCA.", TOPNUM, '.png'), units="cm", width=15, height=12.5, res=250)
  par(bg=NA)
  print(plt.pca)
  dev.off()
  
  loadings                       <- as.data.frame(pca$rotation)
  topX <- 20
  
  pca.1           <-  loadings[ order(loadings$PC1,decreasing=TRUE), ]
  pca.1.topX      <-  pca.1[c(1:topX),]

  pca.2           <-  loadings[ order(loadings$PC2,decreasing=TRUE), ]
  pca.2.topX      <-  pca.2[c(1:topX),]

  print(pca.1.topX)

  
  dd.row     <- as.dendrogram(hclust(dist(t( myNorm.df[select, ]))))
  ddata      <- dendro_data(dd.row)
  labs       <- label(ddata)
  
  labs$group <- labs$label
  labs$group <- gsub("first.*",  "First Trimester", labs$group)
  labs$group <- gsub("second.*", "Second Trimester", labs$group)


  plt.dendro <- ggplot(segment(ddata)) +
                geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
                #geom_text(data=label(ddata),aes(label=label, x=x, y=-0.25, angle=-90, colour=labs$group), hjust=0) +
                scale_colour_manual(name="Sample Group", values = c("control"="red", "first"="purple", "second"="blue"), guide=FALSE) +
                ylim(-3, max(ddata$segments$y)) + ylab("distance") + xlab("") +
                ggtitle(paste0("Dendrogram Top ", TOPNUM, " MV")) +
                theme(axis.text.x = element_blank(),  axis.ticks.x = element_blank()) 
  
  png(paste0(Project, ".QC.dendrogram.", TOPNUM, '.png'), units="cm", width=15, height=12.5, res=250)
  par(bg=NA)
  print(plt.dendro)
  dev.off()
  
  myNorm.df <- ""
  
  return(list(plt.pca, plt.dendro))
}
plt.QC <- customPCADendro(Project, myNorm, 500, sampleTable.champ)

pdf(paste(Project, "_PCA_oxBS.pdf", sep=""),width=5,height=5, onefile=FALSE)
par(bg=NA)
plt.QC[[1]]
dev.off()

plt.QC[[2]]







functionPlotBetaPercentage <- function(myNorm, idx1, sampleTable.champ) {
  myNorm.df           <- as.data.frame(myNorm)
  colnames(myNorm.df) <- sampleTable.champ$Sample_Name
  x                   <-  myNorm.df[,colnames(myNorm.df)[idx1]]
  x                   <- as.data.frame(x)
  xt                  <- colnames(myNorm.df)[idx1]

  plt <- ggplot(x,aes(x=x)) +
         geom_histogram(aes(y=(..count..)/sum(..count..)), binwidth=.2, boundary=0, 
                        colour="black", fill="grey", alpha=0.75, size=0.1) +
         scale_y_continuous(labels = percent, breaks=seq(0,1,0.05)) +
         scale_x_continuous(breaks=seq(0,1,0.2)) +
         xlab(bquote(beta~"(methylation)")) + 
         ylab("Percentage") +
         ggtitle(paste0(xt)) +
         theme(text=element_text(size=6, family="sans"), plot.title=element_text(size=6, family="sans"),
               axis.text.x=element_text(size=6, family="sans"), axis.text.y=element_text(size=6, family="sans"),
               legend.position="none")
 
  png(paste0("SamplePlots/", Project, ".BSoxBS.PercentBeta.", xt, '.png'), 
      units="cm", width=7.5, height=7.5, res=250)
  par(bg=NA)
  print(plt)
  dev.off()


return(plt)
}
functionPlotBetaPercentage(myNorm, 3, sampleTable.champ)
functionPlotBetaPercentage(myNorm, 4, sampleTable.champ)




#
# Calculate DMPs using true mC (oxBS)
#
first_oxBS_vs_second_oxBS_DMP      <- champ.DMP(beta=myNorm, pheno=myLoad$pd$Sample_Group, 
                                                compare.group=c("first_oxBS", "second_oxBS"), arraytype="EPIC")
first_oxBS_vs_second_oxBS_DMP.flt  <- first_oxBS_vs_second_oxBS_DMP[[1]][first_oxBS_vs_second_oxBS_DMP[[1]]$deltaBeta>0,]

head(first_oxBS_vs_second_oxBS_DMP.flt)
nrow(first_oxBS_vs_second_oxBS_DMP.flt)
#write.csv(first_oxBS_vs_second_oxBS_DMP.flt, file=paste0(Project, "_first_vs_second_DMPs_using_oxBS", ".csv"))

#
# Calculate DMRs using true mC (oxBS)
#

# BumphunterDMR
first_oxBS_vs_second_oxBS_DMR       <- champ.DMR(beta=myNorm, pheno=myLoad$pd$Sample_Group, 
                                                 compare.group=c("first_oxBS", "second_oxBS"), arraytype="EPIC")

str(first_oxBS_vs_second_oxBS_DMR)
head(first_oxBS_vs_second_oxBS_DMR$BumphunterDMR)
tail(first_oxBS_vs_second_oxBS_DMR$BumphunterDMR)
nrow(first_oxBS_vs_second_oxBS_DMR$BumphunterDMR)
#write.csv(first_oxBS_vs_second_oxBS_DMR$BumphunterDMR, file=paste0(Project, "_first_vs_second_DMRs_using_oxBS", ".csv"))


# DMRcate
first_oxBS_vs_second_oxBS_DMRcate       <- champ.DMR(beta=myNorm, pheno=myLoad$pd$Sample_Group, method="DMRcate",
                                                     compare.group=c("first_oxBS", "second_oxBS"), arraytype="EPIC")

str(first_oxBS_vs_second_oxBS_DMRcate)
head(first_oxBS_vs_second_oxBS_DMRcate$DMRcateDMR)
nrow(first_oxBS_vs_second_oxBS_DMRcate$DMRcateDMR)

#
# Calculate DMBs using true mC (oxBS)
#

first_oxBS_vs_second_oxBS_DMB <- champ.Block(beta=myNorm,pheno=myLoad$pd$Sample_Group, arraytype="EPIC")

str(first_oxBS_vs_second_oxBS_DMB)
head(first_oxBS_vs_second_oxBS_DMB$Block)
nrow(first_oxBS_vs_second_oxBS_DMB$Block)
#write.csv(first_oxBS_vs_second_oxBS_DMB$Block, file=paste0(Project, "_first_vs_second_DMBs_using_oxBS", ".csv"))









#
# Look at the genomic positions
#
threshold <- 0.8


data(probe.features.epic)
str(probe.features)
probe.features$probe_id <- rownames(probe.features)
head(probe.features)
nrow(probe.features)


#
# Annotate the EPIC samples
#
nrow(myNorm)
myNorm.probes.df                <- as.data.frame(myNorm)
myNorm.probes.df$probe_id       <- rownames(myNorm.probes.df)
myNorm.probes.ann               <- merge(myNorm.probes.df, probe.features, by="probe_id")
myNorm.probes.ann$Strand        <- gsub("F", "1", myNorm.probes.ann$Strand)
myNorm.probes.ann$Strand        <- gsub("R", "-1", myNorm.probes.ann$Strand)

myNorm.probes.ann$promoter      <- myNorm.probes.ann$feature
myNorm.probes.ann$promoter      <- gsub("TSS1500", "promoter", myNorm.probes.ann$promoter)
myNorm.probes.ann$promoter      <- gsub("TSS200", "promoter", myNorm.probes.ann$promoter)

myNorm.probes.ann$first_oxBS_meanBeta  <- rowMeans( myNorm.probes.ann[,grep("first_.*_oxBS",colnames(myNorm.probes.ann))] )
myNorm.probes.ann$second_oxBS_meanBeta <- rowMeans( myNorm.probes.ann[,grep("second_.*_oxBS",colnames(myNorm.probes.ann))] )

head(myNorm.probes.ann)
nrow(myNorm.probes.ann)









functionGetPerSampleFeature <- function(myNorm, Feature, SubFeature, sampleGroup) {
  
 # myNorm           <- melt(myNorm[, c(Feature, "org_meanBeta", "pre_meanBeta", "cntl_meanBeta")])
  #Feature <- paste0("\\b", Feature, "\\b")
  myNorm           <- melt(myNorm[, c("probe_id", Feature, colnames(myNorm)[grep("meanBeta",colnames(myNorm))] )])
  
  colnames(myNorm) <- c("probe_id", "Feature","variable","value")
  myNorm           <- subset(myNorm, Feature==SubFeature)

  #myNorm           <- subset(myNorm, variable!="post_meanBeta")
  
  
  print( head(myNorm))
  
  # myNorm.summary             <- summary(subset(myNorm, myNorm$variable==sampleGroup)$value)
  # myNorm.summary["Std.Dev."] <- round(sd(subset(myNorm, myNorm$variable==sampleGroup)$value),7)
  # 
  # print(myNorm.summary)
  # 
  # myNorm.summary             <- data.frame(sample=sampleGroup, feature=SubFeature, 
  #                                          ymin=(myNorm.summary[[4]]-myNorm.summary[[7]]), lower=myNorm.summary[[2]], middle=myNorm.summary[[4]], 
  #                                          upper=myNorm.summary[[5]], ymax=(myNorm.summary[[4]]+myNorm.summary[[7]]), sd=myNorm.summary[[7]])
  # myNorm.summary$ymin[myNorm.summary$ymin < 0] <- 0
  # 
  # print(myNorm.summary)
  
#  return(myNorm.summary) 
  return(myNorm) 
  
}
head(myNorm.probes.ann)



myNorm.probes.ann.summary           <- myNorm.probes.ann
rownames(myNorm.probes.ann.summary) <- myNorm.probes.ann.summary$probe_id
myNorm.probes.ann.summary           <- myNorm.probes.ann.summary[, c("CHR","MAPINFO", "Strand", 
                                                                     "first_oxBS_meanBeta", "second_oxBS_meanBeta")]
head(myNorm.probes.ann.summary)
nrow(myNorm.probes.ann.summary)
write.table(myNorm.probes.ann.summary, file=paste0(Project, ".myNorm.probes.ann.summary", ".csv"), col.names=TRUE, row.names=TRUE, quote=FALSE, sep=",")


test <- myNorm.probes.ann[, c("probe_id", "cgi", "gene", "feature", "feat.cgi", "promoter", colnames(myNorm.probes.ann)[grep("meanBeta",colnames(myNorm.probes.ann))] )]
head(test)


threshold <- 0.2
ggplot(data=test, aes(x=first_oxBS_meanBeta, y=second_oxBS_meanBeta, colour=abs(first_oxBS_meanBeta-second_oxBS_meanBeta), label=gene)) +
  geom_point(data=subset(test, abs(first_oxBS_meanBeta-second_oxBS_meanBeta) >= 0 ), 
             size=0.25, alpha=0.75) +
  scale_colour_gradient2(name="Methylation\nDifference", low = "grey", mid="grey", high = "darkblue", midpoint = 0.02) +
  geom_abline(intercept = 2*threshold, slope = 1, linetype="dashed") +
  geom_abline(intercept = threshold, slope = 1, linetype="dashed") +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  geom_abline(intercept = -(threshold), slope = 1, linetype="dashed") +
  geom_abline(intercept = -(2*threshold), slope = 1, linetype="dashed") +
  geom_label_repel(data=subset(test, (first_oxBS_meanBeta-second_oxBS_meanBeta) <= -2*threshold), force=10, size=3, colour="black", nudge_x=-0.05, nudge_y=0.05 ) +
  geom_label_repel(data=subset(test, (first_oxBS_meanBeta-second_oxBS_meanBeta) >=  2*threshold), force=10, size=3, colour="black", nudge_x=0.05,  nudge_y=-0.05) +
  scale_x_continuous(name="First Trimester",  limits = c(0,1.01), expand = c(0, 0), breaks=seq(0,1,0.2)) +
  scale_y_continuous(name="Second Trimester", limits = c(0,1.01), expand = c(0, 0), breaks=seq(0,1,0.2)) +
  coord_fixed() +
  theme_bw()


nrow(subset(test, abs(first_oxBS_meanBeta-second_oxBS_meanBeta) >= 2*threshold  ) )


unique(test$feature)
unique(test$cgi)
unique(test$feat.cgi)

test.feature <- myNorm.probes.ann[, c("probe_id", "feature", colnames(myNorm.probes.ann)[grep("meanBeta",colnames(myNorm.probes.ann))] )]
head(test.feature)
nrow(test.feature)
tail(test.feature)

str(test.feature)

test.feature.m <- melt(test.feature, id=c("probe_id","feature"))
head(test.feature.m)

ggplot(data=test.feature.m, aes(x=feature, y=value, fill=variable)) +
  #geom_density_ridges() +
  geom_boxplot(aes(fill=variable), outlier.alpha = 0.0) +
   scale_fill_manual(name="", values = c("second_oxBS_meanBeta"="lightpink1", "first_oxBS_meanBeta"="dodgerblue"), 
                     limits=c("first_oxBS_meanBeta","second_oxBS_meanBeta"), labels=c("First","Second")) +
  ylab( bquote("Methylation ("*beta~"value)")) +
  xlab("") +
  ylim(0,1) +
  theme(text=element_text(size=12,  family="sans"), axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size=12),
        legend.position="top")



cgi.island.first  <- functionGetPerSampleFeature(myNorm.probes.ann, "cgi", "island", "first_meanBeta")  
cgi.island.second <- functionGetPerSampleFeature(myNorm.probes.ann, "cgi", "island", "second_meanBeta")  

feature.body.first  <- functionGetPerSampleFeature(myNorm.probes.ann, "feature", "Body", "first_meanBeta")
feature.body.second <- functionGetPerSampleFeature(myNorm.probes.ann, "feature", "Body", "second_meanBeta")

promoter.first   <- functionGetPerSampleFeature(myNorm.probes.ann, "promoter", "promoter", "first_meanBeta")  
promoter.second  <- functionGetPerSampleFeature(myNorm.probes.ann, "promoter", "promoter", "second_meanBeta")  


feature.summary.tbl          <- rbind(cgi.island.first,   cgi.island.second,   
                                      promoter.first,     promoter.second,    
                                      feature.body.first, feature.body.second )
                                   
head(feature.summary.tbl)
unique(feature.summary.tbl$variable)
unique(feature.summary.tbl$Feature)


pdf(paste(Project, "-GenomicFeatureSummary_First_Vs_Second_Version1.pdf", sep=""), width=6,height=5)
par(bg=NA)
ggplot(data=feature.summary.tbl, aes(x=Feature, y=value)) +
  geom_boxplot(aes(fill=variable), outlier.alpha = 0.0) +
  scale_fill_manual(name="", values = c("second_oxBS_meanBeta"="lightpink1", "first_oxBS_meanBeta"="dodgerblue"), 
                     limits=c("first_oxBS_meanBeta","second_oxBS_meanBeta"), labels=c("First","Second")) +
  scale_x_discrete(name="", labels=c("CpG Islands", "Promoters", "Gene Bodies")) +
  ylab( bquote("Methylation ("*beta~"value)")) +
  ylim(0,1) +
  theme_bw() +
  theme(text=element_text(size=12,  family="sans"), axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size=12),
        legend.position="top")
dev.off()




#
# Calculate some p-values
# K

head(myNorm.probes.ann)

table.for.pearsons.EPIC <- subset(myNorm.probes.ann, cgi=="island")
table.for.pearsons.EPIC <- table.for.pearsons.EPIC[,c("probe_id","first_oxBS_meanBeta", "second_oxBS_meanBeta" )]
nrow(table.for.pearsons.EPIC)
island.first_second.cor  <- cor.test(table.for.pearsons.EPIC[,c("first_oxBS_meanBeta")],table.for.pearsons.EPIC[,c("second_oxBS_meanBeta")])
island.first_second.cor$estimate
island.first_second.cor$p.value

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




feature.summary.tbl2 <- feature.summary.tbl
head(feature.summary.tbl2)
feature.summary.tbl2$Feature <- gsub("island", "CpG Islands", feature.summary.tbl2$Feature)
feature.summary.tbl2$Feature <- gsub("Body", "Gene Bodies", feature.summary.tbl2$Feature)
feature.summary.tbl2$Feature <- gsub("promoter", "Promoters", feature.summary.tbl2$Feature)

feature.summary.tbl2$Feature_ordered = factor(feature.summary.tbl2$Feature, levels=c("CpG Islands","Promoters","Gene Bodies" ))
head(feature.summary.tbl2)


pdf(paste(Project, "-GenomicFeatureSummary_First_Second_DensityVersion.pdf", sep=""), width=10,height=5)
par(bg=NA)
ggplot(feature.summary.tbl2, aes(x = value, y = variable, fill=variable )) + 
  geom_vline(xintercept=0.2, linetype="dashed") +
  geom_vline(xintercept=0.8, linetype="dashed") +
  scale_fill_manual(name="", values = c("first_oxBS_meanBeta"=col_1st, "second_oxBS_meanBeta"=col_2nd), 
                    limits=c("first_oxBS_meanBeta","second_oxBS_meanBeta"), labels=c("First","Second")) +
  scale_x_continuous(breaks=seq(0,1,0.2)) +
  geom_density_ridges(scale=0.75, alpha=0.5) +
  facet_grid(. ~Feature_ordered) +
  xlab( bquote("Methylation ("*beta~"value)")) +
  theme(text=element_text(size=12,  family="sans"),
        axis.text.x = element_text(size=10),
        legend.position="none")
dev.off()








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

pdf(paste(Project, "-GenomicFeatureSummary_Density_and_Box_FreeScale.pdf", sep=""), width=10,height=5)
par(bg=NA)
ggplot(feature.summary.tbl2, aes(y = value, x = variable, fill=variable)) + 
  geom_hline(yintercept=0.2, linetype="dashed") +
  geom_hline(yintercept=0.8, linetype="dashed") +
  geom_flat_violin(position = position_nudge(x = 0.1, y = 0), alpha=0.95, scale="area") +
  geom_boxplot( width=0.1, outlier.shape = NA, alpha=0.95) + 
  coord_flip() +
  scale_fill_manual(name="", values = c("first_oxBS_meanBeta"=col_1st, "second_oxBS_meanBeta"=col_2nd), 
                    limits=c("first_oxBS_meanBeta","second_oxBS_meanBeta"), labels=c("First","Second")) +
  scale_y_continuous(breaks=seq(0,1,0.2)) +
  facet_grid(. ~Feature_ordered) +
  ylab( bquote("Methylation ("*beta~"value)")) +
  theme(text=element_text(size=12,  family="sans"),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10), # angle = 45, hjust = 1),
        legend.position="none")
dev.off()  
  


#------------------------------------------------------------------------------
# END OF SCRIPT
#------------------------------------------------------------------------------