# Transcriptome analysis reveals differences in human placental metabolism, transport and endocrine function across the first-second trimester transition #

**Malwina Prater<sup>a</sup>, Russell S. Hamilton<sup>a,b</sup>, Hong Wa Yung<sup>a</sup>, Andrew M. Sharkey<sup>a,b</sup>, Paul Robson<sup>c</sup>, Eric Jauniaux<sup>d</sup>, D. Stephen Charnock-Jones<sup>a,e,f</sup>, Graham J. Burton<sup>a,§</sup>, Tereza Cindrova-Davies<sup>a,§</sup>**

<sub>
<sup>a</sup> Centre for Trophoblast Research, Department of Physiology, Development and Neuroscience, University of Cambridge, Downing Street, Cambridge, CB2 3EG, UK.
<sup>b</sup> Department of Genetics, University of Cambridge, Downing Street, Cambridge, CB2 3EH
<sup>c</sup> Department of Pathology, University of Cambridge, Tennis Court Road, Cambridge, CB2 1QP
<sup>d</sup> The Jackson Laboratory, The JAX Center for Genetics of Fertility and Reproduction, 10 Discovery Drive Farmington, CT 06032, USA
<sup>e</sup> EGA Institute for Women’s Health, Faculty of Population Health Sciences. University College London, London, WC1E 6BT, UK.
<sup>f</sup> Department of Obstetrics and Gynaecology, University of Cambridge, The Rosie Hospital, Cambridge, CB2 0SW, UK
<sup>g</sup> National Institute for Health Research, Cambridge Biomedical Research Centre.

<sup>§</sup> These authors contributed equally

</sub>

# Citation

Prater, M., Hamilton, R.S., Yung, H-W, Sharkey, A.M., Robson, P, Jauniaux, E., Charnock-Jones, D.S., Burton G.J. & Cindrova-Davies, T. (2020) RNA-Seq reveals changes in human placental metabolism, transport and endocrinology function in the first-second trimester transition. Submitted [DOI]

> Note: RNA-Seq analysis is available on a separate [GitHub](https://github.com/CTR-BFX/2020_Prater_Cindrova) repository.

# Methylation Analysis

## Infinium Methylation EPIC array

Genomic DNA was isolated by QIAamp DNA mini kit (Qiagen, cat. no. 51304) following manufacturer’s instructions. Buffer AL (200&mu;l) was added to the sample, mixed by pulse-vortexing for 15 sec, before incubating at 70&deg;C for 10 min. Absolute Ethanol (200&mu;l) was then added to the sample, and mixed by pulse-vortexing for 15 sec before transferring to the QIAamp Mini spin column and centrifuged at 6000g for 1 min. The Mini spin column was washed once with Buffer AW1 (500&mu;l) following by Buffer AW2 (500&mu;l) before centrifuging at full speed for 1 min. For elution of genomic DNA, DNase-free water (100&mu;l) was added and incubated for 1 min before centrifuging at 6000g for 1 min. The step repeated one more time with another 100&mu;l DNase-free water. DNA concentration of the samples were quantified by NanoDrop and the DNA quality was checked by resolving in 0.8% agarose gel, in which there was a major band visualized at around 10 kbp without obvious smear below, indicating good quality DNA.

Genomic DNA was subjected to oxidative bisulfite (oxBS) conversion using the CEGX TrueMethyl kit (Cambridge Epigenetix / NuGEN,  cat. no. CEGXTMS) and used for microarray-based DNA methylation analysis, performed at GenomeScan (GenomeScan B.V., Leiden, The Netherlands) on the HumanMethylation850 BeadChip (Illumina, Inc., San Diego, CA, U.S.A). The resulting iDAT files were imported and analysed using ChAMP (v2.9.10) [1,2]. The EPIC array assays approximately 865,000 CpG sites  Samples were processed filtering for a probe detection p-value <= 0.01, probes with a beadcount <3 in at least 5% of samples, no CpG and known SNPs[3] at probe starts, probes aligning to multiple locations,  and QC using the on array control probes. In total, 750150 probes on the array passed the filtering and QC steps. The BMIQ[4] method was used to normalise the two probe types (I and II) present on the array. Beta methylation values from the EPIC array range from 0 (unmethylated) to 1 (methylated) and are equivalent to percentage methylation.

DMRs were calculated using the bumpHunter methods in ChAMP, and a methylation difference of 0.2 and adjusted p.value of 0.05 between first and second trimesters was used for filtering DMRs. Probes on the X and Y chromosomes were excluded to minimise sex specific in differential methylation calculations. Pearson’s correlation (R function cor.test) was used to calculate p values between DMR methylation and differential gene expression was calculated Pearson’s. Bedtools was used to determine DMRs overlapping gene bodies and promoters (bedtools closest -D b -d) -a DMRs.bed -b GRCh37.87.gtf.bed and 1.5Kb upstream of the TSS used to define promoters). GO analysis on for DMRs in gene bodies and promoters was performed with ‘goregions’ from missMethyl. Multiple DMRs per gene were manually checked and in all cases the methylation change was in the same direction. Where a DMR overlapped two genes, both were included in the correlation with gene expression. A list of the top 100 DMRs associated with sex-specific methylation in human placentas are from Sung et al. An intersection of sex-specific DMRs with the identified DMRs with associated expression changes above thresholds (log2 FC 1, methylation difference 0.2) show no common genes, suggesting the identified DMRs are not sex-specific. 

EPIC methylation array data have been deposited in the ArrayExpress database at EMBL-EBI under accession number E-MTAB-XXXX (https://www.ebi.ac.uk/arrayexpress/experiments/E-MTABXXXX). Code used to analyse the EPIC array samples is available at https://github.com/darogan/First_Second_Trimester_Methylation

[Rscript 1](Correlate_DMRs_2_RNA.R )
[Rscript 2](CTR_EPIC.PrePost.analysis.R)

### EPIC Methylation Array Sample Table ###

| Sample  | iDAT FileName |
|---------|-----------------|
| first_63_oxBS | 202123800269_R04C01_{Red/Grn}.idat |
| first_64_oxBS | 202123800269_R02C01_{Red/Grn}.idat |
| first_65_oxBS | 202123800167_R06C01_{Red/Grn}.idat |
| first_66_oxBS | 202123800269_R08C01_{Red/Grn}.idat |
| first_67_oxBS | 202128330113_R02C01_{Red/Grn}.idat |
| first_69_oxBS | 202128330113_R04C01_{Red/Grn}.idat |
| first_70_oxBS | 202128330113_R06C01_{Red/Grn}.idat |
| second_71_oxBS | 202139520259_R04C01_{Red/Grn}.idat |
| second_72_oxBS | 202139520259_R06C01_{Red/Grn}.idat |
| second_73_oxBS | 202139520259_R08C01_{Red/Grn}.idat |
| second_74_oxBS | 202123800167_R02C01_{Red/Grn}.idat |
| second_75_oxBS | 202123800269_R06C01_{Red/Grn}.idat |
| second_76_oxBS | 202123800167_R04C01_{Red/Grn}.idat |

[iDAT Sample sheet](SampleSheet_Infinium_MethylationEPIC_103409-001.csv)

##### Sample Clustering #####

<IMG SRC="Figures/CTR_EPIC.First_Second.PCA_oxBS.png" width=400px>
<P>
<IMG SRC="Figures/CTR_EPIC.First_Second.Dendro_oxBS.png" width=400px>

##### Methylation Correlation: First Vs Second trimester #####

|   |   |
|---------|-----------------|
| <IMG SRC="Figures/CTR_EPIC.First_Second.MethylationCorrelation.png" width=300px> | <IMG SRC="Figures/CTR_EPIC.First_Second.AllProbes.Correlations.png" width=300px> |



##### Methylation Vs Genomic Features #####

|   |   |
|---------|-----------------|
| <IMG SRC="Figures/CTR_EPIC.First_Second.GenomicFeatureSelected.DensityOverlay.png" width=400px> | <IMG SRC="Figures/CTR_EPIC.First_Second.All.Feature.Correlations.png" width=400px> |

##### Methylation Vs Know Sex Specific gene expression and methylation  #####

|   |   |
|---------|-----------------|
| Gong et al (2018). Placental polyamine metabolism differs by fetal sex, fetal growth restriction, and preeclampsia. JCI Insight, 3(13). http://doi.org/10.1172/jci.insight.120723 | Gong et al. (2018) Genome-wide oxidative bisulfite sequencing identifies sex-specific methylation differences in the human placenta. Epigenetics. http://doi.org/10.1080/15592294.2018.1429857 |
| <IMG SRC="Figures/CTR_EPIC.First_Second.DMRVsDEGVsSexSpecificDEGs.Venn.png" width=400px> | <IMG SRC="Figures/CTR_EPIC.First_Second.DMRVsDEGVsSexSpecificDMR.Venn.png" width=400px> |



##### Methylation Vs Expression Correlation #####

<IMG SRC="Figures/CTR_EPIC.First_Second.Correlation_DMR_RNA.png" width=400px>



### References ###
1.	Morris, T. J. et al. ChAMP: 450k Chip Analysis Methylation Pipeline. Bioinformatics 30, 428–430 (2014).
2.	Aryee, M. J. et al. Minfi: A flexible and comprehensive Bioconductor package for the analysis of Infinium DNA methylation microarrays. Bioinformatics 30, 1363–1369 (2014).
3.	Zhou, W., Laird, P. W., Shen, H., Infinium, I. & Methylation, D. N. A. Comprehensive characterization, annotation and innovative use of Infinium DNA methylation BeadChip probes. Nucleic Acids Res. 45, e22 (2017).
4.	Teschendorff, A. E. et al. A beta-mixture quantile normalization method for correcting probe design bias in Illumina Infinium 450 k DNA methylation data. Bioinformatics 29, 189–196 (2013).
5. Gong et al (2018). Placental polyamine metabolism differs by fetal sex, fetal growth restriction, and preeclampsia. JCI Insight, 3(13).
6. Gong et al. (2018) Genome-wide oxidative bisulfite sequencing identifies sex-specific methylation differences in the human placenta. Epigenetics.


### Links ###

Description   | URL
------------- | ----------
Publication   | [Journal](https://) and [DOI](https://doi.org/)
Raw Data (EPIC) | ArrayExpress EMBL-EBI [E-MTAB-9312](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-9312)
Raw Data (RNA-Seq) | ArrayExpress EMBL-EBI [E-MTAB-9203](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-9203)



### Contact

Contact Russell S. Hamilton (rsh46@cam.ac.uk) for bioinformatics related queries.
