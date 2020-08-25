#!/bin/bash



cat ../Data/CTR_EPIC.First_Second_first_vs_second_DMRs_using_oxBS.csv | grep -v "seqnames" | awk -F "," '{print $2"\t"$3"\t"$4"\t"$1"\t"$7"\t"$16"\t."}' | sed 's/"chr//g' | sed 's/"//g'  > ../Data/CTR_EPIC.First_Second_first_vs_second_DMRs_using_oxBS.bed

bedtools sort -i ../Data/CTR_EPIC.First_Second_first_vs_second_DMRs_using_oxBS.bed > ../Data/CTR_EPIC.First_Second_first_vs_second_DMRs_using_oxBS.srt.bed

bedtools closest -d  -a ../Data/CTR_EPIC.First_Second_first_vs_second_DMRs_using_oxBS.srt.bed  -b ../Data/GRCh37.87.gtf.bed | awk '{ print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$18"\t"$27}' | sed 's/"//g' | sed 's/;//g' > ../Data/CTR_EPIC.First_Second_first_vs_second_DMRs_using_oxBS.srt.ann.bed

bedtools closest -D b -d  -a ../Data/CTR_EPIC.First_Second_first_vs_second_DMRs_using_oxBS.srt.bed -b ../Data/GRCh37.87.gtf.bed | awk '{ print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$18"\t"$27}' | sed 's/"//g' | sed 's/;//g' > ../Data/CTR_EPIC.First_Second_first_vs_second_DMRs_using_oxBS.srt.ann.TEST.bed
