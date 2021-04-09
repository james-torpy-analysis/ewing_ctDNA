#!/bin/bash

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/ewing_ctDNA"
in_dir="$project_dir/results/bwa/409_016_DBV4V_TAGGCATG-CTCTCTAT_L001"

# call conda env:
snk

# isolate split reads:
samtools view -hf 2048 $in_dir/409_016_DBV4V_TAGGCATG-CTCTCTAT_L001.sorted.bam | \
  samtools view -bh > $in_dir/409_016_DBV4V_TAGGCATG-CTCTCTAT_L001.split.supp.bam

# fetch split read names:
read_names=( $(samtools view $in_dir/409_016_DBV4V_TAGGCATG-CTCTCTAT_L001.split.supp.bam | awk '{print $1}') )
unique_names=($(echo "${read_names[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

# add header to the final split read file:
samtools view -H $in_dir/409_016_DBV4V_TAGGCATG-CTCTCTAT_L001.sam > $in_dir/409_016_DBV4V_TAGGCATG-CTCTCTAT_L001.split.sam

# search for split read primary and supp alignments and add to final split read file:
for n in ${unique_names[@]}; do
  echo $n
  grep $n $in_dir/409_016_DBV4V_TAGGCATG-CTCTCTAT_L001.sam >> $in_dir/409_016_DBV4V_TAGGCATG-CTCTCTAT_L001.split.sam
done;

# remove duplicate rows:
sort $in_dir/409_016_DBV4V_TAGGCATG-CTCTCTAT_L001.split.sam | uniq

# convert to bam and remove sam:
samtools view -bh $in_dir/409_016_DBV4V_TAGGCATG-CTCTCTAT_L001.split.sam > $in_dir/409_016_DBV4V_TAGGCATG-CTCTCTAT_L001.split.bam
rm $in_dir/409_016_DBV4V_TAGGCATG-CTCTCTAT_L001.split.sam

# sort and index:
samtools sort -o $in_dir/409_016_DBV4V_TAGGCATG-CTCTCTAT_L001.split.sorted.bam $in_dir/409_016_DBV4V_TAGGCATG-CTCTCTAT_L001.split.bam
samtools index $in_dir/409_016_DBV4V_TAGGCATG-CTCTCTAT_L001.split.sorted.bam