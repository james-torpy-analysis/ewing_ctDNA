#!/bin/bash

sample_name="409_004_D9YWF_GTAGAGGA-CTCTCTAT_L001"

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/ewing_ctDNA"
in_dir="$project_dir/results/picard/bams/"

# call conda env:
snk

for f in $in_dir/*.bam; do

  printf "\n"
  echo "Isolating split reads from $f..."
  printf "\n"

  prefix=$(echo $f | sed "s/.bam//")

  # isolate split reads:
  samtools view -hf 2048 $f | \
    samtools view -bh > $prefix.split.supp.bam
  
  # fetch split read names:
  read_names=( $(samtools view $prefix.split.supp.bam | awk '{print $1}') )
  unique_names=($(echo "${read_names[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
  
  # add header to the final split read file:
  samtools view -H $in_dir/$sample_name.consensus.sorted.bam > \
    $in_dir/$sample_name.consensus.split.sam
  
  # search for split read primary and supp alignments and add to final split read file:
  for n in ${unique_names[@]}; do
    samtools view $in_dir/$sample_name.consensus.sorted.bam | grep $n  \
      >> $in_dir/$sample_name.split.sam
  done;
  
  # remove duplicate rows:
  sort $in_dir/$sample_name.consensus.split.sam | uniq
  
  # convert to bam and remove sam:
  samtools view -bh $in_dir/$sample_name.consensus.split.sam > \
    $in_dir/$sample_name.consensus.split.bam
  rm $in_dir/$sample_name.consensus.split.sam
  
  # sort and index:
  samtools sort -o $in_dir/$sample_name.consensus.split.sorted.bam \
    $in_dir/$sample_name.consensus.split.bam
  samtools index $in_dir/$sample_name.consensus.split.sorted.bam

done;

