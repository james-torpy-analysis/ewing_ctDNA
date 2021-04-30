#!/bin/bash

### Fetch discordant and split reads from bams ###

# define variables/directories:
project_name="ewing_ctDNA"
sample_name="409_018_DBV4V_AAGAGGCA-CTCTCTAT_L001"

home_dir="/share/ScratchGeneral/jamtor/"
project_dir="$home_dir/projects/$project_name/"
result_dir="$project_dir/results/"

bwa_dir="$result_dir/bwa/"
geneglobe_dir="$result_dir/geneglobe/"

# filter out discordant reads:
#samtools view -H $bwa_dir/$sample_name/$sample_name.sorted.bam > \
#  $bwa_dir/$sample_name/$sample_name.discordant.sam && \
#  samtools view $bwa_dir/$sample_name/$sample_name.sorted.bam | \
#    awk 'BEGIN { OFS = "\t"; } ($3=$7 && $7!="=" && $3!="*" && $4!="0")' | \
#    awk 'BEGIN { OFS = "\t"; } $3="chr"$3' >> \
#      $bwa_dir/$sample_name/$sample_name.discordant.sam && \
#    samtools view -bh $bwa_dir/$sample_name/$sample_name.discordant.sam > \
#      $bwa_dir/$sample_name/$sample_name.discordant.bam && \
#    samtools sort -o $bwa_dir/$sample_name/$sample_name.discordant.sorted.bam \
#      $bwa_dir/$sample_name/$sample_name.discordant.bam && \
#    samtools index $bwa_dir/$sample_name/$sample_name.discordant.sorted.bam && \
#    rm $bwa_dir/$sample_name/$sample_name.discordant.sam

samtools view -h -F 1294 $bwa_dir/$sample_name/$sample_name.sorted.bam | \
  samtools view -bh > $bwa_dir/$sample_name/$sample_name.discordant.bam && \
    samtools sort -o $bwa_dir/$sample_name/$sample_name.discordant.sorted.bam \
      $bwa_dir/$sample_name/$sample_name.discordant.bam && \
    samtools index $bwa_dir/$sample_name/$sample_name.discordant.sorted.bam

samtools view -h -F 1294 $geneglobe_dir/$sample_name/$sample_name.sorted.bam | \
  samtools view -bh > $geneglobe_dir/$sample_name/$sample_name.discordant.bam && \
    samtools sort -o $geneglobe_dir/$sample_name/$sample_name.discordant.sorted.bam \
      $geneglobe_dir/$sample_name/$sample_name.discordant.bam && \
    samtools index $geneglobe_dir/$sample_name/$sample_name.discordant.sorted.bam

# filter out split reads:
samtools view -h -f 2048 $bwa_dir/$sample_name/$sample_name.sorted.bam | \
  samtools view -bh > $bwa_dir/$sample_name/$sample_name.split.bam && \
    samtools sort -o $bwa_dir/$sample_name/$sample_name.split.sorted.bam \
      $bwa_dir/$sample_name/$sample_name.split.bam && \
    samtools index $bwa_dir/$sample_name/$sample_name.split.sorted.bam

samtools view -h -f 2048 $geneglobe_dir/$sample_name/$sample_name.sorted.bam | \
  samtools view -bh > $geneglobe_dir/$sample_name/$sample_name.split.bam && \
    samtools sort -o $geneglobe_dir/$sample_name/$sample_name.split.sorted.bam \
      $geneglobe_dir/$sample_name/$sample_name.split.bam && \
    samtools index $geneglobe_dir/$sample_name/$sample_name.split.sorted.bam