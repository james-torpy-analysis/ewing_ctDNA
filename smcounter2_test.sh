#!/bin/bash

source ~/.bashrc

sample_name="409_021_DBV4V_CTCTCTAC-CTCTCTAT_L001"

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/ewing_ctDNA"
fq_dir="$project_dir/raw_files"
result_dir="$project_dir/results"
in_dir="$result_dir/picard/$sample_name"
int_dir="$result_dir/smcounter2/$sample_name/intermediate_files"
out_dir="$result_dir/smcounter2/$sample_name"
out_pref="$out_dir/$sample_name"

mkdir -p $in_dir
mkdir -p $int_dir

picard_dir="$home_dir/local/bin"

# call snk conda env for samtools/python:
conda activate snkenv

printf "\n\n"
echo "--------------------------------------------------"
echo "converting fastq files to one unmapped bam file..."
echo "--------------------------------------------------"
printf "\n"
java -jar $in_dir/picard.jar FastqToSam \
  FASTQ=$fq_dir/$sample_name\_R1.fastq.gz \
  FASTQ2=$fq_dir/$sample_name\_R2.fastq.gz \
  O=$in_dir/$sample_name.unmapped.bam \
  SM=sample

printf "\n\n"
echo "--------------------------------------------------"
echo "aligning unmapped bam using smCounter2..."
echo "--------------------------------------------------"
printf "\n"

singularity exec \
  -B /share/ScratchGeneral/jamtor/local/lib/qiaseq-smcounter-v2:/rundir \
  -B $int_dir \
  -B $in_dir \
  /share/ScratchGeneral/jamtor/projects/ewing_ctDNA/containers/smcounter-v2/smcounter2_28Apr21v1.sif \
  python /rundir/run.py \
  --runPath $int_dir \
  --bamFile $in_dir/$sample_name.unmapped.bam \
  --outPrefix $out_pref \
  --nCPU 6