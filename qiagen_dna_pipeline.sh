#!/bin/bash

### from qiagen DNA pipeline github, https://github.com/qiaseq/qiaseq-dna ###

# make sure docker engine is running!

# fastqs must be named as 'sample_name_R1_001.fastq.gz' and 
# 'sample_name_R2_001.fastq.gz'

# define and create directories:
sample_name="409_021_DBV4V_CTCTCTAC-CTCTCTAT_L001"
capture_id="CDHS-34925Z-409"

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/ewing_ctDNA"
cont_dir="$project_dir/containers"
working_dir="$project_dir/raw_files/$sample_name"
data_dir="$project_dir/data"
code_dir="$project_dir/code"

mkdir -p $cont_dir
mkdir -p $code_dir
mkdir -p $working_dir

# update params file with sample_name and copy to working_dir:
cp $code_dir/qiaseq-dna/run_sm_counter_v2.params.txt \
  $working_dir/template_params.txt
cat $working_dir/template_params.txt | \
  sed "s/insert_sample_name/$sample_name/g" | \
  sed "s/capture_id/$capture_id/g" \
  > $working_dir/run_sm_counter_v2.params.txt

# shell into container, do not change mounted dir names as this breaks the pipeline!
singularity shell \
  -B $code_dir:/srv/qgen/code/ \
  -B $data_dir:/srv/qgen/data/ \
  -B $working_dir:/srv/qgen/example/ \
  $cont_dir/qiaseq-dna_latest.sif

# change to run directory:
cd /srv/qgen/example
#cp /srv/qgen/code/qiaseq-dna/run_sm_counter_v2.params.txt ./

# edit the bottom of run_consensus.params.txt if you need to change the read set 
# and primer file:

# run the pipeline:
python /srv/qgen/code/qiaseq-dna/run_qiaseq_dna.py \
  run_sm_counter_v2.params.txt \
  v2 \
  single \
  409_021_DBV4V_CTCTCTAC-CTCTCTAT_L001 > run.log 2>&1
exit

