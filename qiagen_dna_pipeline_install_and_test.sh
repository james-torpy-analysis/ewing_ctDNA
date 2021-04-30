#!/bin/bash

### from qiagen DNA pipeline github, https://github.com/qiaseq/qiaseq-dna ###

# make sure docker engine is running!

# define and create directories:
home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/ewing_ctDNA"
cont_dir="$project_dir/containers"
working_dir="$project_dir/raw_files/NEB_S2"
data_dir="$project_dir/data"
code_dir="$project_dir/code"

mkdir -p $cont_dir
mkdir -p $code_dir
mkdir -p $working_dir

# pull the docker image:
qrsh -l mem_requested=16G
cd "/share/ScratchGeneral/jamtor/projects/ewing_ctDNA/containers"
singularity pull docker://qiaseq/qiaseq-dna
exit

# activate conda environment:
source ~/.bashrc
snk

# install gsutil:
pip install gsutil

# get data dependencies:
gsutil -m cp -r gs://qiaseq-dna/data ./

# download test files to raw data dir:
cd $working_dir
wget https://storage.googleapis.com/qiaseq-dna/example/NEB_S2_L001_R1_001.fastq.gz \
https://storage.googleapis.com/qiaseq-dna/example/NEB_S2_L001_R2_001.fastq.gz \
https://storage.googleapis.com/qiaseq-dna/example/DHS-101Z.primers.txt \
https://storage.googleapis.com/qiaseq-dna/example/DHS-101Z.roi.bed ./

# get the latest code from github
cd $code_dir
git clone --recursive https://github.com/qiaseq/qiaseq-dna.git

# change to run directory and copy over parameters file
cd $working_dir
cp $code_dir/qiaseq-dna/run_sm_counter_v2.params.txt ./

## execute python run script using container:
#singularity exec \
#  -B $code_dir:/srv/qgen/code/ \
#  -B $data_dir:/srv/qgen/data/ \
#  -B $working_dir:/srv/qgen/example/ \
#  $cont_dir/qiaseq-dna_latest.sif \
#  python /srv/qgen/code/run_qiaseq_dna.py \
#  /srv/qgen/example/run_sm_counter_v2.params.txt \
#  v2 \
#  single \
#  /srv/qgen/example/NEB_S2 \
#  > run.log 2>&1

# shell into container, do not change mounted dir names as this breaks the pipeline!
singularity shell \
  -B $code_dir:/srv/qgen/code/ \
  -B $data_dir:/srv/qgen/data/ \
  -B $working_dir:/srv/qgen/example/ \
  $cont_dir/qiaseq-dna_latest.sif

# change to run directory and copy over parameters file:
cd /srv/qgen/example
cp /srv/qgen/code/qiaseq-dna/run_sm_counter_v2.params.txt ./

# edit the bottom of run_consensus.params.txt if you need to change the read set 
# and primer file:

# run the pipeline:
python /srv/qgen/code/qiaseq-dna/run_qiaseq_dna.py \
  run_sm_counter_v2.params.txt \
  v2 \
  single \
  NEB_S2 > run.log 2>&1
exit

