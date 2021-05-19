
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

# get the latest code from github
cd $code_dir
git clone --recursive https://github.com/qiaseq/qiaseq-dna.git
