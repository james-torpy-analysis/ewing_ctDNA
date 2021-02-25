project_name="ewing_ctDNA"

sample=$1
bam_dir=$2

home_dir="/share/ScratchGeneral/jamtor/"
project_dir="$home_dir/projects/$project_name/"
manta_dir="results/manta/"
genome_dir="genome/"

source /home/jamtor/.bashrc
conda activate py2.7

mkdir -p $manta_dir/$sample

configManta.py \
--region chr1 --region chr2 --region chr3 \
--region chr4 --region chr5 --region chr6 \
--region chr7 --region chr8 --region chr9 \
--region chr10 --region chr11 --region chr12 \
--region chr13 --region chr14 --region chr15 \
--region chr16 --region chr17 --region chr18 \
--region chr19 --region chr20 --region chr21 \
--region chr22 --region chrX --region chrY \
--tumorBam ../../$bam_dir/$sample.all.sorted.bam \
--referenceFasta ../../$genome_dir/GRCh37.p13.genome.fa \
--runDir $project_dir/$manta_dir/$sample

../../$manta_dir/$sample/runWorkflow.py -m local -j 15
