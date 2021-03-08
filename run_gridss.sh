
sample_name="409_002_D9YW9_GGACTCCT-CTCTCTAT_L001"
ref_file="GRCh37.p13.genome.fa"

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/ewing_ctDNA"
result_dir="$project_dir/results"
in_dir="$result_dir/bwa"
out_dir="$result_dir/gridss/bwa/$sample_name"

genome_dir="$project_dir/genome"
script_dir="$project_dir/scripts"

mkdir -p $out_dir

$script_dir/gridss.sh \
  --reference $genome_dir/$ref_file \
  --output $out_dir/$sample_name.vcf.gz \
  --assembly $out_dir/assembly.bam \
  --threads 6 \
  $in_dir/$sample_name.sorted.bam