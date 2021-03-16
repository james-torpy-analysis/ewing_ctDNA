sample_name="409_002_D9YW9_GGACTCCT-CTCTCTAT_L001"
ref_file="GRCh37.p13.genome.fa"

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/ewing_ctDNA"
result_dir="$project_dir/results"
in_dir="$result_dir/bwa"
out_dir="$result_dir/pindel/bwa/$sample_name"
mkdir -p $out_dir

genome_dir="$project_dir/genome"

samtools view -h $in_dir/$sample_name.sorted.bam \
   | sam2pindel - $sample_name\_pindel_input.txt 150 $sample_name 0 Illumina-PairEnd

pindel \
  -f $genome_dir/$ref_file \
  -p $in_dir/$sample_name\_pindel_input.txt \
  -o $out_dir/$sample_name \
  -T 6 \
  --report_breakpoints \
  --report_interchromosomal_events