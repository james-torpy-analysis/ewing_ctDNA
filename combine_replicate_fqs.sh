home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/ewing_ctDNA"
raw_dir="$project_dir/raw_files"

samplenames=( "409_004" "409_008" "409_009" "409_012" "409_013" "409_014" "409_015" \
  "409_018" "409_023" "409_024" "409_025" "409_028" "409_031" "409_034" "409_035" \
  "409_043" "409_044" "409_045" "409_048" "409_055" "409_056" "409_057" "409_061" \
  "409_062" "409_063" "409_068" "409_069" )

for s in ${samplenames[@]}; do

  replicates=( $(ls $raw_dir | grep $s | sed "s/\///") )
  mkdir -p $raw_dir/$s\_combined

  for r in ${replicates[@]}; do

    gunzip $raw_dir/$r/$r*.fastq.gz

    wc $raw_dir/$r/$r*R1.fastq >> $raw_dir/$s\_merge.log
    wc $raw_dir/$r/$r*R2.fastq >> $raw_dir/$s\_merge.log

    cat $raw_dir/$r/$r*R1.fastq >> $raw_dir/$s\_combined/$s\_combined_R1.fastq
    cat $raw_dir/$r/$r*R2.fastq >> $raw_dir/$s\_combined/$s\_combined_R2.fastq

    wc $raw_dir/$s\_combined/$s\_combined_R1.fastq >> $raw_dir/$s\_merge.log
    wc $raw_dir/$s\_combined/$s\_combined_R2.fastq >> $raw_dir/$s\_merge.log

  done

  gzip $raw_dir/$s\_combined/*

done