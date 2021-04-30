#!/bin/bash

source ~/.bashrc

sample_name="409_021_DBV4V_CTCTCTAC-CTCTCTAT_L001"

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/ewing_ctDNA"
fq_dir="$project_dir/raw_files"
result_dir="$project_dir/results"
int_dir="$result_dir/picard/int_bams/$sample_name"
fq_out_dir="$result_dir/picard/fastqs//$sample_name"
bam_dir="$result_dir/picard/bams/$sample_name"
genome_dir="$project_dir/genome"
stats_dir="$result_dir/read_stats/$sample_name"
script_dir="$project_dir/scripts"

picard_dir="$home_dir/local/bin"
fgbio_dir="$home_dir/local/lib/fgbio/target/scala-2.13"

mkdir -p $int_dir
mkdir -p $fq_out_dir
mkdir -p $bam_dir
mkdir -p $stats_dir

# call snk conda env for samtools/python:
conda activate snkenv

# make function to sort and index bams:
sort_and_index () {
  pref=$(echo $1 | sed "s/\.bam//g")
  samtools sort -o $pref.sorted.bam $pref.bam 
  samtools index $pref.sorted.bam
}

# make function to filter discordant and split reads:
filter_discordant () {
  pref=$(echo $1 | sed "s/\.bam//g")

  sort_and_index $1
  
  # filter and index discordant reads:
  samtools view -h -F 1294 $pref.sorted.bam | samtools view -bh > $pref.discordant.bam
  samtools index $pref.discordant.bam

  # filter and index split reads:
  samtools view -h -f 2048 $pref.sorted.bam | samtools view -bh > $pref.split.bam
  samtools index $pref.split.bam
}


if [ ! -f $fq_out_dir/$sample_name.withUMI.fastq ]; then

  printf "\n\n"
  echo "--------------------------------------------------"
  echo "converting fastq files to one unmapped bam file..."
  echo "--------------------------------------------------"
  printf "\n"
  java -jar $picard_dir/picard.jar FastqToSam \
    FASTQ=$fq_dir/$sample_name\_R1.fastq.gz \
    FASTQ2=$fq_dir/$sample_name\_R2.fastq.gz \
    O=$int_dir/$sample_name.unmapped.bam \
    SM=sample
  
  printf "\n\n"
  echo "--------------------------------------------------"
  echo "removing sequences < 40 bp long..."
  echo "--------------------------------------------------"
  printf "\n"
  samtools view -H $int_dir/$sample_name.unmapped.bam > \
    $int_dir/$sample_name.unmapped.filtered.sam
  samtools view $int_dir/$sample_name.unmapped.bam | \
    awk ' length($10) > 40 ' >> $int_dir/$sample_name.unmapped.filtered.sam
  samtools view -bh $int_dir/$sample_name.unmapped.filtered.sam > \
    $int_dir/$sample_name.unmapped.filtered.bam && \
    rm $int_dir/$sample_name.unmapped.filtered.sam
  
  printf "\n\n"
  echo "--------------------------------------------------"
  echo "extracting UMI and common sequence from R2..."
  echo "--------------------------------------------------"
  printf "\n"
  java -jar $fgbio_dir/fgbio-1.4.0-468a843-SNAPSHOT.jar ExtractUmisFromBam \
    --input=$int_dir/$sample_name.unmapped.filtered.bam \
    --output=$int_dir/$sample_name.unmapped.withUMI.bam \
    --read-structure=1M149T 23M152T \
    --molecular-index-tags=ZA ZB \
    --single-tag=RX
  
  printf "\n\n"
  echo "--------------------------------------------------"
  echo "converting back to single fastq file..."
  echo "--------------------------------------------------"
  printf "\n"
  java -jar $picard_dir/picard.jar SamToFastq I=$int_dir/$sample_name.unmapped.withUMI.bam \
    F=$fq_out_dir/$sample_name.withUMI.fastq \
    INTERLEAVE=true

else

  printf "\n"
  echo "$fq_out_dir/$sample_name.withUMI.fastq already exists, skipping to first alignment step..."
  printf "\n"
  
fi;


##################################################################################################################################
### 1. Uncollapsed bam ###
##################################################################################################################################

printf "\n\n"
echo "--------------------------------------------------"
echo "aligning fastq..."
echo "--------------------------------------------------"
printf "\n"
bwa mem -p -t 5 $genome_dir/GRCh37.p13.genome.fa \
  $fq_out_dir/$sample_name.withUMI.fastq > $int_dir/$sample_name.initial_mapped.bam

# create reference dƒictionary if needed:
#java -jar $picard_dir/picard.jar CreateSequenceDictionary -R $genome_dir/GRCh37.p13.genome.fa

printf "\n\n"
echo "--------------------------------------------------"
echo "merging UMIs and common sequences to bam file..."
echo "--------------------------------------------------"
printf "\n"
java -jar $picard_dir/picard.jar MergeBamAlignment \
  UNMAPPED=$int_dir/$sample_name.unmapped.withUMI.bam \
  ALIGNED=$int_dir/$sample_name.initial_mapped.bam \
  O=$int_dir/$sample_name.initial_mapped_and_UMI.bam \
  R=$genome_dir/GRCh37.p13.genome.fa \
  SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 \
  ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true

printf "\n\n"
echo "--------------------------------------------------"
echo "removing fake read 1 UMI and common sequence..."
echo "--------------------------------------------------"
printf "\n"
# carry over header:
samtools view -H $int_dir/$sample_name.initial_mapped_and_UMI.bam > \
  $int_dir/$sample_name.initial_mapped_and_UMI.formatted.sam

# remove read 1 UMI and common sequence using python script:
python $script_dir/remove_common_sequence.py $int_dir/$sample_name.initial_mapped_and_UMI.bam

# isolate discordant reads:
samtools view -h -F 1294 $bam_dir/$sample_name.uncollapsed.bam | \
  samtools view -bh >  $bam_dir/$sample_name.uncollapsed.discordant.bam

# isolate split reads:
samtools view -h -f 2048 $bam_dir/$sample_name.uncollapsed.bam | \
  samtools view -bh >  $bam_dir/$sample_name.uncollapsed.split.bam


##################################################################################################################################
### 2. Collapsed bam ###
##################################################################################################################################

printf "\n"
echo "skipping UMI error correction as sequences are random..."
printf "\n"

printf "\n\n"
echo "--------------------------------------------------"
echo "grouping read families, keeping reads with min mapq > 10..."
echo "--------------------------------------------------"
printf "\n"
# edits = the allowable number of edits between UMIs
# min-map-q = minimum mapping quality required to keep reads, consider reducing if split reads are removed
java -jar $fgbio_dir/fgbio-1.4.0-468a843-SNAPSHOT.jar GroupReadsByUmi \
  --input=$bam_dir/$sample_name.EWSR1_FLI1_fusion.uncollapsed.bam \
  --output $int_dir/$sample_name.EWSR1_FLI1_fusion.grouped.bam \
  --strategy=Edit \
  --edits=1 \
  --min-map-q=0 \
  --include-non-pf-reads=true \
  --allow-inter-contig=true

filter_discordant $bam_dir/$sample_name.EWSR1_FLI1_fusion.uncollapsed.bam
filter_discordant $int_dir/$sample_name.EWSR1_FLI1_fusion.grouped.bam

printf "\n\n"
echo "--------------------------------------------------"
echo "collapsing read families..."
echo "--------------------------------------------------"
printf "\n"
# --min-reads is the min number of reads required to make a consensus base
java -jar $fgbio_dir/fgbio-1.4.0-468a843-SNAPSHOT.jar CallMolecularConsensusReads \
  --input=$int_dir/$sample_name.grouped.bam \
  --output=$int_dir/$sample_name.consensus.unmapped.bam \
  --error-rate-pre-umi=45 \
  --error-rate-post-umi=30 \
  --min-consensus-base-quality=40 \
  --min-input-base-quality=30 \
  --min-reads=1

printf "\n\n"
echo "--------------------------------------------------"
echo "converting back to single fastq file..."
echo "--------------------------------------------------"
printf "\n"
java -jar $picard_dir/picard.jar SamToFastq I=$int_dir/$sample_name.consensus.unmapped.bam \
  F=$fq_out_dir/$sample_name.consensus.fastq \
  INTERLEAVE=true

printf "\n\n"
echo "--------------------------------------------------"
echo "remapping collapsed reads..."
echo "--------------------------------------------------"
printf "\n"
bwa mem -p -t 5 $genome_dir/GRCh37.p13.genome.fa \
  $fq_out_dir/$sample_name.consensus.fastq > \
  $int_dir/$sample_name.consensus.mapped.bam

printf "\n\n"
echo "--------------------------------------------------"
echo "merging UMIs and common sequences to bam file..."
echo "--------------------------------------------------"
printf "\n"
java -jar $picard_dir/picard.jar MergeBamAlignment \
  UNMAPPED=$int_dir/$sample_name.consensus.unmapped.bam \
  ALIGNED=$int_dir/$sample_name.consensus.mapped.bam \
  O=$bam_dir/$sample_name.consensus.unfiltered.bam \
  R=$genome_dir/GRCh37.p13.genome.fa \
  SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 \
  ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true


##################################################################################################################################
### 3. Filtered bam ###
##################################################################################################################################

printf "\n\n"
echo "--------------------------------------------------"
echo "filtering out low base quality reads..."
echo "--------------------------------------------------"
printf "\n"
java -jar $fgbio_dir/fgbio-1.4.0-468a843-SNAPSHOT.jar FilterConsensusReads \
  --input=$bam_dir/$sample_name.consensus.unfiltered.bam \
  --output=$bam_dir/$sample_name.consensus.filtered.bam \
  --ref=$genome_dir/GRCh37.p13.genome.fa \
  --min-reads=3 \
  --max-read-error-rate=0.05 \
  --max-base-error-rate=0.1 \
  --min-base-quality=40 \
  --max-no-call-fraction=0.1

printf "\n"
echo "pipeline complete, final bams in $bam_dir"


##################################################################################################################################
### 3. Index bams and filter discordant and split reads ###
##################################################################################################################################

for f in $bam_dir/$sample_name*.bam; do

  samtools index $f

  fprefix=$(echo $f | sed "s/\.bam//")
  
  # filter and index discordant reads:
  samtools view -h -F 1294 $f | samtools view -bh > $fprefix.discordant.bam
  samtools index $fprefix.discordant.bam

  # filter and index split reads:
  samtools view -h -f 2048 $f | samtools view -bh > $fprefix.split.bam
  samtools index $fprefix.split.bam

done


##################################################################################################################################
### 4. Index bams and filter discordant and split reads ###
##################################################################################################################################

######
svaba run -t $bam_dir/$sample_name.consensus.unfiltered.bam -G $genome_dir/GRCh37.p13.genome.fa -a $bam_dir/$sample_name -p 6 \
  --override-reference-check

scripts/fix_broken_svaba_vcf.sh $sample_name $bam_dir/$sample_name.svaba.unfiltered.sv.vcf

igvtools index $bam_dir/$sample_name.svaba.unfiltered.sv.vcf


