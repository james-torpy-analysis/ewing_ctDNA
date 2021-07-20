# cell line 1 (A673):

######
# use R:
# load fusions:
fusions <- readRDS("/share/ScratchGeneral/jamtor/projects/ewing_ctDNA/results/fusions/409_001_D9YW9_TCCTGAGC-CTCTCTAT_L001/EWSR1_GOI_fusions.Rdata")
fusions$high_conf_bp$true_positives$fusions$FLI1
######

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/ewing_ctDNA"
out_dir="$project_dir/results/fusions/409_001_D9YW9_TCCTGAGC-CTCTCTAT_L001"
genome_dir="$project_dir/genome"

# isolate chr22 fusion sequences:
grep -A 1 chr22 $genome_dir/GRCh37.p13.genome.fa | tail -1 | \
  cut -c 29684199-29684299 > $out_dir/upstream_of_fusion_chr22_29684299_100bp.txt

grep -A 1 chr22 $genome_dir/GRCh37.p13.genome.fa | tail -1 | \
  cut -c 29684299-29684399 > $out_dir/downstream_of_fusion_chr22_29684299_100bp.txt

# isolate chr11 fusion sequences:
grep -A 1 chr11 $genome_dir/GRCh37.p13.genome.fa | tail -1 | \
  cut -c 128667404-128667505 > $out_dir/upstream_of_fusion_chr11_128667505_100bp.txt

grep -A 1 chr11 $genome_dir/GRCh37.p13.genome.fa | tail -1 | \
  cut -c 128667505-128667605 > $out_dir/downstream_of_fusion_chr11_128667505_100bp.txt


# cell line 2 (ES8):

######
# use R:
# load fusions:
fusions <- readRDS("/share/ScratchGeneral/jamtor/projects/ewing_ctDNA/results/fusions/409_002_D9YW9_GGACTCCT-CTCTCTAT_L001/EWSR1_GOI_fusions.Rdata")
fusions$high_conf_bp$true_positives$fusions$FLI1
######

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/ewing_ctDNA"
out_dir="$project_dir/results/fusions/409_002_D9YW9_GGACTCCT-CTCTCTAT_L001"
genome_dir="$project_dir/genome"

# isolate chr22 fusion sequences:
grep -A 1 chr22 $genome_dir/GRCh37.p13.genome.fa | tail -1 | \
  cut -c 29685836-29685936 > $out_dir/upstream_of_fusion_chr22_29685936_100bp.txt

grep -A 1 chr22 $genome_dir/GRCh37.p13.genome.fa | tail -1 | \
  cut -c 29685936-29686036 > $out_dir/downstream_of_fusion_chr22_29685936_100bp.txt

# isolate chr11 fusion sequences:
grep -A 1 chr11 $genome_dir/GRCh37.p13.genome.fa | tail -1 | \
  cut -c 128645239-128645339 > $out_dir/upstream_of_fusion_chr11_128645339_100bp.txt

grep -A 1 chr11 $genome_dir/GRCh37.p13.genome.fa | tail -1 | \
  cut -c 128645339-128645439 > $out_dir/downstream_of_fusion_chr11_128645339_100bp.txt
