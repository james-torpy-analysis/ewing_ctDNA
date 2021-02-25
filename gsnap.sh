
## build genome index if it does not exist already:
gmap_build -d GRCh37.p13.genome -D /share/ScratchGeneral/jamtor/projects/ewing_ctDNA/genome \
/share/ScratchGeneral/jamtor/projects/ewing_ctDNA/genome/GRCh37.p13.genome.fa

conda activate general

# map reads
gsnap -D /share/ScratchGeneral/jamtor/projects/ewing_ctDNA/genome/ -d GRCh38.primary_assembly.genome --gunzip -t 6 -A sam \
  --find-dna-chimeras 1 -o /share/ScratchGeneral/jamtor/projects/ewing_ctDNA/results/gsnap/409_001_D9YW9_TCCTGAGC-CTCTCTAT_L001_R1.out \
  /share/ScratchGeneral/jamtor/projects/ewing_ctDNA/raw_files/409_001_D9YW9_TCCTGAGC-CTCTCTAT_L001_R1.fastq.gz \
  /share/ScratchGeneral/jamtor/projects/ewing_ctDNA/raw_files/409_001_D9YW9_TCCTGAGC-CTCTCTAT_L001_R2.fastq.gz

  gsnap -D /share/ScratchGeneral/jamtor/projects/ewing_ctDNA/genome/ -d GRCh37.p13.genome --gunzip -t 6 -A sam \
  --find-dna-chimeras 1 -o results/gsnap/ raw_files/409_001_D9YW9_TCCTGAGC-CTCTCTAT_L001_R1.fastq.gz \
  raw_files/409_001_D9YW9_TCCTGAGC-CTCTCTAT_L001_R2.fastq.gz