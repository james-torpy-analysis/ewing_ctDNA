
### build genome index if it does not exist already:
#gmap_build -d GRCh37.p13.genome -D /share/ScratchGeneral/jamtor/projects/ewing_ctDNA/genome \
#/share/ScratchGeneral/jamtor/projects/ewing_ctDNA/genome/GRCh37.p13.genome.fa

snk

# map reads
gsnap -D /share/ScratchGeneral/jamtor/projects/ewing_ctDNA/genome/ -d GRCh37.p13.genome --gunzip -t 6 -A sam \
  --find-dna-chimeras 1 --split-output /share/ScratchGeneral/jamtor/projects/ewing_ctDNA/results/gsnap/409_002_D9YW9_GGACTCCT-CTCTCTAT_L001/409_002_D9YW9_GGACTCCT-CTCTCTAT_L001 \
  /share/ScratchGeneral/jamtor/projects/ewing_ctDNA/raw_files/409_002_D9YW9_GGACTCCT-CTCTCTAT_L001_R1.fastq.gz \
  /share/ScratchGeneral/jamtor/projects/ewing_ctDNA/raw_files/409_002_D9YW9_GGACTCCT-CTCTCTAT_L001_R2.fastq.gz
