
cd "/share/ScratchGeneral/jamtor/projects/ewing_ctDNA/clinSV"

raw_dir="/share/ScratchGeneral/jamtor/projects/ewing_ctDNA/raw_files"
ref_dir="/share/ScratchGeneral/jamtor/projects/ewing_ctDNA/clinSV/refs/refdata-b37"
out_dir="results/test"

mkdir -p $out_dir

singularity run clinsv.sif \
  -i $raw_dir/NA12878_v0.9.bam \
  -ref $ref_dir \
  -p $out_dir