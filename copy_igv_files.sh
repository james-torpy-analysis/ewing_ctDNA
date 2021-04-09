
for f in $(ls results/gsnap); do

  sample=$(echo $f | sed "s/\///")
  echo "copying $sample to igv dir"
  echo ------

  mkdir -p results/igv/$sample

  \cp results/gsnap/$sample/*unpaired*bam* results/igv/$sample
  \cp results/gsnap/$sample/*all.sorted.bam* results/igv/$sample

  \cp results/svaba/gsnap/$sample/$sample.svaba.unfiltered.sv.formatted.vcf \
  results/igv/$sample/$sample.svaba.unfiltered.sv.formatted.gsnap.vcf
  \cp results/svaba/gsnap/$sample/$sample.svaba.unfiltered.sv.formatted.vcf.idx \
  results/igv/$sample/$sample.svaba.unfiltered.sv.formatted.gsnap.vcf.idx
  \cp results/svaba/bwa/$sample/$sample.svaba.unfiltered.sv.formatted.vcf \
  results/igv/$sample/$sample.svaba.unfiltered.sv.formatted.bwa.vcf
  \cp results/svaba/bwa/$sample/$sample.svaba.unfiltered.sv.formatted.vcf.idx \
  results/igv/$sample/$sample.svaba.unfiltered.sv.formatted.bwa.vcf.idx

  \cp results/svaba/gsnap/$sample/$sample.svaba.sv.vcf \
  results/igv/$sample/$sample.svaba.sv.gsnap.vcf
  \cp results/svaba/gsnap/$sample/$sample.svaba.sv.vcf.idx \
  results/igv/$sample/$sample.svaba.sv.gsnap.vcf.idx
  \cp results/svaba/bwa/$sample/$sample.svaba.sv.vcf \
  results/igv/$sample/$sample.svaba.sv.bwa.vcf
  \cp results/svaba/bwa/$sample/$sample.svaba.sv.vcf.idx \
  results/igv/$sample/$sample.svaba.sv.bwa.vcf.idx

done;