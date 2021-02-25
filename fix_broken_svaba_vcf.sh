grep \## 409_003_D9YWF_AGGCAGAA-CTCTCTAT_L001.svaba.unfiltered.sv.vcf > 409_003_D9YWF_AGGCAGAA-CTCTCTAT_L001.svaba.unfiltered.sv.formatted.vcf

awk '{print $1, $2, $3, $4, $5, $6, $7, $8}' OFS='\t' 409_003_D9YWF_AGGCAGAA-CTCTCTAT_L001.svaba.unfiltered.sv.temp.vcf >> 409_003_D9YWF_AGGCAGAA-CTCTCTAT_L001.svaba.unfiltered.sv.formatted.vcf