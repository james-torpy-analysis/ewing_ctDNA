# Run command:
# snakemake --reason --cores 100 --cluster-config cluster.json --cluster 'qsub -pe smp {cluster.cores} -N ewfus.smk -wd '/share/ScratchGeneral/jamtor/projects/ewing_ctDNA/logs' -b y -j y -V -P DSGClinicalGenomics' -j 23

# DAG command:
# snakemake --dag | dot -Tsvg > dag.svg

### This script remaps bams to hg19 and uses SvABA and Manta to identify 
# breakpoints in genomic data ###


# define directories:
project_name = 'ewing_ctDNA'

# define/create directories:
home_dir = '/share/ScratchGeneral/jamtor/'
project_dir = home_dir + 'projects/' + project_name + '/'
results_dir = project_dir + 'results/'
genome_dir = project_dir + 'genome/'
script_dir = project_dir + 'scripts/'

conda_dir = '/share/ClusterShare/thingamajigs/jamtor/local/lib/miniconda3/'
env_dir = conda_dir + 'envs/snkenv/bin/'

fq_dir = 'raw_files/'
align_dir = 'results/BWA_and_picard/bams/'
svaba_dir = 'results/svaba/BWA_and_picard/'
int_dir = 'results/BWA_and_picard/int_bams/'
fusion_dir = 'results/fusions/'
VAF_dir = 'results/VAF_calculation/'

R_dir = "/share/ClusterShare/thingamajigs/jamtor/local/lib/miniconda3/envs/snkenv/bin/"


SAMPLES = list([
    '409_001_D9YW9_TCCTGAGC-CTCTCTAT_L001', '409_002_D9YW9_GGACTCCT-CTCTCTAT_L001', 
    '409_003_D9YWF_AGGCAGAA-CTCTCTAT_L001', '409_004_D9YWF_GTAGAGGA-CTCTCTAT_L001', 
    '409_005_D9YWF_ATCTCAGG-CTCTCTAT_L001', '409_006_DB62M_GGACTCCT-CTCTCTAT_L001', 
    '409_007_DB62M_TAGGCATG-CTCTCTAT_L001', '409_008_DB62M_GTAGAGGA-CTCTCTAT_L001', 
    '409_009_DB62M_GCTCATGA-CTCTCTAT_L001', '409_010_DB62M_ATCTCAGG-CTCTCTAT_L001', 
    '409_011_DBV4V_TAAGGCGA-CTCTCTAT_L001', '409_012_DBV4V_CGTACTAG-CTCTCTAT_L001', 
    '409_013_DBV4V_AGGCAGAA-CTCTCTAT_L001', '409_014_DBV4V_TCCTGAGC-CTCTCTAT_L001', 
    '409_015_DBV4V_GGACTCCT-CTCTCTAT_L001', '409_016_DBV4V_TAGGCATG-CTCTCTAT_L001', 
    '409_017_DBV4V_CGAGGCTG-CTCTCTAT_L001', '409_018_DBV4V_AAGAGGCA-CTCTCTAT_L001', 
    '409_019_DBV4V_GCTCATGA-CTCTCTAT_L001', '409_020_DBV4V_GTAGAGGA-CTCTCTAT_L001', 
    '409_021_DBV4V_CTCTCTAC-CTCTCTAT_L001', '409_022_DCB8V_TAAGGCGA-CTCTCTAT_L001', 
    '409_023_DCB8V_CGTACTAG-CTCTCTAT_L001', '409_024_DCB8V_AGGCAGAA-CTCTCTAT_L001', 
    '409_025_DCB8V_TCCTGAGC-CTCTCTAT_L001', '409_026_DCB8V_GGACTCCT-CTCTCTAT_L001',  
    '409_027_DCKVC_TAGGCATG-CTCTCTAT_L001', '409_027_DCKVC_TAGGCATG-CTCTCTAT_L001', 
    '409_028_DCB8V_CTCTCTAC-CTCTCTAT_L001',  
    '409_029_DCB8V_CGAGGCTG-CTCTCTAT_L001', '409_030_DCB8V_AAGAGGCA-CTCTCTAT_L001', 
    '409_031_DCB8V_GTAGAGGA-CTCTCTAT_L001', '409_032_DCB94_GGACTCCT-CTCTCTAT_L001', 
    '409_033_DCB94_TAGGCATG-CTCTCTAT_L001', '409_034_DCB94_CTCTCTAC-CTCTCTAT_L001', 
    '409_035_DCB94_CGAGGCTG-CTCTCTAT_L001', '409_036_DCB94_AAGAGGCA-CTCTCTAT_L001', 
    '409_037_DCB94_GTAGAGGA-CTCTCTAT_L001', '409_038_DCB8V_GCTCATGA-CTCTCTAT_L001', 
    '409_039_DCB8V_ATCTCAGG-CTCTCTAT_L001', '409_040_DCKVC_GGACTCCT-CTCTCTAT_L001', 
    '409_041_DCCT9_TAGGCATG-CTCTCTAT_L001', '409_042_DCKVC_CTCTCTAC-CTCTCTAT_L001', 
    '409_043_DCKVC_CGAGGCTG-CTCTCTAT_L001', '409_044_DCCT9_AAGAGGCA-CTCTCTAT_L001', 
    '409_045_DCCT9_GTAGAGGA-CTCTCTAT_L001', '409_046_DCCT9_GCTCATGA-CTCTCTAT_L001', 
    '409_047_DCCT9_ATCTCAGG-CTCTCTAT_L001', '409_048_DCB94_TAAGGCGA-CTCTCTAT_L001', 
    '409_049_DCB94_CGTACTAG-CTCTCTAT_L001', '409_050_DCB94_AGGCAGAA-CTCTCTAT_L001', 
    '409_051_DCB94_TCCTGAGC-CTCTCTAT_L001', '409_052_DCB94_GCTCATGA-CTCTCTAT_L001', 
    '409_053_DCB94_ATCTCAGG-CTCTCTAT_L001', '409_054_DCKVC_TAAGGCGA-CTCTCTAT_L001', 
    '409_055_DCCT9_CGTACTAG-CTCTCTAT_L001', '409_056_DCKVC_AGGCAGAA-CTCTCTAT_L001', 
    '409_057_DCKVC_TCCTGAGC-CTCTCTAT_L001', '409_058_DCCT9_GGACTCCT-CTCTCTAT_L001', 
    '409_059_DCCT9_CTCTCTAC-CTCTCTAT_L001', '409_060_DCCT9_TAAGGCGA-CTCTCTAT_L001', 
    '409_061_DCCT9_AGGCAGAA-CTCTCTAT_L001', '409_062_DCCT9_CGAGGCTG-CTCTCTAT_L001', 
    '409_063_DCCT9_TCCTGAGC-CTCTCTAT_L001',
    '409_065_DCKVC_CGTACTAG-CTCTCTAT_L001', '409_066_DCKVC_AAGAGGCA-CTCTCTAT_L001', 
    '409_067_DCKVC_GTAGAGGA-CTCTCTAT_L001', '409_068_DCKVC_GCTCATGA-CTCTCTAT_L001', 
    '409_069_DCKVC_ATCTCAGG-CTCTCTAT_L001'
])

#SAMPLES = list([
#    '409_027_DCKVC_TAGGCATG-CTCTCTAT_L001'
#])

#SAMPLES = list([
#    '409_001_D9YW9_TCCTGAGC-CTCTCTAT_L001', '409_032_DCB94_GGACTCCT-CTCTCTAT_L001',
#    '409_002_D9YW9_GGACTCCT-CTCTCTAT_L001', '409_040_DCKVC_GGACTCCT-CTCTCTAT_L001', 
#    '409_041_DCCT9_TAGGCATG-CTCTCTAT_L001', '409_065_DCKVC_CGTACTAG-CTCTCTAT_L001',
#    '409_066_DCKVC_AAGAGGCA-CTCTCTAT_L001', '409_018_DBV4V_AAGAGGCA-CTCTCTAT_L001', 
#    '409_019_DBV4V_GCTCATGA-CTCTCTAT_L001', '409_033_DCB94_TAGGCATG-CTCTCTAT_L001'
#])

## ANZCHOG abstract:
#SAMPLES = list([
#	'409_016_DBV4V_TAGGCATG-CTCTCTAT_L001', '409_018_DBV4V_AAGAGGCA-CTCTCTAT_L001', 
#	'409_031_DCB8V_GTAGAGGA-CTCTCTAT_L001', '409_014_DBV4V_TCCTGAGC-CTCTCTAT_L001',
#	'409_021_DBV4V_CTCTCTAC-CTCTCTAT_L001', '409_005_D9YWF_ATCTCAGG-CTCTCTAT_L001',
#	'409_019_DBV4V_GCTCATGA-CTCTCTAT_L001', '409_025_DCB8V_TCCTGAGC-CTCTCTAT_L001',
#	'409_007_DB62M_TAGGCATG-CTCTCTAT_L001', '409_006_DB62M_GGACTCCT-CTCTCTAT_L001',
#	'409_012_DBV4V_CGTACTAG-CTCTCTAT_L001', '409_013_DBV4V_AGGCAGAA-CTCTCTAT_L001',
#	'409_008_DB62M_GTAGAGGA-CTCTCTAT_L001','409_020_DBV4V_GTAGAGGA-CTCTCTAT_L001', 
#	'409_015_DBV4V_GGACTCCT-CTCTCTAT_L001',
#	'409_023_DCB8V_CGTACTAG-CTCTCTAT_L001', '409_027_DCB8V_TAGGCATG-CTCTCTAT_L001',
#	'409_062_DCCT9_CGAGGCTG-CTCTCTAT_L001', '409_060_DCCT9_TAAGGCGA-CTCTCTAT_L001',
#	'409_024_DCB8V_AGGCAGAA-CTCTCTAT_L001', '409_026_DCB8V_GGACTCCT-CTCTCTAT_L001',
#	'409_028_DCB8V_CTCTCTAC-CTCTCTAT_L001', '409_055_DCCT9_CGTACTAG-CTCTCTAT_L001',
#	'409_004_D9YWF_GTAGAGGA-CTCTCTAT_L001', '409_009_DB62M_GCTCATGA-CTCTCTAT_L001',
#	'409_011_DBV4V_TAAGGCGA-CTCTCTAT_L001', '409_059_DCCT9_CTCTCTAC-CTCTCTAT_L001',
#	'409_058_DCCT9_GGACTCCT-CTCTCTAT_L001', '409_061_DCCT9_AGGCAGAA-CTCTCTAT_L001',
#	'409_063_DCCT9_TCCTGAGC-CTCTCTAT_L001'
#])

## testing:
#SAMPLES = list([
#    '409_060_DCCT9_TAAGGCGA-CTCTCTAT_L001'
#])

#TYPE = list([
#    'SV', 'all'
#])


#rule all:
#    input:
#        expand(
#            align_dir + '{sample}/{sample}.consensus.bam.bai',
#            sample = SAMPLES
#        )

#rule all:
#    input:
#        expand(
#            svaba_dir + '{sample}/{sample}.svaba.sv.vcf.idx',
#            sample=SAMPLES
#        )

#rule all:
#    input:
#        expand(
#            'logs/completed_jobs/{sample}_complete',
#            sample = SAMPLES
#        )

#rule all:
#    input:
#        expand(
#            fusion_dir + '{sample}/EWSR1_GOI_fusions.Rdata',
#            sample = SAMPLES
#        )

rule all:
    input:
        expand(
            VAF_dir + '{sample}/VAFs.txt',
            sample = SAMPLES
        )


######################################################################################################
### 1. UMI collapse and BWA ###
######################################################################################################

rule BWA_and_umi_collapse:
    input:
        fq1 = fq_dir + '{sample}/{sample}_R1.fastq.gz',
        fq2 = fq_dir + '{sample}/{sample}_R2.fastq.gz'
    output:
        bam = align_dir + '{sample}/{sample}.consensus.bam',
        bai = align_dir + '{sample}/{sample}.consensus.bam.bai',
    threads: 7
    shell:
        'mkdir -p logs/BWA_and_picard; ' +
        'cd logs/BWA_and_picard; ' + 
        script_dir + '1.UMI_collapse.sh' +
            ' {wildcards.sample}' +
            ' 2>&1 {wildcards.sample}.alignment.log'


######################################################################################################
### 2. SvABA ###
######################################################################################################

rule svaba:
   input:
       bam = align_dir + '{sample}/{sample}.consensus.bam',
       bai = align_dir + '{sample}/{sample}.consensus.bam.bai'
   output:
       filt = svaba_dir + '{sample}/{sample}.svaba.sv.vcf',
       unfilt = svaba_dir + '{sample}/{sample}.svaba.unfiltered.sv.vcf'
   threads: 7
   shell:
       'mkdir -p logs/svaba; ' + 
        'cd logs/svaba; ' + 
        'mkdir -p ../../' + svaba_dir + '{wildcards.sample}/; ' +
        'svaba run -t ' + project_dir + '{input.bam} -G ' + 
            genome_dir + 'GRCh37.p13.genome.fa -a ../../' + 
            svaba_dir + 
            '{wildcards.sample}/{wildcards.sample}' + 
            ' -p 6 --override-reference-check' +
            ' --min-overlap 0.13' +
            ' 2> {wildcards.sample}.svaba.errors'

rule format_vcf:
    input:
        svaba_dir + '{sample}/{sample}.svaba.unfiltered.sv.vcf'
    output:
        svaba_dir + '{sample}/{sample}.svaba.unfiltered.sv.formatted.vcf'
    threads: 1
    shell:
        'scripts/fix_broken_svaba_vcf.sh {wildcards.sample} {input}'

rule remove_duds:
    input:
        svaba_dir + '{sample}/{sample}.svaba.unfiltered.sv.formatted.vcf'
    output:
        svaba_dir + '{sample}/{sample}.svaba.semifiltered.sv.formatted.vcf'
    threads: 1
    shell:
        'grep -v NODISC {input} > {output}'

rule vcf_index:
   input:
       filt = svaba_dir + '{sample}/{sample}.svaba.sv.vcf',
       semifilt = svaba_dir + '{sample}/{sample}.svaba.semifiltered.sv.formatted.vcf'
   output:
        filt = svaba_dir + '{sample}/{sample}.svaba.sv.vcf.idx',
        semifilt = svaba_dir + '{sample}/{sample}.svaba.semifiltered.sv.formatted.vcf.idx'
   threads: 2
   shell:
        env_dir + 'igvtools index {input.filt}; ' + 
        env_dir + 'igvtools index {input.semifilt}'

    
######################################################################################################
### 3. Clean up ###
######################################################################################################

rule cleanup:
    input:
        filt = svaba_dir + '{sample}/{sample}.svaba.sv.vcf.idx',
        unfilt = svaba_dir + '{sample}/{sample}.svaba.semifiltered.sv.formatted.vcf.idx'
    output:
        'logs/completed_jobs/{sample}_complete'
    threads: 1
    shell:
        'rm -fr ' + int_dir + '{wildcards.sample}; '
        'touch {output}'


######################################################################################################
### 4. Find fusions ###
######################################################################################################

rule find_fusions:
    input:
        filt = svaba_dir + '{sample}/{sample}.svaba.sv.vcf.idx',
        unfilt = svaba_dir + '{sample}/{sample}.svaba.semifiltered.sv.formatted.vcf.idx'
    output:
        fusion_dir + '{sample}/EWSR1_GOI_fusions.Rdata'
    threads: 7
    shell:
        "mkdir -p logs/find_fusions/{wildcards.sample}/; " + 
        "cd logs/find_fusions/{wildcards.sample}/; " +
        "{R_dir}/R CMD BATCH  --no-save '--args" + 
        " {wildcards.sample}" + 
        "' ../../../scripts/2.find_EWSR1_fusions.R"


######################################################################################################
### 5. Filter bams ###
######################################################################################################

rule filter_bams:
    input:
        fusion_dir + '{sample}/EWSR1_GOI_fusions.Rdata'
    output:
        VAF_dir + '{sample}/Rdata/VAF_calculation_reads.Rdata'
    threads: 8
    shell:
        "mkdir -p logs/filter_bams/{wildcards.sample}/; " + 
        "cd logs/filter_bams/{wildcards.sample}/; " +
        "{R_dir}/R CMD BATCH  --no-save '--args" + 
        " {wildcards.sample}" + 
        "' ../../../scripts/3.filter_bams.R"


######################################################################################################
### 6. Calculate VAFs ###
######################################################################################################

rule calc_VAFs:
    input:
        VAF_dir + '{sample}/Rdata/VAF_calculation_reads.Rdata'
    output:
        VAF = VAF_dir + '{sample}/VAFs.txt',
        read_no = VAF_dir + '{sample}/final_fusion_read_nos.txt'
    threads: 8
    shell:
        "mkdir -p logs/VAF_calculation/{wildcards.sample}/; " + 
        "cd logs/VAF_calculation/{wildcards.sample}/; " +
        "{R_dir}/R CMD BATCH  --no-save '--args" + 
        " {wildcards.sample}" + 
        "' ../../../scripts/4.calculate_VAFs.R"


