# Run command:
# snakemake --reason --cores 100 --cluster-config cluster.json --cluster 'qsub -pe smp {cluster.cores} -N ewfus.smk -wd '/share/ScratchGeneral/jamtor/projects/ewing_ctDNA/logs' -b y -j y -V -P DSGClinicalGenomics' -j 23

# DAG command:
# snakemake --dag | dot -Tsvg > dag.svg

### This script remaps bams to hg19 and uses SvABA and Manta to identify 
# breakpoints in genomic data ###

import pandas as pd

# define variables:
project_name = 'ewing_ctDNA'
capture_id = 'CDHS-34925Z-409'

# define/create directories:
home_dir = '/share/ScratchGeneral/jamtor/'
project_dir = home_dir + 'projects/' + project_name + '/'
results_dir = project_dir + 'results/'
genome_dir = project_dir + 'data/genome/hg19/'
script_dir = project_dir + 'scripts/'
ref_dir = project_dir + 'refs/'

conda_dir = '/share/ClusterShare/thingamajigs/jamtor/local/lib/miniconda3/'
env_dir = conda_dir + 'envs/ewing_ctDNA/bin/'

fq_dir = 'raw_files/'
align_dir = 'results/picard/'
svaba_dir = 'results/svaba/picard/'
int_dir = 'results/picard/int_bams/'
fusion_dir = 'results/fusions/'
VAF_dir = 'results/VAF_calculation/'
variant_dir = 'results/smcounter2/'

# fetch library ids from metadata file:
meta = pd.read_csv(ref_dir + 'metadata.tsv', sep='\\t', engine='python')
SAMPLES = list(meta['Library_id'])
#SAMPLES = list(['409_008_combined', '409_012_combined', '409_031_combined', 
#    '409_036_DCB94_AAGAGGCA-CTCTCTAT_L001', '409_048_combined', '409_063_combined'])

rule all:
    input:
        expand(
            VAF_dir + '{sample}/Rdata/VAF.rds',
            sample = SAMPLES
        )


######################################################################################################
### 1. Trim, align, dedup and find variants ###
######################################################################################################

rule find_variants:
    input:
        fq1 = fq_dir + '{sample}/{sample}_R1.fastq.gz',
        fq2 = fq_dir + '{sample}/{sample}_R2.fastq.gz'
    output:
        fq1 = variant_dir + '{sample}/{sample}.prep.R1.fastq',
        fq2 = variant_dir + '{sample}/{sample}.prep.R2.fastq',
        vcf = variant_dir + '{sample}/{sample}.smCounter.anno.vcf'
    threads: 8
    shell:
        "mkdir -p logs/find_variants/{wildcards.sample}/; " + 
        "cd logs/find_variants/{wildcards.sample}/; " +
        " ../../../scripts/1.smcounter2.sh " +
        "{wildcards.sample} " +
        "{capture_id} 2> {wildcards.sample}.smcounter2.errors"


######################################################################################################
### 2. Dedup ###
######################################################################################################

rule dedup:
    input:
        fq1 = variant_dir + '{sample}/{sample}.prep.R1.fastq',
        fq2 = variant_dir + '{sample}/{sample}.prep.R2.fastq',
        vcf = variant_dir + '{sample}/{sample}.smCounter.anno.vcf'
    output:
        bam = align_dir + '{sample}/{sample}.dedup.sorted.by.coord.bam',
        bai = align_dir + '{sample}/{sample}.dedup.sorted.by.coord.bam.bai',
    threads: 8
    shell:
        'mkdir -p logs/dedup/{wildcards.sample}/; ' +
        'cd logs/dedup/{wildcards.sample}/; ' + 
        ' ../../../scripts/2.UMI_dedup.sh' +
            ' {wildcards.sample}' +
            ' 2>&1 {wildcards.sample}.alignment.log'


######################################################################################################
### 3. Detect fusions and calculate VAFs ###
######################################################################################################

rule detect_and_vaf:
    input:
        align_dir + '{sample}/{sample}.dedup.sorted.by.coord.bam'
    output:
        VAF_dir + '{sample}/Rdata/VAF.rds'
    threads: 8
    shell:
        "mkdir -p logs/detect_and_vaf/{wildcards.sample}/; " + 
        "cd logs/detect_and_vaf/{wildcards.sample}/; " +
        "{env_dir}/R CMD BATCH  --no-save '--args" + 
        " {project_name}" + 
        " {wildcards.sample}" + 
        "' ../../../scripts/3.detect_and_vaf.R"


