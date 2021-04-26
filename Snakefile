# Run command:
# snakemake --reason --use-conda --cores 100 --cluster-config cluster.json --cluster 'qsub -pe smp {cluster.cores} -N ewfus.smk -wd '/share/ScratchGeneral/jamtor/projects/ewing_ctDNA/logs' -b y -j y -V -P DSGClinicalGenomics' -j 23

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


conda_dir = "/share/ClusterShare/thingamajigs/jamtor/local/lib/miniconda3/"
env_dir = conda_dir + "envs/snkenv/bin/"

fq_dir = 'raw_files/'
bwa_dir = 'results/bwa/'
gsnap_dir = 'results/gsnap/'
geneglobe_dir = 'results/geneglobe/'
bwa_svaba_dir = 'results/svaba/bwa/'
gsnap_svaba_dir = 'results/svaba/gsnap/'
geneglobe_svaba_dir = 'results/svaba/geneglobe/'

SAMPLES = list([
    "409_001_D9YW9_TCCTGAGC-CTCTCTAT_L001", "409_002_D9YW9_GGACTCCT-CTCTCTAT_L001", 
    "409_003_D9YWF_AGGCAGAA-CTCTCTAT_L001", "409_004_D9YWF_GTAGAGGA-CTCTCTAT_L001", 
    "409_005_D9YWF_ATCTCAGG-CTCTCTAT_L001", "409_006_DB62M_GGACTCCT-CTCTCTAT_L001", 
    "409_007_DB62M_TAGGCATG-CTCTCTAT_L001", "409_008_DB62M_GTAGAGGA-CTCTCTAT_L001", 
    "409_009_DB62M_GCTCATGA-CTCTCTAT_L001", "409_010_DB62M_ATCTCAGG-CTCTCTAT_L001", 
    "409_011_DBV4V_TAAGGCGA-CTCTCTAT_L001", "409_012_DBV4V_CGTACTAG-CTCTCTAT_L001", 
    "409_013_DBV4V_AGGCAGAA-CTCTCTAT_L001", "409_014_DBV4V_TCCTGAGC-CTCTCTAT_L001", 
    "409_015_DBV4V_GGACTCCT-CTCTCTAT_L001", "409_016_DBV4V_TAGGCATG-CTCTCTAT_L001", 
    "409_017_DBV4V_CGAGGCTG-CTCTCTAT_L001", "409_018_DBV4V_AAGAGGCA-CTCTCTAT_L001", 
    "409_019_DBV4V_GCTCATGA-CTCTCTAT_L001", "409_020_DBV4V_GTAGAGGA-CTCTCTAT_L001", 
    "409_021_DBV4V_CTCTCTAC-CTCTCTAT_L001", "409_022_DCB8V_TAAGGCGA-CTCTCTAT_L001", 
    "409_023_DCB8V_CGTACTAG-CTCTCTAT_L001", "409_024_DCB8V_AGGCAGAA-CTCTCTAT_L001", 
    "409_025_DCB8V_TCCTGAGC-CTCTCTAT_L001", "409_026_DCB8V_GGACTCCT-CTCTCTAT_L001",  
    "409_027_DCB8V_TAGGCATG-CTCTCTAT_L001", "409_027_DCKVC_TAGGCATG-CTCTCTAT_L001", 
    "409_028_DCB8V_CTCTCTAC-CTCTCTAT_L001",  
    "409_029_DCB8V_CGAGGCTG-CTCTCTAT_L001", "409_030_DCB8V_AAGAGGCA-CTCTCTAT_L001", 
    "409_031_DCB8V_GTAGAGGA-CTCTCTAT_L001", "409_032_DCB94_GGACTCCT-CTCTCTAT_L001", 
    "409_033_DCB94_TAGGCATG-CTCTCTAT_L001", "409_034_DCB94_CTCTCTAC-CTCTCTAT_L001", 
    "409_035_DCB94_CGAGGCTG-CTCTCTAT_L001", "409_036_DCB94_AAGAGGCA-CTCTCTAT_L001", 
    "409_037_DCB94_GTAGAGGA-CTCTCTAT_L001", "409_038_DCB8V_GCTCATGA-CTCTCTAT_L001", 
    "409_039_DCB8V_ATCTCAGG-CTCTCTAT_L001", "409_040_DCKVC_GGACTCCT-CTCTCTAT_L001", 
    "409_041_DCCT9_TAGGCATG-CTCTCTAT_L001", "409_042_DCKVC_CTCTCTAC-CTCTCTAT_L001", 
    "409_043_DCKVC_CGAGGCTG-CTCTCTAT_L001", "409_044_DCCT9_AAGAGGCA-CTCTCTAT_L001", 
    "409_045_DCCT9_GTAGAGGA-CTCTCTAT_L001", "409_046_DCCT9_GCTCATGA-CTCTCTAT_L001", 
    "409_047_DCCT9_ATCTCAGG-CTCTCTAT_L001", "409_048_DCB94_TAAGGCGA-CTCTCTAT_L001", 
    "409_049_DCB94_CGTACTAG-CTCTCTAT_L001", "409_050_DCB94_AGGCAGAA-CTCTCTAT_L001", 
    "409_051_DCB94_TCCTGAGC-CTCTCTAT_L001", "409_052_DCB94_GCTCATGA-CTCTCTAT_L001", 
    "409_053_DCB94_ATCTCAGG-CTCTCTAT_L001", "409_054_DCKVC_TAAGGCGA-CTCTCTAT_L001", 
    "409_055_DCCT9_CGTACTAG-CTCTCTAT_L001", "409_056_DCKVC_AGGCAGAA-CTCTCTAT_L001", 
    "409_057_DCKVC_TCCTGAGC-CTCTCTAT_L001", "409_058_DCCT9_GGACTCCT-CTCTCTAT_L001", 
    "409_059_DCCT9_CTCTCTAC-CTCTCTAT_L001", "409_060_DCCT9_TAAGGCGA-CTCTCTAT_L001", 
    "409_061_DCCT9_AGGCAGAA-CTCTCTAT_L001", "409_062_DCCT9_CGAGGCTG-CTCTCTAT_L001", 
    "409_063_DCCT9_TCCTGAGC-CTCTCTAT_L001",
    "409_065_DCKVC_CGTACTAG-CTCTCTAT_L001", "409_066_DCKVC_AAGAGGCA-CTCTCTAT_L001", 
    "409_067_DCKVC_GTAGAGGA-CTCTCTAT_L001", "409_068_DCKVC_GCTCATGA-CTCTCTAT_L001", 
    "409_069_DCKVC_ATCTCAGG-CTCTCTAT_L001"
])

# ANZCHOG abstract:
#SAMPLES = list([
#	"409_016_DBV4V_TAGGCATG-CTCTCTAT_L001", "409_018_DBV4V_AAGAGGCA-CTCTCTAT_L001", 
#	"409_031_DCB8V_GTAGAGGA-CTCTCTAT_L001", "409_014_DBV4V_TCCTGAGC-CTCTCTAT_L001",
#	"409_021_DBV4V_CTCTCTAC-CTCTCTAT_L001", "409_005_D9YWF_ATCTCAGG-CTCTCTAT_L001",
#	"409_019_DBV4V_GCTCATGA-CTCTCTAT_L001", "409_025_DCB8V_TCCTGAGC-CTCTCTAT_L001",
#	"409_007_DB62M_TAGGCATG-CTCTCTAT_L001", "409_006_DB62M_GGACTCCT-CTCTCTAT_L001",
#	"409_012_DBV4V_CGTACTAG-CTCTCTAT_L001", "409_013_DBV4V_AGGCAGAA-CTCTCTAT_L001",
#	"409_008_DB62M_GTAGAGGA-CTCTCTAT_L001","409_020_DBV4V_GTAGAGGA-CTCTCTAT_L001", 
#	"409_015_DBV4V_GGACTCCT-CTCTCTAT_L001",
#	"409_023_DCB8V_CGTACTAG-CTCTCTAT_L001", "409_027_DCB8V_TAGGCATG-CTCTCTAT_L001",
#	"409_062_DCCT9_CGAGGCTG-CTCTCTAT_L001", "409_060_DCCT9_TAAGGCGA-CTCTCTAT_L001",
#	"409_024_DCB8V_AGGCAGAA-CTCTCTAT_L001", "409_026_DCB8V_GGACTCCT-CTCTCTAT_L001",
#	"409_028_DCB8V_CTCTCTAC-CTCTCTAT_L001", "409_055_DCCT9_CGTACTAG-CTCTCTAT_L001",
#	"409_004_D9YWF_GTAGAGGA-CTCTCTAT_L001", "409_009_DB62M_GCTCATGA-CTCTCTAT_L001",
#	"409_011_DBV4V_TAAGGCGA-CTCTCTAT_L001", "409_059_DCCT9_CTCTCTAC-CTCTCTAT_L001",
#	"409_058_DCCT9_GGACTCCT-CTCTCTAT_L001", "409_061_DCCT9_AGGCAGAA-CTCTCTAT_L001",
#	"409_063_DCCT9_TCCTGAGC-CTCTCTAT_L001"
#])

## testing:
#SAMPLES = list([
#    "409_004_D9YWF_GTAGAGGA-CTCTCTAT_L001", "409_018_DBV4V_AAGAGGCA-CTCTCTAT_L001",
#    "409_021_DBV4V_CTCTCTAC-CTCTCTAT_L001"
#])

#TYPE = list([
#    "SV", "all"
#])

#rule all:
#    input:
#        expand(
#            gsnap_dir + '{sample}.SV.sorted.bam.bai',
#            sample=SAMPLES
#        ),
#        expand(
#            svaba_dir + '{sample}/{sample}.discordant.txt.gz',
#            sample=SAMPLES
#        ),
#        expand(
#            manta_dir + '{sample}/results/variants/somaticSV.vcf.gz',
#            sample=SAMPLES
#        )

rule all:
    input:
        expand(
            'logs/completed_jobs/{sample}_complete',
            sample = SAMPLES
        )
#        ,
#        expand(
#            geneglobe_svaba_dir + '{sample}/{sample}.svaba.unfiltered.sv.formatted.vcf.idx',
#            sample = SAMPLES
#        )


######################################################################################################
### 1. GSNAP ###
######################################################################################################

rule gsnap:
    input:
        fq1 = fq_dir + '{sample}_R1.fastq.gz',
        fq2 = fq_dir + '{sample}_R2.fastq.gz'
    output:
        gsnap_dir + '{sample}/{sample}.unpaired_uniq'
    threads: 7
    shell:
        'mkdir -p logs/gsnap; ' + 
        'cd logs/gsnap; ' + 
        env_dir + 'gsnap -D ' + genome_dir + ' -d GRCh37.p13.genome --gunzip -t {threads} -A sam ' +
        '--find-dna-chimeras 1 --split-output ' + project_dir + gsnap_dir + '{wildcards.sample}/{wildcards.sample} ' +
        '../../{input.fq1} ' + 
        '../../{input.fq2} 2> {wildcards.sample}.sort.errors'

rule collate_bams:
    input:
        gsnap_dir + '{sample}/{sample}.unpaired_uniq'
    output:
        gsnap_dir + '{sample}/{sample}.all'
    threads: 1
    shell:
        env_dir + 'samtools view -H {input} > {output}; ' +
        'for f in ' + gsnap_dir + '{wildcards.sample}/*; ' + 
            'do samtools view $f >> {output}; '
        'done;'

rule gsnap_bam:
    input:
        all_bam = gsnap_dir + '{sample}/{sample}.all',
        unpaired_bam = gsnap_dir + '{sample}/{sample}.unpaired_uniq'
    output:
        gsnap_dir + '{sample}/{sample}.all.bam'
    threads: 1
    shell:
        'for f in ' + gsnap_dir + '{wildcards.sample}/*; ' + 
            'do ' + env_dir + 'samtools view -bh $f > ' + '$f.bam; ' + 
        'done;'
        

rule gsnap_sort:
    input:
        gsnap_dir + '{sample}/{sample}.all.bam'
    output:
        gsnap_dir + '{sample}/{sample}.all.sorted.bam'
    threads: 1
    shell:
        'for f in ' + gsnap_dir + '{wildcards.sample}/*.bam; ' + 
            'do id=$(echo $f | sed "s/.bam//"); ' + 
            env_dir + 'samtools sort $f > ' + '$id.sorted.bam; ' + 
        'done'

rule gsnap_index:
    input:
        gsnap_dir + '{sample}/{sample}.all.sorted.bam'
    output:
        all_bai = gsnap_dir + '{sample}/{sample}.all.sorted.bam.bai',
        unpaired_bai = gsnap_dir + '{sample}/{sample}.unpaired_uniq.sorted.bam.bai'
    threads: 1
    shell:
        'for f in ' + gsnap_dir + '{wildcards.sample}/*.sorted.bam; ' + 
            'do ' + env_dir + 'samtools index $f; ' + 
        'done'

rule gsnap_svaba:
   input:
       all_bam = gsnap_dir + '{sample}/{sample}.all.sorted.bam',
       all_bai = gsnap_dir + '{sample}/{sample}.all.sorted.bam.bai'
   output:
       filt = gsnap_svaba_dir + '{sample}/{sample}.svaba.sv.vcf',
       unfilt = gsnap_svaba_dir + '{sample}/{sample}.svaba.unfiltered.sv.vcf'
   threads: 7
   shell:
       'mkdir -p logs/svaba; ' + 
        'cd logs/svaba; ' + 
        'mkdir -p ../../' + gsnap_svaba_dir + '{wildcards.sample}/; ' +
        'svaba run -t ' + project_dir + '{input.all_bam} -G ' + 
            genome_dir + 'GRCh37.p13.genome.fa -a ../../' + 
            gsnap_svaba_dir + 
            '{wildcards.sample}/{wildcards.sample} ' + 
            '-p 6 --override-reference-check' +
            ' 2> {wildcards.sample}.svaba.errors'

rule gsnap_format_vcf:
    input:
        gsnap_svaba_dir + '{sample}/{sample}.svaba.unfiltered.sv.vcf'
    output:
        gsnap_svaba_dir + '{sample}/{sample}.svaba.unfiltered.sv.formatted.vcf'
    threads: 1
    shell:
        'scripts/fix_broken_svaba_vcf.sh {wildcards.sample} {input}'

rule gsnap_vcf_index:
   input:
       filt = gsnap_svaba_dir + '{sample}/{sample}.svaba.sv.vcf',
       unfilt = gsnap_svaba_dir + '{sample}/{sample}.svaba.unfiltered.sv.formatted.vcf'
   output:
        filt = gsnap_svaba_dir + '{sample}/{sample}.svaba.sv.vcf.idx',
        unfilt = gsnap_svaba_dir + '{sample}/{sample}.svaba.unfiltered.sv.formatted.vcf.idx'
   threads: 2
   shell:
        env_dir + 'igvtools index {input.filt}; ' + 
        env_dir + 'igvtools index {input.unfilt}'
        

######################################################################################################
### 2. BWA ###
######################################################################################################

rule bwa:
    input:
        fq1 = fq_dir + '{sample}_R1.fastq.gz',
        fq2 = fq_dir + '{sample}_R2.fastq.gz'
    output:
        bwa_dir + '{sample}/{sample}.sam'
    threads: 7
    shell:
        'mkdir -p logs/bwa; ' + 
        'cd logs/bwa; ' + 
        env_dir + 'bwa mem -t 6 ' + genome_dir + 
            'GRCh37.p13.genome.fa ../../{input.fq1} ' +
            '../../{input.fq2} > ../../{output}; '

rule bwa_bam:
    input:
        bwa_dir + '{sample}/{sample}.sam'
    output:
        bwa_dir + '{sample}/{sample}.bam'
    threads: 1
    shell:
        env_dir + 'samtools view -bh {input} > {output}'

rule bwa_sort:
    input:
        bwa_dir + '{sample}/{sample}.bam'
    output:
        bwa_dir + '{sample}/{sample}.sorted.bam'
    threads: 1
    shell:
        env_dir + 'samtools sort {input} > {output}'

rule bwa_index:
    input:
        bwa_dir + '{sample}/{sample}.sorted.bam'
    output:
        bwa_dir + '{sample}/{sample}.sorted.bam.bai'
    threads: 1
    shell:
        env_dir + 'samtools index {input}'

rule bwa_svaba:
   input:
       bam = bwa_dir + '{sample}/{sample}.sorted.bam',
       bai = bwa_dir + '{sample}/{sample}.sorted.bam.bai'
   output:
       filt = bwa_svaba_dir + '{sample}/{sample}.svaba.sv.vcf',
       unfilt = bwa_svaba_dir + '{sample}/{sample}.svaba.unfiltered.sv.vcf'
   threads: 7
   shell:
       'mkdir -p logs/svaba; ' + 
        'cd logs/svaba; ' + 
        'mkdir -p ../../' + bwa_svaba_dir + '{wildcards.sample}/; ' +
        'svaba run -t ' + project_dir + '{input.bam} -G ' + 
            genome_dir + 'GRCh37.p13.genome.fa -a ../../' + 
            bwa_svaba_dir + 
            '{wildcards.sample}/{wildcards.sample} ' + 
            '-p 6 --override-reference-check' +
            ' 2> {wildcards.sample}.svaba.errors'

rule bwa_format_vcf:
    input:
        bwa_svaba_dir + '{sample}/{sample}.svaba.unfiltered.sv.vcf'
    output:
        bwa_svaba_dir + '{sample}/{sample}.svaba.unfiltered.sv.formatted.vcf'
    threads: 1
    shell:
        'scripts/fix_broken_svaba_vcf.sh {wildcards.sample} {input}'

rule bwa_vcf_index:
   input:
       filt = bwa_svaba_dir + '{sample}/{sample}.svaba.sv.vcf',
       unfilt = bwa_svaba_dir + '{sample}/{sample}.svaba.unfiltered.sv.formatted.vcf'
   output:
        filt = bwa_svaba_dir + '{sample}/{sample}.svaba.sv.vcf.idx',
        unfilt = bwa_svaba_dir + '{sample}/{sample}.svaba.unfiltered.sv.formatted.vcf.idx'
   threads: 2
   shell:
        env_dir + 'igvtools index {input.filt}; ' + 
        env_dir + 'igvtools index {input.unfilt}'
        

#######################################################################################################
#### 4. GeneGlobe ###
#######################################################################################################
#
#rule geneglobe_svaba:
#   input:
#       bam = geneglobe_dir + '{sample}/{sample}.sorted.bam',
#       bai = geneglobe_dir + '{sample}/{sample}.sorted.bam.bai'
#   output:
#       filt = geneglobe_svaba_dir + '{sample}/{sample}.svaba.sv.vcf',
#       unfilt = geneglobe_svaba_dir + '{sample}/{sample}.svaba.unfiltered.sv.vcf'
#   threads: 7
#   shell:
#       'mkdir -p logs/svaba; ' + 
#        'cd logs/svaba; ' + 
#        'mkdir -p ../../' + geneglobe_svaba_dir + '{wildcards.sample}/; ' +
#        'svaba run -t ' + project_dir + '{input.bam} -G ' + 
#            genome_dir + 'GRCh37.p13.genome.fa -a ../../' + 
#            geneglobe_svaba_dir + 
#            '{wildcards.sample}/{wildcards.sample} ' + 
#            '-p 6 --override-reference-check' +
#            ' 2> {wildcards.sample}.svaba.errors'
#
#rule geneglobe_format_vcf:
#    input:
#        geneglobe_svaba_dir + '{sample}/{sample}.svaba.unfiltered.sv.vcf'
#    output:
#        geneglobe_svaba_dir + '{sample}/{sample}.svaba.unfiltered.sv.formatted.vcf'
#    threads: 1
#    shell:
#        'scripts/fix_broken_svaba_vcf.sh {wildcards.sample} {input}'
#
#rule geneglobe_vcf_index:
#   input:
#       filt = geneglobe_svaba_dir + '{sample}/{sample}.svaba.sv.vcf',
#       unfilt = geneglobe_svaba_dir + '{sample}/{sample}.svaba.unfiltered.sv.formatted.vcf'
#   output:
#        filt = geneglobe_svaba_dir + '{sample}/{sample}.svaba.sv.vcf.idx',
#        unfilt = geneglobe_svaba_dir + '{sample}/{sample}.svaba.unfiltered.sv.formatted.vcf.idx'
#   threads: 2
#   shell:
#        env_dir + 'igvtools index {input.filt}; ' + 
#        env_dir + 'igvtools index {input.unfilt}'


######################################################################################################
### 5. Clean up ###
######################################################################################################

rule cleanup:
    input:
        gsnap_svaba_ind = gsnap_svaba_dir + '{sample}/{sample}.svaba.sv.vcf.idx',
        gsnap_unfilt_svaba_ind = gsnap_svaba_dir + '{sample}/{sample}.svaba.unfiltered.sv.formatted.vcf.idx',
        gsnap_unpaired_bai = gsnap_dir + '{sample}/{sample}.unpaired_uniq.sorted.bam.bai',
        bwa_svaba_ind = bwa_svaba_dir + '{sample}/{sample}.svaba.sv.vcf.idx',
        bwa_unfilt_svaba_ind = bwa_svaba_dir + '{sample}/{sample}.svaba.unfiltered.sv.formatted.vcf.idx'
    output:
        'logs/completed_jobs/{sample}_complete'
    threads: 1
    shell:
        'mkdir -p ' + gsnap_dir + '{wildcards.sample}/temp; ' +
        'mv ' + gsnap_dir + '{wildcards.sample}/*sorted* ' + gsnap_dir + '{wildcards.sample}/temp; ' +
        'rm -f ' + gsnap_dir + '{wildcards.sample}/*bam*; ' + 
        'mv ' + gsnap_dir + '{wildcards.sample}/temp/* ' + gsnap_dir + '{wildcards.sample}; ' +
        'mkdir -p ' + bwa_dir + '{wildcards.sample}/temp; ' +
        'mv ' + bwa_dir + '{wildcards.sample}/*sorted* ' + bwa_dir + '{wildcards.sample}/temp; ' +
        'rm -f ' + bwa_dir + '{wildcards.sample}/*bam*; ' + 
        'mv ' + bwa_dir + '{wildcards.sample}/temp/* ' + bwa_dir + '{wildcards.sample}; ' +
        'mkdir -p logs/completed_jobs; '
        'touch {output}'


