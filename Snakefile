# Run command:
# snakemake --reason --use-conda --cores 100 --cluster-config cluster.json --cluster 'qsub -pe smp {cluster.cores} -N ewfus.smk -wd '/share/ScratchGeneral/jamtor/projects/ewing_ctDNA/logs' -b y -j y -V -P TumourProgression' -j 23

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
env_dir = conda_dir + "envs/py37/envs/snk/bin/"

fq_dir = 'raw_files/'
gsnap_dir = 'results/gsnap/'
bwa_dir = 'results/bwa/'
bwa_svaba_dir = 'results/svaba/bwa/'
gsnap_svaba_dir = 'results/svaba/gsnap/'
manta_dir = 'results/manta/'
manta_bin = '/g/data1a/ku3/jt3341/local/lib/manta-1.5.0/bin/'

SAMPLES = list([
    "409_002_D9YW9_GGACTCCT-CTCTCTAT_L001"
#    "409_001_D9YW9_TCCTGAGC-CTCTCTAT_L001", "409_002_D9YW9_GGACTCCT-CTCTCTAT_L001", 
#    "409_003_D9YWF_AGGCAGAA-CTCTCTAT_L001", "409_004_D9YWF_GTAGAGGA-CTCTCTAT_L001", 
#    "409_005_D9YWF_ATCTCAGG-CTCTCTAT_L001", "409_006_DB62M_GGACTCCT-CTCTCTAT_L001", "409_007_DB62M_TAGGCATG-CTCTCTAT_L001", 
#    "409_008_DB62M_GTAGAGGA-CTCTCTAT_L001", "409_009_DB62M_GCTCATGA-CTCTCTAT_L001", "409_010_DB62M_ATCTCAGG-CTCTCTAT_L001", 
#    "409_011_DBV4V_TAAGGCGA-CTCTCTAT_L001", "409_012_DBV4V_CGTACTAG-CTCTCTAT_L001", "409_013_DBV4V_AGGCAGAA-CTCTCTAT_L001", 
#    "409_014_DBV4V_TCCTGAGC-CTCTCTAT_L001", "409_015_DBV4V_GGACTCCT-CTCTCTAT_L001", "409_016_DBV4V_TAGGCATG-CTCTCTAT_L001", 
#    "409_017_DBV4V_CGAGGCTG-CTCTCTAT_L001", "409_018_DBV4V_AAGAGGCA-CTCTCTAT_L001", "409_019_DBV4V_GCTCATGA-CTCTCTAT_L001",  
#    "409_020_DBV4V_GTAGAGGA-CTCTCTAT_L001", "409_021_DBV4V_CTCTCTAC-CTCTCTAT_L001"
])

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
#        expand(
#            gsnap_dir + '{sample}/{sample}.unpaired_uniq.sorted.bam.bai',
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
    threads: 6
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

rule bam:
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
        

rule sort:
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

rule index:
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

rule svaba:
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

rule format_vcf:
    input:
        gsnap_svaba_dir + '{sample}/{sample}.svaba.unfiltered.sv.vcf'
    output:
        gsnap_svaba_dir + '{sample}/{sample}.svaba.unfiltered.sv.formatted.vcf'
    threads: 1
    shell:
        'scripts/fix_broken_svaba_vcf.sh {wildcards.sample} {input}'

rule vcf_index:
   input:
       filt = gsnap_svaba_dir + '{sample}/{sample}.svaba.sv.vcf',
       unfilt = gsnap_svaba_dir + '{sample}/{sample}.svaba.unfiltered.sv.formatted.vcf'
   output:
        filt = gsnap_svaba_dir + '{sample}/{sample}.svaba.sv.vcf.idx',
        unfilt = gsnap_svaba_dir + '{sample}/{sample}.svaba.unfiltered.sv.formatted.vcf.idx'
   threads: 1
   shell:
        env_dir + 'igvtools index {input.filt}; ' + 
        env_dir + 'igvtools index {input.unfilt}'

rule cleanup:
    input:
        svaba_ind = gsnap_svaba_dir + '{sample}/{sample}.svaba.sv.vcf.idx',
        unfilt_svaba_ind = gsnap_svaba_dir + '{sample}/{sample}.svaba.unfiltered.sv.formatted.vcf.idx',
        unpaired_bai = gsnap_dir + '{sample}/{sample}.unpaired_uniq.sorted.bam.bai'
    output:
        'logs/completed_jobs/{sample}_complete'
    threads: 1
    shell:
        'mkdir -p ' + gsnap_dir + '{wildcards.sample}/temp; ' +
        'mv ' + gsnap_dir + '{wildcards.sample}/*sorted* ' + gsnap_dir + '{wildcards.sample}/temp; ' +
        'rm -f ' + gsnap_dir + '{wildcards.sample}/*bam*; ' + 
        'mv ' + gsnap_dir + '{wildcards.sample}/temp/* ' + gsnap_dir + '{wildcards.sample}; ' +
        'mkdir -p logs/completed_jobs; '
        'touch {output}'
        

