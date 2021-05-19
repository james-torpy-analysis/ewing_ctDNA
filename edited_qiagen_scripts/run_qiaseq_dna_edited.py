import ConfigParser
import sys
import multiprocessing
import shutil
# our modules
import core.run_log
import core.run_config
import core.prep
import core.align
import core.umi_filter
import core.umi_mark
import core.umi_merge
import core.primer_clip
import core.samtools
import core.tumor_normal
import core.sm_counter_wrapper
import metrics.sum_specificity
import metrics.sum_uniformity_primer
import metrics.sum_primer_umis
import metrics.sum_all
import metrics.umi_frags
import metrics.umi_depths
import metrics.duplex_summary
import metrics.sum_primer_duplex
import metrics.fraglen_by_rpu
import misc.process_ion
import misc.tvc
import annotate.vcf_complex
import annotate.vcf_annotate
import timeit
sys.argv = [
    'run_qiaseq_dna', 'run_sm_counter_v2.params.txt', 'v2', 'single',\
    '409_004_D9YWF_GTAGAGGA-CTCTCTAT_L001'
]
paramFile = sys.argv[1]
vc = sys.argv[2]
analysis = sys.argv[3]
readSet   = " ".join(sys.argv[4:])
tumorNormal = False
# read run configuration file to memory
cfg = core.run_config.run(readSet,paramFile)
# trim adapters , umi and primers (this module spawns multiple processes)
core.prep.run(cfg)
readFileIn1 = readSet + ".prep.R1.fastq"
readFileIn2 = readSet + ".prep.R2.fastq"
bamFileOut  = readSet + ".align.bam"

if cfg.platform.lower() == "illumina":
    # align trimmed reads to genome using BWA MEM
    core.align.run(cfg, readFileIn1, readFileIn2, bamFileOut)
else: # use tmap for ion torrent reads
    misc.process_ion.alignToGenomeIon(cfg, readFileIn1, bamFileOut)

# call putative unique input molecules using BOTH UMI seq AND genome alignment position on random fragmentation side
bamFileIn  = readSet + ".align.bam"
# the following steps remove split reads:
core.umi_filter.run(cfg, bamFileIn)
core.umi_mark.run(cfg)
metrics.umi_frags.run(cfg)
metrics.umi_depths.run(cfg,vc)
core.umi_merge.run(cfg, bamFileIn)
# soft clip primer regions from read alignments
bamFileIn  = readSet + ".umi_merge.bam"
bamFileOut = readSet + ".primer_clip.bam"
core.primer_clip.run(cfg, bamFileIn, bamFileOut, False)
# additional metrics to generate
metrics.sum_primer_umis.run(cfg) # primer-level umi and read metrics
metrics.sum_specificity.run(cfg) # priming specificity
metrics.sum_uniformity_primer.run(cfg) # primer-level uniformity
if cfg.duplex: # additional metrics for Duplex reads
    metrics.duplex_summary.run(cfg)
    metrics.sum_primer_duplex.run(cfg)
    metrics.fraglen_by_rpu.run(cfg)

# sort the final BAM file, to prepare for downstream variant calling
bamFileIn  = readSet + ".primer_clip.bam"
bamFileOut = readSet + ".bam"
core.samtools.sort(cfg, bamFileIn, bamFileOut)
if cfg.platform.lower() != "illumina": # ion reads
    misc.tvc.run(cfg)

# Run smCounter variant calling
numVariants = core.sm_counter_wrapper.run(cfg, paramFile, vc)
if cfg.platform.lower() != "illumina":
    numVariants = misc.tvc.smCounterFilter(cfg,vc)
    
# create complex variants, and annotate using snpEff
if not tumorNormal:
    post_smcounter_work(numVariants, readSet, cfg, tumorNormal = False)
    # close log file
    core.run_log.close()


    