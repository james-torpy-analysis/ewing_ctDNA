import ConfigParser
import os
import sys
# get read set name
readSet = cfg.readSet
# get standard smCounter parameters from the main run-params.txt file
parser = ConfigParser.SafeConfigParser()
parser.optionxform = str
parser.read(paramFile)
cfgSmCounter = {}
for (paramName, paramVal) in parser.items("smCounter"): 
    cfgSmCounter[paramName] = paramVal

# set up config dictionary to pass to smCounter
cfgSmCounter["outPrefix"] = readSet
cfgSmCounter["bamFile"  ] = readSet + ".bam"
if cfg.platform.lower() != "illumina": # ion reads
    cfgSmCounter["bedTarget"] = readSet + ".tvc_roi.bed"  # subset to tvc variants
else:
    cfgSmCounter["bedTarget"] = cfg.roiBedFile

cfgSmCounter["nCPU"     ] = cfg.numCores
cfgSmCounter["refGenome"] = cfg.genomeFile
cfgSmCounter["isDuplex"]  = cfg.duplex

if cfg.duplex:
    cfgSmCounter["duplexTag"] = cfg.tagNameDuplex

cfgSmCounter["rpu"      ] = cfg.readsPerUmi  # this comes from metrics.umi_frags module
cfgSmCounter["runPath"] = os.getcwd()

import dill
filename = 'before_running_smcounter2.pkl'
dill.dump_session(filename)
# and to load the session again:
#dill.load_session(filename)


sm_counter_v2.run.main(cfgSmCounter)
smCounterThreshold = 6
# need to add the lod quantiles output from smCounter-v2 to umi_depths.summary file
fileoutSummary = open(readSet + ".umi_depths.summary.txt","a")
fileLodQuantiles = readSet + ".umi_depths.variant-calling-lod.bedgraph.quantiles.txt"
if os.path.exists(fileLodQuantiles):
    with open(fileLodQuantiles,"r") as IN:
        for line in IN:
            (metricName, metricVal) = line.strip().split("|")
            metricName = int(metricName.replace("%",""))
            metricVal = float(metricVal)
            thorst = "st" if metricName == 1 else "th"
            fileoutSummary.write("{:6.4f}\t{:2d}{} percentile estimated minimum detectible allele fraction (LOD)\n".format(metricVal, metricName,thorst))
    # remove the temporary file
    os.remove(fileLodQuantiles)
# write smCounter threshold to disk file, for main summary table
fileout = open(readSet + ".smCounter.summary.txt", "w")
fileout.write("{}\tsmCounter variant calling threshold\n".format(smCounterThreshold))
fileout.close()
# return number of primitive variants called
numVariants = -1
cutFile = readSet + ".smCounter.cut.txt"
if os.path.exists(cutFile): # file could be absent - empty bed or bam , or hack for e.g. Duplex
    for line in open(cutFile):
        numVariants += 1
return numVariants
