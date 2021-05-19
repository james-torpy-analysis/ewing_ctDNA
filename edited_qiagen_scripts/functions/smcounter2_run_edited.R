#!/usr/bin/python
import os
import subprocess
import multiprocessing
import datetime
import argparse
import functools
import sys
from collections import defaultdict
import pysam
import datetime
import operator
from collections import defaultdict
import random
import numpy
import string
import logging
import traceback
import scipy.stats
import dill
# define functions:
def argParseInit():  
   global parser
   parser = argparse.ArgumentParser(description = 'smCounter2: variant calling using Unique Molecular Identifiers')
   parser.add_argument('--runPath', default = None, help = 'path to working directory')
   parser.add_argument('--bedTarget', default = None, help = 'BED file')
   parser.add_argument('--bamFile', default = None, help = 'BAM file')
   parser.add_argument('--outPrefix', default = None, help = 'file name prefix')
   parser.add_argument('--nCPU', type = int, default = 1, help = 'number of CPU to use in parallel')
   parser.add_argument('--minBQ', type = int, default = 25, help = 'minimum base quality allowed for analysis')
   parser.add_argument('--minMQ', type = int, default = 50, help = "minimum mapping quality allowed for analysis. If the bam is tagged with its mate's mapq, then the minimum of the R1 and R2 mapq will be used for comparison, if not each read is compared independently.")
   parser.add_argument('--hpLen', type = int, default = 10, help = 'minimum length for homopolymers')
   parser.add_argument('--mismatchThr', type = float, default = 6.0, help = 'average number of mismatches per 100 bases allowed')
   parser.add_argument('--primerDist', type = int, default = 2, help = 'filter variants that are within X bases to primer')
   parser.add_argument('--consThr', type = float, default = 0.8, help = 'threshold on read proportion to determine UMI level consensus')
   parser.add_argument('--rpu', type = float, default = 0.0, help = 'mean read pairs per UMI; default at 0 and will be calculated')
   parser.add_argument('--isRna', action = 'store_true', help = 'RNAseq varinat calling only; default is DNAseq')
   parser.add_argument('--primerSide', type = int, default = 1, help = 'read end that includes the primer; default is 1')
   parser.add_argument('--umiTag', type = str, default = 'Mi', help = 'tag name for normal UMI')
   parser.add_argument('--primerTag', type = str, default = 'pr', help = 'tag name for Primer')
   parser.add_argument('--mqTag', type = str, default = 'MQ', help = 'tag name for MapQ score of mate')
   parser.add_argument('--tagSeparator', type = str, default = '-', help = 'tag seperator for splitting UMI tag') 
   parser.add_argument('--minAltUMI', type = int, default = 3, help = 'minimum requirement of ALT UMIs; default is 3')
   parser.add_argument('--maxAltAllele', type = int, default = 2, help = 'maximum ALT alleles that meet minAltUMI to be reported; default is 2 (tri-allelic variants)')
   parser.add_argument('--refGenome',type = str,help = 'path to the reference fasta file')
   parser.add_argument('--repBed',type = str,help = 'path to the simpleRepeat bed file')
   parser.add_argument('--srBed',type = str,help = 'path to the full repeat bed file')
   parser.add_argument('--ds', type = int, default = 10000, help = 'down sample if number of UMIs greater than this value (RNA only)')
   parser.add_argument('--bamType', type = str, default = 'raw', help = 'raw (default): raw BAM file with UMIs; consensus: consensused BAM file')
   parser.add_argument('--inputVCF', type = str, default = None, help = 'optional input VCF file') 
   parser.add_argument('--sinBkgErrorDist', type = str, default = '/srv/qgen/data/annotation/bkg.error.v2.RData', help = 'background error rate distribution for normal DNA-seq runs') 
   parser.add_argument('--isDuplex', action = 'store_true', help = 'duplex-seq DNA varinat calling only; default is normal DNAseq')
   parser.add_argument('--duplexTag', type = str, default = None, help = 'tag name for duplex UMI')
   parser.add_argument('--minRpu', type = int, default = 2, help = 'minimum read pairs for UMI to be included; default is 1 for normal DNA-seq (before gradually dropping singletons) and 2 for duplex-seq')
   parser.add_argument('--dupBkgErrorDist', type = str, default ='/srv/qgen/data/annotation/duplex.bkg.error.ga.RData', help = 'background G>A error rate distribution for duplex-seq runs')  
#-------------------------------------------------------------------------------------
# find homopolymer sequences
#-------------------------------------------------------------------------------------
def findHp(bedName, outName, minLength, refg, isRna):
   # how much to extend the roi to search for homopolymers
   extensionLen = 0 if isRna else 100
   # loop over roi BED
   outfile = open(outName, 'w')
   for line in open(bedName, 'r'):
      if line.startswith('track name='):
         continue
      lineList = line.strip().split('\t')
      chrom = lineList[0]
      start = int(lineList[1])
      end = int(lineList[2])
      # get reference base
      refSeq = pysam.FastaFile(refg)
      start_coord = start - 1 - extensionLen
      if start_coord < 0:
         start_coord = start
      origRef = refSeq.fetch(reference=chrom, start=start_coord, end=end + extensionLen)
      origRef = origRef.upper()
      hpL = 0
      for i in range(len(origRef))[1:]:
         if origRef[i] == origRef[i-1]:
            continue
         else:
            hpLen = i - hpL 
            realL = hpL - 1 + start - extensionLen
            realR = i - 1  + start - extensionLen
            if hpLen >= minLength and realL <= end and realR >= start:
               outline = '\t'.join([chrom, str(max(realL, start)), str(min(realR, end)), 'HP', str(hpLen), str(realL), str(realR), origRef[hpL]]) + '\n'
               outfile.write(outline)
            hpL = i 
   outfile.close()
#----------------------------------------------------------------------------------------------
# get homopolymer region information
#------------------------------------------------------------------------------------------------
def getHpInfo(bedTarget, refGenome, isRna, hpLen):
   # intersect repeats and target regions
   findHpLen = hpLen if isRna else 6
   findHp(bedTarget, 'hp.roi.bed', findHpLen, refGenome, isRna)
   # gather homopolymer region info
   hpRegion = defaultdict(list)
   with open('hp.roi.bed','r') as IN:
      for line in IN:   
         chrom, regionStart, regionEnd, repType, totalLen, realL, realR, repBase = line.strip().split()
         hpRegion[chrom].append([regionStart, regionEnd, repType, totalLen, realL, realR])
   return(hpRegion)
#----------------------------------------------------------------------------------------------
# get tandem region information
#------------------------------------------------------------------------------------------------
def getTrInfo(bedTarget, repBed, isRna, hpLen):
   # intersect repeats and target regions
   subprocess.check_call('bedtools intersect -a ' + repBed + ' -b ' + bedTarget + ' | bedtools sort -i > rep.roi.bed', shell = True)
   # gather tandem repeat region info
   repRegion = defaultdict(list)
   with open('rep.roi.bed','r') as IN:
      for line in IN:
         chrom, regionStart, regionEnd, repInfo = line.strip().split()[:4]
         unitLen, repLen = repInfo.split("|")[1:3]
         try:
            unitLen_num = float(unitLen)
         except ValueError:
            continue
         try:
            repLen_num = float(repLen)
         except ValueError:
            continue
         if isRna:
            totalLen = int(regionEnd) - int(regionStart)
            if totalLen < hpLen:
               continue
            repLen = str(totalLen / unitLen_num) if unitLen_num > 0 else '0'
            totalLen = str(totalLen)
         else:
            totalLen = str(unitLen_num * repLen_num)
         repBase = repInfo[-1]
         repType = 'RepT'
         repRegion[chrom].append([regionStart, regionEnd, repType, totalLen, unitLen, repLen])
   return(repRegion)
#----------------------------------------------------------------------------------------------
# get other repeats region (simple repeats, low complexity, micro-satelites) information
#------------------------------------------------------------------------------------------------
def getOtherRepInfo(bedTarget, srBed, isRna, hpLen):
   # intersect repeats and target regions
   subprocess.check_call('bedtools intersect -a ' + srBed +  ' -b ' + bedTarget + ' | bedtools sort -i > sr.roi.bed', shell=True)   
   # gather other repeat region info
   srRegion = defaultdict(list)
   with open('sr.roi.bed','r') as IN:
      for line in IN:
         chrom, regionStart, regionEnd, repInfo = line.strip().split()
         repType, totalLen, unitLen, repLen, repBase = repInfo.strip().split("|")
         if repType == 'Simple_repeat':
            repType = 'RepS'
         elif repType == 'Low_complexity':
            repType = 'LowC'
         elif repType == 'Satellite':
            repType = 'SL'
         else:
            repType = 'Other_Repeat'
         if isRna:
            totalLen = int(regionEnd) - int(regionStart)
            if totalLen < hpLen:
               continue
            try:
               unitLen_num = float(unitLen)
               repLen = str(totalLen / unitLen_num) if unitLen_num > 0 else '0'
            except ValueError:
               pass
            totalLen = str(totalLen)
         srRegion[chrom].append([regionStart, regionEnd, repType, totalLen, unitLen, repLen])
   return(srRegion)
#----------------------------------------------------------------------------------------------
# generate locList, where each member is a target site
#------------------------------------------------------------------------------------------------
def getLocList(bedTarget, hpRegion, repRegion, srRegion, isDuplex):
   max_bases_for_interval = 175 if isDuplex else 250
   locList = []
   with open(bedTarget,'r') as IN:
      for line in IN:
         if line.startswith('track name='):
            continue
         lineList = line.strip().split('\t')
         chrom = lineList[0]
         regionStart = int(lineList[1]) + 1   # target region starts from 1-base after 
         regionEnd = lineList[2]
         interval = [] # information for an interval
         nBases = 0 # no. of bases in an interval
         pos = regionStart
         lineEnd = False
         while not lineEnd:
            (hpInfo, srInfo, repInfo) = ('.', '.', '.')
            repTypeSet = set()
            # check if the site is in homopolymer region (not including 1 base before) 
            for (regionStart_tmp, regionEnd_tmp, repType_tmp, totalLen_tmp, realL, realR) in hpRegion[chrom]:
               if pos >= int(regionStart_tmp) - 0 and pos <= int(regionEnd_tmp):
                  repTypeSet.add(repType_tmp)
                  hpInfo = ';'.join([chrom, regionStart_tmp, regionEnd_tmp, totalLen_tmp, realL, realR])
                  break
            # check if the site is in other repeats region (including 1 base before) 
            for (regionStart_tmp, regionEnd_tmp, repType_tmp, totalLen_tmp, unitLen_tmp, repLen_tmp) in srRegion[chrom]:
               if pos >= int(regionStart_tmp) - 1 and pos <= int(regionEnd_tmp):
                  repTypeSet.add(repType_tmp)
                  srInfo = ';'.join([chrom, regionStart_tmp, regionEnd_tmp, totalLen_tmp, unitLen_tmp, repLen_tmp])
                  break
            for [regionStart_tmp, regionEnd_tmp, repType_tmp, totalLen_tmp, unitLen_tmp, repLen_tmp] in repRegion[chrom]:
               if pos >= int(regionStart_tmp) - 1 and pos <= int(regionEnd_tmp):
                  repTypeSet.add(repType_tmp)
                  repInfo = ';'.join([chrom, regionStart_tmp, regionEnd_tmp, totalLen_tmp, unitLen_tmp, repLen_tmp])
                  break
            repType = 'NA' if len(repTypeSet) == 0 else ';'.join(list(repTypeSet))
            interval.append((chrom, str(pos), repType, hpInfo, srInfo, repInfo))
            if nBases == max_bases_for_interval: # restrict interval size
               locList.append(interval)
               interval = []
               nBases = 0
            if str(pos) == regionEnd:
               lineEnd = True
            else:
               nBases += 1
               pos += 1
         if len(interval) > 0:
            locList.append(interval)
   return(locList)
#------------------------------------------------------------------------------------------------
# wrapper function for "vc()" - because Python multiprocessing module does not pass stack trace
#------------------------------------------------------------------------------------------------
def vc_wrapper(general_args, interval):
   timeStart = datetime.datetime.now()
   try:
      output = []
      hqCache = {}
      infoCache = {}
      bamName, minBq, minMq, hpLen, mismatchThr, primerDist, consThr, rpu, primerSide, refg, minAltUmi, maxAltAllele, isRna, ds, bamType, umiTag, primerTag, mqTag, tagSeparator, isDuplex, duplexTag, minRpu = general_args
      chrom = interval[0][0]
      intervalStartPos = interval[0][1]
      intervalEndPos = interval[-1][1]
      # constants and functions to use depending on normal DNA-seq or duplex-seq
      if isDuplex:
         nCols = nColsDup
         defVarFun = dup_defVar
         outBkgFun = dup_outBkg
         outLongFun = dup_outLong
         getUmiFun = dup_getUmi
         dropSingletonFun = dup_dropSingleton
         groupByUmiFun = dup_groupByUmi
         if isRna or bamType != 'raw':
            exit("duplex-seq must be on DNA and with raw BAM input for now.")
      else:
         nCols = nColsSin
         defVarFun = defVar
         outBkgFun = outBkg
         outLongFun = outLong  
         getUmiFun = getUmi
         dropSingletonFun = dropSingleton
         groupByUmiFun = groupByUmi
         if duplexTag != None:
            exit("normal DNA-seq runs don't have duplex tag.")
      refseq = pysam.FastaFile(refg)
      chromLengths = {}
      for idx in range(len(refseq.lengths)):
         chromLengths[refseq.references[idx]] = refseq.lengths[idx]
      i = 0
      for read_pileup in pileup(bamName, chrom, intervalStartPos, intervalEndPos):
         site = interval[i]
         i += 1
         chrom, pos, repType, hpInfo, srInfo, repInfo = site
         if read_pileup is None: # site is not covered at all, pysam simply skips such sites
            origRef = getRef(refseq, chrom, pos)
            outLineLong = '\t'.join([chrom, pos, origRef] + ["0"] * (nCols - 4) + ["LM"]) + "\n"
            out = [outLineLong,""]
            output.append(out)
            continue
         temp = [bamName, chrom, pos, repType, hpInfo, srInfo, repInfo, minBq, minMq, hpLen, mismatchThr, primerDist, consThr, rpu, primerSide, refseq, minAltUmi, maxAltAllele, isRna, ds, bamType, read_pileup, hqCache, infoCache, chromLengths[chrom], umiTag, primerTag, mqTag, tagSeparator, nCols, defVarFun, getUmiFun, groupByUmiFun, dropSingletonFun, outBkgFun, outLongFun, isDuplex, duplexTag, minRpu]
         outLineLong, outLineBkg, hqCache, infoCache = vc(*temp)
         out = [outLineLong, outLineBkg]
         output.append(out)
   except Exception as e:
      out = ("Exception thrown!\n" + traceback.format_exc(), "no_bg")
      output.append(out)
      logger.info("Exception thrown in vc() function at genome location : {pos} in interval : {chrom}:{it1}-{it2}".format(pos = pos, chrom = chrom, it1 = intervalStartPos, it2 = intervalEndPos))
      logger.info(out[0])
      raise Exception(e)   
   refseq.close()
   timeEnd = datetime.datetime.now()   
   logger.info(str(timeEnd - timeStart) + "\t" + "{chrom}:{it1}-{it2}".format(chrom = chrom,it1 = intervalStartPos,it2 = intervalEndPos))
   return output

# define variables:
args = cfgSmCounter

# save session:
filename = 'after_smcounter2_function_import.pkl'
dill.dump_session(filename)
# and to load the session again:
#dill.load_session(filename)

### code begins ###
# global constants
#codePath = os.path.dirname(os.path.abspath(__file__))
pValCode_sin = os.path.join(codePath,'get_pvalue.R')
pValCode_dup = os.path.join(codePath,'get_lr.duplex.R')
locChunkLen = 1000
seed = 10262016
nsim = 5000
parser = None
header_sin = ['CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'sUMT', 'sForUMT', 'sRevUMT', 'sVMT', 'sForVMT', 'sRevVMT', 'sVMF', 'sForVMF', 'sRevVMF', 'VDP', 'VAF', 'RefForPrimer', 'RefRevPrimer', 'primerOR', 'pLowQ', 'hqUmiEff', 'allUmiEff', 'refMeanRpb', 'altMeanRpb', 'rpbEffectSize', 'repType', 'hpInfo', 'simpleRepeatInfo', 'tandemRepeatInfo', 'DP', 'FR', 'MT', 'UFR', 'sUMT_A', 'sUMT_T', 'sUMT_G', 'sUMT_C', 'FILTER']
header_dup = ['CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'sUMT', 'sVMT', 'sVMF', 'dUMT', 'dVMT', 'dVMF', 'DP', 'VDP', 'VAF', 'RefForPrimer', 'RefRevPrimer', 'primerOR', 'pLowQ', 'hqUmiEff', 'allUmiEff', 'refMeanRpb', 'altMeanRpb', 'rpbEffectSize', 'repType', 'hpInfo', 'simpleRepeatInfo', 'tandemRepeatInfo', 'FR', 'MT', 'sForUMT', 'sRevUMT', 'sForVMT', 'sRevVMT', 'UFR', 'sUMT_A', 'sUMT_T', 'sUMT_G', 'sUMT_C', 'dUMT_A', 'dUMT_T', 'dUMT_G', 'dUMT_C', 'FILTER']
header_bkg = ['CHROM', 'POS', 'REF', 'A/G', 'G/A', 'C/T', 'T/C', 'A/C', 'C/A', 'A/T', 'T/A', 'C/G', 'G/C', 'G/T', 'T/G', 'negStrand', 'posStrand', 'AllSMT', 'workflow' ]

# log run start
timeStart = datetime.datetime.now()
print("started " + str(timeStart))

# if argument parser global not assigned yet, initialize it
if parser == None:
   argParseInit()

# get arguments passed in via a lambda object (e.g. from upstream pipeline)
if type(args) is not argparse.Namespace:
   argsList = []
   for argName, argVal in args.iteritems():
      if argName == "isDuplex":
         if argVal: # --isDuplex triggered upstream
            argsList.append("--{0}".format(argName))
      else:
         argsList.append("--{0}={1}".format(argName, argVal))
   args = parser.parse_args(argsList)

for argName, argVal in vars(args).iteritems():
   print(argName, argVal)

# change working directory to runDir and make output directories
if args.runPath != None:
   os.chdir(args.runPath)

# make /intermediate directory to keep the long output
if not os.path.exists('intermediate'):
   os.makedirs('intermediate')

# convert VCF to BED if inputVCF is not 'none'
bedTarget = args.bedTarget if args.inputVCF is None else vcf2bed(args.inputVCF)
# gather repetitive regions information
hpRegion = getHpInfo(bedTarget, args.refGenome, args.isRna, args.hpLen)
repRegion = getTrInfo(bedTarget, args.repBed, args.isRna, args.hpLen)
srRegion = getOtherRepInfo(bedTarget, args.srBed, args.isRna, args.hpLen)
# read in bed file and create a list of positions, annotated with repetitive region
locList = getLocList(bedTarget, hpRegion, repRegion, srRegion, args.isDuplex)

# calculate rpu if args.rpu = 0
if args.rpu == 0.0:
   if args.bamType == 'raw':
      rpu = getMeanRpu(args.bamFile, args.umiTag)
      print("rpu = " + str(round(rpu,1)) + ", computed by smCounter2")
   else:
      rpu = 5.0
      print("rpu = " + str(round(rpu,1)) + ", set by smCounter2 when bamType is consensus and rpu is not given by user")
else:
   rpu = args.rpu
   print("rpu = " + str(round(rpu,1)) + ", given by user")
   
# set primer side
primerSide = 'R1' if args.primerSide == 1 else 'R2'
# set type of input BAM file
bamType = 'raw' if args.bamType == 'raw' else 'consensus'
# select header for normal and duplex-seq runs
header = header_dup if args.isDuplex else header_sin
#----- loop over locs
# prepare to save to disk
outfile_long = open('intermediate/noThres.' + args.outPrefix + '.VariantList.long.txt', 'w')
bkgFileName = 'intermediate/bkg.' + args.outPrefix + '.txt'
outfile_bkg = open(bkgFileName, 'w')
outfile_long.write('\t'.join(header) + '\n')
outfile_bkg.write('\t'.join(header_bkg) + '\n')   
print('runtime' + '\t' + 'interval')
pool = multiprocessing.Pool(args.nCPU)

#import dill
#filename = 'before_processing_bedfile.pkl'
#dill.dump_session(filename)
## and to load the session again:
#dill.load_session(filename)
# process exons/intervals from bed file in parallel

func = functools.partial(
    vc_wrapper, 
    (args.bamFile, args.minBQ, args.minMQ, args.hpLen, args.mismatchThr, 
    args.primerDist, args.consThr, rpu, primerSide, args.refGenome, args.minAltUMI, 
    args.maxAltAllele, args.isRna, args.ds, bamType, args.umiTag, args.primerTag, 
    args.mqTag, args.tagSeparator, args.isDuplex, args.duplexTag, args.minRpu)
)

empty = True
for interval_result in pool.map(func, locList):
   for base_result in interval_result:
      vcOutline,bkgOutline = base_result
      outfile_long.write(vcOutline)
      outfile_bkg.write(bkgOutline)
      empty = False
# clear finished pool
pool.close()
pool.join()
# close output file handles
outfile_long.close()
outfile_bkg.close()

# calculate p-value or likelihood ratio
thres = 'likelihood ratios ' if args.isDuplex else 'p-values '
print("Calculating " + thres + str(datetime.datetime.now()) + "\n")
outfile1 = 'intermediate/noThres.' + args.outPrefix + '.VariantList.long.txt'
outfile2 = 'intermediate/' + args.outPrefix + '.VariantList.long.txt'
outfile_lod = 'intermediate/' + args.outPrefix + '.umi_depths.lod.bedgraph'

if args.isDuplex:
   pValCmd = ' '.join(['Rscript', pValCode_dup, args.dupBkgErrorDist, './', outfile1, outfile2, str(args.minAltUMI)])
else:
   pValCmd = ' '.join(['Rscript', pValCode_sin, args.sinBkgErrorDist, './', outfile1, bkgFileName, str(seed), str(nsim), outfile2, outfile_lod, args.outPrefix, str(rpu), str(args.minAltUMI), str(args.inputVCF).lower()])
if not empty:
   subprocess.check_call(pValCmd, shell=True)
   print("completed p-values " + str(datetime.datetime.now()) + "\n")
   # make VCFs
   vcf.makeVcf('./', outfile2, args.outPrefix, args.refGenome, args.isDuplex)
else:
   print("empty BED or BAM, no variants detected " + str(datetime.datetime.now()) + "\n")

# remove intermediate files
#os.remove('hp.roi.bed')
#os.remove('rep.roi.bed')
#os.remove('sr.roi.bed')
os.remove(outfile1)
# log run completion
timeEnd = datetime.datetime.now()
print("completed running " + str(timeEnd) + "\n")
print("total runtime: "+ str(timeEnd-timeStart) + "\n")  
   