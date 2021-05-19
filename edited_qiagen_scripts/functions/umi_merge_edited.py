import os
import os.path
import subprocess
# get params
readSet          = cfg.readSet
tagNameUmi       = cfg.tagNameUmi
tagNamePrimer    = cfg.tagNamePrimer
tagNameResample  = cfg.tagNameResample
samtoolsDir      = cfg.samtoolsDir
deleteLocalFiles = cfg.deleteLocalFiles
numCores         = cfg.numCores
# sort the original molecules text file by readId (column 15 hard-coded)
fileNameIn = readSet + ".umi_mark.alignments.txt"
cmd = "sort -k15,15 -t\| -T./ --parallel={1} {0} > {0}.tmp.txt 2> \
    {2}.umi_merge.shell.log".format(fileNameIn, numCores,readSet)
subprocess.check_call(cmd, shell=True)
os.rename(fileNameIn + ".tmp.txt", fileNameIn)
print("umi_merge: done sorting molecule-marked alignments text file by read id")
# save BAM header into output file
cmd = samtoolsDir + "samtools view -H " \
+ bamFileIn \
+ " 1> " + readSet + ".umi_merge.temp2.sam" \
+ " 2>>" + readSet + ".umi_merge.shell.log"
subprocess.check_call(cmd, shell=True)
# convert BAM to SAM, rather than use pipe, to enable more parallel Linux gnu 
# sort on very large files (at the expense of more disk I/O)
cmd = samtoolsDir + "samtools view " \
+ " -@ " + numCores  \
+ "    " + bamFileIn \
+ " >  " + readSet + ".umi_merge.temp0.sam" \
+ " 2>>" + readSet + ".umi_merge.shell.log"
subprocess.check_call(cmd, shell=True)
print("umi_merge: done converting input BAM to SAM")
# delete input BAM file if local
if deleteLocalFiles and len(os.path.dirname(bamFileIn)) == 0:
    os.remove(bamFileIn)
    
# sort original BAM file by readId, using Linux sort 
# (NOT samtools sort -n !!!) (could use sambamba instead)
cmd = "sort -k1,1 -T./ --parallel=" + numCores \
+ " "     + readSet + ".umi_merge.temp0.sam" \
+ " 1>> " + readSet + ".umi_merge.temp1.sam" \
+ " 2>> " + readSet + ".umi_merge.shell.log"
subprocess.check_call(cmd, shell=True)
print("umi_merge: done sorting SAM by read id")
os.remove(readSet + ".umi_merge.temp0.sam")
# open SAM files, init read counters
fileout = open(readSet + ".umi_merge.primers.txt","w")
samIn   = open(readSet + ".umi_merge.temp1.sam", "r")
samOut  = open(readSet + ".umi_merge.temp2.sam", "a")  # note this is an append, 
# because header already written
numReadsUmiFile = 0
numReadsSamFile = 0
# read the molecule-marked read alignment text file
with open(fileNameIn) as f:
    line = f.readline()
for line in open(fileNameIn, "r"):
    # parse line
    (pChrom, pStrand, mtLoc, mt, mtReads, mtReads_, mtReadIdx, isResample, fragLen, pLoc5, primer, umiSeq, alignLocR, alignLocP, readId, read1L, read1R, cigar1, read2L, read2R, cigar2) = line.strip().split("|")
    numReadsUmiFile += 1
    # format SAM tag with UMI tag
    tagUmi    = "-".join((pChrom, pStrand, mtLoc, mt))
            
    # spin the sam forward to this read
    numReadsFound = 0
    while True:
        # parse line
        line = samIn.readline()
        if len(line) == 0:
            break
        vals = line.strip().split("\t")
        # get next row if readId not found
        if vals[0] != readId:
            continue

        # put the UMI tag first, in case need to sort the SAM by unique molecule
        outvec = vals[0:11]
        outvec.append(tagNameUmi      + ":Z:" + tagUmi)
        outvec.append(tagNameResample + ":i:" + isResample)
        outvec.extend(vals[11:])
  
        # write output sam
        samOut.write("\t".join(outvec))
        samOut.write("\n")
        numReadsFound   += 1
        numReadsSamFile += 1
        
        # write auxillary file containing the primer at the start of each molecule
        if numReadsFound == 1 and int(mtReadIdx) == 0:
            tagPrimer = "-".join((pChrom, pStrand, pLoc5, str(len(primer))))
            fileout.write("{}|{}\n".format(tagUmi,tagPrimer))
         
        # found both R1 and R2 rows
        if numReadsFound == 2:
            break
         
# done
fileout.close()
samIn.close()
samOut.close()
print("umi_merge: done merging unique molecule tag to sam file")
os.remove(readSet + ".umi_merge.temp1.sam")
   
# debug check - make sure all reads found
if numReadsSamFile != 2 * numReadsUmiFile:
    raise Exception("umi_merge: synchronization error")
    
# convert final file to BAM (not really necessary)
cmd = samtoolsDir + "samtools view -1" \
+ " -@ " + numCores  \
+ "    " + readSet + ".umi_merge.temp2.sam" \
+ " 1> " + readSet + ".umi_merge.bam" \
+ " 2>>" + readSet + ".umi_merge.shell.log" 
subprocess.check_call(cmd, shell=True)
print("umi_merge: done converting SAM to BAM")
os.remove(readSet + ".umi_merge.temp2.sam")
