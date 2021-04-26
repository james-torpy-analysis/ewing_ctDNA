home_dir = '/share/ScratchGeneral/jamtor/'
project_dir = home_dir + 'projects/ewing_ctDNA/'
fq_dir = project_dir + 'raw_files/'
result_dir = project_dir + 'results/'
stats_dir = result_dir + 'stats/'
bam_dir = result_dir + 'picard/bams/'

import pysam
import sys
import re
import os

def minmax(val_list):
    min_val = min(val_list)
    max_val = max(val_list)
    return (min_val, max_val)

orig_fq = fq_dir + "409_004_D9YWF_GTAGAGGA-CTCTCTAT_L001_R1.fastq"








#inbam_name = sys.argv[1]
#stats_name = sys.argv[2]
#discordant_name = sys.argv[3]
#split_name = sys.argv[4]
inbam_name = bam_dir + '409_004_D9YWF_GTAGAGGA-CTCTCTAT_L001.initial_mapped_and_UMI.bam'
stats_name = stats_dir + '409_004_D9YWF_GTAGAGGA-CTCTCTAT_L001.uncollapsed.bam.bc'
discordant_name = bam_dir + '409_004_D9YWF_GTAGAGGA-CTCTCTAT_L001.uncollapsed.discordant.bam'
split_name = bam_dir + '409_004_D9YWF_GTAGAGGA-CTCTCTAT_L001.uncollapsed.split.bam'
out_prefix = stats_dir + re.sub('initial_mapped_and_UMI.bam$', '', os.path.basename(inbam_name))

print('Input bam = ', inbam_name)
print('Stats file = ', stats_name)
print('Discordant bam = ', discordant_name)
print('Split bam = ', split_name)
print('Output prefix = ', out_prefix)

# define filename using pysam:
inbam = pysam.Samfile(inbam_name, 'rb')

# define counters for reads processed and written:
n = 0
w = 0

print('Fetching lengths of sequences...')

lengths=[]
for read in inbam.fetch(until_eof=True):
    n += 1
    # add read seq length to lengths list:
    lengths.append(len(read.seq))

# plot length distribution:




	assert umi1 is not None
	# remove last 11 characters:
	umi1_fix = re.sub('...........$', '', umi1)
	# replace original umi:
	read.set_tag('ZB', umi1_fix, value_type='Z')
	# fetch RX tag:
	umi2 = read.get_tag('RX')
	assert umi2 is not None
	# remove last 11 characters:
	umi2_fix = re.sub('[A-Z]-', '', 
		re.sub('...........$', '', umi1)
	)
	# replace original umi:
	read.set_tag('RX', umi2_fix, value_type='Z')

	# if length = 12, write to file:
	if len(umi1_fix) == 12 & len(umi2_fix) == 12:
		if '-' not in umi2_fix:
			w += 1
			outbam.write(read)
		else:
			print('Hyphen still in RX UMI sequence: ', umi2_fix, ' read=', read)
			break
	else:
		print(
			'UMI sequence wrong length: ZB=', umi1_fix, ' RX=', umi2_fix, ' read=', read
		)

print('All reads processed, checked and written')

inbam.close()
outbam.close()