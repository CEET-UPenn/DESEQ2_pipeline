#import sys
#import numpy as np
#import os
#import pdb
#import gzip
import argparse
import subprocess
#import os.path
import string
from time import sleep

parser = argparse.ArgumentParser()
parser.add_argument("-i", help=".bam filepath",type=str)
parser.add_argument("-n", help="sample name",default='',type=str)
parser.add_argument("-c", help="# seconds to wait before checking",default='120',type=int)
args = parser.parse_args()

Infile = args.i
Name = args.n
Wait = args.c


###################################################
##definitions
###################################################

def runstuff(cmd):
	subprocess.call(cmd, shell=True)
	print("submitted "+cmd)
	
def get_out(cmd):
	stuff = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	#print("returned "+cmd)
	stdout = []
	while True:
		l = stuff.stdout.readline()
		line = l.rstrip()
		if line == b'' and stuff.poll() != None:
			break
		else:
			stdout.append(line)
	return stdout

#checks if previous job name is still running
def wait5(cmd):
	while True:
		stuff = get_out("bjobs -J "+cmd)
		count = stuff[0]
		if "is not found" in str(count):
			print("No unfinished job found, proceeding")
			break
		print("waiting: previous job still running")
		sleep(Wait)


###################################################
##run sample through the pipeline
###################################################


#sort bam
#assemble command
sort = "samtools sort -o "+Name+"_sorted.bam "+ Infile

#submit job
print("sorting input bam...")
runstuff("bsub -J "+Name+" -q voltron_long -o "+Name+"sort_lsf_log.out "+ sort)

sleep(10)
wait5(Name)


#mark duplicates
#assemble mark duplicates command
markdup = "python -u /project/voightlab_01/lorenzk/CEET/DESEQ2_scripts/run_MarkDuplicates.py "+Name+"_sorted.bam "+Name

#submit job
print("marking duplicates...")
runstuff("bsub -J "+Name+" -q voltron_long -o "+Name+"markdup_lsf_log.out "+ markdup)

sleep(10)
wait5(Name)

#rename bam
#no reason to bsub this, just run it. 
rename = "samtools view -H "+Name+"_sorted.md.bam  | sed \"s/SN:/SN:chr/g\" | samtools reheader - "+Name+"_sorted.md.bam > "+Name+"_sorted.id.bam"
print("renaming bam contigs...")
runstuff(rename)
sleep(10)

#run rnaseqc
#assemble rnaseqc command
rnaseqc = "python /project/voightlab_01/lorenzk/CEET/DESEQ2_scripts/run_rnaseqc.py /project/voightlab_01/lorenzk/CEET/DESEQ2_scripts/ref/gencode.v34.GRCh38.genes.collapsed_only.gtf "+Name+"_sorted.id.bam "+Name+" --stranded rf"

#submit job
print("submitting rnaseqc job")
runstuff("bsub -J "+Name+" -q voltron_long -o "+Name+"_rnaseqc_lsf_log.out "+ rnaseqc)


#remove intermediate .bam files
runstuff("rm " +Name+"_sorted.bam ")
runstuff("rm " +Name+"_sorted.md.bam ")
