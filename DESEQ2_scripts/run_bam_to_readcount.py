import sys
import numpy as np
import os
import pdb
import gzip
import argparse
import subprocess
import os.path
import string
import re
from time import sleep

parser = argparse.ArgumentParser()
parser.add_argument("-i", help="list of bamfiles",type=str)
parser.add_argument("-n", help="group name",default='a',type=str)
parser.add_argument("-j", help="number of jobs to run at the same time",default='50',type=int)
parser.add_argument("-c", help="# seconds to wait before checking",default='120',type=int)
args = parser.parse_args()

Infile = args.i
Name = args.n
Njobs = 1+args.j
Wait = args.c

###################################################
##program description
###################################################
#takes in a list of bam files, one per line and runs them through bam_to_readcount_manager.py
#run from the folder the bam files are in. 

###################################################
##definitions
###################################################

def runstuff(cmd):
	subprocess.run(cmd, shell=True)
	print("submitted "+cmd)
	
def get_out(cmd):
	stuff = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	print("returned "+cmd)
	return stuff.stdout

#checks how many jobs are running and waits until less than Njobs are running. 
def wait5(cmd):
	while True:
		stuff = get_out("bjobs -g /lorenzk/"+cmd+" | wc -l")
		count = re.findall(r'\d+', str(stuff))
		if "No unfinished job found" in str(stuff):
			print("returned: No unfinished job found")
			break
		elif (int(count[0]) < Njobs):
			print("count less than njobs")
			break
		sleep(Wait)

###################################################
##read in inputs & make tuples of all outfiles & commands
###################################################
Romp = []

for line in open(Infile,'r'):
	#remove white space at end of line
	line2 = line.rstrip()
	#save as dictionaries to prevent duplication
	Romp.append(line2)

###################################################
##step through list and submit Njobs at a time
###################################################

for bamfile in Romp:
	
	#extract name from .bam file
	
	repname, toss = bamfile.split('.')
	
	#submit job
	wait5(Name)
	runstuff("bsub -g /lorenzk/"+Name+" -o manager_"+repname+".out python /project/voightlab_01/lorenzk/CEET/DESEQ2_scripts/bam_to_readcount_manager.py -n "+repname+" -i "+bamfile)
	sleep(5) #pause before submitting more.
