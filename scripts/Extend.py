#!/usr/bin/python
import argparse
gitRepository='SEDMATCHGITREPO'
version='SEDMATCHGITVERSION'
year=2016
author='Julien Fouret'
contact='julien@fouret.me'

parser = argparse.ArgumentParser(description='Take an region multialignment file and add one speacies from bam file',epilog="Version : "+str(version)+"\n"+str(year)+"\nAuthor : "+author+" for more informations or enquiries please contact "+contact,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-outDir', metavar='/path', required=True, help="Output directory")
parser.add_argument('-aln', metavar='/path', required=True, help="Alignment to extend")
parser.add_argument('-bam', metavar='/path', required=True, help="bam file")
parser.add_argument('-queue', metavar='queue', required=False, help="queue for PBS",default='batch')
parser.add_argument('-ref', metavar='/path', required=True, help="reference fasta file used for alignment")
parser.add_argument('-spec', metavar='name', required=True, help="name of the species to add")
args=parser.parse_args()

# Import libraries
import sys
import os
import re
import time
from jupype import *
from Bio import SeqIO
from Bio.Seq import Seq

# Create output directory structure and logs
rootedDir=RootDir(args.outDir,pbs=True)
rootedDir.logs.writeArgs(args)

# definition of used software
bam2consensus=gitCommand(Git(gitRepository),'bam2consensus.py','bin')
bam2consensus.log()

#Getting absolute path
alnFileName=os.path.abspath(args.aln)
bamFileName=os.path.abspath(args.bam)
refFileName=os.path.abspath(args.ref)

#Define var
reID=re.compile('^(uc\w*)(\.\d*)_hg19_(\d*)_(\d*) (\d*) (\d*) (\d*) (\w+:\d+-\d+)([+-])$')
#Group
#kgID= 1+2
#exon num = 3 
#exon tot = 4
#reg = 8
#strand = 9

# configure batch
Nlim=80
N=0
iterBatch=1

fasta_sequences = SeqIO.parse(open(alnFileName),'fasta')
# step 1 - iterate over the aln File to compute consensus seq
cmdList=list()
for fasta in fasta_sequences:
	name, sequence = fasta.description, str(fasta.seq)
	m=reID.match(name)
	#time.sleep(0.05)#debug
	#print(name)#debug
	if m:
		#time.sleep(0.5)#debug
		#print(name)#debug
		consensusOutPut=rootedDir.results+'/'+m.group(1)+m.group(2)+'/exon'+m.group(3)
		bam2consensusOpt={
			'-outDir':consensusOutPut,
			'-bam':bamFileName,
			'-ref':refFileName,
			'-reg':m.group(8)
		}
		cmdList.append(bam2consensus.create(wd=rootedDir.results,options=bam2consensusOpt))
		N+=1
		if N==Nlim:
			N=0
			batchName='gwAlign_Batch_'+str(iterBatch)
			submitQsubWithPBS(createPBS(cmdList,batchName,queue=args.queue,ppn='2',workdir=rootedDir.results))
			cmdList=list()
			iterBatch+=1
if N!=0:
	batchName='gwAlign_Batch_'+str(iterBatch)
	submitQsubWithPBS(createPBS(cmdList,batchName,queue=args.queue,workdir=rootedDir.results))

# step 2 - iterate over the aln File to add consensus in outFile
fasta_sequences = SeqIO.parse(open(alnFileName),'fasta')
outAln=open(rootedDir.reports+'/algn_extended.fa','w')
for fasta in fasta_sequences:
	name, sequence = fasta.description, str(fasta.seq)
	m=reID.match(name)
	outAln.write('>'+name+"\n")
	outAln.write(sequence+"\n")
	if m:
		strand=m.group(9)
		pos=m.group(8).split(':')[1].split('-')
		pos=[int(pos[0]),int(pos[1])]
		seqLen=max(pos)-min(pos)+1
		consensusOutPut=rootedDir.results+'/'+m.group(1)+m.group(2)+'/exon'+m.group(3)
		while not os.path.exists(consensusOutPut+'/fileRooting.pkl'):
			time.sleep(10)
		header='>'+m.group(1)+m.group(2)+'_'+args.spec+'_'+m.group(3)+'_'+m.group(4)+' '+m.group(5)+' '+m.group(6)+' '+m.group(7)+' '+m.group(8)+m.group(9)+" BAM:BAM-BAM\n"
		outAln.write(header)
		with open(consensusOutPut+'/results/hg19mapped.fa','r') as consensusFile:
			seq=''.join(consensusFile.readlines()[1:]).replace("\n",'')
			if len(seq)!=seqLen:
				print('ERROR SEQ SIZE for '+name)
				sys.exit(1)
		if strand=='-':
			seq = str(Seq(seq).reverse_complement())
		outAln.write(seq+"\n")
outAln.close()
