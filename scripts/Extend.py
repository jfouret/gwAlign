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
parser.add_argument('-queue', metavar='queue', required=False, help="queue for PBS",default='batch')
parser.add_argument('-refs', metavar='/path', required=True, help="references to be used csv file with id;fasta;bam in the order of priority  \n"+
	"id: the id of the species used in the exon alignment file\n"+
	"fasta: path to the fasta file used for short read alignments\n"+
	"bam: path of the sorted and indexed bam file")
parser.add_argument('-spec', metavar='name', required=True, help="name of the species to add")
parser.add_argument('-picard', metavar='/path', required=False, help="picard jar path",default='/export/bin/picard-tools-2.1.0/picard.jar')
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
samtools=Command('samtools',min='1.3')
samtools.versionCtrl()
samtools.log()
java=Command('java')
java.log()
picard_cmd=java.create(options={'-jar':args.picard})
picard=Command(picard_cmd,picard_cmd+' CheckFingerprint --version 2>&1 | sed \'s/(.*$//g\'')
picard.log()
bam2consensus=gitCommand(Git(gitRepository),'bam2consensus.py','bin')
bam2consensus.log()

#Getting absolute path
alnFileName=os.path.abspath(args.aln)
refsFileName=os.path.abspath(args.refs)

#Define var
reID=re.compile('^(uc\w*)(\.\d*)_(\w*)_(\d*)_(\d*) (\d*) (\d*) (\d*)\s?(\w+:\d+-\d+[+-]{1})?$')
#Group
#kgID= 1+2
#spec= 3
#exon num = 4
#exon tot = 5
#reg and strand = 9
reIDhg19=re.compile('^(uc\w*)(\.\d*)_hg19_(\d*)_(\d*) (\d*) (\d*) (\d*)\s?(\w+:\d+-\d+[+-]{1})?$')

# configure batch
Nlim=80
N=0
iterBatch=1
index=0

# extract info from refFile
refDict=dict() # stock the regions + strands for each specie 
specOrder=list() # stock the order of the species
with open(args.refs,'r') as refFile:
	for line in refFile.readlines():
		line=line.rstrip()
		if line!='':
			spec,fasta,bam=line.split(';')
			refDict[spec]=''

			## check for the presence of indexation files

			prefixRef=fasta.replace('.fa','')
			prefixRef=prefixRef.replace('.fna','')
			prefixRef=prefixRef.replace('.fasta','')

			faiRef=fasta+'.fai'
			dictRef=prefixRef+'.dict'

			if not os.path.isfile(dictRef):
				dictOpt={
					'REFERENCE':fasta,
					'OUTPUT':dictRef
				}
				submitOneShell(picard.create(options=dictOpt,sep='=',subprogram='CreateSequenceDictionary'))

			if not os.path.isfile(faiRef):
				faiPos=[fasta]
				submitOneShell(samtools.create(positionals=faiPos,subprogram='faidx'))

			#TODO check if spec names are corect(matching those in the multi aln file)
			#TODO check the presence of hg19

fasta_sequences = SeqIO.parse(open(alnFileName),'fasta')
# step 1 - iterate over the aln File to compute consensus seq
cmdList=list()
NbSpec=0 # 0 species matched yet
NbSpecLim=len(refDict.keys()) # number of species to match before to start searching consensus

for fasta in fasta_sequences:
	name, sequence = fasta.description, str(fasta.seq)
	m=reID.match(name)
	time.sleep(0.05)#debug
	print(name)#debug
	if m and (m.group(3) in refDict.keys()):
		if m.group(9)==None:
			region=''
		else:
			region=m.group(9)
		refDict[m.group(3)]=m.group(9)
		NbSpec=+1
	if NbSpec==NbSpecLim:
		time.sleep(0.05)#debug
		print('GO==>'+name)#debug
		consensusOutPut=rootedDir.results+'/'+m.group(1)+m.group(2)+'/exon'+m.group(4)
		regionOrdList=list()
		for spec in specOrder:
			regionOrdList.append(refDict(spec))
		regions=';'.join(regionOrdList)# regions != region
		bam2consensusOpt={
			'-outDir':consensusOutPut,
			'-reference':refsFileName,
			'-reg':regions
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
	m=reIDhg19.match(name)
	outAln.write('>'+name+"\n")
	outAln.write(sequence+"\n")
	if m:
		pos=m.group(8).split(':')[1][:-1].split('-')
		pos=[int(pos[0]),int(pos[1])]
		seqLen=max(pos)-min(pos)+1
		consensusOutPut=rootedDir.results+'/'+m.group(1)+m.group(2)+'/exon'+m.group(3)
		while not os.path.exists(consensusOutPut+'/fileRooting.pkl'): # infinite loop if bug . . . (sol ? write an error file in error ? through def main: and main())
			time.sleep(10)
		header='>'+m.group(1)+m.group(2)+'_'+args.spec+'_'+m.group(3)+'_'+m.group(4)+' '+m.group(5)+' '+m.group(6)+' '+m.group(7)+' '+" BAM:BAM-BAM\n"
		outAln.write(header)
		with open(consensusOutPut+'/results/hg19mapped.fa','r') as consensusFile:
			seq=''.join(consensusFile.readlines()[1:]).replace("\n",'')
			if len(seq)!=seqLen:
				sys.stderr.write('ERROR SEQ SIZE for '+name)
				sys.exit(1)
		outAln.write(seq+"\n")
outAln.close()

# finish and saveRoot
saveRoot(rootedDir)

sys.exit(0)


reIDhg19=re.compile('^(uc\w*)(\.\d*)_hg19_(\d*)_(\d*) (\d*) (\d*) (\d*)\s?(\w+:\d+-\d+[+-]{1})?$')