#!/usr/bin/python
import argparse
gitRepository='SEDMATCHGITREPO'
version='SEDMATCHGITVERSION'
year=2016
author='Julien Fouret'
contact='julien.fouret12@uniagro.fr'
##parse argument
parser = argparse.ArgumentParser(description='Parse aln and Gblocks file to get positions',epilog="Version : "+str(version)+"\n"+str(year)+"\nAuthor : "+author+" for more informations or enquiries please contact "+contact,formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-gene_dir', metavar='/path', required=True, help="gene ID directory")
parser.add_argument('-ref', metavar='ref_species' , required=False, help="name of the reference specie",default='hg19')
parser.add_argument('-back', metavar='spec1,spec2' , required=True, help="name of the reference specie")
parser.add_argument('-fore', metavar='spec1,spec2' , required=True, help="name of the reference specie")
parser.add_argument('-name', metavar='name', required=True, help="name of phasing-filtering")
parser.add_argument('-v',required=False,action='store_true',help="verbose")

args=parser.parse_args()

import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
import re
import errno

def submitOneShell(cmdString): 
	""" 
	@summary: Submit only one command via a shell and return the stderr and stdout, if stdout or stderr are empty they are written respectively in ./logs/outputs.txt and ./logs/errors.txt 
	@param cmdString: Command to be run on the cluster 
	@type cmdString: str 
	@return: the stdout and the stderr accessible by out and err keys 
	@rtype: dict 
	@requires: subprocess 
	""" 
	import subprocess 
	child=subprocess.Popen(cmdString,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE) 
	out,err=child.communicate() 
	return {'out':out,'err':err} 

os.chdir(args.gene_dir)

try: 
	os.makedirs(args.name) 
except OSError as exc: 
	if exc.errno == errno.EEXIST and os.path.isdir(args.name): 
		pass 
	else: 
		raise 

al_tmp=args.name+'/codon_aln.fa'
al_out=args.name+'/codon_aln_blocks.fa'
pos_out=args.name+'/posDict.tab'

backList=args.back.split(',')
foreList=args.fore.split(',')
refWrite=''
backWrite=''
foreWrite=''

if not os.path.exists('ref_macse_NT.fasta'):
	sys.stderr.write('Error at phase_filter step no ref_macse_NT.fasta file : '+args.name+' '+args.gene_dir+"\n")
	sys.exit()

fasta_sequences = SeqIO.parse(open('ref_macse_NT.fasta'),'fasta')
for fasta in fasta_sequences:
	name, sequence = fasta.description, str(fasta.seq)
	if name in backList:
		backWrite+='>'+name+"\n"+sequence+"\n"
	if name in foreList:
		foreWrite+='>'+name+"\n"+sequence+"\n"
	if name==args.ref:
		refWrite+='>'+name+"\n"+sequence+"\n"

with open(al_tmp,'w') as outAln:
	outAln.write(refWrite+backWrite+foreWrite)

# remove only gap ! ! !
fasta_sequences = SeqIO.parse(open(al_tmp),'fasta')
seqDict=dict()
for fasta in fasta_sequences:
	name, sequence = fasta.description, str(fasta.seq)
	seqDict[name]=sequence
for i in range(0,len(seqDict[args.ref]))[::-1]:
	remove=True
	for species in seqDict.keys():
		if seqDict[species][i]!='-':
			remove=False
	if remove:
		for species in seqDict.keys():
			seqDict[species]=seqDict[species][:i]+seqDict[species][i+1:]

with open(al_tmp,'w') as outAln:
	for species in seqDict.keys():
		outAln.write('>'+species+"\n"+seqDict[species]+"\n")


spec70=str(int(round((0.5*(len(backList)+1)))))

command='cd '+args.name+' ; Gblocks codon_aln.fa -p=Yes -t=c -b1='+spec70+' -b2='+spec70+' -b3=1 -b4=6 -b5=h -k=y'
GblockProc=submitOneShell(command)

if GblockProc['err']!='':
	print(GblockProc['err'])

maskSeq=''
active=False
with open(al_tmp+'-gbMask','r') as maskFile:
	for line in maskFile.readlines():
		line=line.rstrip()
		if active:
			maskSeq+=line.replace(' ','')
		elif line=='>P1;Gblocks':
			active=True

#checkblockRE phasing 
phasedMask=''
count=0
last_letter='#'
for letter in maskSeq:
	count+=1
	if letter=='#':
		if count!=1 and last_letter=='.':
			phasedMask+='.'
		else:
			phasedMask+=letter
	elif letter=='.':
		if count==2 and last_letter=='#':
			phasedMask=phasedMask[:-1]+'..'
		if count==3 and last_letter=='#':
			phasedMask=phasedMask[:-2]+'...'
		else:
			phasedMask+=letter
	if count==3:
		count=0
	last_letter=letter


last_letter='.'
active=False
count=0

intervalls=list()
for letter in phasedMask:
	if active:
		if letter==".":
			end=count
			active=False
			intervalls.append({'start':start,'end':end})
	elif last_letter=='.' and letter=="#":
		active=True
		start=count
	count+=1
	last_letter=letter
if last_letter=='#':
	end=count
	intervalls.append({'start':start,'end':end})



fasta_sequences = SeqIO.parse(open(al_tmp),'fasta')
with open(al_out,'w') as outAln:
	for fasta in fasta_sequences:
		name, sequence = fasta.description, str(fasta.seq)
		if name==args.ref:
			refSequence=sequence
		newSequence=''
		for intervall in intervalls:
			newSequence+=sequence[intervall['start']:intervall['end']]
		outAln.write('>'+name+"\n"+newSequence+"\n")

posDict=dict()

exonList=list()
with open('exons.pos','r') as exonFile:
	for line in exonFile.readlines():
		line=line.rstrip()
		if line !='':
			exonList.append(line)

posDict['ref']=list()
posDict['aln']=list()
posDict['block']=list()
posDict['exon']=list()

count=0
refCount=1
alnCount=1
blockCount=1

if args.v:
	print('ref length: '+str(len(refSequence)))
	print(refSequence)
	print('Mask length: '+str(len(maskSeq)))
	print(maskSeq)
	print('phasedMask length: '+str(len(phasedMask)))
	print(phasedMask)

for letter in refSequence:
	if letter=='-' or letter=='!' or letter=='?':
		posDict['ref'].append('.')
	else:
		posDict['ref'].append(str(refCount))
		refCount+=1
	posDict['aln'].append(str(alnCount))
	alnCount+=1
	if phasedMask[count]=='#':
		posDict['block'].append(str(blockCount))
		blockCount+=1
	else:
		posDict['block'].append('.')
	count+=1
	posDict['exon'].append(exonList[refCount-2])

with open(pos_out,'w') as posFile:
	posFile.write("\t".join(['reference','alignment','exon','block'])+"\n")
	for i in range(0,len(posDict['aln'])):
		posFile.write("\t".join([posDict['ref'][i],posDict['aln'][i],posDict['exon'][i],posDict['block'][i]])+"\n")


