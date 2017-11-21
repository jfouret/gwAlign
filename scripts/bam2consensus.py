#!/usr/bin/env python
import argparse
gitRepository='SEDMATCHGITREPO'
version='SEDMATCHGITVERSION'
year=2017
author='Julien Fouret'
contact='julien@fouret.me'

parser = argparse.ArgumentParser(description='Retrive consensus sequence of a given region based on a sorted bam file',epilog="Version : "+str(version)+"\n"+str(year)+"\nAuthor : "+author+" for more informations or enquiries please contact "+contact,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-outDir', metavar='/path', required=True, help="Output directory")
parser.add_argument('-reference', metavar='/path', required=True, help="csv file with id;fasta;bam in the order of priority  \n"+
	"id: the id of the species used in the exon alignment file\n"+
	"fasta: path to the fasta file used for short read alignments\n"+
	"bam: path of the sorted and indexed bam file")
parser.add_argument('-reg', metavar='chr:start-end+,chr:start-end-,...', required=True, help="region of interest, please respect the format and the order of the ref file")
parser.add_argument('-minCov', metavar='N', required=True, help="min coverage for base calling")
parser.add_argument('-minCovVar', metavar='N', required=True, help="min coverage for base calling")
parser.add_argument('-name', metavar='STR', required=True, help="name fo consensus")
#parser.add_argument('-pval', metavar='N', required=True, help="pval to call the alternative base (H0: reference)")
#parser.add_argument('-gatk', metavar='/path', required=False, help="gatk jar path",default='SEDMATCHGATK')
#parser.add_argument('-picard', metavar='/path', required=False, help="picard jar path",default='SEDMATCHPICARD')
#parser.add_argument('-graphCov', action='store_true', help="Graph the coverage")

args=parser.parse_args()

# Import libraries
import sys
import numpy as np
import os
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import ambiguous_dna
from jupype import *

# Create output directory structure and logs
rootedDir=RootDir(args.outDir)
rootedDir.logs.writeArgs(args)

# Definition of used software
platypus=Command("python /export/source/archive/Platypus_0.8.1/Platypus.py",'echo "(see program log)"')
platypus.log()

minCov=int(args.minCov)
minCovVar=int(args.minCovVar)

iupacDict={
	"GT":"K",
	"TG":"K",
	"AC":"M",
	"CA":"M",
	"CG":"S",
	"GC":"S",
	"AT":"W",
	"TA":"W",
	"AG":"R",
	"GA":"R",
	"CT":"Y",
	"TC":"Y",
}

#get absolute path from inputs
regions=args.reg.split(',')
index=0
refDict=dict()
specOrder=list() # priority of species (closely related ?)
with open(args.reference,'r') as refFile:
	for line in refFile.readlines():
		line=line.rstrip()
		if line!='':
			spec,fasta,bam=line.split(';')
			region=regions[index]
			index+=1
			if region!='' and not (";" in region): # avoid multi locus genomic alignment
				refDict[spec]={
					'fasta':fasta,
					'bam':bam,
					'reg':region
				}
				specOrder.append(spec)
os.chdir(rootedDir.results)

#print(refDict)

## step 1 : Find the species with a maximum of sequencing yield based on mapping

samfile=dict()
	
chosenSpec=None
for spec in refDict.keys():
	samfile[spec]=pysam.AlignmentFile(refDict[spec]['bam'], "rb" )
	reg=refDict[spec]['reg']
	contig,pos=reg[:-1].split(':')
	pos=pos.split('-')
	pos=[int(pos[0]),int(pos[1])]
	start=pos[0]-1 # 0based inclusive
	end=pos[1] # 0based exclusive
	tA,tT,tC,tG=samfile[spec].count_coverage(contig=contig,start=start,stop=end)  #array.arrays
	tA=np.array(tA)
	tT=np.array(tT)
	tC=np.array(tC)
	tG=np.array(tG)
	tcoverage=tA+tT+tC+tG
	tYield=tcoverage.sum()
	#print(spec)#debug
	if chosenSpec==None:
		chosenSpec=spec
		Yield=tYield
		A=tA
		T=tT
		C=tC
		G=tG
		coverage=tcoverage
	elif tYield>Yield:
		chosenSpec=spec
		Yield=tYield
		A=tA
		T=tT
		C=tC
		G=tG
		coverage=tcoverage
	elif tYield==Yield:
		if specOrder.index(chosenSpec)>specOrder.index(spec):
			chosenSpec=spec
			Yield=tYield
			A=tA
			T=tT
			C=tC
			G=tG
			coverage=tcoverage
for spec in refDict.keys():
	if spec!=chosenSpec:
		samfile[spec].close()

## step 2 - Variant calling and consensus

if Yield==0:
	## step 2.1 - create empty consensus if no yield
	with open(rootedDir.results+'/consensus.fa','w') as csFile:
		csFile.write("consensus\n---"+"\n")
else:

	mkdirp(chosenSpec)
	# get mapped reads in case of the iteration is to be done
	with open('mapped.fq','w') 	as mapped:
		for read in samfile[chosenSpec].fetch(region=refDict[chosenSpec]['reg'][:-1]): 
			mapped.write(">%s\n%s\n+\n%s\n"  % (read.query_name, read.query_sequence, "".join(map(chr,[x + 33 for x in read.query_qualities]))))
	samfile[chosenSpec].close()
	## step 2.1 - Variant calling
	#define options and positional args for software, considering PCR free # TODO put this in argument by default
	PlatypusOpt={
		'--regions':refDict[chosenSpec]['reg'][:-1], # [:-1] without strand information yet ...
		'--refFile':refDict[chosenSpec]['fasta'],
		'--bamFiles':refDict[chosenSpec]['bam'],
		'--assemble':'1',
		'--output':rootedDir.results+'/'+chosenSpec+'/call.vcf',
		'--logFileName':rootedDir.results+'/'+chosenSpec+'/call.log',
		'--filterReadPairsWithSmallInserts':'0',
		'--filterReadsWithDistantMates':'0',
		'--filterReadsWithUnmappedMates':'0',
		'--coverageSamplingLevel':'50',
		'--filterDuplicates':'0'
		}

	submitOneShell(platypus.create(options=PlatypusOpt,subprogram='callVariants'))

	#read the vcf file

	vcf=pysam.VariantFile(rootedDir.results+'/'+chosenSpec+'/call.vcf')

	variant=dict()
	for rec in vcf.fetch():
		variant[rec.pos]=dict()
		variant[rec.pos]['ref']=rec.ref
		variant[rec.pos]['alt']=rec.alts
		variant[rec.pos]['info']=rec.info
		variant[rec.pos]['filter']=rec.filter.keys()
		variant[rec.pos]['sample']=rec.samples[rec.samples.keys()[0]]

	spec=chosenSpec
	reg=refDict[chosenSpec]['reg']
	contig,pos=reg[:-1].split(':')
	pos=pos.split('-')
	pos=[int(pos[0]),int(pos[1])]
	start=pos[0]
	end=pos[1]
	seqLen=max(pos)-min(pos)+1
	strand=refDict[chosenSpec]['reg'].split(':')[1][-1]

	# Get left and right sequences to perform iteration
	fastaFile=pysam.FastaFile(refDict[chosenSpec]['fasta'])
	contig_length=fastaFile.get_reference_length(contig)
	if start!=1:
		left=fastaFile.fetch(region=contig+":"+str(max(1,start-200))+"-"+str(start-1))
	else:
		left=''
	if end!=contig_length:
		right=fastaFile.fetch(region=contig+":"+str(end+1)+"-"+str(min(contig_length,end+200)))
	else:
		right=''
	fastaFile.close()

	# code to get consensus
	consensus=""
	for pos in range(start,1+end):
		index=pos-start
		# Have a variant been called with sufficient coverage ?
		if (pos in variant.keys()) and (variant[pos]["sample"]["NR"]>=minCovVar):
			# if unique SNP with no variants
			if (len(variant[pos]['alt'][0]+variant[pos]['ref'][0])==2) and ((variant[pos]["sample"]["GT"]!= (0,1)) or (variant[pos]["sample"]["GT"]!= (1,0))) and (variant[pos]['alt'][0]!='.'):
				consensus+=iupacDict[(variant[pos]['ref'][0]+variant[pos]['alt'][0]).upper()]
			else:
				consensus+=variant[pos]['alt'][0]
		elif A[index] >= max(T[index],C[index],G[index],minCov):
			consensus+="A"
		elif T[index] >= max(C[index],G[index],minCov):
			consensus+="T"
		elif C[index] >= max(G[index],minCov):
			consensus+="C"
		elif G[index] >= minCov:
			consensus+="G"
		else:
			consensus+="-"

	seq_consensus=Seq(consensus,ambiguous_dna)

	#TODO Reverse complement if strand minus 
	if strand=='-':
		consensusRecord=SeqRecord(seq_consensus.reverse_complement(),id=args.name,description="")
	else:
		consensusRecord=SeqRecord(seq_consensus,id=args.name,description="")

	SeqIO.write(consensusRecord,'consensus.fa','fasta')
	#TODO create a bait.fa for iterative mapping

	#TODO name the baits and add this delimition to a bait.gff !!!

# Create output directory structure and logs
saveRoot(rootedDir)
#rootedDir.logs.report()