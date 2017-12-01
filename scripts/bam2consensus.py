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
platypus=Command("python /export/source/git/Platypus/bin/Platypus.py",'echo "(see program log)"')
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
	start=pos[0]-1 # 0-based inclusive
	end=pos[1]     # 0-based exclusive
	tcoverage=[]
	for pos in range(start,end+1):
		al=samfile[spec].fetch(contig=contig,start=start,stop=end)
		cov=0
		for read in al:
			cov+=1
		tcoverage.append(cov)
	tYield=sum(tcoverage)
	#print(spec)#debug
	if chosenSpec==None:
		chosenSpec=spec
		Yield=tYield
		coverage=tcoverage
	elif tYield>Yield:
		chosenSpec=spec
		Yield=tYield
		coverage=tcoverage
	elif tYield==Yield:
		if specOrder.index(chosenSpec)>specOrder.index(spec):
			chosenSpec=spec
			Yield=tYield
			coverage=tcoverage
for spec in refDict.keys():
	if spec!=chosenSpec:
		samfile[spec].close()

reg=refDict[chosenSpec]['reg']
contig,pos=reg[:-1].split(':')
pos=pos.split('-')
pos=[int(pos[0]),int(pos[1])]
start=pos[0] # 1-based inclusive
end=pos[1]     # 1-based inclusive



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
		'--output':rootedDir.results+'/'+chosenSpec+'/call.vcf',
		'--logFileName':rootedDir.results+'/'+chosenSpec+'/call.log',
		'--filterReadPairsWithSmallInserts':'0',
		'--filterReadsWithDistantMates':'0',
		'--filterReadsWithUnmappedMates':'0',
		'--coverageSamplingLevel':'50',
		'--assemblyRegionSize':'800',
		'--assemblerKmerSize':'115',
		'--assembleBrokenPairs':'1',
		'--assembleBadReads':'1',
		'--assembleAll':'1',
		'--filterDuplicates':'0',
		'--maxVariants':"1",
		'--filterVarsByCoverage':'1',
		'--trimSoftClipped':'0',
		'--trimOverlapping':'0',
		'--trimAdapter':'0',
		'--minReads':'1',
		'--trimReadFlank':'0',
		#'--maxVarDist':'1',
		'--mergeClusteredVariants':'0',
		'--minFlank':'0',
		'--minMapQual':'0',
		'--minBaseQual':'0',
		'--minGoodQualBases':'0',
		'--assemble':'1'
		}

	submitOneShell(platypus.create(options=PlatypusOpt,subprogram='callVariants'))

	#read the vcf file

	vcf=pysam.VariantFile(rootedDir.results+'/'+chosenSpec+'/call.vcf')

	variant=dict()
	for rec in vcf.fetch():
		variant[rec.pos]=dict()
		variant[rec.pos]['ref']=rec.ref
		variant[rec.pos]['alt']=rec.alts
		if not isinstance(variant[rec.pos]['alt'], basestring):
			variant[rec.pos]['alt']=variant[rec.pos]['alt'][0]
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
	refSequence=fastaFile.fetch(region=contig+":"+str(start)+"-"+str(end))
	fastaFile.close()
	start_cons=min(start,200)
	end_cons=min(start,200)+end-start
	contig_cons=args.name
	# code to get consensus
	consensus=""
	pos=start
	while pos <= end:
		index=pos-start
		#print("\n\n############START###########\n"+str(pos))
		# Have a variant been called with sufficient coverage ?
		if (pos in variant.keys()) and (variant[pos]["sample"]["NR"]>=minCovVar):
			#if unique SNP with no variants
			if (len(variant[pos]['alt']+variant[pos]['ref'])==2) and ((variant[pos]["sample"]["GT"]!= (0,1)) or (variant[pos]["sample"]["GT"]!= (1,0))) and (variant[pos]['alt'].upper() in ['A','T','C','G']) and (variant[pos]['ref'].upper() in ['A','T','C','G']) :
				consensus+=iupacDict[(variant[pos]['ref'].upper()+variant[pos]['alt']).upper()]
				pos+=1
			elif not ((variant[pos]["sample"]["GT"]!= (0,1)) or (variant[pos]["sample"]["GT"]!= (1,0))):
				consensus+=variant[pos]['alt']
				pos+=len(variant[pos]['ref'])
			else:
				consensus+=variant[pos]['ref']
				pos+=len(variant[pos]['ref'])
		elif coverage[index]>minCov:
			consensus+=refSequence[index]
			pos+=1
		else : 
			consensus+="-"
			pos+=1

	seq_consensus=Seq(consensus,ambiguous_dna)
	bait_reg=contig_cons+":"+str(start_cons)+"-"+str(end_cons)+strand
	bait=left.lower()+consensus.upper()+right.lower()
	with open("bait.fa",'w') as baitFile:
		baitFile.write(">"+args.name+" "+bait_reg+"\n"+bait+"\n")
	#TODO Reverse complement if strand minus 
	if strand=='-':
		consensusRecord=SeqRecord(seq_consensus.reverse_complement(),id=args.name,description="")
	else:
		consensusRecord=SeqRecord(seq_consensus,id=args.name,description="")

	SeqIO.write(consensusRecord,'consensus.fa','fasta')

# Create output directory structure and logs
saveRoot(rootedDir)
#rootedDir.logs.report()