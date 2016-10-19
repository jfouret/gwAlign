#!/usr/bin/python
import argparse
gitRepository='SEDMATCHGITREPO'
version='SEDMATCHGITVERSION'
year=2016
author='Julien Fouret'
contact='julien@fouret.me'

parser = argparse.ArgumentParser(description='Retrive consensus sequence of a given region based on a sorted bam file',epilog="Version : "+str(version)+"\n"+str(year)+"\nAuthor : "+author+" for more informations or enquiries please contact "+contact,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-outDir', metavar='/path', required=True, help="Output directory")
parser.add_argument('-reference', metavar='/path', required=True, help="csv file with id;fasta;bam in the order of priority  \n"+
	"id: the id of the species used in the exon alignment file\n"+
	"fasta: path to the fasta file used for short read alignments\n"+
	"bam: path of the sorted and indexed bam file")
parser.add_argument('-reg', metavar='chr:start-end+;chr:start-end-;...', required=True, help="region of interest, please respect the format and the order of the ref file")
parser.add_argument('-gatk', metavar='/path', required=False, help="gatk jar path",default='/export/bin/source/GenomeAnalysisTK/GenomeAnalysisTK.jar')
parser.add_argument('-picard', metavar='/path', required=False, help="picard jar path",default='/export/bin/picard-tools-2.1.0/picard.jar')
parser.add_argument('-graphCov', action='store_true', help="Graph the coverage")
parser.add_argument('-hg19', metavar='/path', required=False, help="hg19 genome for base mapping",default='/export/data/Genomes/Human/hg19/hg19.fa')
parser.add_argument('-hg19Dict', metavar='/path', required=False, help="hg19 genome for base mapping",default='/export/data/Genomes/Human/hg19/hg19.dict')
args=parser.parse_args()

# Import libraries
import sys
import os
from Bio import SeqIO
from jupype import *

# Create output directory structure and logs
rootedDir=RootDir(args.outDir)
rootedDir.logs.writeArgs(args)

# Definition of used software
samtools=Command('samtools',min='1.3')
samtools.versionCtrl()
samtools.log()
bcftools=Command('bcftools',min='1.3')
bcftools.versionCtrl()
bcftools.log()
java=Command('java')
java.log()
gatk_cmd=java.create(options={'-jar':args.gatk})
gatk=Command(gatk_cmd,gatk_cmd+' -version 2>&1')
gatk.log()
megacc=Command('megacc')
megacc.log()
picard_cmd=java.create(options={'-jar':args.picard})
picard=Command(picard_cmd,picard_cmd+' CheckFingerprint --version 2>&1 | sed \'s/(.*$//g\'')
picard.log()

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
			if region!='':
				refDict[spec]={
					'fasta':fasta,
					'bam':bam,
					'reg':"'"+region+"'"
				}
			specOrder.append(spec)
os.chdir(rootedDir.results)

## step 1 - iteration for all species
NbCallDict=dict()
for spec in refDict.keys():
	## step 1.1 - Creation a vcf-formatted files for mpileup - base calling
	mkdirp(spec)
	#define options and positional args for software
	mpileupOpt={
		'-r':refDict[spec][reg][:-1], # [:-1] without strand information yet ...
		'--reference':efDict[spec][fasta],
		'-v':'',
		'-o':rootedDir.results+'/'+spec+'/call.vcf.gz'
		}
	mpileupPos=[bam]

	submitOneShell(samtools.create(options=mpileupOpt,positionals=mpileupPos,subprogram='mpileup'))

	## step 1.2 - call genotype at vcf format including variant and non-variant

	callOpt={
		'-c':'',
		'-M':'',
		'-O':'v',
		'-o':rootedDir.results+'/'+spec+'/genotype.vcf'
		}
	callPos=[rootedDir.results+'/'+spec+'/call.vcf.gz']

	submitOneShell(bcftools.create(options=callOpt,positionals=callPos,subprogram='call'))

	submitOneShell("rm "+rootedDir.results+'/'+spec+'/genotype.vcf.idx')

	## step 1.3 - check the spe
	NbCallDict[spec]=int(submitOneShell("grep -c -v -P '^#' "+rootedDir.results+'/'+spec+'/genotype.vcf')['out'].rstrip())

## step 2 - choose the to create the consensus based on the species with the most positions called

chosenSpec=None
for spec in NbCallDict.keys():
	if chosenSpec==None:
		chosenSpec=spec
		maxNbCall=NbCallDict[spec]
	elif NbCallDict[spec]>maxNbCall:
		chosenSpec=spec
		maxNbCall=NbCallDict[spec]
	elif NbCallDict[spec]==maxNbCall:
		if specOrder.index(chosenSpec)>specOrder.index(spec):
			chosenSpec=spec
			maxNbCall=NbCallDict[spec]

NbCalls=NbCallDict[chosenSpec]

## step 3.1 - create consensus
reg=refDict[chosenSpec]['reg']
if NbCalls==0:
	import re
	with open(rootedDir.results+'/hg19mapped.fa','w') as csFile:
		csFile.write('>1 '+reg.split('-')[0]+"\n")
		pos=reg[:-1].split(':')[1].split('-')
		pos=[int(pos[0]),int(pos[1])]
		seqLen=max(pos)-min(pos)+1
		csFile.write(re.sub("(.{60})", "\\1\n", '-'*seqLen, 0, re.DOTALL)+"\n")
else:
	gatkOpt={
		'-T':'FastaAlternateReferenceMaker',
		'-R':ref,
		'-o':rootedDir.results+'/'+chosenSpec+'/tmp_consensus.fa',
		'-L':reg[:-1],
		'-V':rootedDir.results+'/'+chosenSpec+'/genotype.vcf',
		'-U':'ALLOW_SEQ_DICT_INCOMPATIBILITY'
	}
	submitOneShell(gatk.create(options=gatkOpt))
	## step 3.2.0 Get the hg19 sequence
	getRefOpt={
		'R':args.hg19,
		'O':rootedDir.results+'/hg19.fa',
		'INTERVAL_LIST':rootedDir.results+'/regions.txt',
	}
	
	submitOneShell('cp '+args.hg19Dict+' '+rootedDir.results+'/regions.txt')
	chromosome,positions=refDict['hg19']['reg'].split(':')
	start,end=positions[:-1].split('-')
	strand=positions[-1]
	#start=str(int(start)+1)
	with open(rootedDir.results+'/regions.txt','a') as regFile:
		regFile.write("\t".join([chromosome,start,end,strand,'.'])+"\n")
	submitOneShell(picard.create(options=getRefOpt,subprogram='ExtractSequences',sep='='))
	if strand=='-':
		#reverse complement
		fasta_sequences=SeqIO.parse(open(rootedDir.results+'/'+chosenSpec+'/tmp_consensus.fa'),'fasta')
		for fasta in fasta_sequences:
			name, sequence = fasta.description, str(fasta.seq.reverse_complement())
			with open(rootedDir.results+'/'+chosenSpec+'/consensus.fa') as consFile:
				consFile.write('>'+name+"\n"+sequence+"\n")
		os.remove(rootedDir.results+'/'+chosenSpec+'/tmp_consensus.fa')
	else:
		# or just rename if strand +
		os.rename(rootedDir.results+'/'+chosenSpec+'/tmp_consensus.fa',rootedDir.results+'/'+chosenSpec+'/consensus.fa')

	submitOneShell('cat '+rootedDir.results+'/hg19.fa '+rootedDir.results+'/'+chosenSpec+'/consensus.fa > '+rootedDir.results+'/'+chosenSpec+'/data.fa')

	## step 3.2.1 Build an alignment consensus-hg19 in aln.fasta
	megaccOpt={
		'-a':gitRepository+'/template/clustal_align_nucleotide.mao',
		'-d':rootedDir.results+'/'+chosenSpec+'/data.fa',
		'-o':rootedDir.results+'/'+chosenSpec+'/aln',
		'-f':'Fasta'
	}
	submitOneShell(megacc.create(options=megaccOpt))

	## step 3.2.2 parse the alignment to prepare to hg19 mapped exon (no gap in hg19)
	seqList=list()
	fasta_sequences=SeqIO.parse(open(rootedDir.results+'/'+chosenSpec+'/aln.fasta'),'fasta')
	for fasta in fasta_sequences:
		name, sequence = fasta.description, str(fasta.seq)
		seqList.append(sequence)
	seqhg19=list(seqList[0])
	seqConsensus=list(seqList[1])
	mappedSeq=''
	for index in range(len(seqhg19)):
		if seqhg19[index]!='-':
			mappedSeq+=seqConsensus[index]
	with open(rootedDir.results+'/hg19mapped.fa','w') as mappedFile:
		mappedFile.write('>hg19mapped '+reg+"\n")
		mappedFile.write(mappedSeq+"\n")

	## step 4.0 - Create coverage.tab file for graph outputs TODO ==> put in a seperate programe just with outDir

	if args.graphCov:
		submitOneShell('echo "Position'+"\t"+'Reference'+"\t"+'Alternative'+"\t"+'Depth'+"\t"+'Genotype"'+' > coverage.tab ; grep -P "^#" -v calling.vcf |grep -e \'INDEL\' -v | awk -F "[\t ;:]" \'{print $2 "\t" $4 "\t" $5 "\t" $8 "\t" $(NF-1) }\' | sed \'s/DP=//g\' >> coverage.tab')
		import pandas as pd
		#import matplotlib
		import plotly.offline as py
		import plotly.graph_objs as go
		#py.init_notebook_mode()
		table=pd.DataFrame.from_csv('coverage.tab',sep="\t",index_col=False)
		tableDict=dict()
		tableDict['Reference Homozygote']=table[(table.Alternative=='.') & (table.Genotype=='0/0')]
		tableDict['Alternative Homozygote']=table[(table.Alternative!='.') & (table.Genotype=='1/1')]
		tableDict['Reference Heterozygote']=table[(table.Alternative=='.') & (table.Genotype!='0/0')]
		tableDict['Alternative Heterozygote']=table[(table.Alternative!='.') & (table.Genotype!='1/1')]
		data=[]
		for key in tableDict.keys():
			if tableDict[key].size!=0:
				trace=go.Bar(
					x=tableDict[key]['Position'],
					y=tableDict[key]['Depth'],
					text="\n"+'Reference: '+tableDict[key].Reference+"\n"+'Alternative: '+tableDict[key].Alternative,
					name=key
				)
				data.append(trace)
				layout=dict(
					title=reg,
					xaxis=dict(
						title='Genomic positions (bases)'
						),
					yaxis=dict(
						title='Genomic depth (bases)'
						)
					)
				fig=go.Figure(data=data,layout=layout)
		py.plot(fig,filename=rootedDir.reports+"/coverage.html",show_link=False,auto_open=False)
		#plot=data.plot.area(title=reg,legend=True)
		#fig=plot.get_figure()
		#fig.savefig(rootedDir.reports+"/coverage.png")
	#Nbcalls=int(submitOneShell('rm '+rootedDir.results+'/tmp.bcf ; grep -P "^#" -v -c calling.vcf')['out'])
	#if Nbcalls>0:
	#	submitOneShell('vcfutils.pl vcf2fq calling.vcf | fastq2fasta.pl > consensus.fa')

# Create output directory structure and logs
saveRoot(rootedDir)
#rootedDir.logs.report()