#!/usr/bin/python
import argparse
gitRepository='SEDMATCHGITREPO'
version='SEDMATCHGITVERSION'
year=2016
author='Julien Fouret'
contact='julien@fouret.me'

parser = argparse.ArgumentParser(description='Retrive consensus sequence of a given region based on a sorted bam file',epilog="Version : "+str(version)+"\n"+str(year)+"\nAuthor : "+author+" for more informations or enquiries please contact "+contact,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-outDir', metavar='/path', required=True, help="Output directory")
parser.add_argument('-bam', metavar='/path', required=True, help="A sorted bam file")
parser.add_argument('-ref', metavar='/path', required=True, help="the reference associated to this bam file")
parser.add_argument('-reg', metavar='chr:start-end', required=True, help="region of interest, please respect the format")
parser.add_argument('-gatk', metavar='/path', required=False, help="gatk jar path",default='/export/bin/source/GenomeAnalysisTK/GenomeAnalysisTK.jar')
parser.add_argument('-picard', metavar='/path', required=False, help="picard jar path",default='/export/bin/picard-tools-2.1.0/picard.jar')
parser.add_argument('-graphCov', action='store_true', help="Graph the coverage")
parser.add_argument('-strict', action='store_true', help="do not consid")
parser.add_argument('-hg19', metavar='/path', required=False, help="hg19 genome for base mapping",default='/export/data/Genomes/Human/hg19/hg19.fa')
parser.add_argument('-hg19Dict', metavar='/path', required=False, help="hg19 genome for base mapping",default='/export/data/Genomes/Human/hg19/hg19.dict')
args=parser.parse_args()

# Import libraries
import sys
import os
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
gatk=Command(gatk_cmd,gatk_cmd+' -version 2>&1',min='3.6')
gatk.versionCtrl()
gatk.log()
megacc=Command('megacc')
picard_cmd=java.create(options={'-jar':args.picard})
picard=Command(picard_cmd,picard_cmd+' CheckFingerprint --version 2>&1 | sed \'s/(.*$//g\'')
picard.log()

#get absolute path from inputs
bam=os.path.abspath(args.bam)
ref=os.path.abspath(args.ref)

os.chdir(rootedDir.results)
## step 1 - Creation a vcf-formatted files for mpileup

#define options and positional args for software
mpileupOpt={
	'-r':args.reg,
	'--reference':ref,
	'-v':'',
	'-o':rootedDir.results+'/tmp.vcf.gz'
	}
mpileupPos=[bam]

submitOneShell(samtools.create(options=mpileupOpt,positionals=mpileupPos,subprogram='mpileup'))

## step 2 - call consensus at vcf format including variant and non-variant

callOpt={
	'-c':'',
	'-M':'',
	'-O':'v',
	'-o':rootedDir.results+'/tmp_calling.vcf'
	}
callPos=[rootedDir.results+'/tmp.vcf.gz']

submitOneShell(bcftools.create(options=callOpt,positionals=callPos,subprogram='call'))
submitOneShell("grep -v -e \"\tINDEL\t\" "+rootedDir.results+'/tmp_calling.vcf > '+rootedDir.results+'/calling.vcf')
submitOneShell("rm "+rootedDir.results+'/tmp_calling.vcf '+rootedDir.results+'/calling.vcf.idx')

## step 3.0 - check for the presence of indexation files

prefixRef=ref.replace('.fa','')
prefixRef=prefixRef.replace('.fna','')
prefixRef=prefixRef.replace('.fasta','')

faiRef=ref+'.fai'
dictRef=prefixRef+'.dict'

if not os.path.isfile(dictRef):

	dictOpt={
		'REFERENCE':ref,
		'OUTPUT':dictRef
	}
	submitOneShell(picard.create(options=dictOpt,sep='=',subprogram='CreateSequenceDictionary'))

if not os.path.isfile(faiRef):
	faiPos=[ref]
	submitOneShell(samtools.create(positionals=faiPos,subprogram='faidx'))

## step 3.0 - check variant calling
NbCalls=submitOneShell("grep -c -v -P '^#' "+rootedDir.results+'/calling.vcf')['out'].rstrip()

## step 3.1 - create consensus
if NbCalls=='0':
	import re
	with open(rootedDir.results+'/hg19mapped.fa','w') as csFile:
		csFile.write('>1 '+args.reg.split('-')[0]+"\n")
		pos=args.reg.split(':')[1].split('-')
		pos=[int(pos[0]),int(pos[1])]
		seqLen=max(pos)-min(pos)+1
		csFile.write(re.sub("(.{60})", "\\1\n", '-'*seqLen, 0, re.DOTALL)+"\n")
else:
	gatkOpt={
		'-T':'FastaAlternateReferenceMaker',
		'-R':ref,
		'-o':rootedDir.results+'/consensus.fa',
		'-L':args.reg,
		'-V':rootedDir.results+'/calling.vcf',
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
	chromosome,positions=args.reg.split(':')
	start,end=positions.split('-')
	#start=str(int(start)+1)
	with open(rootedDir.results+'/regions.txt','a') as regFile:
		regFile.write("\t".join([chromosome,start,end,'+','.'])+"\n")
	submitOneShell(picard.create(options=getRefOpt,subprogram='ExtractSequences',sep='='))
	submitOneShell('cat '+rootedDir.results+'/hg19.fa '+rootedDir.results+'/consensus.fa > '+rootedDir.results+'/data.fa')

	## step 3.2.1 Build an alignment consensus-hg19 in aln.fasta
	megaccOpt={
		'-a':gitRepository+'/template/clustal_align_nucleotide.mao',
		'-d':rootedDir.results+'/data.fa',
		'-o':rootedDir.results+'/aln',
		'-f':'Fasta'
	}
	submitOneShell(megacc.create(options=megaccOpt))

	## step 3.2.2 parse the alignment to prepare to hg19 mapped exon (no gap in hg19)
	from Bio import SeqIO
	seqList=list()
	fasta_sequences=SeqIO.parse(open(rootedDir.results+'/aln.fasta'),'fasta')
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
		mappedFile.write('>hg19mapped '+args.reg+"\n")
		mappedFile.write(mappedSeq+"\n")

	## step 4.0 - Create coverage.tab file for graph outputs

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
					title=args.reg,
					xaxis=dict(
						title='Genomic positions (bases)'
						),
					yaxis=dict(
						title='Genomic depth (bases)'
						)
					)
				fig=go.Figure(data=data,layout=layout)
		py.plot(fig,filename=rootedDir.reports+"/coverage.html",show_link=False,auto_open=False)
		#plot=data.plot.area(title=args.reg,legend=True)
		#fig=plot.get_figure()
		#fig.savefig(rootedDir.reports+"/coverage.png")
	#Nbcalls=int(submitOneShell('rm '+rootedDir.results+'/tmp.bcf ; grep -P "^#" -v -c calling.vcf')['out'])
	#if Nbcalls>0:
	#	submitOneShell('vcfutils.pl vcf2fq calling.vcf | fastq2fasta.pl > consensus.fa')

# Create output directory structure and logs
saveRoot(rootedDir)
#rootedDir.logs.report()