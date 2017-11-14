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
#parser.add_argument('-minCov', metavar='N', required=True, help="min coverage for base calling")
#parser.add_argument('-pval', metavar='N', required=True, help="pval to call the alternative base (H0: reference)")
parser.add_argument('-gatk', metavar='/path', required=False, help="gatk jar path",default='SEDMATCHGATK')
#parser.add_argument('-picard', metavar='/path', required=False, help="picard jar path",default='SEDMATCHPICARD')
#parser.add_argument('-graphCov', action='store_true', help="Graph the coverage")
#parser.add_argument('-hg19', metavar='/path', required=False, help="hg19 genome for base mapping",default='/export/data/Genomes/Human/hg19_ucsc/hg19.fa')
#parser.add_argument('-hg19Dict', metavar='/path', required=False, help="hg19 genome for base mapping",default='/export/data/Genomes/Human/hg19_ucsc/hg19.dict')

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
vcfconsensus=Command('vcf-consensus','vcftools | grep \'VCFtools (\' | sed \'s/VCFtools (v//g\' | sed \'s/)//g\'')
vcfconsensus.log()
java=Command('java')
java.log()
gatk_cmd=java.create(options={'-jar':args.gatk})
gatk=Command(gatk_cmd,gatk_cmd+' -version 2>&1')
gatk.log()
#megacc=Command('megacc')
#megacc.log()
#picard_cmd=java.create(options={'-jar':args.picard})
#picard=Command(picard_cmd,picard_cmd+' CheckFingerprint --version 2>&1 | sed \'s/(.*$//g\'')
#picard.log()
bedtools=Command('bedtools',"bedtools --version | sed -r 's/bedtools v(.*)/\1/g'")
bedtools.log()

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
			if region!='':
				refDict[spec]={
					'fasta':fasta,
					'bam':bam,
					'reg':region
				}
				specOrder.append(spec)
os.chdir(rootedDir.results)

#print(refDict)

## step 1 - iteration for all species
NbCallDict=dict()
for spec in refDict.keys():
	## step 1.1 - Creation a vcf-formatted files with GATK haplotype caller - base calling
	#print("###\n"+str(refDict[spec])+"\n###\n")
	mkdirp(spec)
	#define options and positional args for software, considering PCR free # TODO put this in argument by default
	HaploCallOpt={
		'-L':refDict[spec]['reg'][:-1], # [:-1] without strand information yet ...
		'-T':'HaplotypeCaller',
		'-R':refDict[spec]['fasta'],
		'-I':refDict[spec]['bam'],
		'--pcr_indel_model':'NONE',
		'-o':rootedDir.results+'/'+spec+'/call.vcf'
		}
	mpileupPos=[refDict[spec]['bam']]

	submitOneShell(gatk.create(options=HaploCallOpt))

	## step 1.2 - call genotype at vcf format including variant and non-variant

	#filterOpt={
	#	'-e':"'DP<"+args.minCov+"'"
	#}
	#filterPos=[rootedDir.results+'/'+spec+'/call.vcf.gz']

	#callOpt={
	#	'-O':'v',
	#	'-o':rootedDir.results+'/'+spec+'/genotype.vcf',
	#	'-m':'',
	#	'-M':'',
	#	'--pval-threshold':args.pval
	#}
	#cmdList=list()

	#cmdList.append(bcftools.create(options=filterOpt,positionals=filterPos,subprogram='filter'))
	#cmdList.append(bcftools.create(options=callOpt,subprogram='call'))
	
	#submitShell(cmdList,sep=' | ')

	#submitOneShell("rm "+rootedDir.results+'/'+spec+'/genotype.vcf.idx')

	## step 1.3 - check the spe
	NbCallDict[spec]=submitOneShell("grep -c -v -P '^#' "+rootedDir.results+'/'+spec+'/genotype.vcf')['out'].rstrip()
	submitOneShell("bgzip "+rootedDir.results+'/'+spec+'/call.vcf'+"\n")
	if NbCallDict[spec]=='':
		NbCallDict[spec]=0
	else:
		NbCallDict[spec]=int(NbCallDict[spec])
## step 2 - choose the to create the consensus based on the species with the most positions called

chosenSpec=None
for spec in NbCallDict.keys():
	#print(spec)#debug
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

## step 2.1 - create empty consensus if no calls

if NbCalls==0:
	import re
	reg=refDict["hg19"]['reg']
	chromosome,pos=reg[:-1].split(':')
	pos=pos.split('-')
	pos=[int(pos[0]),int(pos[1])]
	start=pos[0]
	end=pos[1]
	seqLen=max(pos)-min(pos)+1
	with open(rootedDir.results+'/consensus.fa','w') as csFile:
		csFile.write('>1 '+reg.split('-')[0]+"\n")
		csFile.write(re.sub("(.{60})", "\\1\n", '-'*seqLen, 0, re.DOTALL)+"\n")
else:
	spec=chosenSpec
	reg=refDict[chosenSpec]['reg']
	chromosome,pos=reg[:-1].split(':')
	pos=pos.split('-')
	pos=[int(pos[0]),int(pos[1])]
	start=pos[0]
	end=pos[1]
	seqLen=max(pos)-min(pos)+1
	#gatkOpt={
	#	'-T':'FastaAlternateReferenceMaker',
	#	'-R':refDict[chosenSpec]['fasta'],
	#	'-o':rootedDir.results+'/'+chosenSpec+'/tmp_consensus.fa',
	#	'-L':reg[:-1],
	#	'-V':rootedDir.results+'/'+chosenSpec+'/genotype.vcf',
	#	'-U':'ALLOW_SEQ_DICT_INCOMPATIBILITY'
	#}
	#submitOneShell(gatk.create(options=gatkOpt))

	## step 2.2 - create tabix and cov.bed and region.bed and region.bed

	submitOneShell("tabix "+rootedDir.results+'/'+spec+'/call.vcf.gz')
	cmdList=list()
	cmdList.append('zmore '+rootedDir.results+'/'+spec+'/call.vcf.gz')
	cmdList.append('grep -e \'^#\' -v')
	cmdList.append('awk \'{print $1 "\t" ($2 -1) "\t" $2}\'')
	cmdList.append(bedtools.create(subprogram='merge',positionals=[' > '+rootedDir.results+'/'+spec+'/cov.bed']))
	submitShell(cmdList,sep=' | ')
	submitOneShell('echo "'+'\t'.join([chromosome,str(start-1),str(end)])+'" > ' + rootedDir.results+'/'+spec+'/region.bed')
	cmdList=list()
	cmdList.append('zmore '+rootedDir.results+'/'+spec+'/call.vcf.gz')
	cmdList.append('grep \'##contig=<ID='+chromosome+',length=\'')
	cmdList.append('sed -r "s/^##contig=<ID=(.*),length=(.*)>/\\1\t\\2/g" > '+rootedDir.results+'/'+spec+'/region.genome')
	#print('##GENOME LENGTH')
	#print("\n")
	#print(cmdList)
	#print("\n")
	#print('##GENOME LENGTH')
	#print("\n")
 	submitShell(cmdList,sep=' | ')

	## step 2.2 - get the mask.bed

	outRegOpt={
		'-i':rootedDir.results+'/'+spec+'/region.bed',
		'-g':rootedDir.results+'/'+spec+'/region.genome'
	}
	submitOneShell(bedtools.create(subprogram='complement',options=outRegOpt,positionals=[' > '+rootedDir.results+'/'+spec+'/out_region.bed']))

	cmdList=list()
	cmdList.append('cat '+rootedDir.results+'/'+spec+'/out_region.bed '+rootedDir.results+'/'+spec+'/cov.bed')
	cmdList.append(bedtools.create(subprogram='sort'))
	cmdList.append(bedtools.create(subprogram='merge',positionals=[' > '+rootedDir.results+'/'+spec+'/out_and_cov.bed']))
	submitShell(cmdList,sep=' | ')

	outCovOption={
		'-i':rootedDir.results+'/'+spec+'/out_and_cov.bed',
		'-g':rootedDir.results+'/'+spec+'/region.genome'
	}
	submitOneShell(bedtools.create(subprogram='complement',options=outCovOption,positionals=[' > '+rootedDir.results+'/'+spec+'/mask.bed']))

	## step 2.2 - create consensus with mask

	consensusOpt={
		'--mask':rootedDir.results+'/'+spec+'/mask.bed'
	}
	consensusPos=[rootedDir.results+'/'+spec+'/call.vcf.gz','> '+rootedDir.results+'/'+spec+'/tmp_consensus.fa']

	cmdList=list()
	cmdList.append(samtools.create(subprogram='faidx',positionals=[refDict[spec]['fasta'],reg[:-1]]))
	cmdList.append(bcftools.create(options=consensusOpt,positionals=consensusPos,subprogram='consensus'))

	submitShell(cmdList, sep=' | ')

	## step 3.2.0 Get the hg19 sequence
#	getRefOpt={
#		'R':args.hg19,
#		'O':rootedDir.results+'/hg19.fa',
#		'INTERVAL_LIST':rootedDir.results+'/regions.txt',
#	}

	#submitOneShell('cp '+args.hg19Dict+' '+rootedDir.results+'/regions.txt')
	chromosome,positions=refDict[chosenSpec]['reg'].split(':')
	start,end=positions[:-1].split('-')
	strand=positions[-1]
	start=str(int(start)+1)
	#with open(rootedDir.results+'/regions.txt','a') as regFile:
		#regFile.write("\t".join([chromosome,start,end,strand,'.'])+"\n")
	#submitOneShell(picard.create(options=getRefOpt,subprogram='ExtractSequences',sep='='))
	fasta_sequences=SeqIO.parse(open(rootedDir.results+'/'+chosenSpec+'/tmp_consensus.fa'),'fasta')
	for fasta in fasta_sequences:
		if strand=='-': #reverse complement
			name, sequence = fasta.description, str(fasta.seq.reverse_complement())
		else:
			name, sequence = fasta.description, str(fasta.seq)
	sequence=sequence.replace('N','-').replace('n','-')
	with open(rootedDir.results+'/consensus.fa','w') as consFile:
		consFile.write('>'+name+"\n"+sequence+"\n")
	os.remove(rootedDir.results+'/'+chosenSpec+'/tmp_consensus.fa')
		
	#print('cat '+rootedDir.results+'/hg19.fa '+rootedDir.results+'/'+chosenSpec+'/consensus.fa > '+rootedDir.results+'/'+chosenSpec+'/data.fa')
	#submitOneShell('cat '+rootedDir.results+'/hg19.fa '+rootedDir.results+'/'+chosenSpec+'/consensus.fa > '+rootedDir.results+'/'+chosenSpec+'/data.fa')

	## step 3.2.1 Build an alignment consensus-hg19 in aln.fasta
	#megaccOpt={
#		'-a':gitRepository+'/template/clustal_align_nucleotide.mao',
#		'-d':rootedDir.results+'/'+chosenSpec+'/data.fa',
#		'-o':rootedDir.results+'/'+chosenSpec+'/aln',
#		'-f':'Fasta'
#	}
#	submitOneShell(megacc.create(options=megaccOpt))

	## step 3.2.2 parse the alignment to prepare to hg19 mapped exon (no gap in hg19)
#	seqList=list()
#	fasta_sequences=SeqIO.parse(open(rootedDir.results+'/'+chosenSpec+'/aln.fasta'),'fasta')
#	for fasta in fasta_sequences:
#		name, sequence = fasta.description, str(fasta.seq)
#		seqList.append(sequence)
#	seqhg19=list(seqList[0])
#	seqConsensus=list(seqList[1])
#	mappedSeq=''
#	for index in range(len(seqhg19)):
#		if seqhg19[index]!='-':
#			mappedSeq+=seqConsensus[index]
#	with open(rootedDir.results+'/hg19mapped.fa','w') as mappedFile:
#		mappedFile.write('>hg19mapped '+reg+"\n")
#		mappedFile.write(mappedSeq+"\n")

	## step 4.0 - Create coverage.tab file for graph outputs TODO ==> put in a seperate programe just with outDir

#	if args.graphCov:
#		submitOneShell('echo "Position'+"\t"+'Reference'+"\t"+'Alternative'+"\t"+'Depth'+"\t"+'Genotype"'+' > coverage.tab ; grep -P "^#" -v calling.vcf |grep -e \'INDEL\' -v | awk -F "[\t ;:]" \'{print $2 "\t" $4 "\t" $5 "\t" $8 "\t" $(NF-1) }\' | sed \'s/DP=//g\' >> coverage.tab')
#		import pandas as pd
#		#import matplotlib
#		import plotly.offline as py
#		import plotly.graph_objs as go
#		#py.init_notebook_mode()
#		table=pd.DataFrame.from_csv('coverage.tab',sep="\t",index_col=False)
#		tableDict=dict()
#		tableDict['Reference Homozygote']=table[(table.Alternative=='.') & (table.Genotype=='0/0')]
#		tableDict['Alternative Homozygote']=table[(table.Alternative!='.') & (table.Genotype=='1/1')]
#		tableDict['Reference Heterozygote']=table[(table.Alternative=='.') & (table.Genotype!='0/0')]
#		tableDict['Alternative Heterozygote']=table[(table.Alternative!='.') & (table.Genotype!='1/1')]
#		data=[]
#		for key in tableDict.keys():
#			if tableDict[key].size!=0:
#				trace=go.Bar(
#					x=tableDict[key]['Position'],
#					y=tableDict[key]['Depth'],
#					text="\n"+'Reference: '+tableDict[key].Reference+"\n"+'Alternative: '+tableDict[key].Alternative,
#					name=key
#				)
#				data.append(trace)
#				layout=dict(
#					title=reg,
#					xaxis=dict(
#						title='Genomic positions (bases)'
#						),
#					yaxis=dict(
#						title='Genomic depth (bases)'
#						)
#					)
#				fig=go.Figure(data=data,layout=layout)
#		py.plot(fig,filename=rootedDir.reports+"/coverage.html",show_link=False,auto_open=False)
		#plot=data.plot.area(title=reg,legend=True)
		#fig=plot.get_figure()
		#fig.savefig(rootedDir.reports+"/coverage.png")
	#Nbcalls=int(submitOneShell('rm '+rootedDir.results+'/tmp.bcf ; grep -P "^#" -v -c calling.vcf')['out'])
	#if Nbcalls>0:
	#	submitOneShell('vcfutils.pl vcf2fq calling.vcf | fastq2fasta.pl > consensus.fa')

# Create output directory structure and logs
saveRoot(rootedDir)
#rootedDir.logs.report()