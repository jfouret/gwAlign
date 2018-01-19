#!/usr/bin/env python
import argparse
gitRepository='SEDMATCHGITREPO'
version='SEDMATCHGITVERSION'
year=2016
author='Julien Fouret'
contact='julien.fouret12@uniagro.fr'
##parse argument
parser = argparse.ArgumentParser(description='Compute a codon-based alignment with all cds',epilog="Version : "+str(version)+"\n"+str(year)+"\nAuthor : "+author+" for more informations or enquiries please contact "+contact,formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-gene_dir', metavar='/path', required=True, help="Folder correspoding to the gene ID with all species subdirectories")
parser.add_argument('-ref', metavar='ref_species' , required=False, help="name of the reference specie",default='hg19')
parser.add_argument('-macse', metavar='.jar' , required=False, help="jar file path for macse program",default='SEDMATCHMACSE')
parser.add_argument('-virt_mem', metavar='N', required=True,help="allow N Go of memory")
parser.add_argument('-boost_mem', metavar='N',required=True,help="allow N Go of memory within a new job")

args=parser.parse_args()

import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
import re
import glob
from upype import *

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

speciesList=os.walk(args.gene_dir).next()[1]
os.chdir(args.gene_dir)
speciesExons=dict()
reName=re.compile('exon_(\d+)_ext(\d+)')
reNameAlt=re.compile('exon_(\d+)')
refExLen=dict()

###########################################
###### Description of gene directory ######
#
#./[speciesDirectories]/ 
#
#	[speciesN]/
#		cds.fa contains all cds with the following \
#pattern for seq id: exon_(\d+)_ext(\d+)  OR exon_(\d+)    \
#Indeed an initial exon in human can be disrupted in       \
#several parts in another, either because of re-arrangment \
#either resulting from a mis-assembly
#		missingExon[Number]/cds.fa contains the    \
#same pattern for id but a empty sequence '---'. No        \
#coordinates were given in the input file
#		rearrangedExon[number]/cds.fa contains the \
#same pattern for id but the regions coordinates was not   \
#consistent (end<start) so the fasta contains sequence from\
# the input file
#
#	ref.fa To be created the CDS from the reference    \
#sequence.
#	aligned.fa To be created the CDSs from all other   \
#species that have been aligned to the reference either at \
#genomic level or at reads sequencing level.
#

kgID=args.gene_dir.rstrip("/").split("/").pop()

for species in speciesList:
	speciesExons[species]=dict() # Dict to with exon number as key 
	cdsFileList=glob.glob(species+'/cds.fa')+glob.glob(species+'/*/cds.fa')
	for cdsFile in cdsFileList:
		fasta_sequences = SeqIO.parse(open(cdsFile),'fasta')
		for fasta in fasta_sequences:
			name, sequence = fasta.description, str(fasta.seq)
			m=reName.match(name)
			if m:
				exNum=int(m.group(1))
				partNum=int(m.group(2))
				if not exNum in speciesExons[species].keys():
					speciesExons[species][exNum]=dict()
				speciesExons[species][exNum][partNum]=sequence
				if species==args.ref:
					refExLen[exNum]=len(sequence)
			else:
				m=reNameAlt.match(name)
				if m:
					exNum=int(m.group(1))
					partNum=1
					if not exNum in speciesExons[species].keys():
						speciesExons[species][exNum]=dict()
					speciesExons[species][exNum][partNum]=sequence
					if species==args.ref:
						refExLen[exNum]=len(sequence)


speciesCDS=dict()
for species in speciesExons.keys():
	for i in range(1,1+len(speciesExons[species])):
		if i==1:
			speciesCDS[species]=''
		for j in range(1,1+len(speciesExons[species][i])):
			specSeq=speciesExons[species][i][j]
			if len(specSeq)<10*refExLen[i]:
				speciesCDS[species]+=speciesExons[species][i][j]
			else:
				# Warning if the exon is 10 times more long ! 
				sys.stderr.write("Warning 'Removing exon "+str(i)+" 10 time too long for "+species+"' at : "+args.gene_dir+"\n")


with open('ref.fa','w') as refAln:
	with open('aligned.fa','w') as mapAln:
		for species in speciesCDS.keys():
			if species==args.ref:
				refAln.write("\n".join(['>'+species,speciesCDS[species]])+"\n")
			else:
				mapAln.write("\n".join(['>'+species,speciesCDS[species]])+"\n")

with open('exons.pos','w') as exPos:
	for ExNum in range(1,1+len(refExLen)):
		repeat=str(ExNum)+"\n"
		exPos.write(repeat*refExLen[ExNum])

command='java -Xmx'+args.virt_mem+'g -jar '+args.macse+' -prog alignSequences -seq ref.fa -seq_lr aligned.fa -stop 5000 -stop_lr 10000 -fs '+str(100*len(speciesList))+" -fs_lr 10"

alnProc=submit(command)
#TODO relaunch automatically with boosted boost_mem ! 
if alnProc['err']!='':
	command='java -Xmx'+args.boost_mem+'g -jar '+args.macse+' -prog alignSequences -seq ref.fa -seq_lr aligned.fa -stop 5000 -stop_lr 10000 -fs '+str(100*len(speciesList))+" -fs_lr 10"
	# ADD LINE FOR SUBMITTING (change upype) !!!
	schedule(command,"reScheduled_"+kgID,mem=boost_mem+"g")
	with open ('error.txt','a') as errFile:
		errFile.write(alnProc['err'])
		sys.stderr.write('# INFO This error was saved : '+args.gene_dir+"/error.txt\n")


