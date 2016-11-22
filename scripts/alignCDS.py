#!/usr/bin/python
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
parser.add_argument('-only_size' , action='store_true', help="jar file path for macse program")

args=parser.parse_args()

import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
import re
import glob

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

if args.only_size:
	do_align=False
else:
	do_align=True

speciesList=os.walk(args.gene_dir).next()[1]
os.chdir(args.gene_dir)
speciesExons=dict()
reName=re.compile('exon_(\d+)_ext(\d+)')
reNameAlt=re.compile('exon_(\d+)')
refExLen=dict()
for species in speciesList:
	speciesExons[species]=dict()
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
				sys.stderr.write("Waring 'Exon "+str(i)+" too long for "+species+"' at : "+args.gene_dir+"\n")
				do_align=True


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

command='java -Xmx20g -jar '+args.macse+' -prog alignSequences -seq ref.fa -seq_lr aligned.fa -stop 5000 -stop_lr 10000 -fs '+str(100*len(speciesList))

if do_align:
	alnProc=submitOneShell(command)
	if alnProc['err']!='':
		with open ('error.txt','a'):
			errFile.write(alnProc['err'])
			sys.stderr.write('Error with alignments alignment at : '+args.gene_dir+"\n")


