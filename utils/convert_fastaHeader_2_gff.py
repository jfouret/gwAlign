
import argparse

version=1.0
year=2017
author='Julien Fouret'
employer='ViroScan3D'
contact='julien@fouret.me'
Licence=" ".join(["LICENCE V3D-intern",
"The author grants you the permission to use the software under the following conditions:",
"(i) Distribution is limited to ViroScan3D or ProfileXpert, usage via network is considered as a distribution",
"(ii) No modification is allowed without the permissions of the author(s) or ViroScan3D",
"(iii) Commercial use is not allowed without the permission of Viroscan3D",
"(iv) Any use must credit the contribution from the author and ViroScan3D",
"(v) NO WARRANTY is granted",
"(vi) ViroScan3D or the author are not liable for any use of this software"])

desc="convert ucsc fasta header forme exons coordonate in multialignment files in gff format"
dep="Python>=2.7"

##parse argument
parser = argparse.ArgumentParser(description=desc,epilog="Version : "+str(version)+" | "+str(year)+" | Author : "+author+" | Employer : "+employer+" for more informations or enquiries please contact "+contact+" | DEPENDENCES : "+dep+" | "+Licence,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-input', metavar='.fa', required=False, help="input file fasta headers or fasta")
parser.add_argument('-out', metavar='.gff', required=False, help="output gff formatted file")
parser.add_argument('-v',action="store_true",required=False, help="verbose")
args=parser.parse_args()

import re
import sys

reg=re.compile("^>(.*)_.*_([0-9]*)_([0-9]*).* (.*):(.*)-(.*)([+\-]{1})$")

inFile=open(args.input)
outFile=open(args.out,'w')
currentGeneID='None'

for line in inFile.readlines():
	line=line.rstrip()
	m=reg.match(line)
	if m:
		scaff=m.group(4)
		geneID=m.group(1)
		exonsPos=m.group(2)
		exonsTot=m.group(3)
		scaffStart=m.group(5)
		scaffEnd=m.group(6)
		strand=m.group(7)
		if currentGeneID=='None':
			currentGeneID=geneID
			mainScaff=scaff
			mRNA_start=scaffStart
			mRNA_id=geneID
			mRNA_name=geneID
			mRNA_strand=strand
			toWrite=""
		if currentGeneID==geneID and mainScaff==scaff and mRNA_strand==strand :
			toWrite+="\t".join([scaff,"ucsc_alignments","cds",scaffStart,scaffEnd,'.',strand,'0',"ID="+geneID+";Name="+geneID+"_"+exonsPos+"_"+exonsTot+";Parent="+mRNA_id])+"\n"
			mRNA_end=scaffEnd
		elif currentGeneID==geneID and (mainScaff!=scaff or mRNA_strand!=strand):
			outFile.write("\t".join([mainScaff,"ucsc_alignments","mRNA_coding",mRNA_start,mRNA_end,'.',mRNA_strand,'0',"ID="+mRNA_id+";Name="+mRNA_name])+"\n")
			mainScaff=scaff
			mRNA_strand=strand
			mRNA_start=scaffStart
			mRNA_id=mRNA_id+"_"
			currentGeneID=mRNA_id
			outFile.write(toWrite)
			mRNA_end=scaffEnd
			toWrite="\t".join([scaff,"ucsc_alignments","cds",scaffStart,scaffEnd,'.',strand,'0',"ID="+geneID+";Name="+geneID+"_"+exonsPos+"_"+exonsTot+";Parent="+mRNA_id])+"\n"
		else :
			outFile.write("\t".join([mainScaff,"ucsc_alignments","mRNA_coding",mRNA_start,mRNA_end,'.',mRNA_strand,'0',"ID="+mRNA_id+";Name="+mRNA_name])+"\n")
			outFile.write(toWrite)
			mRNA_id=geneID
			currentGeneID=mRNA_id
			mainScaff=scaff
			mRNA_strand=strand
			mRNA_start=scaffStart
			mRNA_end=scaffEnd
			mRNA_name=geneID
			toWrite="\t".join([scaff,"ucsc_alignments","cds",scaffStart,scaffEnd,'.',strand,'0',"ID="+geneID+";Name="+geneID+"_"+exonsPos+"_"+exonsTot+";Parent="+mRNA_id])+"\n"
mRNA_end=scaffEnd
outFile.write("\t".join([mainScaff,"ucsc_alignments","mRNA_coding",mRNA_start,mRNA_end,'.',mRNA_strand,'0',"ID="+geneID+";Name="+mRNA_name])+"\n")
outFile.write(toWrite)

outFile.close()
inFile.close()

sys.exit()

