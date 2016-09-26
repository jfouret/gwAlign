#!/usr/bin/python
import os
import sys
import argparse
import mysql.connector
scriptName='ExonsToProt.py'
version=1.2
year=2016
author='Julien Fouret'
contact='julien.fouret12@uniagro.fr'
##parse argument
parser = argparse.ArgumentParser(description='Sort gene/transcript from an alignement and reassemble all exons. Get the corresponding name(NGNC) of the protein with the accession numer (from ucsc or ncbi) and add kegg, go and uniprot annotations',epilog="Version : "+str(version)+"\n"+str(year)+"\nAuthor : "+author+" for more informations or enquiries please contact "+contact,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-al', metavar='Alignment file', required=True, help="alignment file in fasta format")
parser.add_argument('-out', metavar='Output directory', required=True, help="folder in which to write all computed alignments")
parser.add_argument('-annotate_only',action='store_true', help="only annotation !")
parser.add_argument('-spec', metavar='Species file' , required=True, help="file with all species/assembly name (one per line)")
parser.add_argument('-host', metavar='mysql host' , required=False, help="",default='10.0.0.200')
parser.add_argument('-port', metavar='mysql port' , required=False, help="",default='3306')
parser.add_argument('-db',metavar='Database', required=False, help="name of the database 'ncbi' or 'ucsc'",default='ucsc',choices=['ucsc','ncbi'])
args=parser.parse_args()

sys.path.append('/export/home/jfouret/lib/')
from myfunctions import *

rootedDir=rootDir(args.out)
rootedDir.logs.writeArgs()
##define function
def writedb(query,file_name,header):
	global cnx
	print query
	print("\n")
	cursor=cnx.cursor()
	cursor.execute(query)
	rows=cursor.fetchall()
	bed_file=open(file_name,"w")
	bed_file.write("\t".join(header))
	for row in rows:
		bed_file.write("\n")
		bed_file.write("\t".join(map(str,row)))
	bed_file.close()
	cursor.close()
def Get_annotation(ID):
	global args
	global iterNoName
	global allSymbol
	if args.db=='ucsc':
		field='kgID'
	elif args.db=='ncbi':
		field='refseq'
	#global ID_type
	query_where='WHERE (ref.'+field+'="'+ID+'");'
	query_name=['SELECT hgnc.symbol,ref.kgID,ref.spID,ref.refseq,ref.geneSymbol']
	query_name.append('FROM kgXref ref')
	query_name.append('LEFT JOIN proteinDB.hgncXref hgnc ON ( hgnc.uniProt=ref.spID ) AND (hgnc.refSeq=ref.refseq)')
	query_name.append(query_where)

	query_alias=['SELECT ref.alias']
	query_alias.append('FROM kgAlias ref')
	query_alias.append(query_where)
		
	query_uniprot=['SELECT ref.geneSymbol,feats.start,feats.end,class.val AS `class`,type.val AS `type`,feats.softEndBits AS `OutOfRange`']
	query_uniprot.append('FROM uniProt.feature feats')
	query_uniprot.append('LEFT JOIN kgXref ref ON (ref.spID=feats.acc)')
	query_uniprot.append('LEFT JOIN uniProt.featureClass class ON (class.id=feats.featureClass)')
	query_uniprot.append('LEFT JOIN uniProt.featureType type ON (type.id=feats.featureType)')
	query_uniprot.append(query_where)
		
	query_kegg=['SELECT ref.geneSymbol,kegg.mapID,des.description']
	query_kegg.append('FROM keggPathway kegg')
	query_kegg.append('LEFT JOIN kgXref ref ON (ref.kgID=kegg.kgID)')
	query_kegg.append('LEFT JOIN keggMapDesc des ON (kegg.mapID=des.mapID)')
	query_kegg.append(query_where)

	query_go=['SELECT ref.geneSymbol, goa.goId, goterm.name, goterm.term_type']
	query_go.append('FROM go.goaPart goa')
	query_go.append('LEFT JOIN kgXref ref ON (ref.spID=goa.dbObjectId)')
	query_go.append('LEFT JOIN go.term goterm ON (goa.goId=goterm.acc)')
	query_go.append(query_where)

	writedb(" ".join(query_name),rootedDir.results+"/"+ID+".name.bed",['geneSymbol','kgID','uniprot','refSeq','kgName'])
	name_file=open(rootedDir.results+"/"+ID+".name.bed","r")
	name_file_line1=name_file.readlines()[1]
	gene_symbol=name_file_line1.split("\t")[0]

	hgnc=name_file_line1.split("\t")[0]
	spID=name_file_line1.split("\t")[2]
	refseq=name_file_line1.split("\t")[3]
	kgname=name_file_line1.split("\t")[4]
	gene_symbol=hgnc
	if (hgnc=='None' or hgnc=='') or (spID=='' and refseq=='') :
		gene_symbol='kg_'+kgname
	elif kgname=='None' or kgname=='':
		gene_symbol='sp_'+spID
	elif spID=='None' or spID=='':
		gene_symbol='rs_'+refseq
	elif refseq=='None' or refseq=='':
		gene_symbol='NoID_'+str(iterNoName)
		iterNoName+=1
	gene_symbol=renameFileName(gene_symbol)
	if (gene_symbol in allSymbol.keys()):
		allSymbol[gene_symbol]+=1
		gene_symbol='dup'+str(allSymbol[gene_symbol])+'_'+gene_symbol
	else:
		allSymbol[gene_symbol]=0
	name_file.close()
	os.system('mv '+rootedDir.results+"/"+ID+".name.bed"+' '+rootedDir.results+"/"+gene_symbol+'-'+ID+".name.bed")
	writedb(" ".join(query_uniprot),rootedDir.results+"/"+gene_symbol+'-'+ID+".uniprot.bed",['geneSymbol','start','end','class','type','OutOfRange'])
	writedb(" ".join(query_kegg),rootedDir.results+"/"+gene_symbol+'-'+ID+".kegg.bed",['geneSymbol','mapID','description'])
	writedb(" ".join(query_go),rootedDir.results+"/"+gene_symbol+'-'+ID+".go.bed",['geneSymbol','goId','name','type'])
	writedb(" ".join(query_alias),rootedDir.results+"/"+gene_symbol+'-'+ID+".alias.bed",['alias'])
	return gene_symbol
#start program
allSymbol=dict()

iterNoName=1
cnx = mysql.connector.connect(user='genome',host=args.host,port=args.port,database='hg19')
species_file=open(args.spec,'r')
species=[]
for line in species_file.readlines():
	specie=line.rstrip()
	species.append(specie)
species_file.close()
seq=dict()
exons_file=open(args.al,'r')
line=exons_file.readline()
#gene_id=line.split('_')[0][1:]+"_"+line.split('_')[1]
gene_id=line.split('_')[0][1:]
while line!='':
	if args.annotate_only:
		gene_name=Get_annotation(gene_id)
		line2=exons_file.readline() # pass the sequence
		line=exons_file.readline()
		ref_id=gene_id
		gene_id=line.split('_')[0][1:]
		while (gene_id==ref_id) and (line!=''):
			line2=exons_file.readline() # pass the sequence
                        line=exons_file.readline()
                        gene_id=line.split('_')[0][1:]
	else:
		for specie in species:
			seq[specie]=""
		gene_name=Get_annotation(gene_id)
 		specie=line.split('_')[1]
		seq[specie]=seq[specie]+exons_file.readline().rstrip()
		line=exons_file.readline()
		ref_id=gene_id
		gene_id=line.split('_')[0][1:]
		while (gene_id==ref_id) and (line!=''):
			specie=line.split('_')[1]
			seq[specie]=seq[specie]+exons_file.readline().rstrip()
			line=exons_file.readline()
			gene_id=line.split('_')[0][1:]
		prot_file=open(rootedDir.results+"/"+gene_name+"-"+ref_id+".fa","w")
		for specie in seq:
			prot_file.write(">"+specie+"\n")
			prot_file.write(seq[specie]+"\n")
		prot_file.close()
exons_file.close()
cnx.close()

with open(rootedDir.reports+'/all_withCounts.txt','w') as logFile:
	for key in allSymbol:
		logFile.write(key+"\t"+str(allSymbol[key])+"\n")

with open(rootedDir.reports+'/duplicate.txt','w') as logFile:
	for key in allSymbol:
		if allSymbol[key]>0:
			logFile.write(key+"\n")
sys.exit(0)


