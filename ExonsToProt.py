#!/usr/bin/python

import os
import sys
import argparse
import mysql.connector

version=1.2
year=2016
author='Julien Fouret'
contact='julien.fouret12@uniagro.fr'
##parse argument
parser = argparse.ArgumentParser(description='Sort gene/transcript from an alignement and reassemble all exons. Get the corresponding name of the protein with the accession numer (from ucsc or ncbi)',epilog="Version : "+str(version)+"\n"+str(year)+"\nAuthor : "+author+" for more informations or enquiries please contact "+contact,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-al', metavar='Alignment file', required=True, help="alignment file in fasta format")
parser.add_argument('-out', metavar='Output directory', required=True, help="folder in which to write all computed alignments")
parser.add_argument('-spec', metavar='Species file' , required=True, help="file with all species/assembly name (one per line)")
parser.add_argument('-host', metavar='mysql host' , required=False, help="",default='134.214.189.93')
parser.add_argument('-port', metavar='mysql port' , required=False, help="",default='4040')
parser.add_argument('-db',metavar='Database', required=False, help="name of the database 'ncbi' or 'ucsc'",default='ucsc',choices=['ucsc','ncbi'])
args=parser.parse_args()

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
		bed_file.write("\ŧ".join(map(str,row)))
	bed_file.close()
	cursor.close()

def Get_annotation(ID):
	global args
	
	if args.db=='ucsc':
		field='kgID'
	elif args.db=='ncbi':
		field='refseq'
	
	#global ID_type
	query_where='WHERE ref.'+field+'="'+ID+'";'
		
	query_name=['SELECT ref.geneSymbol,ref.kgID,ref.spID,ref.mRNA,refFlat.chrom,refFlat.strand,refFlat.txStart,refFlat.txEnd,refFlat.cdsStart,refFlat.cdsEnd,CAST(refFlat.exonStarts AS CHAR(10000) CHARACTER SET utf8) AS `exonStarts`,CAST(refFlat.exonEnds AS CHAR(10000) CHARACTER SET utf8) AS `exonEnds`']
	query_name.append('FROM kgXref ref')
	query_name.append('LEFT JOIN refFlat ON ( ref.geneSymbol=refFlat.geneName )')
	query_name.append(query_where)
		
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
	
	writedb(" ".join(query_name),args.out+"/"+ID+".name.bed",['geneSymbol','kgID','uniprot','mRNA','chrom','strand','txStart','txEnd','cdsStart','cdsEnd','exonStart','exonEnd'])
	name_file=open(args.out+"/"+ID+".name.bed","r")
	gene_symbol=name_file.readlines()[1].split("\ŧ")[0]
	name_file.close()
	os.system('mv '+args.out+"/"+ID+".name.bed"+' '+args.out+"/"+gene_symbol+'-'+ID+".name.bed")
	writedb(" ".join(query_uniprot),args.out+"/"+gene_symbol+'-'+ID+".uniprot.bed",['geneSymbol','start','end','class','type','OutOfRange'])
	writedb(" ".join(query_kegg),args.out+"/"+gene_symbol+'-'+ID+".kegg.bed",['geneSymbol','mapID','description'])
	
	
	return gene_symbol

#start program
cnx = mysql.connector.connect(user='genome',host=args.host,port=args.port,database='hg19')

species_file=open(args.spec,'r')
species=[]
for line in species_file.readlines():
	specie=line.rstrip()
	species.append(specie)
species_file.close()

seq=dict()

log_file=open(args.out+"/log","w")
exons_file=open(args.al,'r')

line=exons_file.readline()
#gene_id=line.split('_')[0][1:]+"_"+line.split('_')[1]
gene_id=line.split('_')[0][1:]

while line!='':
	
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


	prot_file=open(args.out+"/"+gene_name+"-"+ref_id+".fa","w")

	for specie in seq:
		prot_file.write(">"+specie+"\n")
		prot_file.write(seq[specie]+"\n")
	prot_file.close()
	
	#TODO
	#if (seq[species[0]].length()%3)!=0:
	#	log_file.write('WANRNING gene '+gene_name' with id '+gene_id+' has a sequence not proportional to 3')
		
exons_file.close()
log_file.close()
cnx.close()
