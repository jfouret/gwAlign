#!/usr/bin/python
import argparse
gitRepository='SEDMATCHGITREPO'
version='SEDMATCHGITVERSION'
year=2016
author='Julien Fouret'
contact='julien@fouret.me'
##parse argument
parser = argparse.ArgumentParser(description='Resume all kegg and Go in one file for functional analysis',epilog="Version : "+str(version)+"\n"+str(year)+"\nAuthor : "+author+" for more informations or enquiries please contact "+contact,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-outDir', metavar='/path', required=True, help="path of the output directory from positive selection analysis")
args=parser.parse_args()

import sys
import os
import re
from jupype import *

rootedDir=loadRoot(args.outDir)


dirList=os.listdir(rootedDir.results)
goDict=dict()
goList=list()
headGO=['goId','name','type','genes']
keggDict=dict()
keggList=list()
headKegg=['mapId','name','genes']
keggFileName=rootedDir.reports+'/kegg.tab'
goFileName=rootedDir.reports+'/go.tab'

for dirName in dirList:
	with open(rootedDir.results+'/'+dirName+'/annotation/consName.txt') as consNameFile:
		gene=consNameFile.readline().rstrip()
	File=open(rootedDir.results+'/'+dirName+'/annotation/go.tab','r')
	File.readline()
	for line in File.readlines():
		line=line.rstrip()
		lineList=line.split("\t")
		goId=lineList[1]
		name=lineList[2]
		gotype=lineList[3]
		if goId in goDict.keys():
			goDict[goId].append(gene)
		else:
			goDict[goId]=[gene]
			goList.append([goId,name,gotype])
	File.close()
	File=open(rootedDir.results+'/'+dirName+'/annotation/kegg.tab','r')
	File.readline()
	for line in File.readlines():
		line=line.rstrip()
		lineList=line.split("\t")
		mapId=lineList[1]
		name=lineList[2]
		if mapId in keggDict.keys():
			keggDict[mapId].append(gene)
		else:
			keggList.append([mapId,name])
			keggDict[mapId]=[gene]
	File.close()
keggFile=open(keggFileName,'w')
keggFile.write("\t".join(headKegg)+"\n")
goFile=open(goFileName,'w')
goFile.write("\t".join(headGO)+"\n")
for lineList in keggList:
	mapId=lineList[0]
	keggFile.write("\t".join(lineList)+"\t"+','.join(keggDict[mapId])+"\n")
for lineList in goList:
	goId=lineList[0]
	goFile.write("\t".join(lineList)+"\t"+','.join(goDict[goId])+"\n")
goFile.close()
keggFile.close()
