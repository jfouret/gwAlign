---
title: 'Manual: gwAlign suite'
author: "Julien Fouret"
mail: "julien@fouret.me"
output:
  html_document:
    toc: true # table of content true
    toc_float: true
    depth: 4  # upto three depths of headings (specified by #, ## and ###)
    theme: united  # many options for theme, this one is my favorite.
    highlight: tango  # specifies the syntax highlighting style
---

# Installation

* Parameters - configuration

Default configuration can be changed in the first part of the makefile. The desired path for installation and the paths of some dependencies have to be given.

* Create executables 

make

* Install 

sudo make install  

(OR to erase previewsly installed)

sudo make force_install

* Uninstall

sudo make uninstall

* Delete executables and docs

make clean

* Documentation html

make doc

# List of softwares

* gwAlign-Extend	

Extend portions of a genome-wide alignment to add a new species without reference genomes, only with DNA-Seq reads aligned vs other references.

![Pipeline for extension of ann alignment to a new specie](./gwAlign-Extend.PNG "Pipeline1")

* gwAlign-Exons2annCDS

DECRIRE ET AJOUTER LA FIGURE

# Usage

For help run executable with -h or --help

# Dependencies

* Python 2.7 (or higher but never tested)
	* jupype
* R
	* markdown
	* knitr
* PBS serveur with qsub
* samtools (from version 1.3)
* bcftools
* vcftools
* bedtools

ADD DEPENDENCIES FOR OTHER ....
