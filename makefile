### CAN BE CHANGED
# install path
INSTALLPATH=/export/bin
# Default picard.jar path
PICARD=/export/bin/picard-tools-2.1.0/picard.jar
# Default marcse.jar path
MACSE=/export/source/script/macse_v1.2.jar
# Default queue name for PBS
QUEUE=batch
#default path for genepythia
GENEPYTHIA=/export/work/batnipah/juTools/funMining/bibliography
#default path for GATK
GATK=/export/bin/source/GenomeAnalysisTK/GenomeAnalysisTK.jar

### DO NOT CHANGE
SHELL=bash
GITREPO=$(shell pwd)
GITREPOSED=$(shell pwd | sed 's/\//\\\//g')
GITVERSION=$(shell git describe --tags | sed 's/^v//g')
SEDGENEPYTHIA=$(shell echo ${GENEPYTHIA} | sed 's/\//\\\//g')
SEDGATK=$(shell echo ${GATK} | sed 's/\//\\\//g')

progs = gwAlign-Extend gwAlign-Exons2annoCDS KeggGo.py bam2consensus.py alignCDS.py phasingBlocks.py
bins=$(addprefix bin/,$(progs))

### makefile core
all : $(bins)

$(bins) : bin 
	sed -e "s/SEDMATCHGITREPO/${GITREPOSED}/g" $(subst bin/,scripts/,$@) | \
	sed -e "s/SEDMATCHGITVERSION/${GITVERSION}/g" | \
	sed -e "s/SEDMATCHGATK/${SEDGATK}/g" | \
	sed -e "s/SEDMATCHPICARD/$(subst /,\\/,${PICARD})/g" | \
	sed -e "s/SEDMATCHGENEPYTHIA/${SEDGENEPYTHIA}/g" | \
	sed -e "s/SEDMATCHQUEUE/${QUEUE}/g" | \
	sed -e "s/SEDMATCHMACSE/$(subst /,\\/,${MACSE})/g" > $@
	chmod 755 $@

bin :
	mkdir bin
	chmod 755 bin

.PHONY : install
install : 
	ln -s ${GITREPO}/bin/gwAlign-Extend ${INSTALLPATH}
	ln -s ${GITREPO}/bin/gwAlign-Exons2annoCDS ${INSTALLPATH}
	ln -s ${GITREPO}/bin/KeggGo.py ${INSTALLPATH}
	ln -s ${GITREPO}/bin/bam2consensus.py ${INSTALLPATH}
	ln -s ${GITREPO}/bin/alignCDS.py ${INSTALLPATH}
	ln -s ${GITREPO}/bin/phasingBlocks.py ${INSTALLPATH}

.PHONY : force_install
force_install : $(bins)
	ln -sf ${GITREPO}/bin/gwAlign-Extend ${INSTALLPATH}
	ln -sf ${GITREPO}/bin/gwAlign-Exons2annoCDS ${INSTALLPATH}
	ln -sf ${GITREPO}/bin/KeggGo.py ${INSTALLPATH}
	ln -sf ${GITREPO}/bin/bam2consensus.py ${INSTALLPATH}
	ln -sf ${GITREPO}/bin/alignCDS.py ${INSTALLPATH}
	ln -sf ${GITREPO}/bin/phasingBlocks.py ${INSTALLPATH}

.PHONY : uninstall
uninstall : 
	rm ${INSTALLPATH}/gwAlign-Extend
	rm ${INSTALLPATH}/gwAlign-Exons2annoCDS

.PHONY : doc
doc : 
	Rscript -e "rmarkdown::render('README.md')"

.PHONY : clean
clean :
	-rm -R bin
	-rm README.html
