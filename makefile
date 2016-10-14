SHELL=sh
GITREPO=$(shell pwd)
GITREPOSED=$(shell pwd | sed 's/\//\\\//g')
GITVERSION=$(shell git describe --tags | sed 's/^v//g')
INSTALLPATH=/usr/local/bin/
all:
	mkdir -p bin
	sed -e "s/SEDMATCHGITREPO/${GITREPOSED}/g" scripts/Extend.py 					| sed -e "s/SEDMATCHGITVERSION/${GITVERSION}/g" > bin/gwAlign-Extend
	sed -e "s/SEDMATCHGITREPO/${GITREPOSED}/g" scripts/Exons2annCDS.py 				| sed -e "s/SEDMATCHGITVERSION/${GITVERSION}/g" > bin/gwAlign-Exons2annCDS
	sed -e "s/SEDMATCHGITREPO/${GITREPOSED}/g" scripts/KeggGo.py 					| sed -e "s/SEDMATCHGITVERSION/${GITVERSION}/g" > bin/KeggGo.py
	sed -e "s/SEDMATCHGITREPO/${GITREPOSED}/g" scripts/bam2consensus.py				| sed -e "s/SEDMATCHGITVERSION/${GITVERSION}/g" > bin/bam2consensus.py
	chmod 755 -R bin
install :
	ln -s ${GITREPO}/bin/gwAlign-Extend ${INSTALLPATH}
	ln -s ${GITREPO}/bin/gwAlign-Exons2annCDS ${INSTALLPATH}
force_install :
	ln -sf ${GITREPO}/bin/gwAlign-Extend ${INSTALLPATH}
	ln -sf ${GITREPO}/bin/gwAlign-Exons2annCDS ${INSTALLPATH}
uninstall :
	rm ${INSTALLPATH}/gwAlign-Extend
	rm ${INSTALLPATH}/gwAlign-Exons2annCDS
doc :
	Rscript -e "rmarkdown::render('README.md')"
clean :
	-rm -R bin
	-rm README.html
