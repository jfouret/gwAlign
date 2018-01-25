# Installation

* Parameters - configuration

```
./configure -h
```

* Create executables 

```
make
```

* Install 

```
sudo make install  
#OR to erase previewsly installed
sudo make force_install
```

* Uninstall

```
sudo make uninstall
```

* Delete executables

```
make clean
```

# List of main softwares

* `gwAlign-Extend`

![Pipeline for extension of ann alignment to a new specie](./gwAlign-Extend.PNG "Pipeline1")

* `gwAlign-Unify`

DECRIRE ET AJOUTER LA FIGURE

# Usage

For help run executable with -h or --help

# Dependencies

* Python 2.7 (or higher but never tested)
	* upype
	* Bio
	* pysam
	* numpy

* HPC with PBS or SLURM with qsub interface
* samtools (from version 1.3)
* Platypus
* PicardTools
* MACSE
