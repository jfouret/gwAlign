setwd('/export/work/batnipah/phylogeny/alignments/gwAlign_ExAnn9/results/uc009xts.3/macrobats')
gencode=read.delim('gencode1.tab',header = F,sep = "\t",col.names = c('codon','aa','aminoacid'))
geneFolder='/export/work/batnipah/phylogeny/alignments/gwAlign_ExAnn9/results/uc009xts.3/macrobats'
alnFile=paste(geneFolder,'codon_aln.fa',sep='/')
posFile=paste(geneFolder,'posDict.tab',sep='/')

library(ggplot2)
library(ggthemes)
library(ggthemr)

read_alignment <- function(file,coding=F){
  
  raw_data <- readLines( file, warn = FALSE ) 
  seq_vector <- c()
  seq_name <- ""
  for (line in raw_data){
    splitter=""
    if (coding==T){splitter=" "}
    # New sequence record? Reset numbering
    if ( grepl("^>", line) ){
      if (seq_name!=''){
        if (coding==T){
          seq_vector[seq_name] <- strsplit(gsub("(.{3})", "\\1 ", gsub(' ','',seq_vector[seq_name])),split=splitter)
        }
        else{
          seq_vector[seq_name] <- strsplit(gsub(' ','',seq_vector[seq_name]),split=splitter)
        }
      }
      seq_name <- sub("^>", "", line)
      seq_vector[seq_name] <- ""     
    }
    else {
      temp_seq <- gsub(" ","",line)
      temp_seq <- gsub("\n","",temp_seq)
      seq_vector[seq_name] <- paste( seq_vector[seq_name], temp_seq, sep=splitter )
    }
    
  }
  if (coding==T){
    seq_vector[seq_name] <- strsplit(gsub("(.{3})", "\\1 ", gsub(' ','',seq_vector[seq_name])),split=splitter)
  }
  else{
    seq_vector[seq_name] <- strsplit(gsub(' ','',seq_vector[seq_name]),split=splitter)
  }
  # Is this an alignment?
  #seq_list <- strsplit(seq_vector, split = splitter)
  #lengths <- sapply(seq_list, length)
  #if ( sum(lengths != lengths[1]) != 0 )
  #stop("Your provided file is not an alignment. Please provide an alignment file in FASTA format to use alignfigR.")
  # Return sequence data parsed into named list
  seq_vector 
}
aln=read_alignment('codon_aln.fa',coding=T)
nbSeq=length(aln)


###### color palette ######
##nucl
#palette=c('red','blue','green','yellow','grey','black')
#names(palette)=c("A","T","C","G","-","!")
#palette Hydrophobe
pH=c('yellow','gold','greenyellow','darkolivegreen3','darkseagreen','seagreen','forestgreen','darkgreen')
names(pH)=c('A','V','I','L','M','F','Y','W')
#palette charged
pc=c('darkred','firebrick','indianred','cornflowerblue','midnightblue')
names(pc)=c('R','H','K','D','E')
#palette polar
pp=c('maroon','violetred','deeppink','turquoise')
names(pp)=c('S','T','N','Q')
#palette 
pstruct=c('cyan','cadetblue','lightcyan','skyblue')
names(pstruct)=c('C','U','P','G')
#paletteAA=palette(rainbow(length(unique(gencode$aa))))
#names(paletteAA)=unique(gencode$aa)
paletteAA=c(pH,pc,pp,pstruct,c('*'='white'))
paletteCodon=paletteAA[gencode$aa]
names(paletteCodon)=gencode$codon
#prot=data.frame(unique(gencode[,c('aa','aminoacid')]))
#prot$aa=factor(prot$aa,levels = names(paletteAA))
#prot$pos=seq(1,21)
#prot$x=1
#legend=ggplot(prot,aes(x=x,y=pos))+
#  #theme_bare+
#  geom_rect(aes(xmin=x,ymin=pos,xmax=x+1,ymax=pos+1,fill=aa))+
#  #geom_text(aes(x=case_w*pos+case_w/2,y=specID+0.5,label=seq),size=sizeLetter)+
#  scale_fill_manual(values=paletteAA)
#id=length(names(aln)):1
#names(id)=names(aln)


### parse the alignment to a ggplot friendly format
init=0
for (spec in names(aln)){
  data_tmp=data.frame(aln[spec])
  names(data_tmp)=c('seq')
  data_tmp$seq=factor(data_tmp$seq)
  data_tmp$specName=spec
  data_tmp$specID=id[spec]
  
  if (init==0){
    seqLen=length(data_tmp$seq)
    data_tmp$pos=1:seqLen
    init=1
    data=data_tmp
  }
  else{
    data_tmp$pos=1:seqLen
    data=rbind(data,data_tmp)
  }
}

library(grid)

# df for positions numerotation
data_tmp$specID=1+data_tmp$specID
data_tmp=data_tmp[c(1,seq(5,seqLen,5)),]
data_tmp$specID=nbSeq+1
# df for names
data_name=data[,c('specName','specID')]
data_name=unique(data_name)
data_name$x=-0.5
data_lim=subset(data_name,specID==1)
names(data_lim)=c('limchar','pos','x')
data_lim$x=-10
data_lim$limchar='.'
data_lim2=subset(data_name,specID==1)
names(data_lim2)=c('limchar','pos','x')
data_lim2$x=seqLen+1
data_lim2$limchar='.'

data_lim=rbind(data_lim,data_lim2)

#graphic options
sizeLetter=28.45276*0.1
theme_bare <- theme(
  axis.line = element_blank(), 
  axis.text.x = element_blank(), 
  axis.text.y = element_blank(),
  axis.ticks = element_blank(), 
  axis.title.x = element_blank(), 
  axis.title.y = element_blank(), 
  legend.position = "none", 
  panel.background = element_blank(), 
  panel.border = element_blank(), 
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(), 
  plot.background = element_blank(),
  legend.margin = unit(0, "cm"),
  plot.margin	=unit(c(0,0,0,0), "cm")
)
case_w=2.15 # largeur per codon in pixels

# extract complementary non-usual codon (ambiguous, frameshift)
compCodon=unique(subset(data,(grepl('-',seq))|(grepl('!',seq)))$seq)
paletteComp=replicate(length(compCodon),'white')
names(paletteComp)=compCodon
paletteFinal=c(paletteCodon,paletteComp)




aln_plot=ggplot(data,aes(x=pos,y=specID))+
  theme_bare+
  geom_rect(aes(xmin=case_w*pos,ymin=specID,xmax=case_w*pos+case_w,ymax=specID+1,fill=seq))+
  geom_text(aes(x=case_w*pos+case_w/2,y=specID+0.5,label=seq),size=sizeLetter)+
  scale_fill_manual(values=paletteFinal)+
  geom_text(data = data_tmp,aes(x=case_w*pos+case_w/2,y=specID+0.5,label=pos),size=sizeLetter)+
  geom_text(data = data_name,aes(x=x,y=specID+0.5,label=specName),size=sizeLetter,hjust=1)+
  geom_text(data = data_lim,aes(x=x,y=pos+0.5,label=limchar),size=sizeLetter,hjust=1)
  


## get exon and block only for codons

data_pos=read.delim('posDict.tab',header=T,sep="\t",na.strings = '.')

data_pos=data_pos[seq(3,dim(data_pos)[1],3),]
data_pos$reference=data_pos$reference/3
data_pos$block=data_pos$block/3
data_pos$alignment=data_pos$alignment/3

data_pos$blockStatus='conserved'
data_pos[is.na(data_pos$block),'blockStatus']='not conserved'
data_pos$pos=1
paletteblock=c('conserved'='darkgreen','not conserved'='white')

data_pos$exonStatus='odd'
data_pos[(data_pos$exon %% 2) == 0,'exonStatus']='even'

paletteExon=c('even'='blue','odd'='black')

block_plot=ggplot(data_pos,aes(x=alignment,y=pos))+
  theme_bare+
  geom_rect(aes(xmin=case_w*alignment,ymin=pos,xmax=case_w*alignment+case_w,ymax=pos+1,fill=blockStatus))+
  scale_fill_manual(values=paletteblock)+
  geom_text(data = data_lim,aes(x=x,y=pos+0.5,label=limchar),size=sizeLetter,hjust=1)

exon_plot=ggplot(data_pos,aes(x=alignment,y=pos))+
  theme_bare+
  geom_rect(aes(xmin=case_w*alignment,ymin=pos,xmax=case_w*alignment+case_w,ymax=pos+1,fill=exonStatus))+
  scale_fill_manual(values=paletteExon)+
  geom_text(data = data_lim,aes(x=x,y=pos+0.5,label=limchar),size=sizeLetter,hjust=1)


### Domain & Motif plot

domFile='/export/work/batnipah/phylogeny/alignments/knownCanonical_Zhang2/results/IFIT2-uc009xts.3.uniprot.bed'
#0-based -->
#start inclusive
#end inclusive

data_dom=read.delim(domFile,header = T,sep = "\t")
data_dom$Nclass=as.integer(data_dom$class)
data_dom$start=data_dom$start+1
data_dom$end=data_dom$end+1
### ATTETIN CHANGE COORDONATES WITH DICTIONARY
data_dom_name=unique(data_dom[,c('class','Nclass')])
data_dom_name$x=-1
nbUniProt=dim(data_dom_name)[1]
data_dom$type=as.character(data_dom$type)
data_dom[data_dom$type=='None','type']=as.character(data_dom[data_dom$type=='None','class'])

dom_plot=ggplot(data_dom,aes(x=start,y=Nclass))+
  theme_bare+
  geom_rect(aes(xmin=start,ymin=Nclass,xmax=end,ymax=Nclass+1))+
  geom_text(aes(x=start+(end-start)/2,y=Nclass+0.5,label=type),size=sizeLetter)+
#  scale_fill_manual(values=paletteFinal)#+
#  geom_text(data = data_tmp,aes(x=case_w*pos+case_w/2,y=specID+0.5,label=pos),size=sizeLetter)+
  geom_text(data = data_dom_name,aes(x=x,y=Nclass+0.5,label=class),size=sizeLetter,hjust=1)+
  geom_text(data = data_lim,aes(x=x,y=pos+0.5,label=limchar),size=sizeLetter,hjust=1)


library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)


plan_uniProt=c(matrix(4,1,nbUniProt))

plan=matrix(c(c(1,1,1,1,1,1,1,1,1,1,1,1,1,2,3),plan_uniProt))
W=(seqLen/3)*case_w
H=(nbSeq+1)/3+2+nbUniProt

test=grid.arrange(aln_plot,block_plot,exon_plot,dom_plot,layout_matrix=plan)

ggsave(filename = 'test.pdf',plot =test ,width = W,height = H,units = 'cm',device = 'pdf',limitsize=F)






