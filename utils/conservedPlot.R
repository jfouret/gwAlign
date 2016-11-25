## Usage : $Rscript conservedPlot.R stat.tab
# stat.tab is produced by gwAlign-Exon2annoCDS -mode stat_phasing
args <- commandArgs(TRUE)

statName=args[1]
pdfName=paste(args[1],'.pdf',sep='')

data=read.delim(statName,na.strings = 'NA')

library(ggplot2)
library(gridExtra)


data_before=data.frame(data[,c('before')])
names(data_before)=c('length')
data_before$status='raw'
data_after=data.frame(data[,c('after')])
names(data_after)=c('length')
data_after$status='conserved'
data_long=rbind(data_before,data_after)

data$delta=data$before-data$after
data$percentage=100*data$after/data$before

both=ggplot(data_long,aes(length,fill=status))+
  theme_bw()+
  geom_density(kernel='gaussian',alpha=0.5,bw=1)+
  scale_fill_manual(values = c('raw'='dodgerblue','conserved'='forestgreen'))+
  coord_cartesian(xlim=c(0,5000))+
  xlab('Size (bp)')
delta=ggplot(data,aes(delta))+
  theme_bw()+
  geom_density(kernel='gaussian',fill='black',bw=1)+
  coord_cartesian(xlim=c(0,2000))+
  xlab('Not conserved bases (bp)')
percentage=ggplot(data,aes(percentage))+
  theme_bw()+
  geom_density(kernel='gaussian',fill='black')+
  xlab('Percentage conserved')

lay <- rbind(c(1,1),
             c(2,3))
assemble=grid.arrange(both,delta,percentage, layout_matrix = lay)
ggsave(pdfName,assemble,width = 6,height = 4)

