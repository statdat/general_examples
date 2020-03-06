#### loading ####
#df
library(reshape2)
library(dplyr)
library(ggplot2)
library(fastDummies)
library(stringr)
library(tidyr)
#plotting
library(ggpubr)
library(cowplot)
library(ggsci)
library(scales)
library(RColorBrewer)
#Tree
#library(ape)
library(ggtree)
library(ggstance)
library(phytools)

setwd('/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/enr_comparison_project/data/fabv_data')
data_o <- '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/enr_comparison_project/data/fabv_data/miscellaneous'
figures_o <- '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/enr_comparison_project/data/fabv_data/miscellaneous/figures'

#### Section 0 ####
#### Making diagnostic plots of checkm and quast outputs and also species composition
# read in checkm output
df_cm_orig <- read.csv('/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/enr_comparison_project/data/hiseq_genome_assemblies/checkm_output_modified.txt',sep='\t')
# read in quast output
df_quality <- read.csv('processed_genomes_metadata.csv')%>%filter(source=='USER')%>%select(id,total_contigs)
# read in taxonomy generation output
df_tax <- read.csv('taxonomic_assignments.csv')%>%rename(id=cleaned_filename)%>%select(id,species)
# making id names equal. This means removing the 205 from all USER genomes
df_tax$id <- str_replace(df_tax$id,'205_','')
df_cm <- df_cm_orig%>%rename(id=Bin.Id)%>%select(id,Completeness,Contamination)
df_cm <- merge(df_cm,df_quality,by='id')
#removing isolates that were not used for pabn testing (all the HS isolates)
df_cm <- df_cm%>%filter(!id%in%c('HS_1','HS_2','HS_3','HS_4','HS_5'))
df_cm <- merge(df_cm,df_tax,by='id',all.x=TRUE)
## Plottig completeness and contamination
# extracting data necessary for histogram
df_hist <- melt(df_cm%>%select(id,Completeness,Contamination))
p1 <- ggplot(df_hist)+
  geom_density(aes(x=value,fill=variable),alpha=1)+
  theme(axis.text.x=element_text(size=45),axis.text.y=element_text(size=45),axis.title.x=element_text(size=45),axis.title.y=element_text(size=45),
        legend.text = element_text(size=60),legend.title=element_blank(),legend.position=c(.6,.9))+
  scale_x_continuous(breaks=seq(0,100,10))+
  scale_y_continuous(breaks=seq(0,.8,.1))+
  scale_fill_jco()+
  xlab(label='Percent')
# saving
save_plot(paste(figures_o,'sec0_0.png',sep='/'),p1, base_height = 18, base_aspect_ratio = 1.4)
## Plotting contig count
df_hist2 <- df_cm%>%select(id,total_contigs)
p2 <- ggplot(df_hist2)+
  geom_histogram(aes(x=total_contigs,fill=),alpha=1)+
  theme(axis.text.x=element_text(size=45),axis.text.y=element_text(size=45),axis.title.x=element_text(size=45),axis.title.y=element_text(size=45),
        legend.text = element_text(size=60),legend.title=element_blank(),legend.position=c(.6,.9))+
  scale_x_continuous(breaks=seq(0,200,25))+
  #scale_y_continuous(breaks=seq(0,.1,.01))+
  scale_fill_jco()+
  xlab(label='Contigs')
# saving
save_plot(paste(figures_o,'sec0_1.png',sep='/'),p2, base_height = 18, base_aspect_ratio = 1.4)
## Plotting species count
# generating count per species and assigning plotting order
df_sp <- df_cm%>%select(id,species)%>%group_by(species)%>%mutate(count=length(id))%>%arrange(desc(count))%>%ungroup()%>%select(species,count)%>%unique()%>%mutate(order=seq(1,nrow(df_sp)))
# making manual colors for plotting
nb.cols <- 13
mycolors <- colorRampPalette(brewer.pal(13, 'Set3'))(nb.cols)
p3 <- ggplot(df_sp,aes(x=reorder(species,order),y=count,fill=species))+coord_flip()+geom_bar(stat='identity')+
  theme(axis.text.x=element_text(size=45),axis.text.y=element_text(size=45),axis.title.x=element_text(size=45),axis.title.y=element_text(size=45),
        legend.position='none')+
  xlab(label='')+
  scale_fill_manual(values = mycolors)
# saving
save_plot(paste(figures_o,'sec0_2.png',sep='/'),p3, base_height = 18, base_aspect_ratio = 1.4)


