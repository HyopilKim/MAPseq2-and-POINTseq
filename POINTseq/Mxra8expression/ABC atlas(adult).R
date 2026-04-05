library(ggplot2)
library(anndata)
library(RColorBrewer)
library(reticulate)



# to find expression of a specific gene. for example, mxra8 (ENSMUSG00000029070)
# to find expression of other genes (Slc11a2=Nramp2, Lrp8=ApoER2)
Mxra8 <- 'ENSMUSG00000029070'
Slc11a2 <- 'ENSMUSG00000023030'
Ldlrad3 <-'ENSMUSG00000048058'
Vldlr <- 'ENSMUSG00000024924'
Lrp8 <- 'ENSMUSG00000028613'
Actb <- 'ENSMUSG00000029580'
Syn1 <- 'ENSMUSG00000037217'
Map2 <- 'ENSMUSG00000015222'

genelist <- c(Mxra8,Slc11a2,Ldlrad3,Vldlr,Lrp8,Actb,Syn1,Map2)
genenamelist <- c('Mxra8','Slc11a2','Ldlrad3','Vldlr','Lrp8','Actb','Syn1','Map2')

#cell cluster data
load('cell_data1.RData')


#read h5ad
import <- c("WMB-10Xv3-CB-log2.h5ad","WMB-10Xv3-CTXsp-log2.h5ad","WMB-10Xv3-HPF-log2.h5ad","WMB-10Xv3-HY-log2.h5ad","WMB-10Xv3-Isocortex-2-log2.h5ad","WMB-10Xv3-Isocortex-1-log2.h5ad","WMB-10Xv3-OLF-log2.h5ad","WMB-10Xv3-P-log2.h5ad","WMB-10Xv3-PAL-log2.h5ad","WMB-10Xv3-MY-log2.h5ad","WMB-10Xv3-MB-log2.h5ad","WMB-10Xv3-STR-log2.h5ad","WMB-10Xv3-TH-log2.h5ad","WMB-10Xv2-CTXsp-log2.h5ad","WMB-10Xv2-HPF-log2.h5ad","WMB-10Xv2-HY-log2.h5ad","WMB-10Xv2-Isocortex-1-log2.h5ad","WMB-10Xv2-Isocortex-2-log2.h5ad","WMB-10Xv2-Isocortex-3-log2.h5ad","WMB-10Xv2-Isocortex-4-log2.h5ad","WMB-10Xv2-MB-log2.h5ad","WMB-10Xv2-OLF-log2.h5ad","WMB-10Xv2-TH-log2.h5ad")

area <-c('CB','CTXsp','HPF','HY','Isocortex','Isocortex','OLF','P','PAL','MY','MB','STR','TH','CTXsp','HPF','HY','Isocortex','Isocortex','Isocortex','Isocortex','MB','OLF','TH')

data_all <- list()

for (j in 1:8) {
  data_inter <- data.frame(cluster_annotation_term_label=character(),gene=numeric(),area=character())
  for (i in 1:23) {
  metadata <- read_h5ad(import[i])
  exp_filt <- data.frame(gene=metadata[,genelist[j]],cell_label=rownames(metadata))
  rm(metadata)
  cell_data2 <- merge(cell_data1,exp_filt)
  rm(exp_filt)
  cell_data2 <- cell_data2[,c('cluster_annotation_term_label',genelist[j])]
  cell_data2$area <- area[i]
  colnames(cell_data2)[2] <-'gene' 
  data_inter <- rbind(data_inter,cell_data2)
  rm(cell_data2)
  gc()
  }
  data_all[[j]] <- data_inter
} 

save(data_all,file='data_all_log.Rdata')


genelist <- c(Mxra8,Slc11a2,Ldlrad3,Vldlr,Lrp8,Actb,Syn1,Map2)



#Coarse cell clusters
data_all_coarse <- list()
for (i in 1:8){
  data_plot <- data_all_coarse[[i]]
  data_plot$cluster_annotation_term_label[data_plot$cluster_annotation_term_label=='Glial'] <- 'Glial'
  data_all_coarse[[i]] <- data_plot
}



save(data_all_coarse,file='data_all_log_coarse.Rdata')
load('data_all_log_coarse.Rdata')


#############for graph



#graph Mxra8, Actb, Syn1, Map2 all cells
data_plot <- data.frame(gene=numeric(),name=character())
for (i in c(1,2,3,4,5)){
  data_tem <- data_all_coarse[[i]]
  data_teem <- data.frame(gene=data_tem$gene[data_tem$cluster_annotation_term_label=='Neuron'],name=genenamelist[i])
  data_plot <- rbind(data_plot,data_teem)
} 

ggplot(data_plot,aes(x=factor(name,levels=genenamelist[c(1,2,3,4,5)]),y=gene))+geom_boxplot(width=0.3,fill='gray',linewidth=1,outlier.shape = NA)+theme_classic()+theme(aspect.ratio=1,axis.text.x=element_text(size=30,angle=45,color='black',hjust=1),axis.text.y=element_text(size=30,color='black'),plot.margin=unit(c(1,1,1,1),'cm'),legend.position = 'none')+labs(y='',x='')+scale_y_continuous(expand=c(0,0),limits = c(0,12),breaks=seq(0,12,4))+labs(fill = NULL)+geom_violin(scale="width",fill=NA,linewidth=1.5,adjust = 5)



#graph according to clusters
data_plot <- data_all_coarse[[1]]

ggplot(data_plot,aes(x=factor(cluster_annotation_term_label,levels=c('Neuron','Glial','Vascular','Immune'),labels=c('Neuronal','Glial','Vascular','Immune')),y=gene),fill=factor(cluster_annotation_term_label))+geom_boxplot(width=0.3,fill='gray',size=1.5,outlier.shape = NA)+theme_classic()+theme(aspect.ratio=1,axis.text.x=element_text(size=30,angle=45,color='black',hjust=1),axis.text.y=element_text(size=30,color='black'),plot.margin=unit(c(1,1,1,1),'cm'),legend.position = 'none')+labs(y='',x='')+scale_y_continuous(expand=c(0,0),limits = c(0,12),breaks=seq(0,12,4))+labs(fill = NULL)+geom_violin(scale="width",fill=NA,linewidth=1.5,adjust = 5)


#graph by area
data_plot <- data_all_coarse[[1]]
data_plot <- data_plot[data_plot$cluster_annotation_term_label=='Neuron',]

ggplot(data_plot,aes(x=factor(area,levels=c('Isocortex','CTXsp','HPF','STR','TH','HY','PAL','MB','OLF','CB','P','MY')),y=gene))+geom_boxplot(width=0.3,fill='gray',size=1.5,outlier.shape = NA)+theme_classic()+theme(aspect.ratio=0.3,axis.text.x=element_text(size=25,angle=45,color='black',hjust=1),axis.text.y=element_text(size=25,color='black'),plot.margin=unit(c(1,1,1,1),'cm'),legend.position = 'none')+labs(y='',x='')+scale_y_continuous(expand=c(0,0),limits = c(0,14),breaks=seq(0,14,4))+geom_violin(scale="width",fill=NA,size=1,adjust = 5)

#t test
Mxra8_neuron <- data_plot$gene[data_plot$name==genenamelist[1]]
pvalues <- list()
for (i in 2:8){
  pvalues[[i]] <- t.test(Mxra8_neuron,data_plot$gene[data_plot$name==genenamelist[i]])
}

#anova
anova_result <- aov(gene ~ name, data = data_plot)
summary(anova_result)

bonferroni <- pairwise.t.test(data_plot$gene, data_plot$name, p.adjust.method = "bonferroni",na.rm = TRUE)


t.test(data_plot$gene[data_plot$name==genenamelist[2]])

save(bonferroni,anova_result,file='Anova and bonferroni among receptors.Rdata')

