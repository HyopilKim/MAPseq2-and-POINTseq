

library(jmuOutlier)
library(MetaNeighbor)
library(R.matlab)
library(SummarizedExperiment)
library(corrplot)


imported <- readMat("forMeta_Corr.mat")

targets26 <- c("MOs-r-Ip","SS-Ip","MOs-r-Con","MOp-Con","SS-Con","Str-r-ip","Str-i-ip","Str-c-ip","Str-r-con","Str-i-con","Str-c-con","TH-mr-Ip","TH-lr-Ip","TH-mc-Ip","TH-lc-Ip","TH-mr-Con","TH-lr-Con","TH-mc-Con","TH-lc-Con","MRN-Ip","PG-Ip","SC-Ip","MY-l-Ip","MY-m-Ip","MY-l-Con","MY-m-Con")


combined_data <- t(as.matrix(rbind(imported$sorteddata.mapseq1,imported$sorteddata.mapseq2)))
labels <-  c(rep(1, dim(imported$sorteddata.mapseq1)[1]), rep(2, dim(imported$sorteddata.mapseq2)[1]))


labels <-  factor(labels,levels=c(1,2),labels = c('MAPseq1','MAPseq2'))

MAPseq1_clusters <- factor(imported$Tsorted.5.mapseq1,levels=c(2,1,3,4,5),labels = c('ITi','IT Str-','IT Str+','ET','CT'))

MAPseq2_clusters <- factor(imported$Tsorted.mapseq2,levels=c(4,3,5,1,2),labels = c('ITi','IT Str-','IT Str+','ET','CT'))

clusters <- c(MAPseq1_clusters,MAPseq2_clusters)


rownames(combined_data) <- targets26
colnames(combined_data) <- 1:length(labels)


se <- SummarizedExperiment(assays=list(counts=combined_data),
                           colData=DataFrame(study_id=labels,
                                             cell_type=clusters))


celltype_NV=MetaNeighborUS(var_genes = targets26,
                             dat = se,
                             study_id = se$study_id,
                             cell_type = se$cell_type)

cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(100))
breaks = seq(0, 1, length=101)
gplots::heatmap.2(celltype_NV,
                  margins=c(10,10),
                  keysize=1,
                  key.xlab="AUROC",
                  key.title="",
                  trace = "none",
                  density.info = "none",
                  col = cols,
                  breaks = breaks,
                  offsetRow=0.1,
                  offsetCol=0.1,
                  cexRow = 1.2,
                  cexCol = 1.2)



#Correlation
data1 <- t(data.frame(imported$mean.MAPseq1))
rownames(data1) <- targets26
colnames(data1) <- c('MAPseq1 ITi','MAPseq1 ITc Str-','MAPseq1 ITc Str+','MAPseq1 ET','MAPseq1 CT')

data2 <- t(data.frame(imported$mean.MAPseq2))
rownames(data2) <- targets26
colnames(data2) <- c('MAPseq2 ITi','MAPseq2 ITc Str-','MAPseq2 ITc Str+','MAPseq2 ET','MAPseq2 CT')

M <- cor(x=data2,y=data1,method='pearson',use='all.obs')

p <- matrix(0, nrow = ncol(data1), ncol = ncol(data2))

for (i in 1:ncol(data2)){
  for (j in 1:ncol(data1)){
    eachp <- perm.cor.test(x=data2[,i],y=data1[,j], method = "pearson", num.sim = 10000)
    p[i,j] <- eachp$p.value}
}


corrplot(M, method = 'color',col = rev(COL2('RdBu', 200)),tl.pos = 'n',cl.pos = 'n')

