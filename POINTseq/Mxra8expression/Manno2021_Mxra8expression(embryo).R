library(Seurat)
install.packages("remotes")
remotes::install_github("mojaveazure/seurat-disk")
library(ggplot2)
library(SeuratDisk)
library(Matrix)
library(tidyr)

loom <- Connect(filename = "your location/dev_all.loom", mode = "r")

# gene names
genes <- loom[["row_attrs"]][["Gene"]][]

# expression matrix: gene x cell
mat <- loom[["matrix"]][, ]
mat_sparse <- Matrix(mat, sparse = TRUE)
saveRDS(mat_sparse,file="mat_sparse.RDS",compress=F)
lib_size <- rowSums(mat_sparse)

# prevent division by zero
lib_size[lib_size == 0] <- 1

# calculate CPM (keep sparse format)
cpm_sparse <- Diagonal(x = 1 / lib_size) %*% mat_sparse * 1e6

class_vec <- loom[["col_attrs"]][["Class"]][]
cell_ids  <- loom[["col_attrs"]][["CellID"]][]

# genes of interest
gene_list <- c("Map2", "Gfap", "Mxra8", "Slc11a2", "Ldlrad3", "Vldlr", "Lrp8")

# find gene indices
gene_idx <- match(gene_list, genes)

# check for missing genes
if (any(is.na(gene_idx))) {
  print(gene_list[is.na(gene_idx)])
  stop("Some genes not found")
}

# select neuron cells only
neuron_cells <- which(class_vec == "Neuron")

# expression matrix (cell x gene)
expr_mat <- cpm_sparse[neuron_cells,gene_idx]
expr_mat2 <- log2(1+expr_mat)

test=as.matrix(expr_mat2)

# assign row/column names
df=data.frame(expr_mat2)

df=data.frame(test)

df_long <- pivot_longer(df, cols = everything(),
                        names_to = "gene",
                        values_to = "value")




ggplot(df_long, aes(x = factor(gene,levels=c("Mxra8", "Slc11a2", "Ldlrad3", "Vldlr", "Lrp8","Gfap","Map2")), y = value))+geom_violin(scale="width",fill=NA,size=1,adjust = 5)+theme_classic()+theme(aspect.ratio=0.6,axis.text.x=element_text(size=40,angle=45,color='black',hjust=1),axis.text.y=element_text(size=35,color='black'),plot.margin=unit(c(1,1,1,1),'cm'),legend.position = 'none')+labs(y='',x='')+scale_y_continuous(expand=c(0,0),limits = c(0,14),breaks=seq(0,14,4))+geom_boxplot(width=0.3,fill='gray',size=1.5,outlier.shape = NA)

