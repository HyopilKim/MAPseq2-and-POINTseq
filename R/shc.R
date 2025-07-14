


library(R.matlab)
library(ggplot2)
library(writexl)

imported <- readMat("your location/forSHC.mat")

data <- as.matrix(imported$scaled.logdata)


#1. using sigclust2 shc (kimes et al., 2017)
shc_result <- shc(data, metric="euclidean", linkage="ward.D2", n_sim = 500, icovest = 2, alpha=0.05)

summary(shc_result)

plot.shc(shc_result)

write.csv(shc_result$p_norm,'p_values.csv')
write.csv(shc_result$nd_type,'node_types.csv')
write.csv(shc_result$hc_dat$height,'Z.csv')
