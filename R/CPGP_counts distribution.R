#data should be a longform dataframe (use tidyr gather function giving key, value variables)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(R.matlab)
library(emmeans)

counts=readMat("you location/cts_distribution.mat")

vta_g <- as.vector(counts$vta.GP)
vta_l <- as.vector(counts$vta.Lat)
vta_dm <- as.vector(counts$vta.dMed)
vta_vm <- as.vector(counts$vta.vMed)
vta_t <- as.vector(counts$vta.t)
vta_r <- as.vector(counts$vta.r)
vta_i <- as.vector(counts$vta.i)
vta_c <- as.vector(counts$vta.c)

snc_g <- as.vector(counts$snc.GP)
snc_l <- as.vector(counts$snc.Lat)
snc_dm <- as.vector(counts$snc.dMed)
snc_vm <- as.vector(counts$snc.vMed)
snc_t <- as.vector(counts$snc.t)
snc_r <- as.vector(counts$snc.r)
snc_i <- as.vector(counts$snc.i)
snc_c <- as.vector(counts$snc.c)

    
# Recounstruction of counts dataframe
vta_g_Counts <- vta_g
vta_g_Region <- rep('GP',length(vta_g_Counts))
vta_g_Source <- rep('VTA',length(vta_g_Counts))

vta_l_Counts <- vta_l
vta_l_Region <- rep('LatCP',length(vta_l_Counts))
vta_l_Source <- rep('VTA',length(vta_l_Counts))

vta_dm_Counts <- vta_dm
vta_dm_Region <- rep('dMedCP',length(vta_dm_Counts))
vta_dm_Source <- rep('VTA',length(vta_dm_Counts))

vta_vm_Counts <- vta_vm
vta_vm_Region <- rep('vMedCP',length(vta_vm_Counts))
vta_vm_Source <- rep('VTA',length(vta_vm_Counts))

vta_t_Counts <- vta_t
vta_t_Region <- rep('CPt',length(vta_t_Counts))
vta_t_Source <- rep('VTA',length(vta_t_Counts))

vta_r_Counts <- vta_r
vta_r_Region <- rep('CPr',length(vta_r_Counts))
vta_r_Source <- rep('VTA',length(vta_r_Counts))

vta_i_Counts <- vta_i
vta_i_Region <- rep('CPi',length(vta_i_Counts))
vta_i_Source <- rep('VTA',length(vta_i_Counts))

vta_c_Counts <- vta_c
vta_c_Region <- rep('CPc',length(vta_c_Counts))
vta_c_Source <- rep('VTA',length(vta_c_Counts))

snc_g_Counts <- snc_g
snc_g_Region <- rep('GP',length(snc_g_Counts))
snc_g_Source <- rep('SNc',length(snc_g_Counts))

snc_l_Counts <- snc_l
snc_l_Region <- rep('LatCP',length(snc_l_Counts))
snc_l_Source <- rep('SNc',length(snc_l_Counts))

snc_dm_Counts <- snc_dm
snc_dm_Region <- rep('dMedCP',length(snc_dm_Counts))
snc_dm_Source <- rep('SNc',length(snc_dm_Counts))

snc_vm_Counts <- snc_vm
snc_vm_Region <- rep('vMedCP',length(snc_vm_Counts))
snc_vm_Source <- rep('SNc',length(snc_vm_Counts))

snc_t_Counts <- snc_t
snc_t_Region <- rep('CPt',length(snc_t_Counts))
snc_t_Source <- rep('SNc',length(snc_t_Counts))

snc_r_Counts <- snc_r
snc_r_Region <- rep('CPr',length(snc_r_Counts))
snc_r_Source <- rep('SNc',length(snc_r_Counts))

snc_i_Counts <- snc_i
snc_i_Region <- rep('CPi',length(snc_i_Counts))
snc_i_Source <- rep('SNc',length(snc_i_Counts))

snc_c_Counts <- snc_c
snc_c_Region <- rep('CPc',length(snc_c_Counts))
snc_c_Source <- rep('SNc',length(snc_c_Counts))


data=data.frame(Counts=c(vta_g_Counts,vta_l_Counts,vta_dm_Counts,vta_vm_Counts,vta_r_Counts,vta_i_Counts,vta_c_Counts,vta_t_Counts,snc_g_Counts,snc_l_Counts,snc_dm_Counts,snc_vm_Counts,snc_r_Counts,snc_i_Counts,snc_c_Counts,snc_t_Counts),Region=c(vta_g_Region,vta_l_Region,vta_dm_Region,vta_vm_Region,vta_r_Region,vta_i_Region,vta_c_Region,vta_t_Region,snc_g_Region,snc_l_Region,snc_dm_Region,snc_vm_Region,snc_r_Region,snc_i_Region,snc_c_Region,snc_t_Region),Source=c(vta_g_Source,vta_l_Source,vta_dm_Source,vta_vm_Source,vta_r_Source,vta_i_Source,vta_c_Source,vta_t_Source,snc_g_Source,snc_l_Source,snc_dm_Source,snc_vm_Source,snc_r_Source,snc_i_Source,snc_c_Source,snc_t_Source))


ggplot(data, aes(x = factor(Region,levels=c('GP','LatCP','dMedCP','vMedCP','CPr','CPi','CPc','CPt')), y = Counts, fill = Source)) +
  geom_boxplot(width = 0.5, position = position_dodge(0.6), outlier.shape = NA) +  # Boxplot 추가
  stat_summary(fun = mean, geom = "point", shape = 19, size = 3, color = "black", position = position_dodge(0.6)) + # Mean 추가
  scale_fill_manual(values = c(rgb(0.9572,0.3062,0.0237), rgb(0.2947,0.3929,0.5251))) + 
  theme_minimal() +
  labs(x = "Region", y = "Number of targets")+scale_y_continuous(expand=c(0,0),breaks = seq(0,18,2))+theme_classic()+theme(plot.margin=unit(c(1,1,1,1),'cm'))+theme(axis.title = element_text(size=30),axis.text = element_text(size=25,color = 'black'),axis.text.x=element_text(angle=45,hjust=1),legend.position='none')+coord_fixed(ratio = 0.3)

result_g <- wilcox.test(snc_g_Counts,vta_g_Counts)
result_l <- wilcox.test(snc_l_Counts,vta_l_Counts)
result_dm <- wilcox.test(snc_dm_Counts,vta_dm_Counts)
result_vm <- wilcox.test(snc_vm_Counts,vta_vm_Counts)
result_r <- wilcox.test(snc_r_Counts,vta_r_Counts)
result_i <- wilcox.test(snc_i_Counts,vta_i_Counts)
result_c <- wilcox.test(snc_c_Counts,vta_c_Counts)
result_t <- wilcox.test(snc_t_Counts,vta_t_Counts)



#Proportion of number of targets across animals
imported=readMat("your location/eachanimal_targetproportion.mat")

rbind(as.numeric(imported[[4]]),as.numeric(imported[[6]]),as.numeric(imported[[2]]),as.numeric(imported[[3]]),as.numeric(imported[[5]]),as.numeric(imported[[1]]))


proportion <- data.frame(source=c('VTA','VTA','VTA','SNc','SNc','SNc'),rbind(as.numeric(imported[[4]]),as.numeric(imported[[6]]),as.numeric(imported[[2]]),as.numeric(imported[[3]]),as.numeric(imported[[5]]),as.numeric(imported[[1]])))

colnames(proportion) <- c('source','1','2','3','4','5','6','>6')

long_data <- proportion %>%
  pivot_longer(cols = starts_with(c('1','2','3','4','5','6','>6')),
               names_to = "xpoint",
               values_to = "value")

summary_data <- long_data %>%
  group_by(source, xpoint) %>%
  summarise(mean = mean(value), 
            se = sd(value)/sqrt(n()), .groups = "drop")

ggplot(summary_data, aes(x = factor(xpoint,levels=c('1','2','3','4','5','6','>6')), y = mean, group = source, color = factor(source))) +
  geom_line(size=2) +
  geom_point() +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2,size=1) +
  labs(x = "X Point", y = "Mean Proportion", color = "Source")+theme(axis.text = element_text(color ='black'))+theme(plot.margin=unit(c(1,1,1,1),'cm'))+labs(x='',y='')+scale_color_manual(values=c(rgb(0.9572,0.3062,0.0237), rgb(0.2947,0.3929,0.5251)))+theme_classic()+theme(legend.position = 'none')+theme(axis.text.y=element_text(size=25),axis.text.x=element_text(size=35),axis.text = element_text(color ='black'))+theme(plot.margin=unit(c(1,1,1,1),'cm'))+coord_fixed(ratio = 8)


anova_results <- aov(value ~ source * xpoint, data=long_data)
summary(anova_results)



emm <- emmeans(anova_results, pairwise ~ source | xpoint)

# Bonferroni
posthoc_result <- summary(emm$contrasts, adjust = "bonferroni")
print(posthoc_result)

