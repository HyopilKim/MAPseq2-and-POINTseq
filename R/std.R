
library(ggplot2)
library(tidyverse)
library(R.matlab)

imported=readMat("your location/std.mat")

snc=as.vector(imported$snc)
vta=as.vector(imported$vta)

snc_labels <- rep("snc", length(snc))
vta_labels <- rep("vta", length(vta))

data1=data.frame(values=c(snc,vta),source=c(snc_labels,vta_labels))

data1

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


data2 <- summarySE(data1,measurevar = 'values',groupvars='source')

ggplot(data2,aes(x=factor(source,levels=c('snc','vta'),labels = c('SNc','VTA')),y=values,fill=source))+geom_bar(stat='identity')+geom_errorbar(aes(ymin=values-se,ymax=values+se),width=0.2)+theme_classic()+theme(aspect.ratio=1.5)+theme(axis.title=element_text(size=40))+theme(axis.text.x=element_text(size=35),axis.text.y=element_text(size=25))+scale_y_continuous(expand=c(0,0))+theme(axis.text = element_text(color ='black'))+theme(plot.margin=unit(c(1,1,1,1),'cm'))+labs(x=NULL,y='SD')+coord_cartesian(ylim = c(0, 0.25))+scale_fill_manual(values=c(rgb(0.9572,0.3062,0.0237), rgb(0.2947,0.3929,0.5251)))+theme(legend.position = "none")


result <- t.test(snc, vta)
summary(result)

