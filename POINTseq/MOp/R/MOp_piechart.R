

library(R.matlab)
library(ggplot2)
library(tidyr)
library(dplyr)

#import forpiechart.mat
imported=readMat("your location/forpiechart.mat")
ident <- imported$identsorted

#Out of 9clusters, select 6clusters of 'ITc STR-','ITi STR-','ITc STR-','ITi STR-','ET','CT'

idx <- imported$Tsorted
idx <- factor(idx, levels=c(1,4,6,2,3,5,7,8,9),labels=c('ITc STR-','ITc STR+','ITc STR+','ITi STR-','ITi STR-','ITi STR+','ET','ET','CT'))

colors <- apply(imported$clusterColors[c(1,4,2,5,7,9),], 1, function(row) rgb(row[1], row[2], row[3]))
  
names = c('ITc STR-','ITc STR+','ITi STR-','ITi STR+','ET','CT')

#MAPseq
MAPseq <- idx[ident==1]

MAPseq_pie <- data.frame(idx=names, value=as.vector(table(MAPseq)),byrow = TRUE)

MAPseq_pie <- MAPseq_pie %>%
  mutate(fraction = value / sum(value),
         ymax = cumsum(fraction),
         ymin = c(0, head(ymax, n = -1)))

MAPseq_inner <- MAPseq_pie[,c('idx','ymax','ymin')]

MAPseq_outer <- data.frame(outeridx=MAPseq_pie$idx,ymax=ifelse(MAPseq_pie$idx%in%names,MAPseq_pie$ymax[MAPseq_pie$idx%in%names],0),ymin=ifelse(MAPseq_pie$idx%in%names,MAPseq_pie$ymin[MAPseq_pie$idx%in%names],0))

ggplot()+
  geom_rect(data=MAPseq_inner, aes(ymax = ymax, ymin = ymin , xmax = 3.4, xmin = 0, fill = factor(idx,levels = names))) +
  coord_polar(theta = "y") +
  xlim(c(0, 4)) +
  theme_void() +
  labs(fill = "") +scale_fill_manual(values = colors)+theme(legend.position = '')


#cux2
cux <- idx[ident==2]
cuxexpec <- c('ITc STR-','ITi STR-','ITc STR+','ITi STR+')
cux_pie <- data.frame(idx=names, value=as.vector(table(cux)),byrow = TRUE)

cux_pie <- cux_pie %>%
  mutate(fraction = value / sum(value),
         ymax = cumsum(fraction),
         ymin = c(0, head(ymax, n = -1)))

cux_inner <- cux_pie[,c('idx','ymax','ymin')]

cux_outer <- data.frame(outeridx=cux_pie$idx,ymax=0,ymin=0)
cux_outer$ymax[cux_pie$idx%in%cuxexpec] <- cux_pie$ymax[cux_pie$idx%in%cuxexpec]
cux_outer$ymin[cux_pie$idx%in%cuxexpec] <- cux_pie$ymin[cux_pie$idx%in%cuxexpec]

ggplot()+
  geom_rect(data=cux_outer, aes(ymax = ymax, ymin = ymin , xmax = 4, xmin = 3.3),fill='black')+
  geom_rect(data=cux_inner, aes(ymax = ymax, ymin = ymin , xmax = 3.3, xmin = 0, fill = factor(idx,levels = names))) +
  coord_polar(theta = "y") +
  xlim(c(0, 4)) +
  theme_void() +
  labs(fill = "") +scale_fill_manual(values = colors)+theme(legend.position = '')



#rbp4
rbp <- idx[ident==3]
rbpexpec <- c('ITc STR-','ITi STR-','ITc STR+','ITi STR+','ET')
rbp_pie <- data.frame(idx=names, value=as.vector(table(rbp)),byrow = TRUE)

rbp_pie <- rbp_pie %>%
    mutate(fraction = value / sum(value),
         ymax = cumsum(fraction),
         ymin = c(0, head(ymax, n = -1)))

rbp_inner <- rbp_pie[,c('idx','ymax','ymin')]

rbp_outer <- data.frame(outeridx=rbp_pie$idx,ymax=0,ymin=0)
rbp_outer$ymax[rbp_pie$idx%in%rbpexpec] <- rbp_pie$ymax[rbp_pie$idx%in%rbpexpec]
rbp_outer$ymin[rbp_pie$idx%in%rbpexpec] <- rbp_pie$ymin[rbp_pie$idx%in%rbpexpec]

ggplot()+
  geom_rect(data=rbp_outer, aes(ymax = ymax, ymin = ymin , xmax = 4, xmin = 3.3),fill='black')+
  geom_rect(data=rbp_inner, aes(ymax = ymax, ymin = ymin , xmax = 3.3, xmin = 0, fill = factor(idx,levels = names))) +
  coord_polar(theta = "y") +
  xlim(c(0, 4)) +
  theme_void() +
  labs(fill = "") +scale_fill_manual(values = colors)+theme(legend.position = '')


#retro
retro <- idx[ident==4]
retroexpec <- c('ITc STR-','ITc STR+')

retro_pie <- data.frame(idx=names, value=as.vector(table(retro)),byrow = TRUE)

retro_pie <- retro_pie %>%
  mutate(fraction = value / sum(value),
         ymax = cumsum(fraction),
         ymin = c(0, head(ymax, n = -1)))

retro_inner <- retro_pie[,c('idx','ymax','ymin')]

retro_outer <- data.frame(outeridx=retro_pie$idx,ymax=0,ymin=0)
retro_outer$ymax[retro_pie$idx%in%retroexpec] <- retro_pie$ymax[retro_pie$idx%in%retroexpec]
retro_outer$ymin[retro_pie$idx%in%retroexpec] <- retro_pie$ymin[retro_pie$idx%in%retroexpec]

ggplot()+
  geom_rect(data=retro_outer, aes(ymax = ymax, ymin = ymin , xmax = 4, xmin = 3.3),fill='black')+
  geom_rect(data=retro_inner, aes(ymax = ymax, ymin = ymin , xmax = 3.3, xmin = 0, fill = factor(idx,levels = names))) +
  coord_polar(theta = "y") +
  xlim(c(0, 4)) +
  theme_void() +
  labs(fill = "") +scale_fill_manual(values = colors)+theme(legend.position = '')

