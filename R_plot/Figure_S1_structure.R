#############################
# Plotting STRUCTURE results 
#############################

### Goal: plot Q tables from STRUCTURE analyses across 
### cross all populations (best K = 2) above separate analysis 
### subsets  within the Gulf of Mexico (K = 3) and the Atlantic Coast (K = 3) 
### populations using 10K random neutral SNPs.

### Set working directory---------------------------------------------------

setwd("~/Dropbox/Mac/Documents/HG/Domestication/Manuscript/Figure_539/Figure_S1")

### Install and load libraries/dependencies---------------------------------


#install.packages('viridis')
library(pophelper)
library(viridis)
library(gridExtra)
library("wesanderson")
### Read in data------------------------------------------------------------

k3 <- readQ('./genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_neutral_pruned_5000_structure_output_K3_run_10_f',filetype='auto')
k4 <- readQ('./genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_neutral_pruned_5000_structure_output_K4_run_1_f',filetype='auto')
k5 <- readQ('./genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_neutral_pruned_5000_structure_output_K5_run_1_f',filetype='auto')
k6 <- readQ('./genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_neutral_pruned_5000_structure_output_K6_run_2_f',filetype='auto')
k7 <- readQ('./genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_neutral_pruned_5000_structure_output_K7_run_1_f',filetype='auto')
k8 <- readQ('./genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_neutral_pruned_5000_structure_output_K8_run_1_f',filetype='auto')
k9 <- readQ('./genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_neutral_pruned_5000_structure_output_K9_run_1_f',filetype='auto')
k10 <- readQ('./genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_neutral_pruned_5000_structure_output_K10_run_10_f',filetype='auto')

# Sample metadata
inds_grps <- read.table('sampleID_pop_n539.txt',header=F,stringsAsFactors = F)
grps <- as.data.frame(inds_grps[,2])
grps[,1] <- as.character(grps[,1])
colnames(grps) <- c('subs')

# Add sample IDs as rownames to qlist
rownames(k3[[1]]) <- inds_grps$V1
rownames(k4[[1]]) <- inds_grps$V1
rownames(k5[[1]]) <- inds_grps$V1
rownames(k6[[1]]) <- inds_grps$V1
rownames(k7[[1]]) <- inds_grps$V1
rownames(k8[[1]]) <- inds_grps$V1
rownames(k9[[1]]) <- inds_grps$V1
rownames(k10[[1]]) <- inds_grps$V1

### Plot CV error------------------------------------------------------------

# par(mfrow=c(1,1))
# plot(CVerror$K,CVerror$CV_error,pch=20,xlab='K',ylab='CV error')
# low <- min(CVerror$CV_error)
# abline(h=low,col='navy',lty=2,lwd=0.5)

### Plot ADMIXTURE plots-----------------------------------------------------

# set population order
pop_order <- c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "DBX1", "DBX2", "DBX3",  "UNC1", "UNC2", "UMFS", "NYH1","NEH1", "NEH2", "MEH2")

# col_gradient <- c( "#0A2C86", "#849cc1",  "#1D92BD", "#8ad5d9", "#93c47d", "#bedbb1", "#a9a9a9", "#dddddd",                    #  DBX1       DBX2      DBX3       UNC1        UNC2       UMFS      NEH1       NEH2       MEH2
#                    "#f9476b", "#fb90a6","#fddae1", "#cf7fbc",  "#e2b2d6", "#fec155", "#e1bb94", "#fbd0a5", "#b58383")
                                  # blue                        # yellow
col_gradient <- c("#FA812F", "#fb90a6",  "#849cc1")
# Groups individuals by subspecies abbreviation

p3 <- plotQ(k3, basesize=10, exportplot=F, returnplot=T,
            grplabspacer = -0.2,  showsp= F,  splabcol="white",
            clustercol=col_gradient,
            barsize = 1, barbordercolour="black",barbordersize=0.01,
            showyaxis=F,showticks=T, 
            showindlab = F, useindlab = T, indlabsize = 5, sharedindlab = F, 
            indlabangle=45,indlabvjust=1,
            showgrplab=F, grplab=grps, grplabjust=0.9, grplabsize=6, grplabpos= 0.2, subsetgrp = pop_order, grplabangle = 45, linesize=1, pointsize=4,
            showdiv=T)
grid.arrange(p3$plot[[1]], nrow=1, widths=c(20,1))

col_gradient <- c("#FA812F", "#fb90a6",  "#849cc1", "#cf7fbc")
# Groups individuals by subspecies abbreviation

p4 <- plotQ(k4, basesize=10, exportplot=F, returnplot=T,
            grplabspacer = -0.2,  showsp= F,  splabcol="white",
            clustercol=col_gradient,
            barsize = 1, barbordercolour="black",barbordersize=0.01,
            showyaxis=F,showticks=T, 
            showindlab = F, useindlab = T, indlabsize = 5, sharedindlab = F, 
            indlabangle=45,indlabvjust=1,
            showgrplab=F, grplab=grps, grplabjust=0.9, grplabsize=6, grplabpos= 0.2, subsetgrp = pop_order, grplabangle = 45, linesize=1, pointsize=4,
            showdiv=T)
grid.arrange(p4$plot[[1]], nrow=1, widths=c(20,1))

col_gradient <- c("#cf7fbc", "#FA812F",  "#fb90a6", "#849cc1","#f9476b")
p5 <- plotQ(k5, basesize=10, exportplot=F, returnplot=T,
            grplabspacer = -0.2,  showsp= F,  splabcol="white",
            clustercol=col_gradient,
            barsize = 1, barbordercolour="black",barbordersize=0.01,
            showyaxis=F,showticks=T, 
            showindlab = F, useindlab = T, indlabsize = 5, sharedindlab = F, 
            indlabangle=45,indlabvjust=1,
            showgrplab=F, grplab=grps, grplabjust=0.9, grplabsize=6, grplabpos= 0.2, subsetgrp = pop_order, grplabangle = 45, linesize=1, pointsize=4,
            showdiv=T)
grid.arrange(p5$plot[[1]], nrow=1, widths=c(20,1))

col_gradient <- c("#fb90a6", "#FA4032",  "#FAAF08", "#FA812F","#cf7fbc", "#849cc1")
p6 <- plotQ(k6, basesize=10, exportplot=F, returnplot=T,
            grplabspacer = -0.2,  showsp= F,  splabcol="white",
            clustercol=col_gradient,
            barsize = 1, barbordercolour="black",barbordersize=0.01,
            showyaxis=F,showticks=T, 
            showindlab = F, useindlab = T, indlabsize = 5, sharedindlab = F, 
            indlabangle=45,indlabvjust=1,
            showgrplab=F, grplab=grps, grplabjust=0.9, grplabsize=6, grplabpos= 0.2, subsetgrp = pop_order, grplabangle = 45, linesize=1, pointsize=4,
            showdiv=T)
grid.arrange(p6$plot[[1]], nrow=1, widths=c(20,1))

col_gradient <- c("#FA812F", "#cf7fbc",  "#FAAF08", "#0A2C86","#fb90a6", "#f9476b", "#849cc1")
p7 <- plotQ(k7, basesize=10, exportplot=F, returnplot=T,
            grplabspacer = -0.2,  showsp= F,  splabcol="white",
            clustercol=col_gradient,
            barsize = 1, barbordercolour="black",barbordersize=0.01,
            showyaxis=F,showticks=T, 
            showindlab = F, useindlab = T, indlabsize = 5, sharedindlab = F, 
            indlabangle=45,indlabvjust=1,
            showgrplab=F, grplab=grps, grplabjust=0.9, grplabsize=6, grplabpos= 0.2, subsetgrp = pop_order, grplabangle = 45, linesize=1, pointsize=4,
            showdiv=T)
grid.arrange(p7$plot[[1]], nrow=1, widths=c(20,1))

col_gradient <- c("#FAAF08", "#FA812F",  "#FA4032", "#0A2C86","#849cc1", "#cf7fbc", "#fb90a6", "#f9476b")
p8 <- plotQ(k8, basesize=10, exportplot=F, returnplot=T,
            grplabspacer = -0.2,  showsp= F,  splabcol="white",
            clustercol=col_gradient,
            barsize = 1, barbordercolour="black",barbordersize=0.01,
            showyaxis=F,showticks=T, 
            showindlab = F, useindlab = T, indlabsize = 5, sharedindlab = F, 
            indlabangle=45,indlabvjust=1,
            showgrplab=F, grplab=grps, grplabjust=0.9, grplabsize=6, grplabpos= 0.2, subsetgrp = pop_order, grplabangle = 45, linesize=1, pointsize=4,
            showdiv=T)
grid.arrange(p8$plot[[1]], nrow=1, widths=c(20,1))

col_gradient <- c(  "#849cc1", "#cf7fbc", "#FA4032", "#8ad5d9","#fb90a6", "#0A2C86", "#f9476b",
                    "#FA812F", "#FAAF08","#fddae1", "#1D92BD",  "#e2b2d6")
show_col(col_gradient)
p9 <- plotQ(k9, basesize=10, exportplot=F, returnplot=T,
            grplabspacer = -0.2,  showsp= F,  splabcol="white",
            clustercol=col_gradient,
            barsize = 1, barbordercolour="black",barbordersize=0.01,
            showyaxis=F,showticks=T, 
            showindlab = F, useindlab = T, indlabsize = 5, sharedindlab = T,
            indlabangle=45,indlabvjust=1,
            showgrplab=F, grplab=grps, grplabjust=0.9, grplabsize=6, grplabpos= 0.2, subsetgrp = pop_order, grplabangle = 45, linesize=1, pointsize=4,
            showdiv=T)
grid.arrange(p9$plot[[1]], nrow=1, widths=c(20,1))

col_gradient <- c(  "#0A2C86", "#FAAF08", "#f9476b", "#8ad5d9","#849cc1", "#FA812F", "#FA4032",
                    "#fb90a6", "#cf7fbc","#fddae1", "#1D92BD",  "#e2b2d6")

p10 <- plotQ(k10, basesize=10, exportplot=F, returnplot=T,
            grplabspacer = -0.2,  showsp= F,  splabcol="white",
            clustercol=col_gradient,
            barsize = 1, barbordercolour="black",barbordersize=0.01,
            showyaxis=F,showticks=T, 
            showindlab = F, useindlab = T, indlabsize = 5, sharedindlab = F, sortind="Cluster9",
            indlabangle=45,indlabvjust=1,
            showgrplab=T, grplab=grps, grplabjust=0.8, grplabsize=8, grplabpos= 0.1, subsetgrp = pop_order, grplabangle = 45, linesize=1, pointsize=4,
            showdiv=T)
grid.arrange(p10$plot[[1]], nrow=1, widths=c(20,1))

library(ggplot2)
library(gridExtra)
library(cowplot)
library(export)

jpeg("Structure_539.jpg", width = 20, height = 10, units = 'in', res = 300)
plot_grid(p3$plot[[1]], p4$plot[[1]], p5$plot[[1]], p6$plot[[1]], p7$plot[[1]], p8$plot[[1]], p9$plot[[1]],p10$plot[[1]], nrow=8, rel_heights = c(1/9, 1/9, 1/9, 1/9, 1/9, 1/9, 1/9, 1/8)) +
  theme(plot.margin = unit(c(1,1,2,1), "cm"))
dev.off()

col_gradient <- c("#fb90a6", "#FA4032",  "#FAAF08", "#FA812F","#cf7fbc", "#849cc1")
p6 <- plotQ(k6, basesize=10, exportplot=F, returnplot=T,
            grplabspacer = -0.8,  showsp= F,  splabcol="white",
            clustercol=col_gradient,
            barsize = 1, barbordercolour=NA,barbordersize=0.0001,
            showyaxis=F,showticks=T, 
            showindlab = F, useindlab = T, indlabsize = 5, sharedindlab = F, 
            indlabangle=45,indlabvjust=1,
            showgrplab=T, grplab=grps, grplabjust=0.5, grplabsize=6, grplabpos= 0.1, subsetgrp = pop_order, grplabangle = 45, linesize=1, pointsize=4,
            showdiv=T)
grid.arrange(p6$plot[[1]], nrow=3, widths=c(20,1))
graph2ppt(file="K6",width=8,height=6)
# 
# k91 <- readQ('./genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_neutral_pruned_5000_structure_output_K9_run_1_f',indlabfromfile=TRUE, filetype='auto')
# k92 <- readQ('./genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_neutral_pruned_5000_structure_output_K9_run_2_f',indlabfromfile=TRUE, filetype='auto')
# 
# # Sample metadata
# inds_grps <- read.table('sampleID_pop_n539.txt',header=F,stringsAsFactors = F)
# grps <- as.data.frame(inds_grps[,2])
# grps[,1] <- as.character(grps[,1])
# colnames(grps) <- c('subs')
# 
# clustercol=c(  "#FAAF08", "#849cc1","#0A2C86", "#f9476b", "#8ad5d9", "#FA812F", "#FA4032",
#                "#fb90a6", "#cf7fbc","#fddae1", "#1D92BD",  "#e2b2d6")
# show_col(col_gradient)
# #p_all <- plotQ(alignK(c(k8, k91, k92)),imgoutput="join",returnplot=T,exportplot=F,basesize=11)
# p_all <- plotQ(alignK(c(k9, k10)),imgoutput="join",
#                clustercol=c(  "#FAAF08", "#849cc1","#0A2C86", "#f9476b", "#8ad5d9", "#FA812F", "#FA4032",
#                               "#fb90a6", "#cf7fbc","#fddae1", "#1D92BD",  "#e2b2d6"), splab=paste0("K=",c(9:10)),
#                basesize=10, exportplot=F, returnplot=T,
#                grplabspacer = -0.2,  showsp= F,  splabcol="white",barsize = 1, barbordercolour="black",barbordersize=0.01,
#                showyaxis=F,showticks=T, 
#                showindlab = F, useindlab = T, indlabsize = 5, sharedindlab = F, 
#                indlabangle=45,indlabvjust=1,
#                showgrplab=T, grplab=grps, grplabjust=0.8, grplabsize=8, grplabpos= 0.1, subsetgrp = pop_order, grplabangle = 45, linesize=1, pointsize=4,
#                showdiv=T)
# grid.arrange(p_all$plot[[1]])

#grid.arrange(p7$plot[[1]],p9$plot[[1]],p14$plot[[1]], nrow=3, rel_heights = c(1/3, 1/3, 1/3))


