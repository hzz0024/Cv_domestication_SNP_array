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

#devtools::install_github('royfrancis/pophelper')
#install.packages('viridis')
library(pophelper)
library(viridis)
library(gridExtra)
library("wesanderson")
### Read in data------------------------------------------------------------

k3 <- readQ('./genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_neutral_pruned_5000_structure_output_K3_run_10_f',filetype='auto')
k9 <- readQ('./genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_neutral_pruned_5000_structure_output_K9_run_10_f',filetype='auto')
k10 <- readQ('./genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_neutral_pruned_5000_structure_output_K10_run_10_f',filetype='auto')

# Sample metadata
inds_grps <- read.table('sampleID_pop_n539.txt',header=F,stringsAsFactors = F)
grps <- as.data.frame(inds_grps[,2])
grps[,1] <- as.character(grps[,1])
colnames(grps) <- c('subs')

# Add sample IDs as rownames to qlist
rownames(k3[[1]]) <- inds_grps$V1
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
            showindlab = F, useindlab = T, indlabsize = 5, sharedindlab = F, sortind="all",
            indlabangle=45,indlabvjust=1,
            showgrplab=F, grplab=grps, grplabjust=0.9, grplabsize=6, grplabpos= 0.2, subsetgrp = pop_order, grplabangle = 45, linesize=1, pointsize=4,
            showdiv=T)
grid.arrange(p3$plot[[1]], nrow=1, widths=c(20,1))

col_gradient <- c(  "#849cc1", "#cf7fbc", "#FA4032", "#8ad5d9","#fb90a6", "#0A2C86", "#f9476b",
                    "#FA812F", "#FAAF08","#fddae1", "#1D92BD",  "#e2b2d6")
show_col(col_gradient)
p9 <- plotQ(k9, basesize=10, exportplot=F, returnplot=T,
            grplabspacer = -0.2,  showsp= F,  splabcol="white",
            clustercol=col_gradient,
            barsize = 1, barbordercolour="black",barbordersize=0.01,
            showyaxis=F,showticks=T, 
            showindlab = F, useindlab = T, indlabsize = 5, sharedindlab = F, sortind="all",
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
            showindlab = F, useindlab = T, indlabsize = 5, sharedindlab = F, sortind="all",
            indlabangle=45,indlabvjust=1,
            showgrplab=T, grplab=grps, grplabjust=0.9, grplabsize=6, grplabpos= 0.5, subsetgrp = pop_order, grplabangle = 45, linesize=1, pointsize=4,
            showdiv=T)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(export)

jpeg("Structure_509.jpg", width = 20, height = 10, units = 'in', res = 300)
plot_grid(p3$plot[[1]],p9$plot[[1]],p10$plot[[1]], nrow=3, rel_heights = c(1/4, 1/4, 1/3))
dev.off()


#grid.arrange(p7$plot[[1]],p9$plot[[1]],p14$plot[[1]], nrow=3, rel_heights = c(1/3, 1/3, 1/3))


