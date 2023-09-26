install.packages("LDcorSV")
library(LDcorSV)
library(tidyverse)

############# example data ###############
data(data.test)
# Allelic doses of 20 markers on a chromosome of 91 Vitis vinifera plants
Geno <- data.test[[1]]
Geno
# Kinship matrix of 91 plants of Vitis vinifera
V.WAIS <- data.test[[2]]
V.WAIS
# Structure population matrix of 91 plants of Vitis vinifera
# in two sub-populations
S.2POP <- data.test[[3]]
S.2POP

data(data.test)
Geno <- data.test[[1]]
info <- apply(Geno, 2, Info.Locus)
info

data(data.test)
V.WAIS <- data.test[[2]]
Inv.V.WAIS <- Inv.proj.matrix.sdp(V.WAIS)
Inv.V.WAIS

data(data.test)
Geno <- data.test[[1]]
V.WAIS <- data.test[[2]]
S.2POP <- data.test[[3]]
LD <- LD.Measures(Geno, V = V.WAIS, S = S.2POP, data = "G", supinfo = TRUE, na.presence = TRUE)
LD

###########################
S.2POP_sup <- 1-S.2POP$Q2pop.1
S.2POP_sup <- as.data.frame(S.2POP_sup)
rownames(S.2POP_sup) <- rownames(S.2POP) 
LD1 <- LD.Measures(Geno, V = V.WAIS, S = S.2POP, data = "G", supinfo = TRUE, na.presence = TRUE)
LD2 <- LD.Measures(Geno, V = V.WAIS, S = S.2POP_sup, data = "G", supinfo = TRUE, na.presence = TRUE)
# it does not matter which column is deleted, only the choice of K matters.


library(dplyr)
library(tidyr)
library(reshape)
setwd("~/Dropbox/Mac/Documents/HG/Domestication/24_LDcorSV")
NYH1 <- read.delim("Empirical.relatedness.NYH1.txt", header = T, sep=' ')
# NYH1_no_sort <- separate(data = NYH1, col = Ind, into = c("left", "right"), sep = "_")
# NYH1_after_sort <- NYH1_no_sort %>% 
#   add_count(left, sort = TRUE)
# NYH1_after_sort <- NYH1_after_sort[,1:3]
# #NYH1_after_sort <- NYH1_after_sort[,1:3][with(NYH1_after_sort[,1:3], order(left, right)), ]
# NYH1_matrix <- cast(NYH1_after_sort, right ~ left)
# xtabs(ritland ~ left + right, data=NYH1_after_sort)

NYH1 %>% 
  #Separate variable Comparison in two columns
  separate(col = Ind,into = c("var1","var2"), sep="_") -> NYH1_sep

x <- unique(NYH1_sep$var1)  

NYH1 %>% 
  #Separate variable Comparison in two columns
  separate(col = Ind,into = c("var1","var2"), sep="_") %>% 
  #Create temporary data.frame 
  {. ->> temp} %>% 
  #Stack temporary data.frame so we have both A1-B1 and B1-A1
  bind_rows(
    temp %>%
      dplyr::rename(var1 = var2,var2 = var1)
  ) %>% 
  #Join with a combination of all levels to make a "complete matrix"
  full_join(
    expand_grid(var1 = x,var2 = x)
  ) %>% 
  mutate(
    #Rounding p-value
    relat = round(ritland,7)) -> NYH1_plot

mat <- cast(NYH1_plot, var1 ~ var2)
V_mat = mat[,-1]
V_mat[is.na(V_mat)] <- 1
rownames(V_mat) = colnames(V_mat)
######################################
## Prepare dataset for Genotype ######
######################################
setwd("~/Dropbox/Mac/Documents/HG/Domestication/24_LDcorSV")
G_file <- read.delim("chr0.NYH1.vcf.012", header = F, sep='\t')
G_df <- G_file[,-1]
snp_file <- read.delim("chr0.NYH1.vcf.012.pos", header = F, sep='\t')
snp_id <- paste0(snp_file$V1, "_", snp_file$V2)
ind_file <- read.delim("chr0.NYH1.vcf.012.indv", header = F, sep='\t')
ind_id <- ind_file$V1
G_df[G_df == "-1"] <- NA
colnames(G_df) = snp_id
rownames(G_df) = ind_id

G_df <- G_df[order(match(rownames(G_df),colnames(V_mat))),]

#######################################
## Prepare dataset for STRUCTURE ######
#######################################
library(dplyr)
setwd("~/Dropbox/Mac/Documents/HG/Domestication/24_LDcorSV")
S_file <- read.delim("S.NYH1.txt", header = F, sep='\t')
S_file %>% arrange(factor(V1, levels = ind_id)) -> S_df
S_df <- S_df[, -1:-2]
rownames(S_df) = ind_id

S_df <- S_df[order(match(rownames(S_df),colnames(V_mat))),]
########################
## LD measurement ######
########################
G_df = G_df[,1:30]
LD <- LD.Measures(G_df, V = V_mat, S = S_df, data = "G", supinfo = TRUE, na.presence = TRUE)
LD

LD <- LD.Measures(Geno, V = V.WAIS, S = S.2POP, data = "G", supinfo = TRUE, na.presence = TRUE)
LD


