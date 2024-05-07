
library(adegenet)
library(vcfR)
library(export)
library(LDheatmap)
require(snpStats)

setwd("~/Dropbox/Mac/Documents/HG/Domestication/05_pairwise_fst/union_1174")

vcftools  = "/Users/HG/Dropbox/Mac/Documents/HG/Github/BioinfoTools/vcftools_0.1.13/bin/vcftools";
plink  = "/Users/HG/Dropbox/Mac/Documents/HG/Domestication/14_ROH/plink";

system(paste(vcftools," --vcf genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe.recode.vcf --snps union_outliers_1174.list --recode --recode-INFO-all --out n_539_union_outliers", sep=""))
pop_info <- read.table("pop_539_sample_list.txt", header=TRUE, sep="\t", stringsAsFactors = TRUE)
pop_info$Pop_correct = factor(pop_info$Pop_correct, levels=c("Wild", "Selected"))

system(paste(vcftools," --vcf genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe.recode.vcf --snps Predictor_loci_classification_best_10perc_unique_37.txt --recode --recode-INFO-all --out n_539_union_outliers", sep=""))
pop_info <- read.table("pop_539_sample_ind.txt", header=TRUE, sep="\t", stringsAsFactors = TRUE)
pop_info$Pop_correct = factor(pop_info$Pop_correct, levels=c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "DBX1", "DBX2", "DBX3",  "UNC1", "UNC2", "UMFS", "NYH1", "NEH1", "NEH2", "MEH2"))

# load vcf file
vcf_file = "n_539_union_outliers.recode.vcf"
#vcf_file = "genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe.recode.vcf"
vcf <- read.vcfR(vcf_file, verbose = FALSE)
Mydata1 <- vcfR2genind(vcf)
Mydata1@pop <- pop_info$Pop_correct
Mydata1

# outlier Fst

library(devtools)
#install_github("jgx65/hierfstat")
library("hierfstat")
library("adegenet")
library(vcfR)
library(MASS)

fst <- genet.dist(Mydata1, method = "WC84") # Pairwise Fst
#mean_fst <- mean(fst[lower.tri(fst, diag = FALSE)], na.rm = TRUE)
write.matrix(fst, file = "pairwise.fst.union1174.pairwise.txt", ,sep = "\t")
write.matrix(fst, file = "pairwise.fst.full.txt", ,sep = "\t")
write.matrix(fst, file = "pairwise.fst.RF37.txt", ,sep = "\t")
write.matrix(fst, file = "pairwise.fst.union1174.txt", ,sep = "\t")


