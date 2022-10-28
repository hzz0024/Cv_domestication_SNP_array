
library(adegenet)
library(adegenet)
library(vcfR)
library(export)
library(LDheatmap)
require(snpStats)
adegenetServer(what = "DAPC")

setwd("~/Dropbox/Mac/Documents/HG/Domestication/10_DAPC_outlier")

vcftools  = "/Users/HG/Dropbox/Mac/Documents/HG/Github/BioinfoTools/vcftools_0.1.13/bin/vcftools";
plink  = "/Users/HG/Dropbox/Mac/Documents/HG/Domestication/14_ROH/plink";

#system(paste(vcftools," --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.recode.vcf --bed pop_n_477_pcadapt_outflank.shared.outlier.igv_Dom_Wild.sliding.zfst.outlier.merged.igv.outlier.merged.igv --recode --recode-INFO-all --out n_509_shared_outliers", sep=""))
system(paste(vcftools," --vcf genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe.recode.vcf --snps Predictor_loci_classification_best_10perc_unique_37.txt --recode --recode-INFO-all --out n_539_shared_outliers", sep=""))

# VCFtools - v0.1.13
# (C) Adam Auton and Anthony Marcketta 2009
# 
# Parameters as interpreted:
#   --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.recode.vcf
# --recode-INFO-all
# --out n_509_shared_outliers
# --recode
# --bed pop_n_477_pcadapt_outflank.shared.10k.outlier.igv_Dom_Wild.sliding.zfst.outlier.merged.igv.outlier.merged.igv
# 
# After filtering, kept 509 out of 509 Individuals
# Outputting VCF file...
# Read 35 BED file entries.
# After filtering, kept 410 out of a possible 141960 Sites
# Run Time = 2.00 seconds

# load the population information
pop_info <- read.table("pop_539_sample_list.txt", header=TRUE, sep="\t", stringsAsFactors = TRUE)
pop_info$Pop_correct = factor(pop_info$Pop_correct, levels=c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "DBX1", "DBX2", "DBX3",  "UNC1", "UNC2", "UMFS", "NYH1", "NEH1", "NEH2", "MEH2"))
# load vcf file
vcf_file = "n_539_shared_outliers.recode.vcf"
#vcf_file = "genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_pruned_neutral_1K.recode.vcf"
vcf <- read.vcfR(vcf_file, verbose = FALSE)
Mydata1 <- vcfR2genind(vcf)
Mydata1@pop <- pop_info$Pop_correct
x <- Mydata1


########### LDheatmap ###########
# vcf <- read.vcfR(vcf_file)
# plink  = "/Users/HG/Dropbox/Mac/Documents/HG/Domestication/14_ROH/plink";
# system(paste(plink, " --vcf n_509_shared_outliers.recode.vcf --allow-extra-chr --make-bed --out n_509_shared_outliers", sep=""))

# eur_gt <- vcf@gt
# eur_snpMat <- t(eur_gt)           
# #define a function to convert the value into 0,1,2
# convertToNumeric <- function(x){
#   gdat <- matrix(NA,nrow = nrow(x), ncol = ncol(x))
#   for (m in 1:nrow(x)){
#     for (n in 1:ncol(x)){
#       a <-as.numeric(unlist(strsplit(x[m,n], "|"))[1]) 
#       
#       b <- as.numeric(unlist(strsplit(x[m,n], "|"))[3])
#       gdat[m,n] <- a+b
#     }
#   }
#   rownames(gdat) <- rownames(x)
#   colnames(gdat) <- colnames(x)
#   return(gdat)
# }
# #convert to snpMatrix - EUR
# gdat_eur <- convertToNumeric(eur_snpMat)
# require(snpStats)
# gdat_eur<-as(x,"SnpMatrix")
# library(LDheatmap)
# #> Warning: package 'LDheatmap' was built under R version 3.4.4
# gdat_eur<-as(gdat_eur,"SnpMatrix")
# #> Warning in asMethod(object): values other than 0, 1 or 2 set to NA
# LDheatmap(gdat_eur,info$filt_snpDist,title='Europeans',add.map=FALSE)

sum(is.na(Mydata1$tab))

# estimate the number of PCs for data analysis
dapc1 <- dapc(x, x@pop)
temp <- optim.a.score(dapc1)
library(export)
graph2ppt(file="Domestication_alpha_score1",width=6,height=5)

# temp$best = 26
x@pop
x@all.names
grp <- find.clusters(Mydata1, max.n.clust=16)
names(grp)
head(grp$grp, 10)
table(pop(x), grp$grp)
table.value(table(pop(x), grp$grp), col.lab=levels(pop(x)),
            row.lab=paste("ori", 1:17))
# formal running
dapc1 <- dapc(x, x@pop, n.pca = 7, n.da=7)
scatter(dapc1)
# adegenetServer(what = "DAPC")
# #reading the genepop file as input
# x <- read.genepop("WAL62_SNP1978.gen",quiet = TRUE)
# x <- read.genepop("WAL605_SNP68.gen",quiet = TRUE)
# x <- read.genepop("WAL547_SNP68.gen",quiet = TRUE)

#dapc1 <- dapc(x, var.contrib=FALSE, scale=FALSE, n.pca=50, n.da=1)

#plot
# cstar=0 will remove the connection lines
# add cell = 0  will remove the ellipses
# cex is the size of the dot
#SMB plot

#           #   MEW1       MEW2       LIW1        LIW2       DBW1      DBW2        NCW1      NCW2         
# col1 <- c("#0A2C86", "#325A98",  "#1D92BD", "#3DB9C1", "#C4E9B3", "#7BD2BF", "#ECF6B9", "#EEE8AA", 
#                #  DBX1       DBX2      DBX3       UNC1        UNC2       UMFS      NEH1       NEH2       MEH2
#                "#F9476B", "#FC709F","#E376B7", "#CF7FBC",  "#A36DC1", "#FEB22B", "#F36616", "#D83B1C", "#FF9117")

         #   MEW1       MEW2       LIW1        LIW2       DBW1      DBW2        NCW1      NCW2       DBX1        
col1 <- c("#325A98", "#325A98",  "#325A98", "#325A98", "#325A98", "#325A98", "#325A98", "#325A98", "#325A98",
          #   DBX2      DBX3       UNC1        UNC2       UMFS      NEH1     NYH1      NEH2       MEH2
           "#F62A00","#F62A00", "#F62A00",  "#F62A00", "#F62A00", "#F62A00", "#F62A00", "#F62A00", "#F62A00")

col_transparent <- adjustcolor(col1,alpha.f = 0.1)

scatter(dapc1,  xlab="Discriminatn Function 1", ylab="Discriminatn Function 2", cex = 0.7, pch=c(rep(17, 9), rep(20,9)), lwd=2, lty=2, solid = .8, cstar=0, scree.da=FALSE,  col = col1, clab = 1 ) #label.inds = list(air = 2, pch = NA))
#legend("right", legend = levels(pop(x)),cex = 1,col = col1, bty = 'n', pch=c(rep(17, 9), rep(20,8)), xpd = TRUE, inset=c(0.2,0.55))
myInset <- function(){
  temp <- dapc1$pca.eig
  temp <- 100* cumsum(temp)/sum(temp)
  plot(temp, col=rep(c("black","lightgrey"),
                     c(dapc1$n.pca,1000)), ylim=c(0,100),
       xlab="PCA axis", ylab="Cumulated variance (%)",
       cex=1, pch=15, type="h", lwd=3)
}
add.scatter(myInset(), posi="topleft",
            inset=c(0.05,0.00), ratio=.14,
            bg=transp("white"))

library(export)
graph2ppt(file="Domestication_outlier_DAPC",width=6,height=6)
####################################################

# outlier Fst

library(devtools)
install_github("jgx65/hierfstat")
library("hierfstat")
library("adegenet")
library(vcfR)
library(MASS)

# see https://popgen.nescent.org/DifferentiationSNP.html for detailed example
# read from a vcf file
setwd("~/Dropbox/Mac/Documents/HG/Domestication/10_DAPC_outlier")

# load the population information
pop_info <- read.table("pop_509_sample_list.txt", header=TRUE, sep="\t", stringsAsFactors = TRUE)
pop_info$Pop_correct = factor(pop_info$Pop_correct, levels=c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "DBX1", "DBX2", "DBX3", "UNC1", "UNC2", "UMFS", "NEH1", "NEH2", "MEH2"))
# load vcf file
vcftools  = "/Users/HG/Dropbox/Mac/Documents/HG/Github/BioinfoTools/vcftools_0.1.13/bin/vcftools";
system(paste(vcftools," --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.recode.vcf --snps Predictor_loci_classification_best_10perc_unique_54.txt --recode --recode-INFO-all --out RF_outlier", sep=""))

vcf_file = "RF_outlier.recode.vcf"
vcf <- read.vcfR(vcf_file, verbose = FALSE)
Mydata1 <- vcfR2genind(vcf)
Mydata1@pop <- pop_info$Pop_correct
Mydata1

# wc(Mydata1) # Weir and Cockerham's estimate
fst <- genet.dist(Mydata1, method = "WC84") # Pairwise Fst
write.matrix(fst, file = "pairwise.fst.txt", ,sep = "\t")

library(ggplot2)
library(ComplexHeatmap)
library(circlize)
dt = read.delim("pairwise.fst.txt", header = TRUE, sep='\t')
head(dt)

plot_dt <- as.matrix(dt)
rownames(plot_dt) <- c("DBX1", "DBX2", "DBX3",  "UNC1", "UNC2", "UMFS", "NEH1", "NEH2", "MEH2")
colnames(plot_dt) <- c("DBX1", "DBX2", "DBX3",  "UNC1", "UNC2", "UMFS", "NEH1", "NEH2", "MEH2")
plot_dt <- as.matrix(plot_dt)
class(plot_dt)

fa = rep(c("Recent", "Old"), times = c(9, 8))
fa_col = c("WILD" = "#aecdc2", "DOM" = "#f0b8b8")
dend1 = cluster_between_groups(plot_dt, fa)

Heatmap(plot_dt, name = "Pairwise Fst",
        col = colorRamp2(c(0, 0.21, 0.42), c("white", "skyblue", "orange")),
        rect_gp = gpar(col = "white", lwd = 2),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        width = unit(14, "cm"), height = unit(8, "cm"),
        row_names_gp = gpar(fontsize = 12),
        column_names_gp = gpar(fontsize = 12),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", plot_dt[i, j]), x, y, gp = gpar(fontsize = 10))
        })
        #top_annotation = HeatmapAnnotation(Group = fa, col = list(Group = fa_col)))

graph2ppt(file="Pairwise_fst",width=10,height=6)
