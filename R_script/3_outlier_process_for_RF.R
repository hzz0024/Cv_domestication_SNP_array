setwd("~/Dropbox/Mac/Documents/HG/Domestication/18_outlier_shared_genotype_heatmap")
rm(list=ls())
#########################
# Check shared outliers #
#########################
library(bigsnpr)
library(stringr)
library(dplyr)
vcftools  = "/Users/HG/Dropbox/Mac/Documents/HG/Github/BioinfoTools/vcftools_0.1.13/bin/vcftools";
plink  = "/Users/HG/Dropbox/Mac/Documents/HG/Domestication/14_ROH/plink";

## all populations
name1 = "pop_n_539_outflank_BP_q05_n_327.bed"
outflank_BP <- read.delim(name1,header = F)
dim(outflank_BP)
name2 = "pop_n_539_PCAdapt_BP_q05_n_918.bed"
PACadapt_BP <- read.delim(name2,header = F)
dim(PACadapt_BP)

# how many intersect snps
length(intersect(outflank_BP$V4, PACadapt_BP$V4))
# 71

### union the outlier candidates

outliers_union <- unique(c(outflank_BP$V4,PACadapt_BP$V4))
# 1174

all_list <- str_split(outliers_union, "_")
chr <- unlist(all_list)[2*(1:length(outliers_union))-1]
pos <- unlist(all_list)[2*(1:length(outliers_union))  ]
df <- data.frame(chr,pos)
df <- as.data.frame(df)
bed_format <- function(df){
  head(df)
  df %>% 
    mutate(chr = str_replace(chr, "1\\b", "NC_035780.1")) %>% 
    mutate(chr = str_replace(chr, "2\\b", "NC_035781.1")) %>% 
    mutate(chr = str_replace(chr, "3\\b", "NC_035782.1")) %>% 
    mutate(chr = str_replace(chr, "4\\b", "NC_035783.1")) %>% 
    mutate(chr = str_replace(chr, "5\\b", "NC_035784.1")) %>% 
    mutate(chr = str_replace(chr, "6\\b", "NC_035785.1")) %>% 
    mutate(chr = str_replace(chr, "7\\b", "NC_035786.1")) %>% 
    mutate(chr = str_replace(chr, "8\\b", "NC_035787.1")) %>%
    mutate(chr = str_replace(chr, "9\\b", "NC_035788.1")) %>% 
    mutate(chr = str_replace(chr, "10\\b", "NC_035789.1"))  -> df
  final_df = data.frame(df$chr, df$pos, df$pos)
  return(final_df)
}
union_list <- bed_format(df)
write.table(outliers_union, file = "union_outliers_1174.list", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(union_list, file = "union_outliers_1174.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

### exclude the union outlier candidates
setwd("~/Dropbox/Mac/Documents/HG/Domestication/00_vcf")
system(paste(vcftools," --vcf genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe.recode.vcf --exclude-bed union_outliers_1174.bed --recode --recode-INFO-all --out genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_neutral", sep=""))
# VCFtools - v0.1.13
# (C) Adam Auton and Anthony Marcketta 2009
# 
# Parameters as interpreted:
#   --vcf genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe.recode.vcf
# --recode-INFO-all
# --out genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_neutral
# --recode
# --exclude-bed union_outliers_1174.bed
# 
# After filtering, kept 539 out of 539 Individuals
# Outputting VCF file...
# Read 1175 BED file entries.
# After filtering, kept 140502 out of a possible 141676 Sites
# Run Time = 22.00 seconds

system(paste(plink, " --vcf genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_neutral.recode.vcf --allow-extra-chr --make-bed --out genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_neutral", sep=""))

### produce LD-clumping dataset from neutral data 140,502 SNPs

sub_name="genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_neutral"
  
f_bk = paste0(sub_name, ".bk")
if (file.exists(f_bk)) {
  #Delete file if it exists
  file.remove(f_bk)
}

# part 1 SNP clumpping and data preparation
snp_readBed(paste0(sub_name, ".bed"))
# this will create a .rds file
obj.bigSNP <- snp_attach(paste0(sub_name, ".rds"))
G <- obj.bigSNP$genotypes
SNPs <- obj.bigSNP$map$marker.ID
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
# check if there is any missing values as NA
#big_counts(G, ind.col = 1:dim(G)[1]) # normally the data include missing values
# genotype imputation
G <- snp_fastImputeSimple(G, method = c("mean0"), ncores = 8) # mean0 is based on rounded mean
#big_counts(G, ind.col = 1:dim(G)[1]) # check if NAs are 0
# LD clumping using r2 = 0.2
newpc <- snp_autoSVD(G, infos.chr = CHR, infos.pos = POS, thr.r2 = 0.2, size = 10) # size is the window size of 10K
# extract SNPs after clumpping
which_pruned = attr(newpc, 'subset')
keep_snp_ids = SNPs[which_pruned]
write.table(keep_snp_ids, file = paste0(sub_name, "_pruned_SNP.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
print(paste0("SNPs after clumpping is: ", length(keep_snp_ids), " out of ", dim(obj.bigSNP$map)[1]))

# generate thinned vcf file
system(paste(vcftools," --vcf ",sub_name,".recode.vcf", " --snps ", sub_name, "_pruned_SNP.txt", " --recode --recode-INFO-all --out ", sub_name, "_pruned", sep=""))

# VCFtools - v0.1.13
# (C) Adam Auton and Anthony Marcketta 2009
# 
# Parameters as interpreted:
#   --vcf genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_neutral.recode.vcf
# --recode-INFO-all
# --out genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_neutral_pruned
# --recode
# --snps genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_neutral_pruned_SNP.txt
# 
# After filtering, kept 539 out of 539 Individuals
# Outputting VCF file...
# After filtering, kept 105672 out of a possible 140502 Sites
# Run Time = 17.00 seconds

system(paste(plink, " --vcf ",sub_name,"_pruned.recode.vcf", " --allow-extra-chr --make-bed --out ", sub_name,"_pruned", sep=""))

### shared outliers among the candidatess

x <- list(
  A = outflank_BP$V4, 
  B = PACadapt_BP$V4
)


library("ggvenn")
names(x) <- c("Outflank","PCAdapt")
fill_color = c("#0073C2FF", "#EFC000FF")
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

# # extract the outliers and format for deep-learning
# system(paste(vcftools," --vcf genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe.recode.vcf --snps shared_snp_pcadapt_outflank.txt --recode --recode-INFO-all --out genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_outlier_n_233", sep=""))
# system(paste(plink, " --vcf genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_outlier_n_233.recode.vcf --allow-extra-chr --make-bed --out genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_outlier_n_233", sep=""))
# 
# system(paste(vcftools," --vcf genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe.recode.vcf --snps union_outliers_1575.list --recode --recode-INFO-all --out genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_outlier_n_1575", sep=""))
# system(paste(plink, " --vcf genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_outlier_n_1575.recode.vcf --allow-extra-chr --make-bed --out genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_outlier_n_1575", sep=""))
# 
# f_bk = paste0(sub_name, ".bk")
# if (file.exists(f_bk)) {
#   #Delete file if it exists
#   file.remove(f_bk)
# }
# 
# sub_name = "genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_outlier_n_1575"
# snp_readBed(paste0(sub_name, ".bed"))
# # this will create a .rds file
# 
# toMatrix <- function(G){
#   Gm = matrix(nrow=length(G[,1]), ncol=length(G[1,]))
#   for(i in seq(length(G[,1]))){
#     Gm[i,]=G[i,]
#   }
#   return(Gm)
# }
# 
# obj.bigSNP <- snp_attach(paste0(sub_name, ".rds"))
# G <- obj.bigSNP$genotypes
# SNPs <- obj.bigSNP$map$marker.ID
# CHR <- obj.bigSNP$map$chromosome
# POS <- obj.bigSNP$map$physical.pos
# # check if there is any missing values as NA
# #big_counts(G, ind.col = 1:dim(G)[1]) # normally the data include missing values
# # genotype imputation
# G <- snp_fastImputeSimple(G, method = c("mean0"), ncores = 8) # mean0 is based on rounded mean
# Gm = toMatrix(G)
# G_coded = add_code256(big_copy(Gm, type="raw"), code=bigsnpr:::CODE_012)
# write.table(G_coded[1:539,], file = paste0(sub_name, ".csv"), sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)


####################
# genotype heatmap #
####################
setwd("~/Dropbox/Mac/Documents/HG/Domestication/18_outlier_shared_genotype_heatmap")
#remotes::install_github("JimWhiting91/genotype_plot")
library(devtools)
#install('/Users/HG/Dropbox/Mac/Documents/HG/Github/BioinfoTools/genotype_plot-master')
#install.packages("wesanderson")
library("dplyr")
library("ggplot2")
library(GenotypePlot)
library(R.utils)
library(gghighlight)
library(ggtext)
library(patchwork)
library(plotrix)
library(plyr)
library(qqman)
library(qvalue)
library(GenotypePlot)
library(reshape2)
library(tidyr)
library(zoo)
library(infer)
options(dplyr.summarise.inform = FALSE)
library(bigsnpr)
library("wesanderson")
library("directlabels")
library(OutFLANK)
library(pcadapt)
library(adegenet)
library(poppr)
library(vcfR)
library(maps)
library(mapdata)
library("plyr")
library(RColorBrewer)
library(matrixStats)
library(data.table)
library(ggdendro)
library(ggridges)

setwd("~/Dropbox/Mac/Documents/HG/Domestication/18_outlier_shared_genotype_heatmap")
##################################################################################################
vcftools  = "/Users/HG/Dropbox/Mac/Documents/HG/Github/BioinfoTools/vcftools_0.1.13/bin/vcftools";
plink  = "/Users/HG/Dropbox/Mac/Documents/HG/Domestication/14_ROH/plink";

# for population targeted of contrasts
system(paste(vcftools," --vcf genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe.recode.vcf --bed pop_n_477_pcadapt_outflank.shared.outlier.igv_Dom_Wild.sliding.zfst.outlier.merged.igv.outlier.merged.igv --recode --recode-INFO-all --out pop_n_477_outliers", sep=""))
system(paste(plink, " --vcf pop_n_477_outliers.recode.vcf --allow-extra-chr --make-bed --out pop_n_477_outliers", sep=""))

################### for across chromosomes ###############
# population map
popmap <- read.table("./pop_n_477_outliers.pop.list", header = TRUE)
popmap$Type <- factor(popmap$Type, levels=c("Selected","Wild"))
popmap$POP <- factor(popmap$POP, levels=c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "DBX1", "DBX2", "DBX3",  "UNC1", "UNC2", "UMFS", "NEH1", "NEH2", "MEH2"))
popmap.ni <- popmap

gp.popmap <-setNames(data.frame(matrix(ncol = 2, nrow = length(popmap.ni$POP))), c("ind", "pop"))
gp.popmap$ind <- popmap.ni$Sample
gp.popmap$pop <- popmap.ni$POP

target <- c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "DBX1", "DBX2", "DBX3",  "UNC1", "UNC2", "UMFS", "NEH1", "NEH2", "MEH2")
#target <- rev(target)
gp.popmap <-gp.popmap %>% arrange(factor(pop, levels = target))

gp.label <- read.table("pop_n_477_outliers.popmap.txt", header = TRUE)
gp.label$Type <- factor(gp.label$Type, levels=c("Wild", "Selected"))
gp.label$POP <- factor(gp.label$POP, levels=c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "DBX1", "DBX2", "DBX3",  "UNC1", "UNC2", "UMFS", "NEH1", "NEH2", "MEH2"))

bgzip = "/Users/HG/Dropbox/Mac/Documents/Backup/Ryan_workplace/packages/ROHan/rohan/lib/htslib/bgzip"
tabix = "/Users/HG/Dropbox/Mac/Documents/Backup/Ryan_workplace/packages/ROHan/rohan/lib/htslib/tabix"
bcftools = "/Users/HG/Dropbox/Mac/Documents/HG/Domestication/04_pcadapt/bcftools"
system(paste(paste(bgzip, " -c pop_n_477_outliers.recode.vcf > pop_n_477_outliers.recode.vcf.gz"), sep=""))
system(paste(paste(tabix, " -f -p vcf pop_n_477_outliers.recode.vcf.gz"), sep=""))

# color plate
pal1 <- wes_palette("Zissou1")
pal2 <- wes_palette("GrandBudapest1")
pal3 <- wes_palette("GrandBudapest2")
pal4 <- wes_palette("Royal2")
col_pal <- c(pal1[1], pal1[2], pal2[3], pal3[3], pal3[4], pal3[2],pal1[2],pal4[2],pal4[1],pal4[3])
#   MEW1       MEW2       LIW1        LIW2       DBW1      DBW2        NCW1      NCW2         
cbPalette <- c("#0A2C86", "#325A98",  "#1D92BD", "#3DB9C1", "#7BD2BF", "#C4E9B3", "#ECF6B9", "#EEE8AA", 
               #  DBX1       DBX2      DBX3       UNC1        UNC2       UMFS      NEH1       NEH2       MEH2
               "#F9476B", "#FC709F","#E376B7", "#CF7FBC",  "#A36DC1", "#FEB22B", "#F36616", "#D83B1C", "#FF9117")

chrom_outlier_GP_plot <- function(chrom, vcf, n_ind){
  plot_pal <- c(col_pal[as.integer(chrom)],wes_palette("GrandBudapest2"))
  #plot_pal <- c("#3B9AB2", "#E6A0C4" ,"#C6CDF7" ,"#D8A499" ,"#7294D4")
  outlier.gp <- genotype_plot(vcf=vcf,chr=chrom, popmap = gp.popmap, start=1, end=1000000000,cluster=F,colour_scheme=plot_pal) 
  
  dfx <- min(ggplot_build(outlier.gp$genotypes)$data[[1]]$x) 
  dfmax <- max(ggplot_build(outlier.gp$genotypes)$data[[1]]$x)
  dflen <- length(ggplot_build(outlier.gp$genotypes)$data[[1]]$x)/n_ind # this is represent the number of SNPs
  dloc <- (dfmax-dfx)/dflen
  print(dloc)
  print(dflen)
  
  if(dloc > 50000 & dflen > 50){
    dloc <- dloc/3
  }
  if(dloc > 50000 & dflen < 50){
    dloc <- dloc/2
  }
  #ddelta <- dloc*dflen/20
  if(dflen > 500){
    #ddelta <- dloc*dflen/20
    ddelta <- dloc*50
    dddelta <- ddelta/10000
  }else if(dflen < 500 & dflen > 200){
    ddelta <- dloc*15
    dddelta <- ddelta/1000
  }else if(dflen < 200 & dflen > 100){
    ddelta <- dloc*10
    dddelta <- ddelta/100
  }else if(dflen < 100 & dflen > 20){
    ddelta <- dloc*5
    dddelta <- ddelta/10
  } else if(dflen == 1){
    ddelta <- 10000
    dddelta <- ddelta/2
  } else{
    ddelta <- dloc*5
    dddelta <- ddelta/1
  }
  dfx <- dfx -ddelta/2
  
  #gp.label$X <- min(ggplot_build(outlier.gp$genotypes)$data[[1]]$x-3000000)
  
  gp.label$X <- c(dfx-0,dfx-ddelta/3,dfx-ddelta/2,dfx-ddelta*3/4,dfx-ddelta,dfx-ddelta*3/4,dfx-ddelta/2,dfx-ddelta/3,dfx-0,
                  dfx-0,dfx-ddelta/2,dfx-ddelta*3/4,dfx-ddelta,dfx-ddelta,dfx-ddelta*3/4,dfx-ddelta/2,dfx-0
                  )
  
  y_axis <- ggplot_build(outlier.gp$genotypes)$data[[2]]$yintercept
  gp.label$Y <- rollmean(y_axis, 2)
  
  fig <-  outlier.gp$genotypes +
    geom_point(data=gp.label,aes(x=X,y=Y,color=POP,shape=Type),size=6, alpha=0.95)  +
    scale_x_continuous(expand = expand_scale(mult = c(0.1, 0)))+
    #scale_x_continuous(limits = c(min(gp.label$X) - ddelta/2 , max(ggplot_build(outlier.gp$genotypes)$data[[1]]$x))) + 
    scale_shape_manual(values=c(18, 20), name="Origin") + 
    scale_color_manual(values=cbPalette , name="Population/Line") +
    theme(legend.position="right") +
    theme(legend.title = element_text(size = 12),
          legend.text=element_text(size=12)) + 
    # guides(fill = "none")+ 
    # guides(color="none")+
    # guides(shape = "none") +
    theme(text=element_text(family="Times New Roman", face="bold", size=12, colour="black"))
  
  assign(paste("gp.plot",chrom,"out", sep="."),fig, envir = globalenv())
}

for(i in c(1)){
  chrom_outlier_GP_plot(i,"pop_n_477_outliers.recode.vcf.gz", 539)
}

gp.plot.1.out<- gp.plot.1.out+
  guides(fill = guide_legend(override.aes = list(size=6)),shape=guide_legend(override.aes = list(size=6, shape=c(23,21)) ))+ 
  guides(color=guide_legend(override.aes = list(size=6)))

chrom_outlier_GP_plot_no_label <- function(chrom, vcf, n_ind){
  plot_pal <- c(col_pal[as.integer(chrom)],wes_palette("GrandBudapest2"))
  #plot_pal <- c("#3B9AB2", "#E6A0C4" ,"#C6CDF7" ,"#D8A499" ,"#7294D4")
  outlier.gp <- genotype_plot(vcf=vcf,chr=chrom, popmap = gp.popmap, start=1, end=1000000000,cluster=F,colour_scheme=plot_pal) 
  
  fig <-  outlier.gp$genotypes +
    scale_shape_manual(values=c(18, 20), name="Origin") + 
    scale_color_manual(values=cbPalette , name="Population/Line") +
    theme(legend.position="right") +
    theme(legend.title = element_text(size = 12),
          legend.text=element_text(size=12)) + 
    guides(fill = "none")+
    guides(color="none")+
    guides(shape = "none") +
    theme(text=element_text(family="Times New Roman", face="bold", size=12, colour="black"))
  
  
  assign(paste("gp.plot",chrom,"out", sep="."),fig, envir = globalenv())
}

#chrom_outlier_GP_plot_no_label(2,"pop_n_477_outliers.recode.vcf.gz", 539)

for(i in c(5)){
  chrom_outlier_GP_plot_no_label(i,"pop_n_477_outliers.recode.vcf.gz", 539)
}

layout <-" 
AAAAAAAAAAAAAA
AAAAAAAAAAAAAA
AAAAAAAAAAAAAA
AAAAAAAAAAAAAA
AAAAAAAAAAAAAA
BBBBCCDDDDDDDD
BBBBCCDDDDDDDD
BBBBCCDDDDDDDD
BBBBCCDDDDDDDD
BBBBCCDDDDDDDD
BBBBCCDDDDDDDD
"

tiff("figure2.tiff", units="in", width=12, height=8, res=300)
p/gp.plot.1.out  + gp.plot.5.out + plot_layout(design=layout, guides="collect")
dev.off()

layout <-" 
BBBDDDDDDDDDDD
BBBDDDDDDDDDDD
BBBDDDDDDDDDDD
BBBDDDDDDDDDDD
BBBDDDDDDDDDDD
BBBDDDDDDDDDDD
"

tiff("figure2_part3.tiff", units="in", width=12, height=10, res=300)
gp.plot.1.out + gp.plot.5.out + plot_layout(design=layout, guides="collect")
dev.off()
