library(devtools)
install_github("jgx65/hierfstat")
library("hierfstat")
library("adegenet")
library(vcfR)
library(MASS)

# see https://popgen.nescent.org/DifferentiationSNP.html for detailed example
# read from a vcf file
setwd("~/Dropbox/Mac/Documents/HG/Domestication/13_diversity_Fis_hierfstat")

# bash scrip to subset the 1K SNPs from neutral pruned SNPs
# 15_relatedness.sh
# bcftools query -f '%CHROM\t%POS\t%ID\n' genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_neutral_pruned.recode.vcf  > neutral_SNPs.list
# shuf -n 1000 neutral_SNPs.list | cut -f3 -d$'\t' > neutral_SNPs_n_1K_list.txt
# vcftools --vcf genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_neutral_pruned.recode.vcf --snps neutral_SNPs_n_1K_list.txt --recode --recode-INFO-all --out genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_neutral_pruned_1K

# load the population information
pop_info <- read.table("pop_539_sample_list.txt", header=TRUE, sep="\t", stringsAsFactors = TRUE)
pop_info$Pop_correct = factor(pop_info$Pop_correct, levels=c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "DBX1", "DBX2", "DBX3", "UNC1", "UNC2", "UMFS", "NYH1", "NEH1", "NEH2", "MEH2"))
pop_info$Pop_correct = factor(pop_info$Pop_correct.1, levels=c("Native", "Selected"))
# load vcf file
vcftools  = "/Users/HG/Dropbox/Mac/Documents/HG/Github/BioinfoTools/vcftools_0.1.13/bin/vcftools";
#system(paste(vcftools," --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_pruned.recode.vcf --remove UMFS_2_outlier.txt --recode --recode-INFO-all --out genetyped_data_n_507_maf05_maxmiss095_popmiss095_hwe_pruned", sep=""))

#vcf_file = "genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_pruned_no_outlier.recode.vcf"
vcf_file = "genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_neutral_pruned_10K.recode.vcf"
vcf <- read.vcfR(vcf_file, verbose = FALSE)
Mydata1 <- vcfR2genind(vcf)
Mydata1@pop <- pop_info$Pop_correct
Mydata1

# vcf_file = "genetyped_data_n_507_maf05_maxmiss095_popmiss095_hwe_pruned.recode.vcf"
# #vcf_file = "genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_pruned_neutral_1K.recode.vcf"
# vcf <- read.vcfR(vcf_file, verbose = FALSE)
# Mydata2 <- vcfR2genind(vcf)
# Mydata2@pop <- pop_info$Pop_correct
# Mydata2

# Estimates allelic richness, the rarefied allelic counts, per locus and population
#Arich <- allelic.richness(Mydata1,min.n=NULL,diploid=TRUE)
#colMeans(x=Arich$Ar, na.rm = TRUE)

Arich <- allelic.richness(Mydata1,min.n=NULL,diploid=TRUE)
ind_mean <- colMeans(x=Arich$Ar, na.rm = TRUE)
ind_mean
# MEW1     MEW2     LIW1     LIW2     DBW1     DBW2     NCW1     NCW2     DBX1     
# 1.835990 1.845247 1.864615 1.866921 1.870363 1.873105 1.869544 1.841741 1.803380 
# DBX2     DBX3     UNC1     UNC2     UMFS     NYH1     NEH1     NEH2     MEH2 
# 1.766539 1.842716 1.722316 1.701973 1.711456 1.453512 1.809908 1.781295 1.664748 
wilcox.test(ind_mean[1:9],ind_mean[10:18], alternative = "greater")
mean(ind_mean[1:9])
mean(ind_mean[10:18])
# 
# Wilcoxon rank sum exact test
# 
# data:  ind_mean[1:9] and ind_mean[10:18]
# W = 77, p-value = 0.0002468
# alternative hypothesis: true location shift is greater than 0
write.table(ind_mean, file = "539.individual.allelic.richness.txt", sep = "\t", quote = FALSE,
            row.names = T, col.names = F)

# These statistics come from the package hierfstat Fst following Nei (1987) on genind object
basicstat <- basic.stats(Mydata1, diploid = TRUE, digits = 2) 
names(basicstat)
# $overall
# Ho   Hs   Ht  Dst  Htp Dstp  Fst Fstp  Fis Dest 
# 0.24 0.26 0.28 0.02 0.28 0.02 0.07 0.07 0.08 0.03 
# mean Ho per population
Ho <- colMeans(x=basicstat$Ho, na.rm = TRUE)
# MEW1     MEW2     LIW1     LIW2     DBW1     DBW2     NCW1     NCW2     DBX1     
# 0.241115 0.240736 0.245775 0.248049 0.242163 0.241222 0.239758 0.231845 0.242839 
# DBX2     DBX3     UNC1     UNC2     UMFS     NYH1     NEH1     NEH2     MEH2 
# 0.238546 0.251046 0.247326 0.226499 0.231323 0.185497 0.246155 0.247423 0.250486 
# mean He per population
write.table(Ho, file = "pop.Ho.txt", sep = "\t", quote = FALSE,
            row.names = T, col.names = F)

wilcox.test(Ho[1:9],Ho[10:18], alternative = "greater")

He <- colMeans(x=basicstat$Hs, na.rm = TRUE)
write.table(He, file = "pop.He.txt", sep = "\t", quote = FALSE,
            row.names = T, col.names = F)

# MEW1     MEW2     LIW1     LIW2     DBW1     DBW2     NCW1     NCW2     DBX1     
# 0.270641 0.271958 0.271036 0.272392 0.275185 0.276896 0.276797 0.267170 0.262844 
# DBX2     DBX3     UNC1     UNC2     UMFS     NYH1     NEH1     NEH2     MEH2 
# 0.256842 0.275769 0.244617 0.238100 0.238581 0.170907 0.267275 0.260375 0.233630 

wilcox.test(He[1:9],He[10:18], alternative = "greater")

# Wilcoxon rank sum exact test
# 
# data:  He[1:9] and He[10:18]
# W = 72, p-value = 0.001995
# alternative hypothesis: true location shift is greater than 0


# div <- summary(Mydata1)
# names(div)
# 
wc(Mydata1) # Weir and Cockerham's estimate
fst <- genet.dist(Mydata1, method = "WC84") # Pairwise Fst
write.matrix(fst, file = "pairwise.fst.txt", ,sep = "\t")

# MEW1          MEW2          LIW1          LIW2          DBW1          DBW2          NCW1          NCW2          DBX1          DBX2          DBX3          UNC1          UNC2          UMFS          NYH1          NEH1          NEH2
# MEW2 -0.0009593305                                                                                                                                                                                                                                
# LIW1  0.0346639148  0.0305511685                                                                                                                                                                                                                  
# LIW2  0.0339778291  0.0302968939  0.0010724035                                                                                                                                                                                                    
# DBW1  0.0418246994  0.0375860041  0.0098443592  0.0093000441                                                                                                                                                                                      
# DBW2  0.0423146753  0.0379062865  0.0103591098  0.0101741504  0.0005883800                                                                                                                                                                        
# NCW1  0.0527838439  0.0482495241  0.0222212044  0.0213112054  0.0115309962  0.0118031009                                                                                                                                                          
# NCW2  0.0653728290  0.0606550983  0.0347848789  0.0346985940  0.0249523770  0.0256331618  0.0301914010                                                                                                                                            
# DBX1  0.0683318649  0.0638312812  0.0381871443  0.0370123790  0.0288022386  0.0281400654  0.0396447751  0.0538945015                                                                                                                              
# DBX2  0.0696225364  0.0686706711  0.0656557000  0.0655516763  0.0652489138  0.0648144498  0.0735280461  0.0872928873  0.0912087668                                                                                                                
# DBX3  0.0438672386  0.0405141282  0.0424130348  0.0417346448  0.0422651668  0.0421497505  0.0485189181  0.0661932846  0.0677659146  0.0665671191                                                                                                  
# UNC1  0.1113811671  0.1064435456  0.0805505569  0.0790577116  0.0700659624  0.0702496086  0.0623160974  0.0886295359  0.0993809539  0.1314432782  0.1060913346                                                                                    
# UNC2  0.1164045047  0.1130352259  0.0863312778  0.0873234019  0.0773625597  0.0779444912  0.0816738565  0.0570957521  0.1077941300  0.1395347981  0.1164003183  0.1434559961                                                                      
# UMFS  0.0676687152  0.0675870283  0.0800422934  0.0801382864  0.0866077730  0.0874753143  0.0964398639  0.1099992358  0.1131196753  0.1335211869  0.1101340550  0.1591749079  0.1616785919                                                        
# NYH1  0.2141561642  0.2118194076  0.1888917338  0.1830966801  0.1985912309  0.1988881983  0.2060322351  0.2194704327  0.2236693378  0.2433496270  0.2146087653  0.2743396169  0.2787651641  0.2615076279                                          
# NEH1  0.0487987771  0.0443718318  0.0650506527  0.0641664653  0.0702813248  0.0691875629  0.0768031934  0.0918576494  0.0966679938  0.0800571437  0.0354324508  0.1330171790  0.1419255920  0.1350096083  0.2345979982                            
# NEH2  0.0519283042  0.0468007475  0.0712864735  0.0706000373  0.0763904929  0.0756587170  0.0845780269  0.0984912019  0.1043557389  0.0735429630  0.0463601714  0.1415966116  0.1493538720  0.1379843670  0.2424746086  0.0138342741              
# MEH2  0.0858078262  0.0838242073  0.1014988189  0.1014738736  0.1081389012  0.1079199851  0.1197938346  0.1324937576  0.1362969421  0.1328129197  0.1052473456  0.1793881003  0.1863947225  0.1601039611  0.2797840519  0.1073001323  0.1063059950

# using Kmeans and DAPC in adegenet 
# set.seed(5); dapc_a_score <- dapc(Mydata1,Mydata1$pop, n.pca = 20,n.da=10)
# temp_score <- optim.a.score(dapc_a_score)
# 
# dapc1 <- dapc(Mydata1, Mydata1$pop, n.pca = 9, n.da = 2) 
# scatter(dapc1, legend=TRUE, solid=.5) # plot of the group
# graph2ppt(file="DAPC",width=8,height=5)
# 
# percent= dapc1$eig/sum(dapc1$eig)*100
# barplot(percent, ylab="Percent of genetic variance explained by eigenvectors", names.arg=round(percent,2))


################
#### ggplot ####
################
library(cowplot)
df <- read.delim("He_Ar_summary.csv", header = TRUE, sep=',')
df$Source <- factor(df$Source, levels=c("Wild populations", "Selected line"))
wilcox.test(df$Ne[1:8],df$Ne[9:17])

Ar <- df %>%
  ggplot(aes(Source, Ar, fill=Source)) +
  geom_violin(trim = FALSE, alpha=0.4, show.legend = FALSE) +
  geom_boxplot(width = 0.1, aes(fill=Source))  +
  scale_fill_manual(values = c("#049DD9", "#F25C05"))+
  xlab(NULL) + ylab("Allelic richness (Ar)")+
  theme_classic()+
  theme(text = element_text(size=20),
        legend.position = "none")
  #geom_point(position = "jitter", alpha = 0.7, size = 3)
He <- df %>%
  ggplot(aes(Source, He, fill=Source)) +
  geom_violin(trim = FALSE, alpha=0.4, show.legend = FALSE) +
  geom_boxplot(width = 0.1, aes(fill=Source))  +
  scale_fill_manual(values = c("#049DD9", "#F25C05"))+
  xlab(NULL) + ylab("Expected heterozygosity (He)")+
  theme_classic()+
  theme(text = element_text(size=20),
        legend.position = "none")

#   MEW1        MEW2      LIW1        LIW2       DBW1      DBW2        NCW1      NCW2         
col <- c( "#0A2C86", "#325A98",  "#1D92BD", "#3DB9C1", "#C4E9B3", "#7BD2BF", "#ECF6B9", "#EEE8AA", 
          #  DBX1       DBX2      DBX3       UNC1        UNC2       UMFS      NEH1       NEH2       MEH2
          "#F9476B", "#FC709F","#E376B7", "#CF7FBC",  "#A36DC1", "#FEB22B", "#F36616", "#D83B1C", "#FF9117")
df$Sites = factor(df$Sites, levels=c("MEW1", "MEW2","LIW1","LIW2","DBW1","DBW2","NCW1","NCW2","DBX1","DBX2","DBX3","UNC1","UNC2","UMFS","NEH1","NEH2", "MEH2"))

Ne <- ggplot(df, aes(x=Sites, y=Ne, fill=Sites)) +
  geom_bar(stat="identity", fill=col) +
  geom_errorbar(aes(ymin=Down,ymax=Up), width=0.5)+
  theme_classic() + 
  theme(axis.text.x = element_text(angle=45, hjust=1),
        text = element_text(size=20),) +
  labs(x=NULL, y = "Effective population size (Ne)")  +  
  scale_fill_manual(values=col) 

N1 <- Ne + coord_cartesian(ylim = c(0, 700))
N2 <- Ne + coord_cartesian(ylim = c(10000, 50000))
up_row <- plot_grid(He, Ar)
plot_grid(
  up_row,
  N2,
  ncol = 1
)
graph2ppt(file="Diversity1",width=10,height=6)

