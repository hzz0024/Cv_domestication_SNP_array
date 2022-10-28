library(randomForest)
library(RcppCNPy)
library(export)

# Format for random forest ###################################################
# set up directory
setwd("~/Dropbox/Mac/Documents/HG/Domestication/19_random_forest/format/")
vcftools  = "/Users/HG/Dropbox/Mac/Documents/HG/Github/BioinfoTools/vcftools_0.1.13/bin/vcftools";
plink  = "/Users/HG/Dropbox/Mac/Documents/HG/Domestication/14_ROH/plink";

# for the union set of outliers from PCAdapt and Outflank
system(paste(vcftools," --vcf genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe.recode.vcf --bed union_outliers_1174.bed --recode --recode-INFO-all --out pop_n_539_outlier_n_1174", sep=""))
system(paste(plink, " --vcf pop_n_539_outlier_n_1174.recode.vcf --allow-extra-chr --make-bed --out pop_n_539_outlier_n_1174", sep=""))

sub_name = "pop_n_539_outlier_n_1174"
f_bk = paste0(sub_name, ".bk")
if (file.exists(f_bk)) {
  #Delete file if it exists
  file.remove(f_bk)
}

snp_readBed(paste0(sub_name, ".bed"))
# this will create a .rds file

toMatrix <- function(G){
  Gm = matrix(nrow=length(G[,1]), ncol=length(G[1,]))
  for(i in seq(length(G[,1]))){
    Gm[i,]=G[i,]
  }
  return(Gm)
}

obj.bigSNP <- snp_attach(paste0(sub_name, ".rds"))
G <- obj.bigSNP$genotypes
SNPs <- obj.bigSNP$map$marker.ID
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
# check if there is any missing values as NA
#big_counts(G, ind.col = 1:dim(G)[1]) # normally the data include missing values
# genotype imputation
G <- snp_fastImputeSimple(G, method = c("mean0"), ncores = 8) # mean0 is based on rounded mean
Gm = toMatrix(G)
G_coded = add_code256(big_copy(Gm, type="raw"), code=bigsnpr:::CODE_012) # the genotype will be covert to 0,0.5,1 using excel

snp_id = read.delim("pop_n_539_outlier_n_1174.bim", header = FALSE, sep = "\t")
popmap_fam = read.delim("pop_n_539_outlier_n_1174.fam", header = FALSE, sep = " ")
df<-data.frame(id=popmap_fam$V1,pop=substr(popmap_fam$V1, 1,4), geno=G_coded[1:539,])
colnames(df) = c("Sample", "Pop", snp_id$V2)
write.table(df, file = paste0(sub_name, ".csv"), sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE) 

# Start to run random forest ###################################################
setwd("~/Dropbox/Mac/Documents/HG/Domestication/19_random_forest")
# Import the data set, which comprises 509 individuals genotyped at 606 biallelic loci, where 0=homozygote 1, 1=heterozygote, 2=homozygote 2
# Each individual also has a binary phenotype - selected or wild - where 1=wild and 2=selected lines
# The objective is to identify loci associated with convergent response to captive conditions
class_data <- read.csv("pop_n_539_outlier_n_1174.csv",row.names=1) #This is the simulated data from Table S8

# First, explore the overall distribution of the phenotype
hist(class_data$Class) 

# linear transfer the input #####################
# # wild/selected individuals appears to differ between the populations, and we should correct
# # for this stratification before conducting RF to minimize the risk of false positive associations.
# # Specifically, we will correct the genotypes using the approach of Zhao et al.(2012). 
# # We will not correct the phenotype, as it is binary and we want to maintain its categorical distribution.
# class_data_corrected <- class_data # create another data frame for the corrected genotypes
# class_data_corrected[,3:1176] <- NA  # Keep columns 1-2 with population ID and phenotype, but then replace with NA's over which you can write the residuals.
# 
# # Now correct the genotypes using the regression/residual method. We're using a standard linear regression because Zhao et al. 2012 found that the correction procedure is robust to selection of the link function
# for (i in 3:ncol(class_data)){
#   print(i)
#   LM_SNP_i <- lm(class_data[,i] ~ factor(class_data$Pop)) # apply linear model to all loci and the response
#   class_data_corrected[,i] <- LM_SNP_i$residuals
#   colnames(class_data_corrected)[i]<-colnames(class_data)[i] 
#   if(i%%50==0) print(i)
# }
# 
# # Verify that the residuals have been written to the data frame properly, using the last column as an example
# class_data_corrected[,1176]-LM_SNP_i$residuals  #Should all be zero if correct
# 
# # Export a copy of the corrected data for future reference (This corresponds to the data in Table S9)
# write.csv(class_data_corrected,file="pop_n_539_outlier_n_1174_corrected.csv",row.names=FALSE)

# count the sample size for random forest running ###########
# Before running Random Forest, let's also check for an imbalance in the response variable 
# because - as discussed in the manuscript - any imbalances can bias the results.
length(which(class_data$Class==1)) # 279 susceptible individuals total
length(which(class_data$Class==2)) # 260 resistant individuals total

sample_cnt <- as.integer(min(length(which(class_data$Class==1)), length(which(class_data$Class==2)) )*(2/3)) # 153 susceptible individuals total
sample_size <- c(sample_cnt,sample_cnt)
###########################################################################################################################################
###########################################################################################################################################

# Now run Random Forest analysis. Since this is a binary trait, we need to conduct a classification RF

# First, we need to optimize mtry by running different values of mtry at different values of ntree. 

# We will run mtry values of sqrt(p), 2*sqrt(p), 0.1(p), 0.2(p), p/3, and p, where p is the number of loci
# We will initially run each of these mtry values at ntree=100 to 1000 (by increments of 100). 
# We are looking for a plateau where the out-of-bag error rate (OOB-ER) stops decreasing with larger values of ntree
# Once we reach the plateau, we will choose the mtry value that minimizes the OOB-ER.
p = 1174


results_optimization <- matrix(data=NA , nrow = 0, ncol = 3)
for (i in seq(from = 1, to = 1000 , by = 1)){  # values of ntree
  print(i)
  for (j in c(as.integer(sqrt(p)), as.integer(2*sqrt(p)), as.integer(0.1*p), as.integer(0.2*p), as.integer(p/3), p)){    #values of mtry based on 1000 total loci
    rf_ij <- randomForest(x = class_data[,3:1176], y = as.factor(class_data$Class), importance=TRUE ,proximity=TRUE, ntree=i, mtry=j, strata=as.factor(class_data$Class), sampsize=sample_size)
    results_optimization <- rbind(results_optimization, c(i,j,tail(rf_ij$err.rate,1)[1]))
  }
}

# Clean up the file format
results_optimization<-as.data.frame(results_optimization)
colnames(results_optimization)<-c("ntree", "mtry","OOB_ER")

# Now plot results to see if there's a plateau
plot(results_optimization$ntree[results_optimization$mtry == 34],results_optimization$OOB_ER[results_optimization$mtry == 34], type="l", col="black", xlab="ntree",ylab="OOB-ER",ylim=c(0,0.3))
lines(results_optimization$ntree[results_optimization$mtry == 68],results_optimization$OOB_ER[results_optimization$mtry == 68], col="blue")
lines(results_optimization$ntree[results_optimization$mtry == 117],results_optimization$OOB_ER[results_optimization$mtry == 117], col="green")
lines(results_optimization$ntree[results_optimization$mtry == 234],results_optimization$OOB_ER[results_optimization$mtry == 234], col="purple")
lines(results_optimization$ntree[results_optimization$mtry == 391],results_optimization$OOB_ER[results_optimization$mtry == 391], col="orange")
lines(results_optimization$ntree[results_optimization$mtry == 1174],results_optimization$OOB_ER[results_optimization$mtry == 1174], col="red")

mytheme = list(
  theme_classic()+
    theme(panel.background = element_blank(),strip.background = element_rect(colour="black", fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.position="right", strip.text = element_text(face="plain", size=12),
          axis.text=element_text(face="plain",size=12),axis.title = element_text(face="plain",size=14),plot.title = element_text(face = "bold", hjust = 0.5,size=12))
)

p1 <- ggplot(results_optimization, aes(x=ntree, y=OOB_ER, color=as.factor(mtry)))+
  #geom_point(shape=1, alpha=0.05) +
  #facet_grid(mtry~.)+
  geom_smooth(linetype="dashed")+
  theme_bw()+
  mytheme+
  labs(color='Number of predictors')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Number of trees")+ylab("Out-of-bag error rate (OOB-ER)")
  #labs(caption = "34, 68, 117, 234, 391, 1174 represent sqrt(p), , 0.1(p), 2*sqrt(p), 0.2(p), p/3, and p, where p is the number of loci")
p1
graph2ppt(file="1_OOB_ER_vs_ntree", width=8, height=6) 

# This plot shows that mtry=p is the best in terms of OOB-ER (although mtry=0.2p and p/3 are very similar), and that the OOB-ER has reached a plateau. 
# Therefore, we will use mtry=p for our Random Forest analyses. 
# Note that this plot differs from Figure 3 in the manuscript and thus demonstrates that optimal parameter values will vary based on each data set.

###########################################################################################################################################
###########################################################################################################################################
# Now begin the full Random Forest analyses

ntree = 1000

# Recall that even though we optimized mtry, we must now run a larger number of trees in order to achieve convergence of importance values between forests.
# As a starting point, we will grow 25,000 trees and increase if necessary. We do not need to worry about this increase in ntree affecting our mtry optimization,
# since the OOB-ER reached a plateau for a given mtry value after about 400 trees.

rf_all_1 = randomForest(x = class_data[,3:1176], y = as.factor(class_data$Class), importance=TRUE ,proximity=TRUE, mtry=391, ntree=ntree, strata=as.factor(class_data$Class), sampsize=sample_size)
#save(rf_all_1,file="rf_all_1.Rdata")

rf_all_2 = randomForest(x = class_data[,3:1176], y = as.factor(class_data$Class), importance=TRUE ,proximity=TRUE, mtry=391, ntree=ntree, strata=as.factor(class_data$Class), sampsize=sample_size)
#save(rf_all_2,file="rf_all_2.Rdata")

#Check correlation of locus importance values between forests 
importance_rf_all_1<-data.frame(importance(rf_all_1,type=1)) #type=1 is mean decrease in accuracy for classification, so a large, positive value means that permuting the variable led to a big decrease in prediction accuracy (which is indicative of an important locus)
colnames(importance_rf_all_1)<-c("importance")
importance_rf_all_2<-data.frame(importance(rf_all_2,type=1))
colnames(importance_rf_all_2)<-c("importance")

cor(importance_rf_all_1,importance_rf_all_2) # 0.8846539 # A correlation of 0.98 for locus importance values between forests is extremely good, so we'll use ntree trees for the remaining forests

rf_all_3 = randomForest(x = class_data[,3:1176], y = as.factor(class_data$Class), importance=TRUE ,proximity=TRUE, mtry=391, ntree=ntree, strata=as.factor(class_data$Class), sampsize=sample_size)
#save(rf_all_3,file="rf_all_3.Rdata")
importance_rf_all_3<-data.frame(importance(rf_all_3,type=1))
colnames(importance_rf_all_3)<-c("importance")

############################################################################################################################################
############################################################################################################################################

# The predictive ability of classification trees is measured by the out-of-bag error rate.  An error rate is calculated for each tree within a forest. 
# We will use the error rate from the last tree in the forest, which takes all previous trees into account and thus represents the error rate after the model stabilizes/converges

rf_all_1_err.rate <- rf_all_1$err.rate[ntree]
rf_all_2_err.rate <- rf_all_2$err.rate[ntree]
rf_all_3_err.rate <- rf_all_3$err.rate[ntree]

#Combine importance (mean decrease in accuracy) values of each locus across the three forests
importance_rf_all <-cbind(rownames(importance_rf_all_1),importance_rf_all_1,importance_rf_all_2, importance_rf_all_3)
colnames(importance_rf_all)<-c("Variable","Importance1","Importance2", "Importance3")

# Export importance values for future reference
write.csv(importance_rf_all,file="rf_importance_values_outlier_n_1174.csv",row.names=FALSE)
read.csv("rf_importance_values_outlier_n_1174.csv")

############################################################################################################################################
############################################################################################################################################

# Now conduct RF on subsets of the data to identify a group of loci that may be predictive of disease Class. 
# For each subset, we will use mtry=p since that is the optimal setting that we previously found.

##### Best 2% 

names_best_2perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.98))]
names_best_2perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.98))]
names_best_2perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.98))]
names_best_2perc_unique<-unique(c(names_best_2perc_1,names_best_2perc_2,names_best_2perc_3))

# Extract genotypes 
genotypes_2perc<-class_data[,colnames(class_data) %in% names_best_2perc_unique]

# Now conduct RF on this subset
rf_2perc_1 = randomForest(x = genotypes_2perc, y = as.factor(class_data$Class), importance=TRUE ,proximity=TRUE, mtry=length(names_best_2perc_unique), ntree=ntree, strata=as.factor(class_data$Class), sampsize=sample_size)
#save(rf_2perc_1,file="rf_2perc_1.Rdata")

rf_2perc_2 = randomForest(x = genotypes_2perc, y = as.factor(class_data$Class), importance=TRUE ,proximity=TRUE, mtry=length(names_best_2perc_unique), ntree=ntree, strata=as.factor(class_data$Class), sampsize=sample_size)
#save(rf_2perc_2,file="rf_2perc_2.Rdata")

rf_2perc_3 = randomForest(x = genotypes_2perc, y = as.factor(class_data$Class), importance=TRUE ,proximity=TRUE, mtry=length(names_best_2perc_unique), ntree=ntree, strata=as.factor(class_data$Class), sampsize=sample_size)
#save(rf_2perc_3,file="rf_2perc_3.Rdata")

# Extract and save the out of bag error rate for this subset of loci
rf_2perc_1_err.rate <- rf_2perc_1$err.rate[ntree]
rf_2perc_2_err.rate <- rf_2perc_2$err.rate[ntree]
rf_2perc_3_err.rate <- rf_2perc_3$err.rate[ntree]

rm(rf_2perc_1,rf_2perc_2,rf_2perc_3) # remove the objects to save memory in R 
##### Best 3% 

names_best_3perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.97))]
names_best_3perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.97))]
names_best_3perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.97))]
names_best_3perc_unique<-unique(c(names_best_3perc_1,names_best_3perc_2,names_best_3perc_3))

# Extract genotypes
genotypes_3perc<-class_data[,colnames(class_data) %in% names_best_3perc_unique]

# Now conduct RF on this subset
rf_3perc_1 = randomForest(x = genotypes_3perc, y = as.factor(class_data$Class), importance=TRUE ,proximity=TRUE, mtry=length(names_best_3perc_unique), ntree=ntree, strata=as.factor(class_data$Class), sampsize=sample_size)
#save(rf_3perc_1,file="rf_3perc_1.Rdata")

rf_3perc_2 = randomForest(x = genotypes_3perc, y = as.factor(class_data$Class), importance=TRUE ,proximity=TRUE, mtry=length(names_best_3perc_unique), ntree=ntree, strata=as.factor(class_data$Class), sampsize=sample_size)
#save(rf_3perc_2,file="rf_3perc_2.Rdata")

rf_3perc_3 = randomForest(x = genotypes_3perc, y = as.factor(class_data$Class), importance=TRUE ,proximity=TRUE, mtry=length(names_best_3perc_unique), ntree=ntree, strata=as.factor(class_data$Class), sampsize=sample_size)
#save(rf_3perc_3,file="rf_3perc_3.Rdata")

# Extract and save the out of bag error rate for this subset of loci
rf_3perc_1_err.rate <- rf_3perc_1$err.rate[ntree]
rf_3perc_2_err.rate <- rf_3perc_2$err.rate[ntree]
rf_3perc_3_err.rate <- rf_3perc_3$err.rate[ntree]

rm(rf_3perc_1,rf_3perc_2,rf_3perc_3) # remove the objects to save memory in R 

##### Best 4% 

names_best_4perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.96))]
names_best_4perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.96))]
names_best_4perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.96))]
names_best_4perc_unique<-unique(c(names_best_4perc_1,names_best_4perc_2,names_best_4perc_3))

# Extract genotypes
genotypes_4perc<-class_data[,colnames(class_data) %in% names_best_4perc_unique]

# Now conduct RF on this subset
rf_4perc_1 = randomForest(x = genotypes_4perc, y = as.factor(class_data$Class), importance=TRUE ,proximity=TRUE, mtry=length(names_best_4perc_unique), ntree=ntree, strata=as.factor(class_data$Class), sampsize=sample_size)
#save(rf_4perc_1,file="rf_4perc_1.Rdata")

rf_4perc_2 = randomForest(x = genotypes_4perc, y = as.factor(class_data$Class), importance=TRUE ,proximity=TRUE, mtry=length(names_best_4perc_unique), ntree=ntree, strata=as.factor(class_data$Class), sampsize=sample_size)
#save(rf_4perc_2,file="rf_4perc_2.Rdata")

rf_4perc_3 = randomForest(x = genotypes_4perc, y = as.factor(class_data$Class), importance=TRUE ,proximity=TRUE, mtry=length(names_best_4perc_unique), ntree=ntree, strata=as.factor(class_data$Class), sampsize=sample_size)
#save(rf_4perc_3,file="rf_4perc_3.Rdata")

# Extract and save the out of bag error rate for this subset of loci
rf_4perc_1_err.rate <- rf_4perc_1$err.rate[ntree]
rf_4perc_2_err.rate <- rf_4perc_2$err.rate[ntree]
rf_4perc_3_err.rate <- rf_4perc_3$err.rate[ntree]

rm(rf_4perc_1,rf_4perc_2,rf_4perc_3) # remove the objects to save memory in R 

##### Best 5% 

names_best_5perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.95))]
names_best_5perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.95))]
names_best_5perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.95))]
names_best_5perc_unique<-unique(c(names_best_5perc_1,names_best_5perc_2,names_best_5perc_3))

# Extract genotypes
genotypes_5perc<-class_data[,colnames(class_data) %in% names_best_5perc_unique]

# Now conduct RF on this subset
rf_5perc_1 = randomForest(x = genotypes_5perc, y = as.factor(class_data$Class), importance=TRUE ,proximity=TRUE, mtry=length(names_best_5perc_unique), ntree=ntree, strata=as.factor(class_data$Class), sampsize=sample_size)
#save(rf_5perc_1,file="rf_5perc_1.Rdata")

rf_5perc_2 = randomForest(x = genotypes_5perc, y = as.factor(class_data$Class), importance=TRUE ,proximity=TRUE, mtry=length(names_best_5perc_unique), ntree=ntree, strata=as.factor(class_data$Class), sampsize=sample_size)
#save(rf_5perc_2,file="rf_5perc_2.Rdata")

rf_5perc_3 = randomForest(x = genotypes_5perc, y = as.factor(class_data$Class), importance=TRUE ,proximity=TRUE, mtry=length(names_best_5perc_unique), ntree=ntree, strata=as.factor(class_data$Class), sampsize=sample_size)
#save(rf_5perc_3,file="rf_5perc_3.Rdata")

# Extract and save the out of bag error rate for this subset of loci
rf_5perc_1_err.rate <- rf_5perc_1$err.rate[ntree]
rf_5perc_2_err.rate <- rf_5perc_2$err.rate[ntree]
rf_5perc_3_err.rate <- rf_5perc_3$err.rate[ntree]

rm(rf_5perc_1,rf_5perc_2,rf_5perc_3) # remove the objects to save memory in R 

##### Best 10% 

names_best_10perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.90))]
names_best_10perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.90))]
names_best_10perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.90))]
names_best_10perc_unique<-unique(c(names_best_10perc_1,names_best_10perc_2,names_best_10perc_3))

# Extract genotypes
genotypes_10perc<-class_data[,colnames(class_data) %in% names_best_10perc_unique]

# Now conduct RF on this subset
rf_10perc_1 = randomForest(x = genotypes_10perc, y = as.factor(class_data$Class), importance=TRUE ,proximity=TRUE, mtry=length(names_best_10perc_unique), ntree=ntree, strata=as.factor(class_data$Class), sampsize=sample_size)
#save(rf_10perc_1,file="rf_10perc_1.Rdata")

rf_10perc_2 = randomForest(x = genotypes_10perc, y = as.factor(class_data$Class), importance=TRUE ,proximity=TRUE, mtry=length(names_best_10perc_unique), ntree=ntree, strata=as.factor(class_data$Class), sampsize=sample_size)
#save(rf_10perc_2,file="rf_10perc_2.Rdata")

rf_10perc_3 = randomForest(x = genotypes_10perc, y = as.factor(class_data$Class), importance=TRUE ,proximity=TRUE, mtry=length(names_best_10perc_unique), ntree=ntree, strata=as.factor(class_data$Class), sampsize=sample_size)
#save(rf_10perc_3,file="rf_10perc_3.Rdata")

# Extract and save the out of bag error rate for this subset of loci
rf_10perc_1_err.rate <- rf_10perc_1$err.rate[ntree]
rf_10perc_2_err.rate <- rf_10perc_2$err.rate[ntree]
rf_10perc_3_err.rate <- rf_10perc_3$err.rate[ntree]

rm(rf_10perc_1,rf_10perc_2,rf_10perc_3) # remove the objects to save memory in R 

##### Best 20% 

names_best_20perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.80))]
names_best_20perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.80))]
names_best_20perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.80))]
names_best_20perc_unique<-unique(c(names_best_20perc_1,names_best_20perc_2,names_best_20perc_3))

# Extract genotypes
genotypes_20perc<-class_data[,colnames(class_data) %in% names_best_20perc_unique]

# Now conduct RF on this subset
rf_20perc_1 = randomForest(x = genotypes_20perc, y = as.factor(class_data$Class), importance=TRUE ,proximity=TRUE, mtry=length(names_best_20perc_unique), ntree=ntree, strata=as.factor(class_data$Class), sampsize=sample_size)
#save(rf_20perc_1,file="rf_20perc_1.Rdata")

rf_20perc_2 = randomForest(x = genotypes_20perc, y = as.factor(class_data$Class), importance=TRUE ,proximity=TRUE, mtry=length(names_best_20perc_unique), ntree=ntree, strata=as.factor(class_data$Class), sampsize=sample_size)
#save(rf_20perc_2,file="rf_20perc_2.Rdata")

rf_20perc_3 = randomForest(x = genotypes_20perc, y = as.factor(class_data$Class), importance=TRUE ,proximity=TRUE, mtry=length(names_best_20perc_unique), ntree=ntree, strata=as.factor(class_data$Class), sampsize=sample_size)
#save(rf_20perc_3,file="rf_20perc_3.Rdata")

# Extract and save the out of bag error rate for this subset of loci
rf_20perc_1_err.rate <- rf_20perc_1$err.rate[ntree]
rf_20perc_2_err.rate <- rf_20perc_2$err.rate[ntree]
rf_20perc_3_err.rate <- rf_20perc_3$err.rate[ntree]

rm(rf_20perc_1,rf_20perc_2,rf_20perc_3) # remove the objects to save memory in R 

##### Best 30% 

names_best_30perc_1<-rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs=0.70))]
names_best_30perc_2<-rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs=0.70))]
names_best_30perc_3<-rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs=0.70))]
names_best_30perc_unique<-unique(c(names_best_30perc_1,names_best_30perc_2,names_best_30perc_3))

# Extract genotypes
genotypes_30perc<-class_data[,colnames(class_data) %in% names_best_30perc_unique]

# Now conduct RF on this subset
rf_30perc_1 = randomForest(x = genotypes_30perc, y = as.factor(class_data$Class), importance=TRUE ,proximity=TRUE, mtry=length(names_best_30perc_unique), ntree=ntree, strata=as.factor(class_data$Class), sampsize=sample_size)
#save(rf_30perc_1,file="rf_30perc_1.Rdata")

rf_30perc_2 = randomForest(x = genotypes_30perc, y = as.factor(class_data$Class), importance=TRUE ,proximity=TRUE, mtry=length(names_best_30perc_unique), ntree=ntree, strata=as.factor(class_data$Class), sampsize=sample_size)
#save(rf_30perc_2,file="rf_30perc_2.Rdata")

rf_30perc_3 = randomForest(x = genotypes_30perc, y = as.factor(class_data$Class), importance=TRUE ,proximity=TRUE, mtry=length(names_best_30perc_unique), ntree=ntree, strata=as.factor(class_data$Class), sampsize=sample_size)
#save(rf_30perc_3,file="rf_30perc_3.Rdata")

# Extract and save the out of bag error rate for this subset of loci
rf_30perc_1_err.rate <- rf_30perc_1$err.rate[ntree]
rf_30perc_2_err.rate <- rf_30perc_2$err.rate[ntree]
rf_30perc_3_err.rate <- rf_30perc_3$err.rate[ntree]

rm(rf_30perc_1,rf_30perc_2,rf_30perc_3) # remove the objects to save memory in R 


# Now combine all of the error rates from the subsets and for all loci to identify a group for the backward purging approach

All_initial_err.rate <- rbind(cbind(rf_all_1_err.rate,rf_all_2_err.rate,rf_all_3_err.rate),
                              cbind(rf_2perc_1_err.rate,rf_2perc_2_err.rate,rf_2perc_3_err.rate),
                              cbind(rf_3perc_1_err.rate,rf_3perc_2_err.rate,rf_3perc_3_err.rate),
                              cbind(rf_4perc_1_err.rate,rf_4perc_2_err.rate,rf_4perc_3_err.rate),
                              cbind(rf_5perc_1_err.rate,rf_5perc_2_err.rate,rf_5perc_3_err.rate),
                              cbind(rf_10perc_1_err.rate,rf_10perc_2_err.rate,rf_10perc_3_err.rate),
                              cbind(rf_20perc_1_err.rate,rf_20perc_2_err.rate,rf_20perc_3_err.rate),
                              cbind(rf_30perc_1_err.rate,rf_30perc_2_err.rate,rf_30perc_3_err.rate))

# Plot error rates for the various subsets
All_initial_err.rate<-data.frame(All_initial_err.rate)
All_initial_err.rate$Number_loci<-c(1174,length(names_best_2perc_unique),length(names_best_3perc_unique),length(names_best_4perc_unique),length(names_best_5perc_unique),length(names_best_10perc_unique),length(names_best_20perc_unique),length(names_best_30perc_unique))
rownames(All_initial_err.rate)<-c("All","Best2%","Best3%","Best4%","Best5%","Best10%","Best20%","Best30%")
All_initial_err.rate$Average<-apply(All_initial_err.rate[,1:3],1,mean)

# Write error rates to file for future reference
write.csv(All_initial_err.rate,file="All_initial_err_rate_classification_1174.csv")

# Plot error rates as well
par(mar=c(5,6,3,3))
plot(All_initial_err.rate$Number_loci,All_initial_err.rate$Average,log="x", pch=19,xlab="Number of Loci", ylab="OOB Error Rate",cex.lab=1.5,cex.axis=1.5)
graph2ppt(file="2_OOB_ER_vs_topSNP", width=8, height=6) 
# Based on this table and plot, the best 4% of loci have the lowest error rate
# As a conservative measure, I'll run backward purging RF with the best 5% loci

#################### Backward purging approach
names_purging <- names_best_5perc_unique

genotypes_purging<-class_data[,colnames(class_data) %in% names_purging]

rf_purging_1 = randomForest(x=genotypes_purging, y = as.factor(class_data$Class), importance=TRUE ,proximity=TRUE, mtry=length(names_purging), ntree=ntree, strata=as.factor(class_data$Class), sampsize=sample_size)
#save(rf_purging_1,file="rf_purging_1.Rdata")
rf_purging_2 = randomForest(x=genotypes_purging, y = as.factor(class_data$Class), importance=TRUE ,proximity=TRUE, mtry=length(names_purging), ntree=ntree, strata=as.factor(class_data$Class), sampsize=sample_size)
#save(rf_purging_2,file="rf_purging_2.Rdata")
rf_purging_3 = randomForest(x=genotypes_purging, y = as.factor(class_data$Class), importance=TRUE ,proximity=TRUE, mtry=length(names_purging), ntree=ntree, strata=as.factor(class_data$Class), sampsize=sample_size)
#save(rf_purging_3,file="rf_purging_3.Rdata")

names_all_iterations<-list()
names_all_iterations[[length(names_purging)]]<-names_purging
error_rate_best<-data.frame(V1=1:length(names_purging),V2=1:length(names_purging),V3=1:length(names_purging))
rownames(error_rate_best)<-1:length(names_purging)
error_rate_best[length(names_purging),] <- c(rf_purging_1$err.rate[ntree],rf_purging_2$err.rate[ntree],rf_purging_3$err.rate[ntree])

for (i in 1:(length(names_purging)-2)){  # RF cannot be conducted with 1 locus, which is why the loop is from 1:length(names_purging)-2
  print(i)
  imp_purging_1<-data.frame(importance(rf_purging_1,type=1))
  imp_purging_2<-data.frame(importance(rf_purging_2,type=1))
  imp_purging_3<-data.frame(importance(rf_purging_3,type=1))
  rownames(imp_purging_1)<-colnames(genotypes_purging)
  colnames(imp_purging_1)<-"Mean_Decrease_Accuracy1"
  rownames(imp_purging_2)<-colnames(genotypes_purging)
  colnames(imp_purging_2)<-"Mean_Decrease_Accuracy2"
  rownames(imp_purging_3)<-colnames(genotypes_purging)
  colnames(imp_purging_3)<-"Mean_Decrease_Accuracy3"
  all_imp<-cbind(imp_purging_1,imp_purging_2,imp_purging_3)
  all_imp$average<-apply(all_imp[,1:3],1,mean)
  dont_keep<-which(all_imp[,'average']==min(all_imp[,'average']))
  if (length(dont_keep)==1) {
    table_keep<-all_imp[-dont_keep,]
  } else {
    table_keep<-all_imp[-dont_keep[sample(length(dont_keep),1)],]
  }
  names_keep<-rownames(table_keep)
  names_all_iterations[[length(names_purging)-i]]<-names_keep
  genotypes_purging<-class_data[,colnames(class_data) %in% names_keep]
  rf_purging_1 = randomForest(x=genotypes_purging, y = as.factor(class_data$Class), importance=TRUE ,proximity=TRUE, mtry=length(names_keep), ntree=ntree, strata=as.factor(class_data$Class), sampsize=sample_size)
  rf_purging_2 = randomForest(x=genotypes_purging, y = as.factor(class_data$Class), importance=TRUE ,proximity=TRUE, mtry=length(names_keep), ntree=ntree, strata=as.factor(class_data$Class), sampsize=sample_size)
  rf_purging_3 = randomForest(x=genotypes_purging, y = as.factor(class_data$Class), importance=TRUE ,proximity=TRUE, mtry=length(names_keep), ntree=ntree, strata=as.factor(class_data$Class), sampsize=sample_size)
  error_rate_best[length(names_purging)-i,] <- c(rf_purging_1$err.rate[ntree],rf_purging_2$err.rate[ntree],rf_purging_3$err.rate[ntree])
}

error_rate_best$Average<-apply(error_rate_best,1,mean)
write.csv(error_rate_best, file="Backward_purging_OOB-ER_classification_best_5perc_unique.csv") # Save the error rates

# Now plot the backward purging results. Omit error rates from one locus since RF cannot be conducted with just one locus
plot(seq(2,nrow(error_rate_best),1),error_rate_best$Average[-c(1)],xlab="Number of Loci", ylab="OOB Error Rate",cex.lab=1.5,cex.axis=1.5,pch=16)
graph2ppt(file="3_OOB_ER_vs_5perc_unique", width=8, height=6) 
# Which group of loci yields the lowest error rate?
which(error_rate_best$Average==min(error_rate_best$Average[-c(1)])) #37 loci have the lowest OOB-ER
best_loci <- which(error_rate_best$Average==min(error_rate_best$Average[-c(1)]))
# Export the names of the predictor loci
#write.csv(names_all_iterations[[best_loci]],file="Predictor_loci_classification_tutorial.csv")
write.csv(names_all_iterations[[37]],file="Predictor_loci_classification_best_10perc_unique_37.csv")

#### Final test and plot
library(dplyr)
library(caTools)
library(stringr)
set.seed(123)
class_data <- read.csv("pop_n_539_outlier_n_1174.csv",row.names=1) #This is the simulated data from Table S8
class_data$Class[class_data$Class == 1] <- 0
class_data$Class[class_data$Class == 2] <- 1
# extrac the dataset for 54 top SNPs
class_data_extrat<-data.frame(class_data$Class, class_data[,colnames(class_data) %in% names_all_iterations[[37]]])
# split for evalucation
split = sample.split(class_data_extrat$class_data.Class, SplitRatio = 2/3)
training_set = subset(class_data_extrat, split == TRUE)
test_set = subset(class_data_extrat, split == FALSE)

classifier = randomForest(x = training_set[-1],
                          y = as.factor(training_set$class_data.Class),
                          mtry=length(names_all_iterations[[37]]), ntree=ntree, random_state = 0, strata=as.factor(training_set$class_data.Class))

y_pred = predict(classifier, newdata = test_set[-1])
cm = table(test_set$class_data.Class, y_pred)
cm 
colnames(cm) <- c("Wild", "Selected")
rownames(cm) <- c("Wild", "Selected")
cm_plot <- data.frame(cm)
colnames(cm_plot) = c("Prediction", "Reference", "Freq")

mytheme = list(
  theme_classic()+
    theme(panel.background = element_blank(),strip.background = element_rect(colour="black", fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),legend.position="", strip.text = element_text(face="plain", size=12),
          axis.text=element_text(face="plain",size=14),axis.title = element_text(face="plain",size=16),plot.title = element_text(face = "bold", hjust = 0.5,size=12))
)

p2 <- cm_plot %>% 
  mutate(Prediction = factor(Prediction, levels = c("Wild", "Selected"))) %>%
  group_by(Reference) %>% 
  mutate(
    total = sum(Freq),
    frac_fill = if_else(Prediction == Reference, Freq / total, Freq / total),
  ) %>%
  ggplot(aes(Prediction, Reference, fill = frac_fill)) +
  geom_tile() +
  geom_text(aes(label = str_c(Freq, ", ", round(frac_fill * 100), "%")), size = 6) +
  scale_fill_gradient(low = "white", high = "#badb33") +
  scale_x_discrete(position = "bottom") +
  geom_tile(color = "black", fill = "black", alpha = 0)+
  mytheme+
  xlab("Predication class")+ylab("True class")
p2

graph2ppt(file="4_cm_topSNP", width=4, height=3) 

comp = data.frame(rownames(test_set), substr(rownames(test_set), 1,4), test_set$class_data.Class, y_pred)
colnames(comp) = c("Ind", "Pop", "True", "Pred")
comp_mismatch <- comp[which(comp$True != comp$Pred),]
library(dplyr)
comp_final <- comp_mismatch %>% 
  dplyr::group_by(Pop) %>% 
  dplyr::summarise(count = n())
comp_final
# cnt_UNC
(1)/sum(comp_final$count) # 0.33 or 33%
# # A tibble: 5 × 2
# Pop   count
# <chr> <int>
# 1 DBX1      1
# 2 DBX2      1
# 3 MEW1      2
# 4 MEW2      1
# 5 UNC2      1


