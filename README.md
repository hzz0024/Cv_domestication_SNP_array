# Consequences of domestication in Eastern oyster

#### Repository for Zhao et al. 2023 "Consequences of domestication in Eastern oyster: insights from whole genome analyses"

Folders below lists the code and raw vcf used for domestication study:

|Folder Name| Contents|
|-----------|---------|
|[R_script](/R_script)| Contains R code for each of the data analyses step |
|[R_plot](/R_plot)| Contains R code and output for Figures in the manuscript|
|[Raw_vcf](/Raw_vcf)| Contains raw vcf for domestication study | 

[1_pcadapt.R](https://github.com/hzz0024/Cv_domestication_SNP_array/blob/main/R_script/1_pcadapt.R)

R code for PCAdapt outlier identification. The “best practice” approach used the naïve clumped SNPs for neutral parameterization and performed the statistical test on the full SNPs

[2_outflank.R](https://github.com/hzz0024/Cv_domestication_SNP_array/blob/main/R_script/2_outflank.R)

R code for OutFLANK outlier identification. The “best practice” approach used the naïve clumped SNPs for neutral parameterization and performed the statistical test on the full SNPs

[3_outlier_process_for_RF.R](https://github.com/hzz0024/Cv_domestication_SNP_array/blob/main/R_script/3_outlier_process_for_RF.R)

Data formatting for the random forest analysis

[4_RF_domestication.R](https://github.com/hzz0024/Cv_domestication_SNP_array/blob/main/R_script/4_RF_domestication.R)

R code for Random forest outlier detection. Random forest (RF) analysis was performed to identify SNPs informatics for wild-selected assignments using the R package randomForest. We followed a two-step backward purging approach to develop a group of SNPs that best predict the response variable (see Brieuc et al. 2015)

[5_diversity_hierfstat.R](https://github.com/hzz0024/Cv_domestication_SNP_array/blob/main/R_script/5_diversity_hierfstat.R)

R code for population diversity indices, including observed (Ho) and expected (He) heterozygosity, and allelic richness (Ar) using hierfstat v 0.5-10 (Goudet 2005). Population differentiation was estimated for all pairs of populations using the Weir & Cockerham estimator of FST (Weir & Clark Cockerham 1984) implemented in hierfstat.

[6_relatedness.R](https://github.com/hzz0024/Cv_domestication_SNP_array/blob/main/R_script/6_relatedness.R)

R code for relatedness measurement using the R package Demerelate v 0.9-3 (Kraemer & Gerlach 2017). The Ritland estimator was used because it has been shown to have the least bias with small sample sizes (Ritland 1996). 

[7_outlier_DAPC_fst.R](https://github.com/hzz0024/Cv_domestication_SNP_array/blob/main/R_script/7_outlier_DAPC_fst.R)

Genetic cluster patterns that were inferred from RF outliers using a discriminant analysis of principal components (DAPC) function implemented in the R package Adegenet (Jombart & Ahmed 2011)

[8_ROH.R](https://github.com/hzz0024/Cv_domestication_SNP_array/blob/main/R_script/8_ROH.R)

Runs of homozygosity (ROH) measurement following a standardized protocol  (Gorssen et al. 2021) using Plink v 1.9 

[9_SROH.R](https://github.com/hzz0024/Cv_domestication_SNP_array/blob/main/R_script/9_SROH.R)

Sum of Runs of homozygosity (SROH) measurement following a standardized protocol  (Gorssen et al. 2021) using Plink v 1.9 

[10_ROH_plot.R](https://github.com/hzz0024/Cv_domestication_SNP_array/blob/main/R_script/10_ROH_plot.R)

ROH plot R script

[11_LDdecay.R](https://github.com/hzz0024/Cv_domestication_SNP_array/blob/main/R_script/11_LDdecay.R)

The LD for each SNP pair estimated using PopLDdecay v. 3.41 (Zhang et al. 2019). The full dataset was used for the LD decay analyses. 

[12_GO_enrichment.R](https://github.com/hzz0024/Cv_domestication_SNP_array/blob/main/R_script/12_GO_enrichment.R)

Gene annotation and GO enrichment test
