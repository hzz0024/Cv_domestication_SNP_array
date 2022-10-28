devtools::install_github("whitlock/OutFLANK")
library(OutFLANK)  # outflank package
library(vcfR)
library(bigsnpr)
library(ggplot2)
library(plyr)
library(stringr)
library(related)
library(bigstatsr)
library(vcfR)
library(RColorBrewer)
library(fields)
library(pcadapt)
library(devtools)
source("manhattan.R")
# http://rstudio-pubs-static.s3.amazonaws.com/305384_9aee1c1046394fb9bd8e449453d72847.html
# set up the working directory
setwd("~/Dropbox/Mac/Documents/HG/Domestication/11_outflank")

################################
#  process the input vcf input #
################################

vcftools  = "/Users/HG/Dropbox/Mac/Documents/HG/Github/BioinfoTools/vcftools_0.1.13/bin/vcftools";
plink  = "/Users/HG/Dropbox/Mac/Documents/HG/Domestication/14_ROH/plink";

system(paste(vcftools," --vcf genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe.recode.vcf --keep pop_n_539.txt --maf 0.05 --recode --recode-INFO-all --out pop_n_539", sep=""))
system(paste(plink, " --vcf pop_n_539.recode.vcf --allow-extra-chr --make-bed --out pop_n_539", sep=""))

# define the Matrix function
toMatrix <- function(G){
  Gm = matrix(nrow=length(G[,1]), ncol=length(G[1,]))
  for(i in seq(length(G[,1]))){
    Gm[i,]=G[i,]
  }
  return(Gm)
}

##########################
# outflank best practice #
##########################

outflank_BP <- function(sub_name, k, ylim){
  #sub_name = "no_DBX1_UNC_MEH_MEW_MEH_n_340"
  # remove .bk file is already exists
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
  
  # part 2 calculate fst 
  Gm = toMatrix(G)
  G_coded = add_code256(big_copy(Gm, type="raw"), code=bigsnpr:::CODE_012)
  # load the population information and calculate Fst
  pop <- read.delim(paste0(sub_name,"_list.txt"),header = T)$Pop
  my_fst <- MakeDiploidFSTMat(toMatrix(G_coded), locusNames = paste0(CHR,'_',POS), popNames = pop)
  
  # Data checks: FST vs. FSTNoCorr
  ggplot(my_fst, aes(x=FST, y=FSTNoCorr))+
    geom_point(alpha=0.1, size=0.5) +
    geom_abline(intercept = 0, slope = 1, col = "red",linetype = "dashed")+
    theme_bw() +
    theme( 
      legend.position="none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    ) 
  
  # Check if FST distribution follows chi-square distribution 
  ggplot(my_fst, aes(x=FSTNoCorr)) + geom_histogram(color="black", fill="white", bins = length(seq(0,0.35, by=0.001)))+
    theme_classic() +
    theme( 
      legend.position="none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    )
  
  # Removing low Heterozygosity variants results in a more chi-square looking FST distribution
  ggplot(my_fst[my_fst$He>0.1,], aes(x=FSTNoCorr)) + geom_histogram(color="black", fill="white", bins = length(seq(0,0.35, by=0.001)))+
    theme_classic() +
    theme( 
      legend.position="none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    )
  
  ### part 3 running outFlank
  # k <- 2 ## Number of pops only focus on wild vs selected lines
  out_ini <- OutFLANK(my_fst, NumberOfSamples=k,
                      RightTrimFraction = 0.05, LeftTrimFraction = 0.35,
                      qthreshold = 0.05, Hmin = 0.1) 
  # Plot results to compare chi-squared distribution vs. actual FST distribution
  OutFLANKResultsPlotter(out_ini, withOutliers = TRUE,
                         NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                           FALSE, RightZoomFraction = 0.05, titletext = NULL)
  # zoom-in (zoom=TRUE) to compare chi-squared distribution vs. actual FST distribution
  OutFLANKResultsPlotter(out_ini, withOutliers = TRUE,
                         NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                           TRUE, RightZoomFraction = 0.15, titletext = NULL)
  
  ### part 4 LD pruning
  # Evaluating OutFLANK with pruned data
  # plot(my_fst$He[which_pruned], my_fst$FST[which_pruned])
  Fst_LD <- my_fst[which_pruned,] 
  #dim(Fst_LD)
  
  # after LD and He filtering
  Fst_LD_He <- Fst_LD[Fst_LD$He>0.1,]
  #plot(Fst_LD_He$He, Fst_LD_He$FST)
  
  # # Trimming without He trimming
  # out_trim1 <- OutFLANK(Fst_LD, NumberOfSamples=k)
  # OutFLANKResultsPlotter(out_trim1, withOutliers = TRUE,
  #                        NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
  #                          TRUE, RightZoomFraction = 0.15, titletext = NULL)
  # Trimming loci with low He
  out_trim <- OutFLANK(Fst_LD_He, NumberOfSamples=k)
  head(out_trim$results)
  plot(out_trim$results$He, out_trim$results$FST)
  OutFLANKResultsPlotter(out_trim, withOutliers = TRUE,
                         NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                           FALSE, RightZoomFraction = 0.1, titletext = NULL)
  OutFLANKResultsPlotter(out_trim, withOutliers = TRUE,
                         NoCorr = TRUE, Hmin = 0.1, binwidth = 0.005, Zoom =
                           TRUE, RightZoomFraction = 0.15, titletext = NULL)
  
  # part 5 best practice
  out_BP <- pOutlierFinderChiSqNoCorr(my_fst, Fstbar = out_trim$FSTNoCorrbar, 
                                      dfInferred = out_trim$dfInferred, qthreshold = 0.05, Hmin=0.1)
  head(out_BP)
  # retain the SNPs with He > Hmin (0.1)
  Hmin = 0.1
  keepers = which((out_BP$FSTNoCorr > 0) & (out_BP$He >= Hmin))
  out_BP = out_BP[keepers,]
  # count how many outlier candidates with q-value < 0.05
  length(out_BP$OutlierFlag[out_BP$OutlierFlag==TRUE])
  
  # notice how the output is ordered differently
  my_out <- out_BP$OutlierFlag==TRUE
  # plot(out_BP$He, out_BP$FST, pch=19, col=rgb(0,0,0,0.1))
  # points(out_BP$He[my_out], out_BP$FST[my_out], col="blue")
  
  # plot the p-values
  hist(out_BP$pvaluesRightTail, xlab='Raw p-values')
  # plot the q-values
  hist(out_BP$qvalues)
  
  # Histogram of P-values 
  ggplot(out_BP, aes(x=pvaluesRightTail)) + geom_histogram(color="black", fill="white", bins = length(seq(0,0.35, by=0.001)))+
    theme_classic() +
    theme( 
      legend.position="none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    )
  
  #check how many outliers with q-value < 0.05
  cnt <- sum(out_BP$qvalues<0.05, na.rm=TRUE)
  print(paste0("number of outliers from best practise in ", sub_name, " is: ", cnt))
  
  # plot(out_ini$results$He, out_ini$results$FST, pch=20, col="grey")
  # points(out_ini$results$He[out_ini$results$qvalues<0.05], out_ini$results$FST[out_ini$results$qvalues<0.05], pch=21, col="blue")
  # ### Note how OutFLANK identifies potential outliers at He < 0.1, even though
  # ### these loci were excluded in the trimming algorithm
  # top_candidates <- out_ini$results$qvalues<0.05 & out_ini$results$He>0.1
  # plot(out_ini$results$He, out_ini$results$FST, pch=20, col="grey")
  # points(out_ini$results$He[top_candidates], out_ini$results$FST[top_candidates], pch=21, col="blue")
  
  # start to make Mahattan plot
  CHR = sapply(strsplit(out_BP$LocusName, "_"), "[[", 1)
  POS = sapply(strsplit(out_BP$LocusName, "_"), "[[", 2)
  daf = data.frame(CHR=CHR, POS=POS, SNP=paste0(CHR,"_",POS), Ps=out_BP$pvalues, qvalue=out_BP$qvalues)
  daf = daf[complete.cases(daf), ]
  daf$CHR = as.numeric(daf$CHR)
  daf$POS = as.numeric(daf$POS)
  outlier_daf <- daf[which(daf$qvalue < 0.05),]
  outlier_SNP <- outlier_daf$SNP
  outlier <- data.frame(CHR=outlier_daf$CHR, POS=outlier_daf$POS, POS=outlier_daf$POS, snp_id=outlier_daf$SNP)
  write.table(outlier, file = paste0(sub_name, "_outflank_BP_q05_n_", dim(outlier)[1], ".bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  don <- daf %>% 
    dplyr::group_by(CHR) %>% 
    dplyr::summarise(chr_len=max(POS)) %>% 
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    dplyr::select(-chr_len) %>% 
    left_join(daf, ., by=c("CHR"="CHR")) %>%
    arrange(CHR, POS) %>%
    mutate(BPcum=POS+tot) %>%
    mutate(is_highlight=ifelse(SNP %in% outlier_SNP, "yes", "no"))
  axisdf = don %>% dplyr::group_by(CHR) %>% dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  write.table(don, file = paste0(sub_name, "_qs.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  outlier_plot <- ggplot(don, aes(x=BPcum, y=-log10(qvalue))) +
    geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=0.5) +
    scale_color_manual(values = rep(c("grey", "#4F4F4F"), 22 )) +
    scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0.2)) + # remove space between plot area and x axis
    geom_point(data=subset(don, is_highlight=="yes"), color="red", size=0.8) +
    xlab("Chromosome") + 
    ylab("-log10(qvalue)")+
    theme_bw() +
    theme( 
      legend.position="none",
      #panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    ) +
    ylim(0, ylim) +  # for better visualization, need to remove this y scale limit for complete plot
    labs(caption=paste0(dim(don[-log10(don$qvalue) > ylim,])[1], " out of ", dim(outlier)[1], " outliers not shown due to Y-axis setting < ", ylim))
  ggsave(paste0("figs/", sub_name, ".BP.jpg"), outlier_plot, width = 9, height = 7)
}

outflank_BP("pop_n_539", 2, 2)

          