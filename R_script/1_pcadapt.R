setwd("~/Dropbox/Mac/Documents/HG/Domestication/04_pcadapt")
library(OutFLANK)
library(vcfR)
library(bigsnpr)
library(ggplot2)
library(dplyr)
# to install export, one need to isntall X11(https://www.xquartz.org/), and then type devtools::install_github("tomwenseleers/export")
library(export)
library(qvalue)
library(stringr)
library(pcadapt)
library(export)
library(ggrepel)
source("manhattan.R")

#############################################
## Prepare dataset for PCAdapt running ######
#############################################
vcftools  = "/Users/HG/Dropbox/Mac/Documents/HG/Github/BioinfoTools/vcftools_0.1.13/bin/vcftools";
plink  = "/Users/HG/Dropbox/Mac/Documents/HG/Domestication/14_ROH/plink";

system(paste(vcftools," --vcf genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe.recode.vcf --keep pop_n_539.txt --maf 0.05 --recode --recode-INFO-all --out pop_n_539", sep=""))
system(paste(plink, " --vcf pop_n_539.recode.vcf --allow-extra-chr --make-bed --out pop_n_539", sep=""))

#############################################
## Function for PCadapt best practise (BP) ##
#############################################

toMatrix <- function(G){
  Gm = matrix(nrow=length(G[,1]), ncol=length(G[1,]))
  for(i in seq(length(G[,1]))){
    Gm[i,]=G[i,]
  }
  return(Gm)
}

sub_name = "pop_n_539"
# colours   DOM       WILD     
col <- c("#ec9c9d", "#8ea4bf")
alpha = 0.05

pcadapt_BP <- function(sub_name, alpha){
  
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
  system(paste(vcftools," --vcf ",sub_name,".recode.vcf", " --snps ", sub_name, "_pruned_SNP.txt", " --remove remove.txt --recode --recode-INFO-all --out ", sub_name, "_pruned", sep=""))
  system(paste(plink, " --vcf ",sub_name,"_pruned.recode.vcf", " --allow-extra-chr --make-bed --out ", sub_name,"_pruned", sep=""))
  
  # part 2 Population structure patterns from thinned data
  file <- read.pcadapt(paste0(sub_name,"_pruned.bed"), type = "bed")
  bim_file = read.delim(paste0(sub_name,"_pruned.bim"), header = FALSE, sep='\t')
  # Scree plot
  x <- pcadapt(input = file, K = 10)
  # The eigenvalues that correspond to random variation lie on a straight line whereas the ones that correspond to population structure lie on a steep curve. It is recommended to keep PCs that correspond to eigenvalues to the left of the straight line (Cattell’s rule).
  scree_plot <- plot(x, option = "screeplot")
  # Score plot: another option to choose the number of PCs is based on the ‘score plot’ that displays population structure.
  sample_file = read.delim(paste0(sub_name,"_list.txt"), header = TRUE, sep='\t')
  poplist.names <- sample_file$Pop
  # PCA plot
  PCA_plot <- plot(x, option = "scores", i=1, j=2, pop=poplist.names, col=col)
  
  # part 3 Best practise for outlier detection
  Gm = toMatrix(G)
  G_coded = add_code256(big_copy(Gm, type="raw"), code=bigsnpr:::CODE_012)
  Gm = Gm[,snp_MAF(G_coded) != 0] # eliminate monomorphic SNPs
  G_coded = add_code256(big_copy(Gm, type="raw"), code=bigsnpr:::CODE_012)
  #res_pcadapt <- snp_gc(snp_pcadapt(G_coded, U.row = newpc$u[,1:2])) # use PC1 and PC2 for outlier identification
  res_pcadapt <- snp_gc(snp_pcadapt(G_coded, U.row = newpc$u[,1:6])) # use defined PCs for outlier identification
  res_p <- predict(res_pcadapt,log10 = F)
  hist(res_p, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
  # Choosing a cutoff for outlier detection
  qval <- qvalue(res_p)$qvalues
  outliers <- which(qval < alpha)
  print(paste0("number of outliers from best practise in ", sub_name, " is: ", length(outliers)))
  # start to make Mahattan plot
  CHR = CHR[which(snp_MAF(G_coded) != 0)]
  POS = POS[which(snp_MAF(G_coded) != 0)]
  daf = data.frame(CHR=CHR, POS=POS, SNP=paste0(CHR,"_",POS), Ps=res_p, qvalue=qval)
  outlier_SNP <- paste0(CHR[outliers],'_',POS[outliers])
  outlier_daf <- data.frame(CHR=CHR[outliers], POS=POS[outliers], POS=POS[outliers], snp_id=outlier_SNP)
  write.table(outlier_daf, file = paste0(sub_name, "_PCAdapt_BP_q05_n_", length(outliers), ".bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  #jpeg("Mahattan_Ind509_best_practice_outlier_PC1-2_10K.jpg", width = 16, height = 9, units = 'in', res = 300)
  # customize mahattan plot see https://www.r-graph-gallery.com/101_Manhattan_plot.html
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
    ylim(0, 9) +  # for better visualization, need to remove this y scale limit for complete plot
    labs(caption=paste0(dim(don[-log10(don$qvalue) > 15,])[1], " out of ", length(outliers), " outliers not shown due to Y-axis setting < 15"))
  
  bottom_row <- plot_grid(scree_plot, PCA_plot, labels = c('A', 'B'), align = 'h', rel_widths = c(1, 1.3))
  final_plot <- plot_grid(bottom_row, outlier_plot, labels = c('', 'C'), ncol = 1, rel_heights = c(1, 1.2))
  ggsave(paste0("figs/", sub_name, ".jpg"), final_plot, width = 9, height = 7)
}

pcadapt_BP("pop_n_539", 0.05)



