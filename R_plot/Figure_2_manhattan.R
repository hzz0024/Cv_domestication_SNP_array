rm(list=ls())
library(ggplot2)
library(dplyr)
library(grid)
library(gtable)
library(grDevices)
library(patchwork)
library(wesanderson)

setwd("~/Dropbox/Mac/Documents/HG/Domestication/Manuscript/Figure_539/Figure_2")

manhattan_plot_up <- function(pname, outlier_color, x_label){
  #pname="pop_n_539_qs_pcadapt.txt"
  #outlier_color = "red"
  #x_label = "pcadapt"
  don = read.delim(pname, header = TRUE, sep='\t')
  axisdf = don %>% dplyr::group_by(CHR) %>% dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  sig_data <- don %>% 
    subset(is_highlight == "yes")
  notsig_data <- don %>% 
    subset(is_highlight == "no") %>%
    group_by(CHR) %>% 
    sample_frac(1)
  
  plot_dat <- bind_rows(sig_data, notsig_data)
  
  ggplot(plot_dat, aes(x=BPcum, y=-log10(qvalue))) +
    # Show all points
    # Show all points
    geom_point(aes(color=as.factor(CHR)), alpha=0.5, size=0.2) +
    scale_color_manual(values = rep(c("grey90", "grey50"), 22 )) +
    
    #geom_point(aes(colour=factor(seq(1, length(points.col)))), alpha=0.9, size=0.8) +
    #scale_color_manual(values = points.col) +
    # custom X axis:
    scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
    # Custom the theme:
    # Add highlighted points
    geom_point(data=subset(plot_dat, is_highlight=="yes"), color=outlier_color, alpha=0.3, size=0.3) +
    facet_grid(cols = vars(CHR), scales = "free_x", space="free_x") +
    ylab(as.expression(bquote(atop(.(x_label), -log[10]~qvalue))))+
    xlab("Chromosome")+
    #ylim(-0.5, 0.5)+
    theme_bw() +
    theme(
      axis.title=element_text(size=18),
      text=element_text(color ="black",size = 16),
      legend.position="none",
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      axis.title.x=element_blank(), # remove the x axis title
      panel.border = element_blank(),
      axis.line.x = element_line(color="black", size = 0.5),
      axis.line.y = element_line(color="black", size = 0.5),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    )
}

manhattan_plot_pcadapt <- function(pname, outlier_color, x_label){
  #pname="pop_n_477_outflank_BP_q05_n_276_outflank.bed"
  #outlier_color = "red"
  #x_label = "pcadapt"
  don = read.delim(pname, header = TRUE, sep='\t')
  axisdf = don %>% dplyr::group_by(CHR) %>% dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  sig_data <- don %>% 
    subset(is_highlight == "yes")
  notsig_data <- don %>% 
    subset(is_highlight == "no") %>%
    group_by(CHR) %>% 
    sample_frac(1)
  
  plot_dat <- bind_rows(sig_data, notsig_data)
  
  ggplot(plot_dat, aes(x=BPcum, y=-log10(qvalue))) +
    # Show all points
    # Show all points
    geom_point(aes(color=as.factor(CHR)), alpha=0.5, size=0.2) +
    scale_color_manual(values = rep(c("grey90", "grey50"), 22 )) +
    
    #geom_point(aes(colour=factor(seq(1, length(points.col)))), alpha=0.9, size=0.8) +
    #scale_color_manual(values = points.col) +
    # custom X axis:
    scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
    # Custom the theme:
    # Add highlighted points
    geom_point(data=subset(plot_dat, is_highlight=="yes"), color=outlier_color, alpha=0.3, size=0.3) +
    facet_grid(cols = vars(CHR), scales = "free_x", space="free_x") +
    ylab(as.expression(bquote(atop(.(x_label), -log[10]~qvalue))))+
    #xlab("Chromosome")+
    ylim(0, 5)+
    theme_bw() +
    theme( 
      text=element_text(color ="black",size = 12),
      legend.position="none",
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      axis.title.x=element_blank(), # remove the x axis title
      panel.border = element_blank(),
      axis.line.x = element_line(color="black", size = 0.5),
      axis.line.y = element_line(color="black", size = 0.5),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    )
}

manhattan_plot_bottom <- function(pname, outlier_color, x_label){
  #pname="Dom_Wild.qs.txt"
  #outlier_color = "red"
  #x_label = "Window scan"
  don = read.delim(pname, header = TRUE, sep='\t')
  axisdf = don %>% dplyr::group_by(CHR) %>% dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  sig_data <- don %>% 
    subset(is_highlight == "yes")
  notsig_data <- don %>% 
    subset(is_highlight == "no") %>%
    group_by(CHR) %>% 
    sample_frac(1)
  
  plot_dat <- bind_rows(sig_data, notsig_data)
  
  ggplot(plot_dat, aes(x=BPcum, y=-log10(qvalue))) +
    # Show all points
    # Show all points
    geom_point(aes(color=as.factor(CHR)), alpha=0.5, size=0.2) +
    scale_color_manual(values = rep(c("grey90", "grey50"), 22 )) +
    
    #geom_point(aes(colour=factor(seq(1, length(points.col)))), alpha=0.9, size=0.8) +
    #scale_color_manual(values = points.col) +
    # custom X axis:
    scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
    # Custom the theme:
    # Add highlighted points
    geom_point(data=subset(plot_dat, is_highlight=="yes"), color=outlier_color, alpha=0.3, size=0.3) +
    facet_grid(cols = vars(CHR), scales = "free_x", space="free_x") +
    ylab(as.expression(bquote(atop(.(x_label), -log[10]~qvalue))))+
    xlab("Chromosome")+
    #ylim(-0.5, 0.5)+
    theme_bw() +
    theme( 
      axis.title=element_text(size=18),
      text=element_text(color ="black",size = 16),
      legend.position="none",
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      #axis.title.x=element_blank(), # remove the x axis title
      panel.border = element_blank(),
      axis.line.x = element_line(color="black", size = 0.5),
      axis.line.y = element_line(color="black", size = 0.5),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    )
}

pal1 <- wes_palette("Zissou1")
pal2 <- wes_palette("GrandBudapest1")
pal3 <- wes_palette("GrandBudapest2")
pal4 <- wes_palette("Royal2")
col_pal <- c(pal1[1], pal1[2], pal3[4])

pname="pop_n_539_qs_pcadapt.txt"
don = read.delim(pname, header = TRUE, sep='\t')
unique(don$tot)
# pcadapt tot
# 0
# 65642311
# 127393718
# 204426069
# 263178747
# 361721604
# 412815720
# 470600906
# 546541227
# 650669346
data_vline <- read.delim("pcadapt_RF_outlier_coordinates.txt", header = TRUE, sep='\t')
data_vline$Color <- as.character(data_vline$Color)
p1 <- manhattan_plot_up(pname="pop_n_539_qs_pcadapt.txt", outlier_color = "red", x_label = "pcadapt")+
  geom_vline(data = data_vline, mapping=aes(xintercept = BPcum), alpha=0.4, color=data_vline$Color)

pname="pop_n_539_qs_outflank.txt"
don = read.delim(pname, header = TRUE, sep='\t')
unique(don$tot)
# tot
# 0
# 65642311
# 127393718
# 204426069
# 263178747
# 361721604
# 412815720
# 470600906
# 546541227
# 650669346

data_vline <- read.delim("outflank_RF_outlier_coordinates.txt", header = TRUE, sep='\t')
p2 <- manhattan_plot_bottom(pname="pop_n_539_qs_outflank.txt", outlier_color = "red", x_label = "OutFLANK")+
  geom_vline(data = data_vline, mapping=aes(xintercept = BPcum), alpha=0.4, color=data_vline$Color)

# tot
# 0
# 65639000
# 127386000
# 203565000
# 262314000
# 359529000
# 410620000
# 468401000
# 544326000
# 648389000
#data_vline <- read.delim("Zfst_common_outlier_coordinates.txt", header = TRUE, sep='\t')
#p3 <- manhattan_plot_bottom_Zfst("Dom_Wild.qs.txt", outlier_color = "red", x_label = "Window scan")+
#  geom_vline(data = data_vline, mapping=aes(xintercept = BPcum), alpha=0.2, color=data_vline$Color)

layout <-" 
AAAAAA
BBBBBB
"

top.mar=0.2
right.mar=0.2
bottom.mar=0.2
left.mar=0.2

p <- p1+p2 + plot_layout(design=layout, guides='collect') +
  theme(plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),units="inches"))
  
#grid::grid.draw(grid::textGrob("Chromosome", y=0.015, gp = gpar(col = "black", fontsize = 16)))
#grid::grid.draw(grid::textGrob(ylab, x=0.02, gp = gpar(col = "black", fontsize = 16)))

tiff("figure2_part1.tiff", units="in", width=16, height=5, res=300)
p
dev.off()

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
system(paste(vcftools," --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.recode.vcf --bed pop_n_477_pcadapt_outflank.shared.10k.outlier.igv_Dom_Wild.sliding.zfst.outlier.merged.igv.outlier.merged.igv --recode --recode-INFO-all --out pop_n_477_outliers", sep=""))
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
  
  gp.label$X <- c(dfx-0,dfx-ddelta/2,dfx-ddelta*3/4,dfx-ddelta,dfx-ddelta,dfx-ddelta*3/4,dfx-ddelta/2,dfx-0,
                  dfx-0,dfx-ddelta/3,dfx-ddelta/2,dfx-ddelta*3/4,dfx-ddelta,dfx-ddelta*3/4,dfx-ddelta/2,dfx-ddelta/3,dfx-0)
  
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
  chrom_outlier_GP_plot(i,"pop_n_477_outliers.recode.vcf.gz", 509)
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

#chrom_outlier_GP_plot_no_label(2,"pop_n_477_outliers.recode.vcf.gz", 509)

for(i in c(2,5)){
  chrom_outlier_GP_plot_no_label(i,"pop_n_477_outliers.recode.vcf.gz", 509)
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
p/gp.plot.1.out + gp.plot.2.out + gp.plot.5.out + plot_layout(design=layout, guides="collect")
dev.off()

layout <-" 
BBBBCCDDDDDDDD
BBBBCCDDDDDDDD
BBBBCCDDDDDDDD
BBBBCCDDDDDDDD
BBBBCCDDDDDDDD
BBBBCCDDDDDDDD
"

tiff("figure2_part3.tiff", units="in", width=12, height=10, res=300)
gp.plot.1.out + gp.plot.2.out + gp.plot.5.out + plot_layout(design=layout, guides="collect")
dev.off()


