install.packages("dartR")
install.packages("Demerelate")
library(Demerelate)
library(dartR)
library(adegenet)
library(vcfR)

setwd("~/Dropbox/Mac/Documents/HG/Domestication/15_relatedness/Demerelate_1k_539")
# load the population information
pop_info <- read.table("pop_539_sample_list.txt", header=TRUE, sep="\t", stringsAsFactors = TRUE)
pop_info$Pop_correct = factor(pop_info$Pop_correct, levels=c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "DBX1", "DBX2", "DBX3", "UNC1", "UNC2", "UMFS", "NYH1", "NEH1", "NEH2", "MEH2"))
# load vcf file
vcf_file = "genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_neutral_pruned_1K.recode.vcf"
#vcf_file = "genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe_pruned_neutral_1K.recode.vcf"
vcf <- read.vcfR(vcf_file, verbose = FALSE)
Mydata1 <- vcfR2genind(vcf)
Mydata1@pop <- pop_info$Pop_correct
Mydata1

Mydata1.gl <- (gi2gl(Mydata1))
Mydata1.de <- gl2demerelate(Mydata1.gl, verbose = 1) # https://rdrr.io/cran/dartR/man/gl2demerelate.html
dim(Mydata1.de)

#Loci.test(Mydata1.de, bt=1000, ref.pop=NA, object=TRUE, value="ritland", file.output=TRUE)
#pop.relate <- Emp.calc(Mydata1.de, value="ritland", ref.pop="NA")
Demerelate(Mydata1.de, ref.pop="NA",
           tab.dist="NA", Fis=FALSE, NA.rm=FALSE,
           object=TRUE, pairs=1000, iteration=100,
           value="ritland", file.output=TRUE)

#################################
######## Relatdeness plot #######
#################################
library(export)
library(beanplot)
library(ggplot2)
library(hrbrthemes)
setwd("~/Dropbox/Mac/Documents/HG/Domestication/15_relatedness/Demerelate_1k_539/Demerelate-0.4453/")

data_process <- function(headname, pop){
  #headname = "DBW1"
  name = paste0("Empirical.relatedness.",headname, ".txt")
  DT = read.delim(name, header = TRUE, sep=' ')
  dat <- data.frame(relat=DT$ritland, Population=pop) # @change
  dat$relat <- as.numeric(dat$relat)
  return(dat)
}
# build data.frame
plotdat = NULL
for(name in c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "DBX1", "DBX2", "DBX3", "UNC1", "UNC2", "UMFS", "NYH1", "NEH1", "NEH2", "MEH2")){
  dat <- data_process(name, name)
  if(is.null(plotdat)){
    plotdat = dat}
  else{
    plotdat = rbind(plotdat, dat)}
}

# summarize mean and sd for relatedness
library(dplyr)
plotdat_group  <- dplyr::group_by(plotdat, Population) # group data by party *and* congress
mean_sd <- dplyr::summarize(plotdat_group,
          mean = mean(relat),
          sd = sd(relat))
write.table(mean_sd, file = "relatedness_1K_mean_sd.txt", sep = "\t", quote = FALSE,
            row.names = T, col.names = F)

# Mann Whitney U / Wilcoxon Rank-Sum Test 
wild = NULL
for(name in c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "DBX1")){
  dat <- data_process(name, name)
  if(is.null(wild)){
    wild = dat}
  else{
    wild = rbind(wild, dat)}
}

wild_mean <- wild %>%
  group_by(Population) %>%
  dplyr::summarize(Mean = mean(relat, na.rm=TRUE))

sel = NULL
for(name in c("DBX2", "DBX3", "UNC1", "UNC2", "UMFS", "NYH1", "NEH1", "NEH2", "MEH2")){
  dat <- data_process(name, name)
  if(is.null(sel)){
    sel = dat}
  else{
    sel = rbind(sel, dat)}
}

sel_mean  <- sel %>%
  group_by(Population) %>%
  dplyr::summarize(Mean = mean(relat, na.rm=TRUE))

wilcox.test(wild_mean$Mean,sel_mean$Mean)

plotdat$Population = factor(plotdat$Population, levels=c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "DBX1", "DBX2", "DBX3", "UNC1", "UNC2", "UMFS", "NYH1", "NEH1", "NEH2", "MEH2"))
jpeg("Relatedness_plot.jpg", width = 10, height = 5, units = 'in', res = 150)
p <- ggplot(plotdat, aes(x=Population, y=relat, fill=Population)) +
  geom_boxplot(width = 0.8, outlier.size = 0.1) + 
  #geom_jitter(alpha = 0.03, width = 0.2)+
  theme_classic() + 
  theme(axis.text.x = element_text(angle=45,hjust=1)) +
  #stat_summary(fun = mean, geom='point', shape=20, size=2, color='red', fill='red') +
  labs(x=NULL, y = "Ritland estimate of relatedness") 
#p +  scale_fill_manual(values=col) + 
  #scale_y_continuous(limits=c(0, 0.25)) +
#  theme(text = element_text(size=20),
#        legend.position = "none")
p
graph2ppt(file="relatdeness",width=10,height=6)

# Ne vs Nb
Ne_Nb_mean <- read.table("Ne_summary.txt", header=TRUE, sep="\t", stringsAsFactors = TRUE)
wilcox.test(Ne_Nb_mean$Ne,Ne_Nb_mean$Nb)
