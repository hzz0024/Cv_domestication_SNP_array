####################################
### format the PopLDdecay output ###
####################################
setwd("~/Dropbox/Mac/Documents/HG/Domestication/17_popLDdecay/")

#PopLDdecay outputs a file with LD results for all pairs of sites, with following columns
#Dist	Mean_r^2	Mean_D'	Sum_r^2	Sum_D'	NumberPairs
#Dist is distance between the pairs of SNPs, Mean_r^2 is the mean squared correlation coefficient (r2) between SNPs
#Sum_r^2 is the sum of squared correlation coefficient (r2) between SNPs, NumberPairs is the number of LD comparisons in each window distance 

rm(list=ls())

format<-function(ngsLD_output, win){
  #ngsLD_output = "chr1.DBW1"
  dat = read.delim(ngsLD_output, header = TRUE, sep='\t')
  plot(dat$X.Dist, y=dat$Mean_r.2, xlab="Distance between SNPs (bp)", ylab=expression(R^2), pch=20, col=rgb(0,0,0,alpha=0.2), main="Decay of Linkage Disequilibrium")
  ### Add a smoothed moving average line
  midpt=rep(0, length(dat$X.Dist))# Set up an empty vector to hold the mean r-squared for each "bin"
  LD.averages=data.frame(dat$X.Dist, dat$Mean_r.2, midpt) # Create the results table (to hold calculations)
  colnames(LD.averages) = c("bins", "r2", "midpt")
  ### Go through all of the "bin" values,
  ### For each one, find the subset of data that
  ### falls in that bin, and get the mean r-squared
  for (i in 1:length(midpt)) {
    if (i < 10) {
      LD.averages$midpt[i]=(LD.averages$bins[i]+LD.averages$bins[i]+10)/2
    }
    else{
      LD.averages$midpt[i]=(LD.averages$bins[i]+LD.averages$bins[i]+win)/2
    }
  }
  return(LD.averages)
}

####################################
### save as rdata format (list)  ###
####################################
out = list()
chr <- seq(0,9)
pop <- c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "DBX1", "DBX2", "DBX3",  "UNC1", "UNC2", "UMFS", "NEH1","NYH1", "NEH2", "MEH2") # "NYH1"
for (c in chr) {
  for (p in pop){
    name0 = paste0("chr", c,".", p)
    out[[name0]] = format(paste0("chr", c,".", p), 100)
  }
}
save(out, file = "./output/PopLDdecay.RData")

####################################
###   start plot for r2 pattern  ###
####################################
library(export)
library(ggplot2) # cut_interval()
setwd("~/Dropbox/Mac/Documents/HG/Domestication/17_popLDdecay/")
# We generated mean r2 within 100 bp bins of distances using r-code.
# 3 columns: bins, r2, midpt
load("./output/PopLDdecay.RData")
pop <- substr(names(out), 6, 10)
meta <- read.csv('./output/meta_pops_V1.csv', stringsAsFactors=T) # added by @HG replace with meta_pops_V1.csv for all pops
site <- meta$Type[match(pop,meta$Site.Abb)]
#site <- meta$Site.Abb[match(pop,meta$Site.Abb)]
#site <- factor(site, levels = c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "UMFS", "MEH2", "NEH1", "NEH2", "DBX1", "DBX2", "DBX3", "UNC1", "UNC2"))
#tiff("output/LDanalysis.makeCurves.500bpBINS.jpg", units="in", width=16, height=12, res=300)
#quartz(width=8,height=5)
par(mfrow=c(1,2),mar=c(2, 2, 2,2)) # increase the first two will show x and y-axis labels, last one is controling the width
#par(mar=c(1,1,1,1))
stats <- c()

#                    MEW1         MEW2       LIW1      LIW2        DBW1       DBW2      NCW1        NCW2
col_gradient <- c( "#0A2C86", "#849cc1",  "#1D92BD", "#8ad5d9", "#93c47d", "#bedbb1", "#a9a9a9", "#dddddd", 
                   #  DBX1       DBX2      DBX3       UNC1        UNC2       UMFS       NYH1       NEH1       NEH2       MEH2
                   "#f9476b", "#fb90a6","#fddae1", "#cf7fbc",  "#e2b2d6", "#fec155", "#C05805" , "#e1bb94", "#fbd0a5", "#b58383")
#########         Dom        wild      
cbPalette <- c( "orange", "skyblue") 

xlim_max <- 500000
ylim_max <- 0.6

#for (i in 1:17) # uncomment this line and comment line blew to plot the decay patter for each population
for (i in 1:2)
{
  plot(1,1,xlim=c(0,xlim_max),ylim=c(0,ylim_max),xlab="distance (bp)", ylab="LD (r^2)",type="n")
  #par(new=TRUE)
  tmp <- out[site==levels(site)[i]]
  
  for(j in 1:length(tmp))
  {
    xbar <- tmp[[j]]
    Cstart <- c(C=0.1)
    CDist <- function(n,C,distance)
    {
      ((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance)))
    }
    n=30 ###### change the value to corresponding sample size
    xbar = xbar[which(xbar$r2 != "Inf"),] ###### to get rid of inf r2 values
    modelC = try(nls(r2 ~ CDist(n,C,midpt), data=xbar, start=Cstart, control=nls.control(maxiter=200)))#
    #error = function(e) modelC=print("oops"))
    if(length(modelC)>1)
    {
      xbar$prd <- predict(modelC, newdata = xbar) #@ newdata = xbar added by HG following the https://stackoverflow.com/questions/33309792/r-predict-function-returning-too-many-values 
      halfdecay = (max(xbar$prd))*0.5
      halfdecaydist <- xbar$midpt[which.min(abs(xbar$prd-halfdecay))]
      dist.LD.is.10perc <- xbar$midpt[which.min(abs(xbar$prd-0.1))] 
      stats <- rbind(stats,data.frame(pop=substr(names(tmp)[j],6,10),
                                      chr=substr(names(tmp)[j],1,4),
                                      halfdecay,
                                      halfdecaydist,
                                      dist.LD.is.10perc))
      lines(xbar$midpt, xbar$prd, col=alpha(cbPalette[i],.5), lwd=1)
    }
    else {stats <- rbind(stats,data.frame(pop=substr(names(tmp)[j],6,10),
                                          chr=substr(names(tmp)[j],1,4),
                                          halfdecaydist=NA,
                                          halfdecaydist=NA,
                                          dist.LD.is.10perc=NA))}
  }
}

write.table(stats,"output/PopLDdecay.stats.500K.100bpBINS.csv",sep=",",row.names=F,quote = F)

graph2ppt(file="output/LDdecay_500K.100bpBINS.pptx", width=12, height=8)
#dev.off()

###########################################
# final plot with NYH1 highlighted as red #
###########################################
load("./output/PopLDdecay.RData")
pop <- substr(names(out), 6, 10)
meta <- read.csv('./output/meta_pops_V2.csv', stringsAsFactors=T) # added by @HG replace with meta_pops_V1.csv for all pops
site <- meta$Type[match(pop,meta$Site.Abb)]
#site <- meta$Site.Abb[match(pop,meta$Site.Abb)]
#site <- factor(site, levels = c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "UMFS", "MEH2", "NEH1", "NEH2", "DBX1", "DBX2", "DBX3", "UNC1", "UNC2"))
#tiff("output/LDanalysis.makeCurves.500bpBINS.jpg", units="in", width=16, height=12, res=300)
#quartz(width=8,height=5)
par(mfrow=c(2,2),mar=c(2, 2, 2,2)) # increase the first two will show x and y-axis labels, last one is controling the width
#par(mar=c(1,1,1,1))
stats <- c()

#                    MEW1         MEW2       LIW1      LIW2        DBW1       DBW2      NCW1        NCW2
col_gradient <- c( "#0A2C86", "#849cc1",  "#1D92BD", "#8ad5d9", "#93c47d", "#bedbb1", "#a9a9a9", "#dddddd", 
                   #  DBX1       DBX2      DBX3       UNC1        UNC2       UMFS       NYH1       NEH1       NEH2       MEH2
                   "#f9476b", "#fb90a6","#fddae1", "#cf7fbc",  "#e2b2d6", "#fec155", "#C05805" , "#e1bb94", "#fbd0a5", "#b58383")
#########         Dom        wild     NYH1 
cbPalette <- c( "red", "orange", "skyblue") 

xlim_max <- 500000
ylim_max <- 0.6

#for (i in 1:17) # uncomment this line and comment line blew to plot the decay patter for each population
for (i in 1:3)
{
  plot(1,1,xlim=c(0,xlim_max),ylim=c(0,ylim_max),xlab="distance (bp)", ylab="LD (r^2)",type="n")
  #par(new=TRUE)
  tmp <- out[site==levels(site)[i]]
  
  for(j in 1:length(tmp))
  {
    xbar <- tmp[[j]]
    Cstart <- c(C=0.1)
    CDist <- function(n,C,distance)
    {
      ((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance)))
    }
    n=30 ###### change the value to corresponding sample size
    xbar = xbar[which(xbar$r2 != "Inf"),] ###### to get rid of inf r2 values
    modelC = try(nls(r2 ~ CDist(n,C,midpt), data=xbar, start=Cstart, control=nls.control(maxiter=200)))#
    #error = function(e) modelC=print("oops"))
    if(length(modelC)>1)
    {
      xbar$prd <- predict(modelC, newdata = xbar) #@ newdata = xbar added by HG following the https://stackoverflow.com/questions/33309792/r-predict-function-returning-too-many-values 
      halfdecay = (max(xbar$prd))*0.5
      halfdecaydist <- xbar$midpt[which.min(abs(xbar$prd-halfdecay))]
      dist.LD.is.10perc <- xbar$midpt[which.min(abs(xbar$prd-0.1))] 
      stats <- rbind(stats,data.frame(pop=substr(names(tmp)[j],6,10),
                                      chr=substr(names(tmp)[j],1,4),
                                      halfdecay,
                                      halfdecaydist,
                                      dist.LD.is.10perc))
      lines(xbar$midpt, xbar$prd, col=alpha(cbPalette[i],.5), lwd=1)
    }
    else {stats <- rbind(stats,data.frame(pop=substr(names(tmp)[j],6,10),
                                          chr=substr(names(tmp)[j],1,4),
                                          halfdecaydist=NA,
                                          halfdecaydist=NA,
                                          dist.LD.is.10perc=NA))}
  }
}

write.table(stats,"output/PopLDdecay.stats.500K.100bpBINS.csv",sep=",",row.names=F,quote = F)

graph2ppt(file="output/LDdecay_500K.100bpBINS.NYH1.pptx", width=10, height=14)
#dev.off()

#############################
###     plot halfdecay    ###
#############################

rm(list=ls())
library(scales) # transparency
library(ggplot2)
library(lmerTest)
library(dplyr)
library(plyr)
meta <- read.csv('./output/meta_pops_V1.csv')
stats <- read.csv('./output/PopLDdecay.stats.500K.100bpBINS.csv')
colnames(stats) <- c("pop","chr","halfdecay","halfdecaydist","dist.LD.is.10perc")
stats$chr <- paste("sca",substr(stats$chr,4,5),sep="")
stats$natnon <- meta$Type[match(stats$pop,meta$Site.Abb)]
stats$pop2 <- meta$Site.Abb[match(stats$pop,meta$Site.Abb)] ### corrected name
stats$pop2 <- factor(stats$pop2, levels = c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "DBX1", "DBX2", "DBX3",  "UNC1", "UNC2", "UMFS", "NEH1","NYH1", "NEH2", "MEH2")) #
# Wilcoxon rank sum test
wilcox.test(stats$halfdecaydist[stats$natnon == "Selected"],stats$halfdecaydist[stats$natnon == "Wild"], alternative = "greater")
# median values
median(stats$halfdecaydist[stats$natnon == "Selected"])
median(stats$halfdecaydist[stats$natnon == "Wild"])

ddply(stats,~pop,summarise,mean=mean(halfdecaydist),sd=sd(halfdecaydist))

m <- lmer(log(halfdecaydist)~chr*natnon+(1|pop),data=stats)
print(anova(m))

#stats$halfdecaydist = stats$halfdecaydist/1000

gg <- ggplot(stats, aes(x=pop2, y = log(halfdecaydist), fill=natnon)) +
  geom_boxplot(width = 0.8, outlier.shape = NA) + 
  geom_jitter(alpha = 0.01, width = 0.2)+
  ylim(min(log(stats$halfdecaydist)), max(log(stats$halfdecaydist))) +
  #scale_y_continuous(breaks = seq(0, 50000, by = 10000))+
  xlab("") +
  ylab("log(Half-decay distance) (bp)") +
  scale_fill_manual(values=alpha(c("orange","skyblue"),.5)) +
  theme(legend.position = "none")+
  theme_classic() + 
  theme(axis.text.x = element_text(angle=45,hjust=1))+
  theme(text = element_text(size = 20)) 
  
gg
graph2ppt(file="output/halfdecay_dis_500K.100bpBINS.pptx", width=7, height=6)

################################
###  halfdecay by chrLength  ###
################################
library(lsmeans)
chrs <- read.delim('output/Table-chr.key.txt')
meta <- read.csv('./output/meta_pops_V1.csv')
stats <- read.csv('output/PopLDdecay.stats.500K.100bpBINS.csv')
colnames(stats) <- c("pop","chr","halfdecay","halfdecaydist","dist.LD.is.10perc")
stats$natnon <- meta$Type[match(stats$pop,meta$Site.Abb)]
stats$natnon <- factor(stats$natnon)
stats$size.kbp <- chrs$size.bp[match(stats$chr,chrs$name1)]/1000

m <- lmer(log(halfdecaydist)~size.kbp*natnon+(1|pop),data=stats)
print(anova(m))

m.interaction <- lm(log(halfdecaydist)~size.kbp*natnon,data=stats)
anova(m.interaction)
m.interaction$coefficients
m.lst <- lstrends(m.interaction, "natnon", var="size.kbp")
m.lst
pairs(m.lst)

m2 <- lmer(log(halfdecaydist)~size.kbp+(1|pop),data=stats[stats$natnon=="Selected",])
print(anova(m2))
m3 <- lmer(log(halfdecaydist)~size.kbp+(1|pop),data=stats[stats$natnon=="Wild",])
print(anova(m3))

#pdf('output/halfdecay.by.chrLength.500bpBINS.pdf',width=8,height=5)
par(mfrow=c(1,2),mar=c(2,2,2,2))
for (i in 1:2)
{
  plot(1,1,xlim=c(min(stats$size.kbp),max(stats$size.kbp)),ylim=c(0,max(stats$halfdecaydist)),xlab="Chromosome length (kbp)", ylab="Half-decay distance (bp)",type="n")
  tmp <- stats[stats$natnon==unique(stats$natnon)[i],]
  points(tmp$size.kbp, log(tmp$halfdecaydist), col=alpha(c("orange","skyblue")[i],0.8),pch=19,cex=1.5)
  abline(lm(halfdecaydist~size.kbp,data=tmp),col=c("orange","skyblue")[i],lwd=2)
}  

graph2ppt(file="output/halfdecay_CHRlength_500K.100bpBINS.pptx", width=10, height=6)


