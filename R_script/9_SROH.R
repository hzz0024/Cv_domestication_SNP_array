#====================================

#====================================
#Rscript for analysis of Runs Of Homozygosity in PLINK (Gorssen et al. (submitted at GSE)
#====================================

#====================================

#This script is distributed under the license of Creative Commons CC-BY 4.0
#More information can be found at https://creativecommons.org/licenses/by/4.0/deed

#If you use this script in your analysis, cite the according paper:
#Gorssen, W., Meyermans, R., Janssens, S., & Buys, N. (2020). An ROH island repository reveals signatures of selection for 442 populations in 8 different livestock and pet species

#                          (Full citation will be added when available)

#More information on the used quality control and plink settings can be found in:
#Meyermans, R., Gorssen, W., Buys, N., & Janssens, S. (2020). How to study runs of homozygosity using PLINK? A guide for analyzing medium density SNP data in livestock and pet species. BMC genomics, 21(1), 1-14.

#For any questions concerning this script, please contact the authors at: wim.gorssen@kuleuven.be or roel.meyermans@kuleuven.be

#======
#====
#Before you run this script
#====
#======

#The following script requires a .bim, .bed and .fam file of a set of individuals from the same species with only 
#autosomal chromosomes (no unmapped, nor sex chromosomes). The script handles every FID in your dataset as different breeds/populations.
#The script loops for every FID and performs Quality control and Runs of Homozygosity analysis population specific
#First, QC is performed per breed/population/FID. Hereafter, ROH analyses is performed per breed/population/FID and
#results saved in the specified directory.
#Directories and names of your dataset need to match the format used in this script, if you want it to run automatically and write  results.
#
#Your dataset needs to be in a folder with following extension: ./Data/species/  (fill in name of species)
#and your dataset has to be given the following standard name: speciesdatasetnumber_beforeQC.bim
#for example:"./Data/pig/pig1_beforeQC.bim"
#Both fam, bim and bed files should be present in this file
#The output will appear in 2 standard folders which will be created in this script: 
#1.  ./Output/species  (fill in name of species): for example: ./output/pig
#2.  ./Output/all_chr_ROH/species  (fill in name of species): for example: ./output/all_chr_ROH/pig



# Indicate where the executable of plinkv1.9 is on your device 
plink  = "/Users/HG/Dropbox/Mac/Documents/HG/Domestication/14_ROH/plink"; ## Needs to be the full path to the program ./plink

#set working directory: this is the path with the directory specified above

setwd("~/Dropbox/Mac/Documents/HG/Domestication/14_ROH/All_n_539_formal_run_SROH//")

#===
#Load and install following functions and packages
#===

#install.packages('data.table')
library(data.table)
#install.packages('plyr')
library(plyr)
#install.packages('xlsx') # xlsx may not be loaded for macOS m1, need to install the java aarch64 https://download.oracle.com/java/17/archive/jdk-17.0.1_macos-aarch64_bin.dmg
library(xlsx)
#install.packages('qqman')
library(qqman)
#install.packages('ggplot2')
library(ggplot2)

#===
#Function to visualize pairs plots: panel.cor, used in the script
#=== 

panel.cor <- function(x, y, digits=2, prefix="", cex.cor) 
{
  usr <- par("usr"); on.exit(par(usr)) 
  par(usr = c(0, 1, 0, 1)) 
  r <- (cor(x, y,use = "complete.obs")) 
  txt <- format(c(r, 0.123456789), digits=digits)[1] 
  txt <- paste(prefix, txt, sep="") 
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt) 
  
  test <- cor.test(x,y) 
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = anyNA , 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " ")) 
  
  text(0.5, 0.5, txt, cex = cex) 
  text(.8, .8, Signif, cex=cex, col=2) 
}  

#====================================

#####################################
#1. Specify Parameter settings
#####################################

#====================================



#===
#quality control settings
#===


mind<-0.10  #Mind: tolerated call rate of individuals => a mind of 0.10 corresponds to minimum 90% call-rate
maf<-0 #Minor allele frequency: Choose 0 if you don't want to perform maf (recommended for ROH analyses)
hwe<-0 #Hardy weinberg equilibrium: choose  0 if you don't want to perform hwe (recommended for ROH analyses)
geno<-0.05 #geno: call rate of SNPs => a geno of 0.05 corresponds to minimum 95% call-rate
#No LD pruning in this script (as recommended by Meyermans et al. 2020)

#===
#Specify ROH settings: divided in 2 classes (more info: https://www.cog-genomics.org/plink/1.9/ibd#homozyg)
#===

#=
#1. Window settings:
#=

#these settings are used to define window length to scan genome and
#assign 'scores' to individual snps. These scores determine if a snp is called in a 
#segment or not via threshold parameter

het<-1 #--homozyg-window-het AND --homozyg-het (it is important to specify both in PLINK; in this script het is used equal for both flags)
mis<-5 #--homozyg-window-missing
windowsnp<-0 ##--homozyg-window-snp => this parameter will be adjusted within script to parameter l!
windowthreshold<-0.05 #--homozyg-window-threshold => this parameter will be adjusted within script (t=floor(Nout+1L,3)!
Nout<-2 #The number of outer SNPs you don't want to allow for ROH-analysis (calculation of windowthreshold)

#=
#2.Segment settings:
#=

#After segment calling via window settings: extra requirements for segments!
kb<-200 #--homozyg-kb: the minimum segment length in kb
l<-0    #--homozyg-snp => this parameter will be adjusted within script!
gap<-1000 #--homozyg-gap: the maximal gap in kb between two SNPs within an ROH
density<-50 #--homozyg-density #the minimal average density of an ROH (expressed as 1SNP/xxx kb, so in our example: minimal 1SNP/150kb)
bin <- 1000000



#==
#Specify minimum number of animals that have to be present in population to perform ROH-analysis for that population (set to 1 if you want to do it for all populations)
#==

minimum_number_animals_population<-15


#====================================

#####################################
#2. Specify species to use and dataset
#####################################

#====================================



#===  
#choose species you are using: We optimized this script for the following species: "cat","chicken","dog","horse","cow","goat","sheep","pig", "water buffalo"
#===  

species<-"oyster"


#===  
#choose dataset number for this species you want to use
#===  

dataset<-1


#=====================
#If you did everything right, from this point on, you should not have to change any code anymore. 
#you can now just run the script
#=====================


#==
#Load fam file
#==

fam<-read.table(paste("./Data/",species,"/",species,dataset,"_beforeQC.fam",sep=""),h=F)
fam<-fam[1:2]#only keep first two columns
colnames(fam)<-c("FID","IID")#rename columns

#===
#Make list of all breeds/populations/FIDs in your dataset (FID's) and amount of animals 
#===

pop<-data.frame(table(fam$FID))
colnames(pop)<-c("FID","N_animals")
pop <- with(pop,  pop[order(pop$FID,decreasing = FALSE) , ])#sort alphabetically

#Check the 'pop' dataframe if all populations (FID) are named correctly and if the number of individuals (N_animals) is correct

#===
#Load map file, determine amount of chromosomes (no sex chromosomes!), and length of genome (in kb) covered by SNPs (sum of distance between first and last SNP per chromosome)
#===

#load bim file: determine number of chromosomes
map<-read.table(paste("./Data/",species,"/",species,dataset,"_beforeQC.bim",sep=""),h=F)
colnames(map)<-c("CHR","SNP_name","POS","BP","A1","A2")    #adapt colnames
map_chr<-max(map$CHR)    #Determine number of chromosomes without sex chromosomes

#Calculate length in bp per chromosome and length of total genome covered (in kb)
genome_length=0
#for loop per chromosome
for (ks in 1:map_chr){
  #map file of chromosome i
  tmp_map<-map[map$CHR==ks,]
  #determine genome length of this chromosome and add to total genome length
  genome_length<-genome_length+((max(tmp_map$BP)-min(tmp_map$BP))/1000) #divide by thousand to get kb
}

#===
#Create directories to write  output
#===

#=
#1. create output directory and directory to write  summary results of all breeds within species
#=

dir.create("./output")
dir.create("./output/all_chr_ROH")

#=
#2. create map 'species' within ./output directory and ./output/all_chr_ROH directory
#=

dir.create(paste("./output/",species,sep=""))
dir.create(paste("./output/all_chr_ROH/",species,sep=""))


#====================================

#####################################
#Run loop per FID: population specific Quality control and ROH analyses using PLINK
#####################################

#====================================

#remove possible summary information with Sroh from previous run
rm(summary_ROH_breeds)

#start loop for each population
for (k in 1:nrow(pop)){
  
  #Print the current population
  print(paste(species,dataset,pop[k,1]),sep=" ")
  
  #At beginning of run: remove files in directory (from possible previous runs)! 
  #=>if error in next run, do not proceed with previous
  unlink(x="./Data/ROH_analyse_population*", recursive = FALSE)
  unlink(x="./Data/DummyGenotypes_ROH*", recursive = FALSE)
  
  #Choose population to use for ROH-analysis:
  #Give this population a standard name and format: ROH_analyse_x
  chosen_pop<-as.data.frame(as.character(pop[k,1]))
  write.table(chosen_pop,paste("./Data/",species,"/Keep_FID.txt",sep=""),row.names = F,col.names = F,quote=FALSE)
  
  #Load bed and bim beforeQC
  system(paste(plink," --chr-set ",map_chr," --bfile ./Data/",species,"/",species,dataset,"_beforeQC"," --keep-fam ./Data/",species,"/Keep_FID.txt --allow-extra-chr --make-bed --out ./Data/beforeQC",sep=""))
  
  
  #If there are less animals before QC than you specified there may be: skip to next population
  fam_pop<-read.table(paste("./Data/beforeQC.fam",sep=""),h=F)
  if(nrow(fam_pop)<minimum_number_animals_population){next}
  
  
  #====================================
  
  #####################################
  #Start of Quality Control
  #####################################
  
  #====================================
  
  #===
  #1. mind
  #===
  system(paste(plink," --chr-set ",map_chr," --bfile ./Data/beforeQC --mind ",mind," --allow-extra-chr --make-bed --out ./Data/beforeQC_mind",sep=""))
  
  #if all animals fail QC: go to next population
  if (file.exists(paste("./Data/beforeQC_mind.bed",sep=""))=="FALSE") {
    next
  }
  
  #===
  #2. Check for duplicates
  #===
  
  ## Calculate IBS of all pairwise combinations 
  system(paste(plink," --chr-set ",map_chr," --bfile ./Data/beforeQC_mind --genome --allow-extra-chr --make-bed --out ./Data/beforeQC",sep=""))
  
  
  #if there are less than 2 animals left: no IBS possible: skip this population!:
  if (file.exists(paste("./Data/beforeQC.genome",sep=""))=="FALSE") {
    next
  }
  
  
  #load IBS-matrix (comparable qua design to relationshipmatrix: nxn)
  genome<-read.table("./Data/beforeQC.genome",h=T)
  #Which animals have relatedness higher than 95%?
  tmp<-genome[as.numeric(genome$PI_HAT)>=0.95,]
  #If there are more than 0 animals with more than 95% relatedness: remove these animals
  if (nrow(tmp)>0){
    #Rbind iid1 and iid2
    test1<-tmp[c("FID1","IID1")]
    colnames(test1)<-c("FID2","IID2")
    test<-rbind(test1,tmp[c("FID2","IID2")])
    #remove these animals
  } else {test<-c()}
  #write  file with animals to remove
  write.table(test,"./Data/failIBD.txt",row.names = F,col.names = F,quote = F)
  
  #remove these animals
  system(paste(plink," --chr-set ",map_chr," --bfile ./Data/beforeQC_mind --remove ./Data/failIBD.txt  --allow-extra-chr --make-bed --out ./Data/beforeQC_IBD",sep=""))
  
  #===
  #3. Maf
  #===
  
  #If you have specified maf to be =0: do not perform maf!
  if(maf==0){
    system(paste(plink," --chr-set ",map_chr," --bfile ./Data/beforeQC_IBD --allow-extra-chr --make-bed --out ./Data/beforeQC_maf",sep=""))
  }
  if(maf!=0){
    system(paste(plink," --chr-set ",map_chr," --bfile ./Data/beforeQC_IBD --maf ",maf," --allow-extra-chr --make-bed --out ./Data/beforeQC_maf",sep=""))
  }
  
  #===
  #4. hwe
  #===
  
  #If you have specified hwe to be =0: do not perform maf!
  if(hwe==0){
    system(paste(plink," --chr-set ",map_chr," --bfile ./Data/beforeQC_maf --allow-extra-chr --make-bed --out ./Data/beforeQC_hwe",sep=""))
  }
  if(maf!=0){
    system(paste(plink," --chr-set ",map_chr," --bfile ./Data/beforeQC_maf --hwe ",hwe," --allow-extra-chr --make-bed --out ./Data/beforeQC_hwe",sep=""))
  }
  
  
  #===
  #5. geno
  #===
  
  system(paste(plink," --chr-set ",map_chr," --bfile ./Data/beforeQC_hwe --geno ",geno," --allow-extra-chr --make-bed --out ./Data/beforeQC_final",sep=""))
  
  #Load QC-controlled datafile and make ped and map
  system(paste(plink," --chr-set ",map_chr," --bfile ./Data/beforeQC_final --allow-extra-chr --make-bed --recode --out ./Data/ROH_analyse_population_x",sep=""))
  
  #if final bed file is not formed: skip this population!:
  if (file.exists(paste("./Data/ROH_analyse_population_x.bed",sep=""))=="FALSE") {
    next
  }
  
  
  #===
  #At the end of QC: remove files in directory! =>if error in next run, do not proceed with previous
  #===
  
  unlink(x="./Data/beforeQC*", recursive = FALSE)
  
  
  
  
  
  #====================================
  
  #####################################
  #Start of ROH analysis
  #####################################
  
  #====================================
  
  
  
  #Choose one of the two options below: do a for loop per chromosome (results per chromosome per FID) or for all chromosomes
  
  #Option 1: Per chromosome  loop => if you want to loop per chromosome and do ROH analysis on a chromosome basis
  #for (chr in c(paste("1-",map_chr,sep=""),1:map_chr)){
  
  #Option 2: use all chromosomes loop
  for (chr in c(paste("1-",map_chr,sep=""))){
    
    #Print species, dataset and population and chromosomes that are under evaluation
    print(paste("-species-",species,"-dataset-",dataset,"-population-",as.character(pop[k,1]),"-chromosome-",chr,sep=" "))
    
    
    #Make directory for outpur of breed/population/FID 1
    population=as.character(pop[k,1])
    
    
    #Create population directory
    dir.create(paste("./output/",species,"/",population,sep=""))
    
    #Store location of directory to write  results
    directory<-paste("./output/",species,"/",population,sep="")
    
    
    
    #======
    #1.2.2 Define dataset
    #======
    
    #===
    #define data frame to fill in during for loop: general statistics populational level
    #===
    
    #Create empty matrix to later fill in with results
    Sroh<-data.frame(matrix(ncol=21,nrow=1))
    #Colnames
    names(Sroh)<-c("Population","run","Het","Nanimals","NROH","Nsnps","KBROH","dummy_length","SROH", "SROH02_05", "SROH05_1", "SROH1_2","SROH2_4","SROH4_8","SROH8_16","SROH16_","SROH1_","Fgrm","Fhom","Funi","l")
    
    #Write  details of number allowed heterozygotes and missing snps in ROH
    Sroh[1,c("run")]<-paste(het,"_het_",mis,"_mis",sep="")
    
    #Write  number of animals in population
    Sroh[1,c("Nanimals")]<-pop[k,2]
    
    #Load specific fam file for chosen population
    fam<-read.table(paste("./Data/ROH_analyse_population_x.fam",sep=""),h=F)
    colnames(fam)<-c("FID","IID")
    
    #If there are less than x animals in fam file, skip population
    if (nrow(fam)<minimum_number_animals_population) {
      next
    }
    
    #===
    #define data frame to fill in during for loop: per animal statistics
    #===
    
    #Make empty data frame to fill in:
    ROH_animal<-data.frame(matrix(ncol=2,nrow=nrow(fam)))
    names(ROH_animal)<-c("FID","IID")
    
    #FID and IID information
    ROH_animal[,1:2]<-fam[,1:2]
    
    
    #======
    #1.2.3. ROH-analyse on whole dataset without pruning:
    #======
    
    Sroh[c(1),c("Population")]<-population
    
    #===
    #1. calculate Mean heterozygosity population using --hardy in plink
    #===
    
    system(paste(plink," --chr-set ",map_chr," --allow-extra-chr --bfile ./Data/ROH_analyse_population_x --chr ",chr," --hardy --out ./Data/ROH_analyse_population_x",sep=""))
    
    #Load results
    plink.hwe <- read.table(file="./Data/ROH_analyse_population_x.hwe",fill=TRUE,header=TRUE)
    #write  observed heterozygosity
    heter<-mean(plink.hwe$O.HET.,na.rm = T) #
    
    #Write  mean heterozygosity value in data frame
    Sroh[c(1),c("Het")]<-mean(plink.hwe$O.HET.,na.rm = T)
    
    #===
    #2. Mean heterozygosity per animal: --het in plink
    #===
    
    system(paste(plink," --chr-set ",map_chr," --allow-extra-chr --bfile ./Data/ROH_analyse_population_x --chr ",chr," --het --out ./Data/ROH_analyse_population_x",sep=""))
    
    #load results
    plink.het <- read.table(file="./Data/ROH_analyse_population_x.het",fill=TRUE,header=TRUE)
    #Calculate observed heterozygosity per animal
    plink.het$O.Het<-(plink.het$N.NM.-plink.het$O.HOM.)/plink.het$N.NM.
    
    #Write  results
    ROH_animal$O_het<-plink.het$O.Het
    
    #===
    # 3 Number of SNPs per window
    #===
    
    #l=log(0.05/(number_of_SNPs*number of animals))/log(1-mean_heterozygosity)
    l <- (log(0.05/(nrow(map)*nrow(fam)))/log(1-heter)) #Lencz et al. 2007 Lencz, T., Lambert, C., DeRosse, P., Burdick, K. E., Morgan, T. V., Kane, J. M., ... & Malhotra, A. K. (2007). Runs of homozygosity reveal highly penetrant recessive loci in schizophrenia. Proceedings of the National Academy of Sciences, 104(50), 19942-19947.
    #or Purfield et al.2012 (Purfield, D. C., Berry, D. P., McParland, S., & Bradley, D. G. (2012). Runs of homozygosity and population history in cattle. BMC genetics, 13(1), 70.)
    #Round l parameter to an integer number
    l<-round(l)
    if (is.na(l)){l<-50} #In a rare occasion l is set to NA, set l = 50 in this case
    
    #Change window-snp to l
    windowsnp<-l ##--homozyg-window-snp
    
    #Write  results
    Sroh[c(1),c("l")]<-l
    
    #Calculate windowthreshold based on formula from Meyermans&Gorssen et al, 2020:
    windowthreshold<-floor(1000*((Nout+1)/l))/1000
    
    
    #===
    #4. ROH analysis in Plink (--homozyg)
    #===
    
    #ROH command in PLINK with all parameters as specified at start of script
    system(paste(plink," --chr-set ",map_chr," --allow-extra-chr --bfile ./Data/ROH_analyse_population_x --chr ",chr," --homozyg group --homozyg-window-het ",het," --homozyg-het ",het," --homozyg-window-missing ",mis," --homozyg-snp ",l," --homozyg-kb ",kb," --homozyg-window-snp ",windowsnp,"  --homozyg-window-threshold ",windowthreshold,"  --homozyg-density ",density," --homozyg-gap ",gap," --out ./Data/ROH_analyse_population_x",sep=""))
    
    #Load ROH data individual animals
    ROH.indiv <- fread(file="./Data/ROH_analyse_population_x.hom.indiv",fill=TRUE,header=TRUE)
    #Load roh data per ROH
    ROH <- fread(file="./Data/ROH_analyse_population_x.hom",fill=TRUE,header=TRUE)
    #Load roh summary per SNP
    populationx.ROH.summary<- fread(file="./Data/ROH_analyse_population_x.hom.summary",fill=TRUE,header=TRUE)
    
    #If next file does not exist: do not write 
    if(file.exists("./Data/ROH_analyse_population_x.hom.overlap")){
      ROH.overlap<- fread(file="./Data/ROH_analyse_population_x.hom.overlap",fill=TRUE,header=TRUE)
      write.table(as.matrix(ROH.overlap),paste(directory,"/ROH.overlap_",population,"_.txt",sep=""),row.names = F,col.names = F,quote=FALSE)
    }
    
    #Write  ROH.indiv and other files
    write.table(as.matrix(ROH),paste(directory,"/ROH_",population,".txt",sep=""),row.names = F,col.names = F,quote=FALSE)
    write.table(as.matrix(populationx.ROH.summary),paste(directory,"/populationx.ROH.summary_",population,"_.txt",sep=""),row.names = F,col.names = F,quote=FALSE)
    
    
    #===
    #5 Calculating the genome coverage length of the analysis
    #===
    
    #Create a dummy animal: one individual with the same map file as in the analysis, but with completely homozygous genotype
    
    #Load map and ped file
    dummy.map.basis <- fread(file="./Data/ROH_analyse_population_x.map")
    dummy.ped.basis <- fread(file="./Data/ROH_analyse_population_x.ped")
    
    #Write  number of snps
    Sroh[1,c("Nsnps")]<-nrow(dummy.map.basis)
    
    #make large file with only A's to construct homozygous dummy animal
    dummy_basis<-t(as.matrix((rep("A",30000000))))
    
    #make dummy pedigree: 6 columns + 2*NR of SNPs (all SNPs A)
    #First six columns are always equal to fam file!
    dummy_ped<-c((as.matrix(fam[1,1:6])),dummy_basis[c(7:ncol(dummy.ped.basis))])
    dummy_ped[1:20]
    dummy.ped.basis<-as.matrix(t(dummy_ped))
    
    #Write  dummy pedigree and map
    write.table(x=as.matrix(dummy.ped.basis),file="./Data/DummyGenotypes_ROH_basis.ped",row.names=FALSE,col.names=FALSE,sep=" ",quote=FALSE,na="-9")
    write.table(x=as.matrix(dummy.map.basis),file="./Data/DummyGenotypes_ROH_basis.map",row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE,na="-9")
    
    #ROH-analysis on dummy pedigree with same parameter settings
    system(paste(plink," --chr-set ",map_chr," --allow-extra-chr --file ./Data/DummyGenotypes_ROH_basis --chr ",chr," --homozyg group --homozyg-window-het ",het," --homozyg-het ",het," --homozyg-window-missing ",mis," --homozyg-snp ",l," --homozyg-kb ",kb," --homozyg-window-snp ",windowsnp,"  --homozyg-window-threshold ",windowthreshold,"  --homozyg-density ",density," --homozyg-gap ",gap,"  --out ./Data/DummyGenotypes_ROH_basis",sep="")) 
    
    #read files from dummy analysis
    dummy.hom_VPF <- fread(file="./Data/DummyGenotypes_ROH_basis.hom.indiv",fill=TRUE,header=TRUE)
    dummy_length<-dummy.hom_VPF$KB[1] #Kb 2.051.210   (median total length (kb): 2.457.910 is the length of the porcine genome (https://www.ncbi.nlm.nih.gov/genome/?term=sus+scrofa))
    
    #Write  coverage of dummy in data frame
    Sroh[1,c("dummy_length")]<-dummy.hom_VPF$KB[1]
    
    #Write  total number of ROHs
    Sroh[1,c("NROH")]<-sum(ROH.indiv$NSEG)
    
    #Write  total length of ROHs
    Sroh[1,c("KBROH")]<-mean(ROH.indiv$KB)
    
    
    
    #===
    #6. Calculate ROH per length class: "SROH","SROH02_05","SROH05_1", "SROH1_2","SROH2_4","SROH4_8","SROH8_16","SROH16_","SROH1_"
    #===
    
    #Store total Sroh as total kb of ROH divided by dummy ROH length (maximal detectable ROH using current settings)
    ROH.indiv$F_ROHall<-ROH.indiv$KB/1
    
    #ROH per individual when only considering runs Smaller than 500 kb
    for (j in 1:nrow(fam)) {
      iid<-as.character(ROH.indiv$IID)[j]
      tmp<-ROH[ROH$IID==iid & ROH$KB>=200 & ROH$KB<500,]
      if (nrow(tmp)>0) {
        ROH.indiv$SROH02_05[ROH.indiv$IID==iid]<-sum(tmp$KB)/1
        ROH.indiv$NSEG02_05[ROH.indiv$IID==iid]<-nrow(tmp)
      }
      else {
        ROH.indiv$SROH02_05[ROH.indiv$IID==iid]<-0
        ROH.indiv$NSEG02_05[ROH.indiv$IID==iid]<-0
      }
    }
    
    #ROH per individual when only considering runs Smaller than 2000 kb
    for (j in 1:nrow(fam)) {
      iid<-as.character(ROH.indiv$IID)[j]
      tmp<-ROH[ROH$IID==iid & ROH$KB>=500 & ROH$KB<1000,]
      if (nrow(tmp)>0) {
        ROH.indiv$SROH05_1[ROH.indiv$IID==iid]<-sum(tmp$KB)/1
        ROH.indiv$NSEG05_1[ROH.indiv$IID==iid]<-nrow(tmp)
      }
      else {
        ROH.indiv$SROH05_1[ROH.indiv$IID==iid]<-0
        ROH.indiv$NSEG05_1[ROH.indiv$IID==iid]<-0
      }
    }
    
    #ROH per individual when only considering runs Smaller than 2000 kb
    for (j in 1:nrow(fam)) {
      iid<-as.character(ROH.indiv$IID)[j]
      tmp<-ROH[ROH$IID==iid & ROH$KB>=1000 & ROH$KB<2000,]
      if (nrow(tmp)>0) {
        ROH.indiv$SROH1_2[ROH.indiv$IID==iid]<-sum(tmp$KB)/1
        ROH.indiv$NSEG1_2[ROH.indiv$IID==iid]<-nrow(tmp)
      }
      else {
        ROH.indiv$SROH1_2[ROH.indiv$IID==iid]<-0
        ROH.indiv$NSEG1_2[ROH.indiv$IID==iid]<-0
      }
    }
    
    #ROH per individual when only considering runs between 2000 and 4000 kb
    for (j in 1:nrow(fam)) {
      iid<-as.character(ROH.indiv$IID)[j]
      tmp<-ROH[ROH$IID==iid & ROH$KB>=2000 & ROH$KB<4000,]
      if (nrow(tmp)>0) {
        ROH.indiv$SROH2_4[ROH.indiv$IID==iid]<-sum(tmp$KB)/1
        ROH.indiv$NSEG2_4[ROH.indiv$IID==iid]<-nrow(tmp)
      }
      else {
        ROH.indiv$SROH2_4[ROH.indiv$IID==iid]<-0
        ROH.indiv$NSEG2_4[ROH.indiv$IID==iid]<-0
      }
    }
    
    #ROH per individual when only considering runs between 4000 and 8000 kb
    for (j in 1:nrow(fam)) {
      iid<-as.character(ROH.indiv$IID)[j]
      tmp<-ROH[ROH$IID==iid & ROH$KB>=4000 & ROH$KB<8000,]
      if (nrow(tmp)>0) {
        ROH.indiv$SROH4_8[ROH.indiv$IID==iid]<-sum(tmp$KB)/1
        ROH.indiv$NSEG4_8[ROH.indiv$IID==iid]<-nrow(tmp)
      }
      else {
        ROH.indiv$SROH4_8[ROH.indiv$IID==iid]<-0
        ROH.indiv$NSEG4_8[ROH.indiv$IID==iid]<-0
        
      }
    }
    
    #ROH per individual when only considering runs between 8 and 16 kb
    for (j in 1:nrow(fam)) {
      iid<-as.character(ROH.indiv$IID)[j]
      tmp<-ROH[ROH$IID==iid & ROH$KB>=8000 & ROH$KB<16000,]
      if (nrow(tmp)>0) {
        ROH.indiv$SROH8_16[ROH.indiv$IID==iid]<-sum(tmp$KB)/1
        ROH.indiv$NSEG8_16[ROH.indiv$IID==iid]<-nrow(tmp)
      }
      else {
        ROH.indiv$SROH8_16[ROH.indiv$IID==iid]<-0
        ROH.indiv$NSEG8_16[ROH.indiv$IID==iid]<-0
      }
    }
    
    #ROH per individual when only considering runs larger than 16000kb
    
    for (j in 1:nrow(fam)) {
      iid<-as.character(ROH.indiv$IID)[j]
      tmp<-ROH[ROH$IID==iid & ROH$KB>=16000,]
      if (nrow(tmp)>0) {
        ROH.indiv$SROH16_[ROH.indiv$IID==iid]<-sum(tmp$KB)/1
        ROH.indiv$NSEG16_[ROH.indiv$IID==iid]<-nrow(tmp)
      }
      else {
        ROH.indiv$SROH16_[ROH.indiv$IID==iid]<-0
        ROH.indiv$NSEG16_[ROH.indiv$IID==iid]<-nrow(tmp)
        
      }
    }
    
    #ROH per individual when only considering runs larger than 1000kb
    
    for (j in 1:nrow(fam)) {
      iid<-as.character(ROH.indiv$IID)[j]
      tmp<-ROH[ROH$IID==iid & ROH$KB>=1000,]
      if (nrow(tmp)>0) {
        ROH.indiv$SROH1_[ROH.indiv$IID==iid]<-sum(tmp$KB)/1
        ROH.indiv$NSEG1_[ROH.indiv$IID==iid]<-nrow(tmp)
      }
      else {
        ROH.indiv$SROH1_[ROH.indiv$IID==iid]<-0
        ROH.indiv$NSEG1_[ROH.indiv$IID==iid]<-nrow(tmp)
        
      }
    }
    
    #Join  in one file with summary statistics per animal
    ROH_animal<-join(ROH_animal,ROH.indiv)
    ROH_animal$PHE<-NULL #Remove column PHE
    
    #Write  SROH
    Sroh[1,c("SROH")]<-mean(ROH.indiv$F_ROHall)
    
    #Write  SROH02_05
    Sroh[1,c("SROH02_05")]<-mean(ROH.indiv$SROH02_05)
    
    #Write  SROH05_1
    Sroh[1,c("SROH05_1")]<-mean(ROH.indiv$SROH05_1)
    
    #Write  SROH1_2
    Sroh[1,c("SROH1_2")]<-mean(ROH.indiv$SROH1_2)
    
    #Write  SROH2_4
    Sroh[1,c("SROH2_4")]<-mean(ROH.indiv$SROH2_4)
    
    #Write  SROH4_8
    Sroh[1,c("SROH4_8")]<-mean(ROH.indiv$SROH4_8)
    
    #Write  SROH8_16
    Sroh[1,c("SROH8_16")]<-mean(ROH.indiv$SROH8_16)
    
    #Write  SROH16
    Sroh[1,c("SROH16_")]<-mean(ROH.indiv$SROH16_)
    
    #Write  SROH5
    Sroh[1,c("SROH1_")]<-mean(ROH.indiv$SROH1_)
    
    #Write  ROH.indiv and other files
    write.table(ROH.indiv,paste(directory,"/ROH.indiv_",population,"_.txt",sep=""),row.names = F,col.names = F,quote=FALSE)
    
    
    #================
    ## Inbreeding as Fhat1: --ibc command in plink (~variance-standardized relationship minus 1, Fhat2 (~Plink) and Fhat3 (~Yang) 
    #================
    
    #Run command for different inbreeding values via PLINK
    system(paste(plink," --chr-set ",map_chr," --allow-extra-chr --bfile ./Data/ROH_analyse_population_x --chr ",chr," --ibc --out ./Data/ROH_analyse_population_x",sep=""))
    
    
    # Fhat1 = Diagonal of GRM (based on variance in additive genetic values): not completely true
    # Fhat2 = Plink inbreeding based on expected heterozygosity (excess heterozygotes)
    # Fhat3 = Correlation between uniting gametes
    
    #Load results
    grm.ibc <- read.table("./Data/ROH_analyse_population_x.ibc",header=TRUE)
    
    #Join these results to summary file
    ROH_animal<-join(ROH_animal,grm.ibc)
    #Remove column nomiss
    ROH_animal$NOMISS<-NULL
    
    #Write  ROH.indiv and other files
    write.table(ROH.indiv,paste(directory,"/grm.ibc_",population,"_.txt",sep=""),row.names = F,col.names = F,quote=FALSE)
    
    #Rearange order of columns
    ROH_animal<-ROH_animal[c("FID","IID","O_het","F_ROHall","SROH02_05", "SROH05_1","SROH1_2","SROH2_4","SROH4_8","SROH8_16","SROH1_",
                             "SROH16_","NSEG", "NSEG02_05", "NSEG05_1", "NSEG1_2","NSEG2_4","NSEG4_8","NSEG8_16","NSEG1_",
                             "NSEG16_","KB","KBAVG","Fhat1","Fhat2","Fhat3")]
    
    #round following columns to 2 digits
    ROH_animal[,c(4:12,22:24)]<-round(ROH_animal[,c(4:12,22:24)]*100,2)
    
    #Write  summary file ROH_animal
    write.xlsx(ROH_animal,file=paste(directory,"/ROH_animal_Results_",population,".xls",sep=""),row.names = F,col.names = T,showNA = F)
    
    #Write  Fgrm
    Sroh[1,c("Fgrm")]<-mean(grm.ibc$Fhat1)
    
    #Write  Fhom
    Sroh[1,c("Fhom")]<-mean(grm.ibc$Fhat2)
    
    #Write  Funi
    Sroh[1,c("Funi")]<-mean(grm.ibc$Fhat3)
    
    #Correlational plots between inbreeding levels: only if there are more than 10 observations and observed ROH
    if(nrow(ROH_animal)>10 & sum(ROH_animal$F_ROHall,na.rm=T)>0.01){
      df<-ROH_animal
      #Remove na columns
      df <- df[,colSums(is.na(df))<nrow(df)]
      png(filename=paste(directory,"/Pairs_plot_inbreeding_methods_",population,"_.png",sep=""),width = 1300,height = 1000)
      pairs(df[c("F_ROHall","SROH16_","SROH1_","Fhat1","Fhat2","Fhat3")], lower.panel=panel.smooth, upper.panel=panel.cor)
      dev.off()
    }
    
    #======
    #1.2.5. Write  results
    #======
    
    
    #store results under new name
    SROH_file<-Sroh
    #calculate mean number of ROH per animal
    SROH_file$Mean_NROH<-SROH_file$NROH/SROH_file$Nanimals
    #calculate total kb of ROH population
    SROH_file$KBROH_total<-SROH_file$KBROH*SROH_file$Nanimals
    #calculate percentage of coverage (dummy length over genome length)
    SROH_file$perc_coverage<-SROH_file$dummy_length/genome_length
    #calculate Sroh using genome length
    SROH_file$SROH_totalgenome<-SROH_file$KBROH/genome_length
    
    
    #Round results
    SROH_file[c(3,9:18,22:23)]<-round(SROH_file[c(3,9:18,22:23)]*100,2)
    
    #Write  summary file
    write.xlsx(SROH_file,file=paste(directory,"/ROH_Results_",population,".xls",sep=""),row.names = F,col.names = T,showNA = F)
    
    #==
    #Make one big summary file for ALL breeds
    #==
    #if this is first round (summary_ROH_breeds does not exist yet)
    if(exists("summary_ROH_breeds")=="FALSE"){
      summary_ROH_breeds<-SROH_file}
    #else: rbind output from every breed
    else{summary_ROH_breeds<-rbind(summary_ROH_breeds,SROH_file)}
    
    #Write  txt file
    write.table(summary_ROH_breeds, file=paste("./output/all_chr_ROH/",species,"/Summary_ROH_per_breed_",species,".txt",sep=""),col.names = TRUE,row.names = FALSE,quote = FALSE,sep="|")
    
    
    #Write  log file with all specifications
    
    log.txt<-paste("species:",species,"\n",
                   "dataset:",dataset,"\n",
                   "population:",population,"\n",
                   "Number of individuals before QC:",nrow(fam_pop),"\n",
                   "Number of SNPs before QC:",nrow(map),"\n","\n",
                   "QC-settings","\n",
                   "Call rate per individual =",1-mind,"\n",
                   "maf =",maf,"\n",
                   "hwe =",hwe,"\n",
                   "geno =",geno,"\n","\n",
                   "ROH-settings","\n",
                   "Number of snps in window =",windowsnp,"\n",
                   "Number of allowed heterozygotes =",het,"\n",
                   "Number of allowed missing =",mis,"\n",
                   "Number of outer SNPs not allowed for threshold calculation =",Nout,"\n",
                   "Threshold to call segments from window output =",windowthreshold,"\n",
                   "Density parameter for segments =",density,"\n",
                   "gap parameter for segments =",gap,"\n",
                   "Minimum Length parameter for segments =",kb,"\n",
                   "Minimum number of snps for segments =",l,"\n","\n",
                   "Number of individuals after QC:",nrow(fam),"\n",
                   "Number of SNPs after QC:",nrow(dummy.map.basis),"\n",
                   
                   sep=" ")
    cat(log.txt)
    
    #To write  file without blank space at the end: capture output!
    capture.output(cat(log.txt), file=paste(directory,"/Log_ROH_settings_",population,"_",dataset,".txt",sep=""))
    
    
    #If you have studied all chromosomes: make manhattan plots of ROH and look at density distribution snps
    if (chr==paste("1-",map_chr,sep="")) {
      
      #====================================
      
      #####################################
      #4. Distribution and density of chromosomes per animal: manhattan plots
      #####################################
      
      #====================================
      
      
      #Load results of every ROH detected
      ROH <- fread(file="./Data/ROH_analyse_population_x.hom",fill=TRUE,header=TRUE)
      populationx.ROH.summary<- fread(file="./Data/ROH_analyse_population_x.hom.summary",fill=TRUE,header=TRUE)
      
      
      #==============
      #4.1. SNPs in ROH per animal (1= appeared in roh; 0 is did not appear in roh)
      #==============
      
      
      #==============
      #4.2. manhattan plot
      #==============
      
      
      #Calculate percentage of animals with snp in ROH:
      populationx.ROH.summary$percentage_ROH<-round(populationx.ROH.summary$UNAFF/nrow(fam)*100,2)
      
      #Write  table with % of animals in roh per SNP
      SNP_results<-populationx.ROH.summary[,c(1:3,5,6)]
      colnames(SNP_results)<-c("CHR","SNP","BP","N_animals_ROH","Perc_animals_ROH")
      SNP_results$population<-population #also write  name of population
      write.table(SNP_results,paste(directory,"/SNP_results_ROH_",population,".txt",sep=""),row.names = F,col.names = T,quote=FALSE)
      write.table(SNP_results,paste("./output/all_chr_ROH/",species,"/SNP_results_ROH_",population,".txt",sep=""),row.names = F,col.names = T,quote=FALSE)
      
      
      
      
      #use jitter function to improve figures
      populationx.ROH.summary$jitter_mean<-jitter(populationx.ROH.summary$UNAFF)
      
      #Calculate percentage of animals with snp in ROH:
      populationx.ROH.summary$percentage_ROH<-populationx.ROH.summary$jitter_mean/nrow(fam)*100
      
      
      #===
      #Manhattan plot expressed as % of animals
      #===
      
      
      #Save manhattan plot as png
      png(filename=paste(directory,"/",population,"_N=",nrow(fam),"Incidence_SNP_ROH_allChr_density",density,".png",sep=""),width = 1300,height = 1000, pointsize=30)
      
      manhattan(populationx.ROH.summary,chr="CHR",bp="BP",snp="SNP",p="percentage_ROH",logp=F,
                cex = 1, cex.axis = 0.8,
                ylim=c(0,100),ylab="Percentage of SNP incidence in ROH",suggestiveline=F,genomewideline = F,
                col=c("#99A8B4FF","darkblue") )
      mytitle = population
      mysubtitle = paste("N=",nrow(fam))
      mtext(side=3, line=2, at=0.00, adj=0, cex=1.5, mytitle)
      mtext(side=3, line=1, at=0.00, adj=0, cex=1, mysubtitle)
      
      
      dev.off()
      
      
      #Write  version to seperate map: allows to scan quickly for patterns!
      png(filename=paste("./output/all_chr_ROH/",species,"/",population,"_N=",nrow(fam),"Incidence_SNP_ROH_allChr_density",density,".png",sep=""),width = 1300,height = 1000, pointsize=30)
      
      manhattan(populationx.ROH.summary,chr="CHR",bp="BP",snp="SNP",p="percentage_ROH",logp=F,
                cex = 1, cex.axis = 0.8,
                ylim=c(0,100),ylab="Percentage of SNP incidence in ROH",suggestiveline=F,genomewideline = F,
                col=c("#99A8B4FF","darkblue") )
      mytitle = population
      mysubtitle = paste("N=",nrow(fam))
      mtext(side=3, line=2, at=0.00, adj=0, cex=1.5, mytitle)
      mtext(side=3, line=1, at=0.00, adj=0, cex=1, mysubtitle)
      dev.off()
      
      #===
      #Manhattan plot expressed as number of animals
      #===                  
      
      png(filename=paste(directory,"/",population,"N=",nrow(fam),"Incidence_Nranimals_SNP_ROH_allChr_density",density,".png",sep=""),width = 1300,height = 1000, pointsize=30)
      
      manhattan(populationx.ROH.summary,chr="CHR",bp="BP",snp="SNP",p="UNAFF",logp=F,
                cex = 1, cex.axis = 0.8,
                ylim=c(0,nrow(fam)),ylab="Number of animals with SNP in ROH",suggestiveline=F,genomewideline = F,
                col=c("#99A8B4FF","darkblue") )
      mytitle = population
      mysubtitle = paste("N=",nrow(fam))
      mtext(side=3, line=2, at=0.00, adj=0, cex=1.5, mytitle)
      mtext(side=3, line=1, at=0.00, adj=0, cex=1, mysubtitle)
      
      dev.off()
      
      
      #===
      #Manhattan plot expressed as percentage of animals per chromosome
      #===    
      
      #in percentages
      pdf(paste(directory,"/",population,"perc_SNP_ROH_perChr_density",density,"_.pdf",sep=""),width = 10)
      
      #for loop per chromosome
      for (i in 1:map_chr){
        
        #get chromosome i
        tmp<-populationx.ROH.summary[populationx.ROH.summary$CHR==i,]
        
        if (nrow (tmp)>0){
          
          manhattan(tmp,chr="CHR",bp="BP",snp="SNP",p="percentage_ROH",logp=F,
                    cex = 1, cex.axis = 0.8,
                    ylim=c(-0.5,100),ylab="Percentage of SNP incidence in ROH",suggestiveline=F,genomewideline = F,
                    col=c("darkblue") )
          mytitle = population
          mysubtitle = paste("N=",nrow(fam))
          mtext(side=3, line=2, at=0.00, adj=0, cex=1.5, mytitle)
          mtext(side=3, line=1, at=0.00, adj=0, cex=1, mysubtitle)
          
        }
        
      }
      
      dev.off()
      
      
      
      #==============
      #4.3. Visualisation of SNP density per chromomsome: detection of regions were it is hard to find ROH! (low density)
      #==============
      
      
      ggplot(map) +
        geom_histogram(aes(x=BP),binwidth = 1e6) +
        facet_wrap(~CHR,ncol=4) +
        ggtitle(paste("Density of the pruned SNPs for each chromosome",population)) +
        xlab("Position in the genome") +
        ylab("SNP density (SNPs/Mb)") +
        theme_bw()
      
      ggsave(paste(directory,"/SNP_density_",population,".png",sep=""),height = 9,width=10)
      
      
      
      
      #======================================================
      #======================================================
      #ROH island definition: based on Gorssen et al. (2019) and Purfield et al. (2017)
      #======================================================
      #======================================================
      
      
      #Rename dataset 
      roh_inc_all_snps<-populationx.ROH.summary
      
      #===
      #  #Determine dataset with only ROH-island snps based on criterium of Gorssen et al; (2019)
      #===
      
      #Only keep top SNP per bin (1mb): this avoids different weighting (1000 snps per bin vs 10 snps per bin)
      
      #Remove certain files
      rm(fill_vector,bins,bin_fill,All_chr)
      rm(list=ls(pattern="CHR_"))
      
      #For loop per chromosome
      for (z in 1:map_chr){
        #Make empty vector to fill in bins of 1Mb per chromsome
        fill_vector<-data.frame(matrix(nrow=0,ncol=9))
        
        #store chromosome number
        bin_fill<-roh_inc_all_snps[roh_inc_all_snps$CHR==z,]
        #How many bins of 1 MB are there?
        bins<-max(bin_fill$BP)%/%(bin)
        #Bin for loop
        for (r in 0:bins){#bin loop
          #Print chromosome and bin
          print(paste("CHR ",z," bin from ",r*0.5," MB to ",r*0.5+0.5," MB",sep=""))
          xmin<-r*bin
          xmax<-r*bin+bin
          #Look at highest percentage of ROH incidence within specific bin
          SNP_highest_incidence_in_bin<-roh_inc_all_snps[roh_inc_all_snps$CHR==z & roh_inc_all_snps$BP>xmin & roh_inc_all_snps$BP<xmax,]
          #if there are less than 2 animals: skip
          if(nrow(SNP_highest_incidence_in_bin)<2){
            next
          }
          
          #Order SNP_highest_incidence_in_bin based on ROH_incidence
          SNP_highest_incidence_in_bin<-SNP_highest_incidence_in_bin[order(SNP_highest_incidence_in_bin$percentage_ROH,decreasing = T),]
          SNP_highest_incidence_in_bin$bin_begin<-xmin/bin
          SNP_highest_incidence_in_bin$bin_end<-xmax/bin
          #Write  SNP with highest ROH-incidence in vector
          fill_vector[r+1,]<-SNP_highest_incidence_in_bin[1,]
        }
        assign(paste("CHR_",z,sep=""),fill_vector)
      }
      
      #Bind all chromosomes together
      All_chr<-CHR_1#Start with chromosome 1
      #for loop to bind other chromosomes
      for (z in 2:map_chr){
        All_chr<-rbind(All_chr,get(paste("CHR_",z,sep="")))
      }
      
      
      #adapt colnames
      colnames(All_chr)<-c(colnames(roh_inc_all_snps),"Bin_begin","Bin_end")
      
      #Remove na rows
      All_chr<-na.omit(All_chr)
      
      #Calculate ROH island threshold
      threshold<-All_chr
      
      
      #Make columns character
      threshold$CHR<-as.character(threshold$CHR)
      threshold$SNP<-as.character(threshold$SNP)
      #add population name
      threshold$population<-population
      
      #===
      #Calculate z-scores (z-score= (observation-mean)/sd)
      #===
      
      threshold$Z_score <- (threshold$percentage_ROH - mean(threshold$percentage_ROH)) / sd(threshold$percentage_ROH)
      
      #Calculate p-value based on z-score table
      threshold$p_value<-pnorm(threshold$Z_score)
      
      #Determine p-value to put thresthold: what is the % of SNPs in ROH above which we call a SNP in an ROH-island
      p<-min(threshold$percentage_ROH[threshold$p_value>=0.999])
      
      #if threshold is below 30%: make 30%
      if(p<30){p<-30}
      #if threshold is above 80%: make 80%
      if(p>80){p<-80}
      
      #Extract SNPs in ROH island: SNPs with % of ROH larger or equal than p
      snps_in_ROH_island<-threshold[threshold$percentage_ROH>=p,]
      
      #Write  summary file (if there are any snps in roh island)
      if(nrow(snps_in_ROH_island)>0){
        write.xlsx(snps_in_ROH_island,file=paste(directory,"/SNPs_in_ROH_island_",population,".xls",sep=""),row.names = F,col.names = T,showNA = F)
      }
      
      #===============
      #Make manhattan plot with threshold visualized
      #===============
      
      #make chromosome back numeric
      threshold$CHR<-as.numeric(threshold$CHR)
      roh_inc_all_snps$CHR<-as.numeric(roh_inc_all_snps$CHR)
      
      
      #Save manhattan plot as png
      png(filename=paste(directory,"/",population,"_N=",nrow(fam),"Incidence_SNP_ROH_allChr_1snp_per_mb_ROH_island_threshold.png",sep=""),width = 1300,height = 1000, pointsize=30)
      
      #Manhattan plot
      manhattan(populationx.ROH.summary,chr="CHR",bp="BP",snp="SNP",p="percentage_ROH",logp=F,
                cex = 1, cex.axis = 0.8,
                ylim=c(0,100),ylab="Percentage of SNP incidence in ROH",suggestiveline=F,genomewideline = p,
                col=c("#99A8B4FF","darkblue") )
      mytitle = population
      mysubtitle = paste("N=",nrow(fam))
      mtext(side=3, line=2, at=0.00, adj=0, cex=1.5, mytitle)
      mtext(side=3, line=1, at=0.00, adj=0, cex=1, mysubtitle)
      
      
      dev.off()
      
      
      
    }#End of ROH_plots  loop
    
  }#End of chromosomes loop
  
}#End of loop populations per dataset



