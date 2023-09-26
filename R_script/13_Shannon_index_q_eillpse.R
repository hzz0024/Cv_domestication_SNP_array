setwd("~/Dropbox/Mac/Documents/HG/Domestication/Manuscript/Figure_539/Figure_S1")

### Install and load libraries/dependencies---------------------------------

#devtools::install_github('royfrancis/pophelper')
#install.packages('viridis')
library(pophelper)
library(viridis)
library(gridExtra)
library(vegan)
library(dplyr)
### Read in data------------------------------------------------------------

k3 <- readQ('./genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_neutral_pruned_5000_structure_output_K3_run_10_f',filetype='auto')
k4 <- readQ('./genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_neutral_pruned_5000_structure_output_K4_run_1_f',filetype='auto')
k5 <- readQ('./genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_neutral_pruned_5000_structure_output_K5_run_1_f',filetype='auto')
k6 <- readQ('./genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_neutral_pruned_5000_structure_output_K6_run_2_f',filetype='auto')
k7 <- readQ('./genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_neutral_pruned_5000_structure_output_K7_run_1_f',filetype='auto')
k8 <- readQ('./genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_neutral_pruned_5000_structure_output_K8_run_1_f',filetype='auto')
k9 <- readQ('./genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_neutral_pruned_5000_structure_output_K9_run_10_f',filetype='auto')
k10 <- readQ('./genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_neutral_pruned_5000_structure_output_K10_run_10_f',filetype='auto')

# https://www.flutterbys.com.au/stats/tut/tut13.2.html for different diversity index
#S <- apply(k3[[1]] >0.05,1,sum)
k3[[1]][k3[[1]] < 0.05]<- 0
k3_shan <- as.data.frame(diversity(k3[[1]], index = "shannon"))
k3_shan$pop <- inds_grps$V2
k3_shan$ind <- inds_grps$V1
colnames (k3_shan) <- c("div", "pop", "ind")
k3_admix <- k3_shan %>% 
  group_by(as.factor(pop)) %>%
  summarize(admix=mean(div))
colnames (k3_admix) <- c("pop", "admix")
k3_admix <- k3_admix[order(k3_admix$pop),]
#k3_admix <- k3_admix[which(k3_admix$pop != "NYH1"),]

k6[[1]][k6[[1]] < 0.05]<- 0
k6_shan <- as.data.frame(diversity(k6[[1]], index = "shannon"))
k6_shan$pop <- inds_grps$V2
k6_shan$ind <- inds_grps$V1
colnames (k6_shan) <- c("div", "pop", "ind")
k6_admix <- k6_shan %>% 
  group_by(as.factor(pop)) %>%
  summarize(admix=mean(div))
colnames (k6_admix) <- c("pop", "admix")
k6_admix <- k6_admix[order(k6_admix$pop),]

k9[[1]][k9[[1]] < 0.05]<- 0
k9_shan <- as.data.frame(diversity(k9[[1]], index = "shannon"))
k9_shan$pop <- inds_grps$V2
k9_shan$ind <- inds_grps$V1
colnames (k9_shan) <- c("div", "pop", "ind")
k9_admix <- k9_shan %>% 
  group_by(as.factor(pop)) %>%
  summarize(admix=mean(div))
colnames (k9_admix) <- c("pop", "admix")
k9_admix <- k9_admix[order(k9_admix$pop),]
#k9_admix <- k9_admix[which(k9_admix$pop != "NYH1"),]

k10[[1]][k10[[1]] < 0.05]<- 0
k10_shan <- as.data.frame(diversity(k10[[1]], index = "shannon"))
k10_shan$pop <- inds_grps$V2
k10_shan$ind <- inds_grps$V1
colnames (k10_shan) <- c("div", "pop", "ind")
k10_admix <- k10_shan %>% 
  group_by(as.factor(pop)) %>%
  summarize(admix=mean(div))
colnames (k10_admix) <- c("pop", "admix")
k10_admix <- k10_admix[order(k10_admix$pop),]
#k10_admix <- k10_admix[which(k10_admix$pop != "NYH1"),]

# Sample metadata
inds_grps <- read.table('sampleID_pop_n539.txt',header=F,stringsAsFactors = F)
# grps <- as.data.frame(inds_grps[,2])
# grps[,1] <- as.character(grps[,1])
# colnames(grps) <- c('subs')

# Add sample IDs as rownames to qlist

#######################
##### formal run ######
#######################
# discrete color
library(RColorBrewer)
# gradient color
library(viridis)
library(colorspace)
library(ggplot2)
library(ggforce)
library(dplyr) 
library(tibble)
library(patchwork)
library(export)
library(gdsfmt)
library(SNPRelate) # if there is something wrong with gfortran see link here https://thecoatlessprofessor.com/programming/cpp/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/

# for subsets n=509
# vcf
setwd("~/Dropbox/Mac/Documents/HG/Domestication/01_pca")
# 509
vcf.fn <- "genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_neutral_pruned.recode.vcf"
# VCF => GDS
snpgdsVCF2GDS(vcf.fn, "genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_neutral_pruned.recode.gds", method="biallelic.only")
# summary
snpgdsSummary("genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_neutral_pruned.recode.gds")
# Open the GDS file
genofile <- snpgdsOpen("genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe_neutral_pruned.recode.gds")

pca <- snpgdsPCA(genofile,autosome.only=FALSE)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],    # the second eigenvector
                  EV4 = pca$eigenvect[,4],    # the second eigenvector
                  stringsAsFactors = FALSE)
print(tab)
tab_pop = read.delim("sample_eigen_539_edit.txt", header = TRUE, sep='\t')

tab_pop$population_nonum = stringr::str_remove(tab_pop$pop, "[0-9]+")
# count how many individuals in each population
tab_pop %>% count(pop) 

n <- 18
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#col_divergent = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#                    MEW1         MEW2       LIW1      LIW2        DBW1       DBW2      NCW1        NCW2      DBX1     
col_gradient <- c( "#0A2C86", "#849cc1",  "#1D92BD", "#8ad5d9", "#93c47d", "#bedbb1", "#a9a9a9", "#dddddd","#f9476b", 
                   #   DBX2      DBX3       UNC1        UNC2       UMFS       NYH1       NEH1       NEH2       MEH2
                   "#fb90a6","#fddae1", "#cf7fbc",  "#e2b2d6", "#fec155", "#C05805" , "#e1bb94", "#fbd0a5", "#b58383")

#col_gradient = viridis_pal(option = "C")(17)  # n = number of colors seeked
col_gradient = rep(c("#049DD9", "#F25C05"), c(9, 9))
# PC1-2 for individual populations
order1 = c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "DBX1", "DBX2", "DBX3",  "UNC1", "UNC2", "UMFS", "NYH1", "NEH1", "NEH2", "MEH2")
tab_pop$pop <-factor(tab_pop$pop, levels=order1)
# for 2 clusters
order2 = c("Wild", "Selected")
tab_pop$pop1 <-factor(tab_pop$pop1, levels=order2)

p1 <- ggplot(tab_pop, aes(x = EV1, y = EV2)) + 
  #geom_mark_ellipse(aes(fill = pop, label = pop, color=pop), show.legend = FALSE, alpha = 0.2, label.fontface = c("plain"), con.cap=0, label.buffer = unit(0, 'mm'))+
  geom_point(size=3, aes(color=pop, shape=pop1), alpha = 0.6)+
  scale_color_manual(values=col_gradient , name="Population/Line") +
  scale_shape_manual(values=c(15, 17), name="Origin") +
  guides(fill = guide_legend(override.aes=list(shape=17)))+
  #scale_color_manual(values = c("#F25C05", "#049DD9"))+
  #scale_shape_manual(values=c(rep(c(0:2, 5, 15:18), 2), 6), breaks=order1)+
  #scale_shape_manual(values=c(rep(15, 9), rep(17,8)))+
  scale_x_continuous(paste("PC1 (",round(pc.percent[1],3),"%", ")",sep="")) + 
  scale_y_continuous(paste("PC2 (",round(pc.percent[2],3),"%",")",sep=""))+
  theme(legend.position="right",
        legend.title = element_text(size = 12),
        legend.text=element_text(size=12)) + 
  theme(text=element_text(family="Times New Roman", size=12, colour="black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  theme(panel.background = element_rect(fill = 'white', colour = "Black"))+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))

p2 <- p1 + 
  stat_ellipse(aes(fill = pop), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = F) +
  #scale_fill_manual(values = c("#F25C05", "#049DD9"))+
  scale_fill_manual(values = col_gradient) +
  guides(colour = guide_legend(order = 2), 
         shape = guide_legend(override.aes = list(shape = c(0,2)),
                              order = 1))
p2

# Calculate the ellipse size for each population 
# Plot object
# p = ggplot (data, aes (x = x, y = y))+
#   geom_point()+
#   stat_ellipse(segments=201) # Default is 51. We use a finer grid for more accurate area.
# https://stackoverflow.com/questions/38782051/how-to-calculate-the-area-of-ellipse-drawn-by-ggplot2

pb = ggplot_build(p2)

size=c()
for (i in seq(1,18)) {
  # Get ellipse coordinates from plot for each group
  el = pb$data[[2]][which(pb$data[[2]]$group == i),][c("x","y")]
  # Center of ellipse
  ctr = MASS::cov.trob(el)$center  # Per @Roland's comment
  
  # Calculate distance to center from each point on the ellipse
  dist2center <- sqrt(rowSums((t(t(el)-ctr))^2))
  
  # Calculate area of ellipse from semi-major and semi-minor axes. 
  # These are, respectively, the largest and smallest values of dist2center. 
  size_ = pi*min(dist2center)*max(dist2center)
  size = c(size, size_)
}

ellipse_size = data.frame(size, levels(tab_pop$pop))
colnames(ellipse_size) = c("size", "pop")
ellipse_size <- ellipse_size[order(ellipse_size$pop),]

cor.test(ellipse_size$size, k3_admix$admix)
cor.test(ellipse_size$size, k6_admix$admix)
cor.test(ellipse_size$size, k9_admix$admix)
cor.test(ellipse_size$size, k10_admix$admix)

dt_k3 <- cbind(ellipse_size, k3_admix, "k_3") 
colnames (dt_k3) = c("size", "pop", "pop1", "admix", "K_value")
dt_k6 <- cbind(ellipse_size, k6_admix, "k_6") 
colnames (dt_k6) = c("size", "pop", "pop1", "admix", "K_value")
dt_k9 <- cbind(ellipse_size, k9_admix, "k_9") 
colnames (dt_k9) = c("size", "pop", "pop1", "admix", "K_value")
dt_k10 <- cbind(ellipse_size, k10_admix, "k_10") 
colnames (dt_k10) = c("size", "pop", "pop1", "admix", "K_value")

dt_all <- rbind(dt_k3,dt_k6, dt_k9,dt_k10)
dt_all$K_value <- factor(dt_all$K_value, levels = c("k_3", "k_6", "k_9", "k_10"))

mytable <- dt_all %>% 
  group_by(as.factor(K_value)) %>%
  summarize(cor=cor(admix, size),
            p = cor.test(admix, size)$p.value)%>%
  mutate_if(is.numeric, ~sprintf("%.3f",.))

mytheme = list(
  theme_classic()+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),legend.position="", strip.text = element_text(face="plain", size=12),
          axis.text=element_text(face="plain",size=12),axis.title = element_text(face="plain",size=14),plot.title = element_text(face = "bold", hjust = 0.5,size=12))
)

ggplot (data=dt_all, aes(x= size, y = admix, color=as.factor(K_value)))+
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values=c("#83D1C4", "#86AC41", "#F17950", "#78517C"))+
  geom_smooth(method=lm , level=0.95, fill="lightgrey") +
  theme_classic()+
  #annotation_custom(tableGrob(mytable), xmin=0, xmax=0.005, ymin=1, ymax=1.5)+
  #facet_grid(.~K_value)+
  ylab("Mean admixture level")+
  xlab("Ellipse size in PCA")+
  labs(color='K value') +
  theme(strip.text = element_text(face="plain", size=12),
        legend.title=element_text(size=14), legend.text=element_text(size=14),
        axis.text=element_text(face="plain",size=12),axis.title = element_text(face="plain",size=14),plot.title = element_text(face = "bold", hjust = 0.5,size=12),
        panel.border = element_rect(fill = NA, color = "black"),
        text=element_text(face="bold", size=8, colour="black"),
        axis.text.x = element_text(angle = 75, vjust = 1, hjust=1))+
  facet_grid(~K_value,scales = 'free_x', space = 'free_x')

setwd("~/Dropbox/Mac/Documents/HG/Domestication/Manuscript/Figure_539/Figure_S4")
graph2ppt(file="Figure_S4_q005.pptx", width=10, height=6)
# ANOVA

anova <- aov(size ~ admix*as.factor(K_value), data=dt_all)
summary(anova)


###############
# relatedness #
###############
setwd("~/Dropbox/Mac/Documents/HG/Domestication/Manuscript/Figure_539/Figure_S4")
relat_df = read.delim("Relatedness.txt", header = TRUE, sep='\t')
cor.test(ellipse_size$size, relat_df$Relatedness)
relat_size <- data.frame(relat_df$Pop, ellipse_size$size, relat_df$Relatedness)
colnames(relat_size) = c("pop", "size", "Relatedness")

mytable <- relat_size %>% 
  summarize(cor=cor(Relatedness, size),
            p = cor.test(Relatedness, size)$p.value)%>%
  mutate_if(is.numeric, ~sprintf("%.3f",.))

ggplot (data=relat_size, aes(x= size, y = Relatedness, color=pop))+
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method=lm , level=0.95, fill="lightgrey") +
  theme_classic()
graph2ppt(file="Figure_S4_related.pptx", width=6, height=6)

#########
# ANOVA #
#########

df = data.frame(relat_df$Pop, ellipse_size$size, relat_df$Relatedness, dt_all$admix[dt_all$K_value == "k_6"])
colnames(df) = c("pop", "size", "relatedness", "admix")
anova <- aov(size ~ admix*relatedness, data=df)
summary(anova)
