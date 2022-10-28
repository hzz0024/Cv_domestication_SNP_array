# Patterns of ROH analysis
# This script contains the analyses for the first part of the paper,
# exploring ROH variation across the genome, ROH islands and deserts,
# and the association between ROH prevalence and recombination rate variation.

library(data.table)
library(tidyverse)
source("theme_simple.R")
library(windowscanr)
library(cowplot)
library(gt)
library(grid)
library(ggplotify)
library(patchwork)
library(viridis)
library(scales)
library(ggridges)
library(ggpubfigs)
library(gt)
library(lme4)
library(ggpubfigs)
library(broom)
library(dplyr)
library(export)
library(zoo)
options(scipen=9999)

#################
# plot for FROH #
#################
setwd("~/Dropbox/Mac/Documents/HG/Domestication/14_ROH/All_n_539_formal_plot")
library(cowplot)
library(reshape2)

df <- read.delim("Summary_ROH_per_breed_oyster_FROH.csv", header = TRUE, sep=',')
df$Source <- factor(df$Source, levels=c("Wild", "Selected"))
#                    MEW1         MEW2       LIW1      LIW2        DBW1       DBW2      NCW1        NCW2
col_gradient <- c( "#0A2C86", "#849cc1",  "#1D92BD", "#8ad5d9", "#93c47d", "#bedbb1", "#a9a9a9", "#dddddd", 
                   #  DBX1       DBX2      DBX3       UNC1        UNC2       UMFS    NYH1     NEH1       NEH2       MEH2
                   "#f9476b", "#fb90a6","#fddae1", "#cf7fbc",  "#e2b2d6", "#fec155", "#C05805" , "#e1bb94", "#fbd0a5", "#b58383")
order1 = c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "DBX1", "DBX2", "DBX3",  "UNC1", "UNC2", "UMFS", "NYH1", "NEH1", "NEH2", "MEH2")
df$FID <-factor(df$FID, levels=order1)
df$FROH1_=df$FROH1_/100

df %>%  
  dplyr::group_by(FID) %>% 
  dplyr::summarise(Mean = mean(FROH1_)) -> mean_FROH

mean(mean_FROH$Mean[1:9])
mean(mean_FROH$Mean[10:18])

FROH <- df %>%
  ggplot(aes(FID, FROH1_, fill=FID)) +
  #geom_violin(trim = FALSE, alpha=0.4, show.legend = FALSE, width=6) +
  geom_boxplot(width = 0.4, alpha=0.8, aes(fill=FID), outlier.shape = NA)  +
  geom_jitter(alpha=0.2, width = 0.1)+
  scale_fill_manual(values = col_gradient) +
  scale_color_manual(values = col_gradient) +
  #scale_fill_manual(values = c("#049DD9", "#F25C05"))+
  xlab(NULL) + ylab(expression("Genomic"~"Inbreeding"~"Coefficient"~"("~italic(F)[ROH]~")"))+
  theme_classic()+
  theme(#=strip.placement = 'outside',
    panel.spacing = unit(0, "lines"),
    legend.position="none",
    axis.title.y = element_text(margin=margin(r=-50))) + # control the label position y is r 
  theme(text=element_text(family="Times New Roman", size=12, colour="black"),
        axis.text.x = element_text(angle = 75, vjust = 1, hjust=1))
 
FROH <- FROH + ggtitle("(a)")
FROH

kruskal.test(FROH1_ ~ Source, data = df)

#wilcox.test(df$FROH1_[df$Source == "Selected"], df$FROH1_[df$Source == "Wild"])

# DT = read.delim("ROH_analyse_population_x.hom.txt", header = TRUE, sep='\t')
# plotdat <- data.frame(ROH=DT$NSEG, Population=DT$FID)
# 
# wilcox_dt = read.delim("ROH_analyse_population_x.hom.wildsel.txt", header = TRUE, sep='\t')
# wilcox.test(wilcox_dt$NSEG[wilcox_dt$FID == "WILD"],wilcox_dt$NSEG[wilcox_dt$FID == "SEL"])

#########################
# plot for SROH vs NROH #
#########################
library(ggrepel)
setwd("~/Dropbox/Mac/Documents/HG/Domestication/14_ROH/All_n_539_formal_plot/")

ROH_plot <- read_delim("TableS3_SROH_NROH_summary.txt", delim = "\t")

order1 = c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "DBX1", "DBX2", "DBX3",  "UNC1", "UNC2", "UMFS", "NYH1", "NEH1", "NEH2", "MEH2")
ROH_plot$Pop <-factor(ROH_plot$Pop, levels=order1)
col <- c( "#F62A00", "#325A98")

SROH_NSEG <- ggplot(ROH_plot, aes(x = SROH1_, y = NSEG1_, color = Origin, shape=Origin)) +
  geom_point(size=2.5)+
  scale_color_manual(values = col) +
  geom_label_repel(aes(label = Pop),size = 2.4, fill = "white", box.padding = 0.7, max.overlaps = 20, show.legend = FALSE)+
  xlab("Sum total length of ROH (Mb)") + ylab("Total number of ROH (NROH)") + 
  theme_classic()+
  theme(legend.title =element_text(size = 12),
        legend.text=element_text(size=12), 
        legend.position = c("top"))+
  theme(legend.position=c(.2,.85), 
        text=element_text(family="Times New Roman", size=12, colour="black"))

SROH_NSEG <- SROH_NSEG + ggtitle("(b)")
SROH_NSEG

#############################
# plot for ROH distribution #
#############################
#setwd("~/Dropbox/Mac/Documents/HG/Github/References/ROH/sheep_ID-master")
setwd("~/Dropbox/Mac/Documents/HG/Domestication/14_ROH/All_n_539_formal_plot/")
# Chr lengths
chr_data <- read_delim("./chromosome_info_Cv30.txt", delim = "\t") %>% 
        dplyr::rename(size_BP = Length, CHR = Part) %>% 
        mutate(size_KB = size_BP / 1000)

autosomal_genome_size <- chr_data %>% 
        .[2:11, ] %>% 
        summarise(sum_KB = sum(size_KB)) %>% 
        as.numeric()

#~~ ROH for survival data subset
file_path <- "./ROH/All_ROH_distribution.txt"
roh_lengths <- fread(file_path)

# max ROH length
roh_lengths[which.max(roh_lengths$KB), ]

# ROH overview -----------------------------------------------------------------
# descriptive ROH statistics
num_roh_per_ind <- roh_lengths %>% group_by(IID) %>% tally() 
summary(num_roh_per_ind$n)
sd(num_roh_per_ind$n)

# inbreeding coefficients
froh <- roh_lengths %>%
        dplyr::group_by(IID) %>%
        dplyr::summarise(KBAVG = mean(KB), KBSUM = sum(KB)) %>%
        mutate(FROH = KBSUM/autosomal_genome_size) %>% 
        mutate(FROH_cent = FROH - mean(FROH))
mean(froh$FROH)
range(froh$FROH)

# longest ROH
roh_lengths %>% arrange(-KB)

# longest ROH proportional to chr size
chr_sizes <- chr_data %>% .[-1, ] %>% 
        mutate(CHR = str_replace(CHR, "Chromosome ", "")) %>% 
        mutate(CHR = as.integer(CHR))

roh_lengths %>% 
        arrange(desc(KB)) %>% 
        left_join(chr_sizes, by = "CHR") %>% 
        mutate(prop_chr = KB / size_KB) %>% 
        arrange(desc(prop_chr))
#plot(froh$KBSUM, num_roh_per_ind$n)

# ROH length and abundance in the 1% least and most inbred individuals 
num_roh_per_ind %>% 
        left_join(froh) %>% 
       # top_frac(-0.01, FROH) %>%    # top 1% least inbred individuals
        top_n(-7, FROH) %>%      # top 1% most inbred individuals
        summarise(mean(n), mean(KBAVG))

# Supplementary Figure FROH / ROH across individuals ---------------------------
p_froh <- ggplot(froh, aes(FROH)) +
        geom_histogram(bins = 100,  color = "white", size = 0.1, position = "identity",
                       alpha = 1, fill = "#4c566a") +
        ylab("individuals") +
        xlab(expression(inbreeding~coefficient~F[ROH])) +
        scale_y_continuous(expand = c(0, 0)) 
        #theme_simple(grid_lines = FALSE, axis_lines = TRUE, base_size = 12) 
p_froh

p_roh <- ggplot(num_roh_per_ind, aes(n)) +
        #geom_histogram(binwidth = 1,  fill = "#E5E9F0", color = "black",  size = 0.1) +
        geom_histogram(binwidth = 1, color = "white", size = 0.1, position = "identity",
                       alpha = 1, fill = "#4c566a") +
        ylab("individuals") +
        xlab("ROH per genome") +
        scale_y_continuous(expand = c(0, 0)) 
        theme_simple(grid_lines = FALSE, axis_lines = TRUE, base_size = 12, 
                     base_family = "Lato") 
p_roh

p_roh_dist <- p_froh + p_roh + plot_annotation(tag_levels = 'a') &
                theme(plot.tag = element_text(face = "bold"))
p_roh_dist
ggsave("figs/Sup_ROH_dist.jpg", p_roh_dist, width = 7, height = 2.5)

#############################
# plot for ROH heatmap      #
#############################
setwd("~/Dropbox/Mac/Documents/HG/Domestication/14_ROH/All_n_539_formal_plot/")
all_roh <- roh_lengths %>% 
         dplyr::group_by(IID) %>% 
         dplyr::summarise(sum_roh = sum(KB)) %>% 
         ungroup() %>% 
         arrange(desc(sum_roh))

# edit sample_id_to_modify.txt file to have a vector following your target list
write_delim(all_roh, path = "./sample_id_to_modify.txt")
# after editing, now reload the sample list for next step
sample_list <- read.delim("sample_id_to_modify_sorted.csv", header = TRUE, sep=',')

# for 2 clusters

df <- roh_lengths %>%
        mutate(POS1 = POS1 / 1e+6,
               POS2 = POS2 / 1e+6,
               MB = KB / 1000)

df <- df %>% 
        mutate(IID = factor(IID, levels = sample_list$IID))

yax <- data.frame(IID = fct_inorder(levels(df$IID))) %>%
        mutate(yax = seq(from = 2,
                         to = 2*length(unique(df$IID)),
                         by = 2)) 
num_ind <- dim(yax)[1]

df <- left_join(df, yax, by = "IID")

order1 = c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "DBX1", "DBX2", "DBX3",  "UNC1", "UNC2", "UMFS", "NYH1", "NEH1", "NEH2", "MEH2")
df$FID <-factor(df$FID, levels=order1)

shade <- df %>%
        dplyr::group_by(CHR) %>%
        dplyr::summarise(min = min(POS1), max = max(POS2)) %>%
        mutate(min = case_when(CHR == 2 | CHR == 4 | CHR == 6 | CHR == 8 | CHR == 10  ~ 0,
                               TRUE ~ min)) %>%
        mutate(max = case_when(CHR == 2 | CHR == 4 | CHR == 6 | CHR == 8 | CHR == 10  ~ 0,
                               TRUE ~ max))

#col <- c("#2b2c2e", "#FF0000")
#col <- c("#1E3231", "#9aadbf")
col <- c( "#f67321", "#000000")
chr_names <- as.character(1:10)
names(chr_names) <- as.character(1:10)
#chr_names[c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)] <- ""

# df is the freqment plot without labels
df %>% 
        filter(MB > 0) %>% 
        filter(CHR %in% 1:10) %>% 
        ggplot() +
        geom_rect(data=shade, aes(xmin=min, xmax=max, ymin=0, ymax=num_ind*2 + 1), 
                  alpha=0.5, fill = "#eceff4") + # "#f7f7f7" "#eceff4"
        geom_hline(data = yax, aes(yintercept = yax), color = "#ffffff", size = 0.4) + 
        #geom_hline(data = yax, aes(yintercept = yax), color = "black", size = 0.4) +
        geom_rect(aes(xmin = POS1, xmax = POS2, ymin = yax - 1.6, ymax = yax + 1.8, 
                      fill = as.factor(CHR)),  col = "grey", size = 0, alpha = 1, inherit.aes=FALSE) + 
        scale_fill_manual(values = rep(col, 10)) + 
        scale_color_manual(values = rep(col, 10)) +
        scale_y_reverse(expand = c(0, 0)) +
        theme_simple(axis_lines = TRUE, grid_lines = FALSE, base_size = 13, base_family = "Helvetica") +
        facet_grid(~CHR,scales = 'free_x', space = 'free_x', switch = 'x',
                   labeller = as_labeller(chr_names)) +
        theme(#=strip.placement = 'outside',
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_blank(),
                panel.spacing = unit(0, "lines"),
                plot.margin = margin(r = 0.5, l = 0.1, b = 0.1, t = 0.1, unit = "cm"),
                axis.line.x = element_blank(),
                legend.position="right",
                axis.title.x = element_text(margin=margin(t=-25)), # control the label position x is t
                axis.title.y = element_text(margin=margin(r=-30)), # control the label position y is r
                axis.text.y = element_text(colour = "white"),
                axis.line.y = element_blank()) +
        coord_cartesian(clip = 'off') +
        xlab("Chromosome") +
        ylab("Individuals") -> ROH_per_ind
ROH_per_ind
# start to add the popultion labels
# first is the number of samples per populations
FID_cnt <- sample_list %>%
  dplyr::group_by(FID) %>% 
  dplyr::summarise(n=n())
target <- c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "DBX1", "DBX2", "DBX3",  "UNC1", "UNC2", "UMFS", "NYH1", "NEH1", "NEH2", "MEH2")
FID_cnt <- FID_cnt[match(target, FID_cnt$FID),]
FID_cnt

gp.label <- data.frame(FID_cnt$FID, FID_cnt$n, c(rep("Wild", 9), rep("Selected", 9)))
colnames(gp.label) <- c("POP", "N", "Type")
gp.label$Y <- rollmean(c(0, cumsum(gp.label$N*2)),2)
gp.label$X <- x_axis <- min(ggplot_build(ROH_per_ind)$data[[1]]$xmin)-10
gp.label$POP <- factor(gp.label$POP, levels=c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "DBX1", "DBX2", "DBX3",  "UNC1", "UNC2", "UMFS", "NYH1", "NEH1", "NEH2", "MEH2"))
gp.label$CHR <- 1 # this step is important to ensure the points are added to CHR 1 panel!
#                    MEW1         MEW2       LIW1      LIW2        DBW1       DBW2      NCW1        NCW2
cbPalette <- c( "#0A2C86", "#849cc1",  "#1D92BD", "#8ad5d9", "#93c47d", "#bedbb1", "#a9a9a9", "#dddddd", 
                   #  DBX1       DBX2      DBX3       UNC1        UNC2       UMFS       NYH1     NEH1       NEH2       MEH2
                   "#f9476b", "#fb90a6","#fddae1", "#cf7fbc",  "#e2b2d6", "#fec155", "#C05805" , "#e1bb94", "#fbd0a5", "#b58383")

gp.line <- data.frame(seq(1:length(c(0, cumsum(gp.label$N*2)))) , c(0, cumsum(gp.label$N*2)))
colnames(gp.line) <- c("Number", "Pos")

DROH <-  ROH_per_ind +
  geom_point(data = gp.label,aes(x=X,y=Y,color=POP,shape=Type),size=4, alpha=0.95) +
  geom_hline(data = gp.line, aes(yintercept = Pos), color = "grey", size=0.8, alpha=0.5) + 
  scale_x_continuous(expand = expand_scale(mult = c(0.1, 0)))+
  #scale_x_continuous(limits = c(min(gp.label$X) - ddelta/2 , max(ggplot_build(outlier.gp$genotypes)$data[[1]]$x))) + 
  scale_shape_manual(values=c(15, 17), name="Origin") + 
  scale_color_manual(values=cbPalette, name="Population/Line") +
  theme(legend.position="right") +
  theme(legend.title = element_text(size = 10),
        legend.text=element_text(size=10)) + 
  guides(fill = "none")+ 
  # guides(color="none")+
  # guides(shape = "none") +
  theme(text=element_text(family="Times New Roman", size=12, colour="black"))

DROH <- DROH+ ggtitle("(c)")
DROH
layout <-" 
AAAAAAAAAAAAAA
AAAAAAAAAAAAAA
AAAAAAAAAAAAAA
BBBBBCCCCCCCCC
BBBBBCCCCCCCCC
BBBBBCCCCCCCCC
BBBBBCCCCCCCCC
"

final_ROH <- FROH + SROH_NSEG + DROH + plot_layout(design=layout, guides="keep")
final_ROH
ggsave("figs/final_ROH.jpg", final_ROH, width = 12, height = 10, dpi = 300)
graph2ppt(file="figs/final_ROH", width=12, height=10) 



