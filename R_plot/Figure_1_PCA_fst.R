
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
# see https://github.com/clairemerot/Intro2PopGenomics/blob/master/3.2.3/PCA/script_PCA_from_vcf.R for original code
library(gdsfmt)
library(SNPRelate) # if there is something wrong with gfortran see link here https://thecoatlessprofessor.com/programming/cpp/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/

# this is the function for ellipse
geom_enterotype <- function(mapping = NULL, data = NULL, stat = "identity",  position = "identity", 
                            alpha = 0.3, prop = 0.5, ..., lineend = "butt", linejoin = "round", 
                            linemitre = 1, arrow = NULL, na.rm = FALSE, parse = FALSE, 
                            nudge_x = 0, nudge_y = 0, label.padding = unit(0.15, "lines"), 
                            label.r = unit(0.15, "lines"), label.size = 0.1, 
                            show.legend = TRUE, inherit.aes = TRUE) {
  library(ggplot2)
  # create new stat and geom for PCA scatterplot with ellipses
  StatEllipse <- ggproto("StatEllipse", Stat, 
                         required_aes = c("x", "y"), 
                         compute_group = function(., data, scales, level = 0.75, segments = 51, ...) {
                           library(MASS)
                           dfn <- 2
                           dfd <- length(data$x) - 1
                           if (dfd < 3) {
                             ellipse <- rbind(c(NA, NA))
                           } else {
                             v <- cov.trob(cbind(data$x, data$y))
                             shape <- v$cov
                             center <- v$center
                             radius <- sqrt(dfn * qf(level, dfn, dfd))
                             angles <- (0:segments) * 2 * pi/segments
                             unit.circle <- cbind(cos(angles), sin(angles))
                             ellipse <- t(center + radius * t(unit.circle %*% chol(shape)))
                           }
                           ellipse <- as.data.frame(ellipse)
                           colnames(ellipse) <- c("x", "y")
                           return(ellipse)
                         })
  
  # write new ggproto 
  GeomEllipse <- ggproto("GeomEllipse", Geom, 
                         draw_group = function(data, panel_scales, coord) {
                           n <- nrow(data)
                           if (n == 1) 
                             return(zeroGrob())
                           munched <- coord_munch(coord, data, panel_scales)
                           munched <- munched[order(munched$group), ]
                           first_idx <- !duplicated(munched$group)
                           first_rows <- munched[first_idx, ]
                           grid::pathGrob(munched$x, munched$y, default.units = "native", 
                                          id = munched$group, 
                                          gp = grid::gpar(col = first_rows$colour, 
                                                          fill = alpha(first_rows$fill, first_rows$alpha), lwd = first_rows$size * .pt, lty = first_rows$linetype))
                         }, 
                         default_aes = aes(colour = "NA", fill = "grey20", size = 0.5, linetype = 1, alpha = NA, prop = 0.5), 
                         handle_na = function(data, params) {
                           data
                         }, 
                         required_aes = c("x", "y"), 
                         draw_key = draw_key_path
  )
  
  # create a new stat for PCA scatterplot with lines which totally directs to the center
  StatConline <- ggproto("StatConline", Stat, 
                         compute_group = function(data, scales) {
                           library(miscTools)
                           library(MASS)
                           df <- data.frame(data$x,data$y)
                           mat <- as.matrix(df)
                           center <- cov.trob(df)$center
                           names(center)<- NULL 
                           mat_insert <- insertRow(mat, 2, center )
                           for(i in 1:nrow(mat)) {
                             mat_insert <- insertRow( mat_insert, 2*i, center )
                             next
                           }
                           mat_insert <- mat_insert[-c(2:3),]
                           rownames(mat_insert) <- NULL
                           mat_insert <- as.data.frame(mat_insert,center)
                           colnames(mat_insert) =c("x","y")
                           return(mat_insert)
                         },
                         required_aes = c("x", "y")
                         
  )
  
  # create a new stat for PCA scatterplot with center labels
  StatLabel <- ggproto("StatLabel" ,Stat,
                       compute_group = function(data, scales) {
                         library(MASS)
                         df <- data.frame(data$x,data$y)
                         center <- cov.trob(df)$center
                         names(center)<- NULL 
                         center <- t(as.data.frame(center))
                         center <- as.data.frame(cbind(center))
                         colnames(center) <- c("x","y")
                         rownames(center) <- NULL
                         return(center)
                       },
                       required_aes = c("x", "y")
  )
  
  
  layer1 <- layer(data = data, mapping = mapping, stat = stat, geom = GeomPoint, 
                  position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
                  params = list(na.rm = na.rm, ...))
  layer2 <- layer(stat = StatEllipse, data = data, mapping = mapping, geom = GeomEllipse, position = position, show.legend = FALSE, 
                  inherit.aes = inherit.aes, params = list(na.rm = na.rm, prop = prop, alpha = alpha, ...))
  layer3 <- layer(data = data, mapping = mapping, stat =  StatConline, geom = GeomPath, 
                  position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
                  params = list(lineend = lineend, linejoin = linejoin, 
                                linemitre = linemitre, arrow = arrow, na.rm = na.rm, ...))
  if (!missing(nudge_x) || !missing(nudge_y)) {
    if (!missing(position)) {
      stop("Specify either `position` or `nudge_x`/`nudge_y`", 
           call. = FALSE)
    }
    position <- position_nudge(nudge_x, nudge_y)
  }
  layer4 <- layer(data = data, mapping = mapping, stat = StatLabel, geom = GeomLabel, 
                  position = position, show.legend = FALSE, inherit.aes = inherit.aes, 
                  params = list(parse = parse, label.padding = label.padding, 
                                label.r = label.r, label.size = label.size, na.rm = na.rm, ...))
  return(list(layer1,layer2,layer3,layer4))
}

source("individual_pca_functions.R")
# for subsets n=509
# vcf
setwd("~/Dropbox/Mac/Documents/HG/Domestication/01_pca")
vcftools  = "/Users/HG/Dropbox/Mac/Documents/HG/Github/BioinfoTools/vcftools_0.1.13/bin/vcftools";
#system(paste(vcftools," --vcf genetyped_data_n_509_maf05_maxmiss095_popmiss095_hwe.recode.vcf --remove UMFS_2_outlier.txt --recode --recode-INFO-all --out genetyped_data_n_507_maf05_maxmiss095_popmiss095_hwe", sep=""))

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
# output the tab contents for modification
write.table(tab, "sample_eigen_539.txt", row.names=F, sep="\t", quote=F,col.names=T)
tab_pop = read.delim("sample_eigen_539_edit.txt", header = TRUE, sep='\t')
tab_pop = read.delim("sample_eigen_509_edit.txt", header = TRUE, sep='\t')

tab_pop$population_nonum = stringr::str_remove(tab_pop$pop, "[0-9]+")
# count how many individuals in each population
tab_pop %>% count(pop) 

n <- 18
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#col_divergent = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#                    MEW1         MEW2       LIW1      LIW2        DBW1       DBW2      NCW1        NCW2
col_gradient <- c( "#0A2C86", "#849cc1",  "#1D92BD", "#8ad5d9", "#93c47d", "#bedbb1", "#a9a9a9", "#dddddd", 
                   #  DBX1       DBX2      DBX3       UNC1        UNC2       UMFS       NYH1       NEH1       NEH2       MEH2
                   "#f9476b", "#fb90a6","#fddae1", "#cf7fbc",  "#e2b2d6", "#fec155", "#C05805" , "#e1bb94", "#fbd0a5", "#b58383")

#col_gradient = viridis_pal(option = "C")(17)  # n = number of colors seeked
#col_gradient = rep(c("#049DD9", "#F25C05"), c(9, 8))
# PC1-2 for individual populations
order1 = c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "DBX1", "DBX2", "DBX3",  "UNC1", "UNC2", "UMFS", "NYH1", "NEH1", "NEH2", "MEH2")
tab_pop$pop <-factor(tab_pop$pop, levels=order1)
# for 2 clusters
order2 = c("Wild", "Selected")
tab_pop$pop1 <-factor(tab_pop$pop1, levels=order2)

p1 <- ggplot(tab_pop, aes(x = EV1, y = EV2)) + 
  geom_mark_ellipse(aes(fill = pop, label = pop, color=pop), show.legend = FALSE, alpha = 0.2, label.fontface = c("plain"), con.cap=0, label.buffer = unit(0, 'mm'))+
  geom_point(size=2, aes(color=pop, shape=pop1), alpha = 0.8)+
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
  #stat_ellipse(aes(fill = pop), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = F) +
  #scale_fill_manual(values = c("#F25C05", "#049DD9"))
  scale_fill_manual(values = col_gradient) +
  guides(colour = guide_legend(order = 2), 
         shape = guide_legend(override.aes = list(shape = c(0,2)),
                              order = 1))
p2
graph2ppt(file="PCA_539_1",width=4,height=3)
graph2ppt(file="PCA_509_1",width=8,height=6)

jpeg("PCA_509.jpg", width = 8, height = 6, units = 'in', res = 300)
p2
dev.off()

tab_pop = read.delim("sample_eigen_509_edit.txt", header = TRUE, sep='\t')

tab_pop$population_nonum = stringr::str_remove(tab_pop$pop, "[0-9]+")
# count how many individuals in each population
tab_pop %>% count(pop) 

n <- 17
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#col_divergent = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#                    MEW1         MEW2       LIW1      LIW2        DBW1       DBW2      NCW1        NCW2
col_gradient <- c( "#0A2C86", "#849cc1",  "#1D92BD", "#8ad5d9", "#93c47d", "#bedbb1", "#a9a9a9", "#dddddd", 
                   #  DBX1       DBX2      DBX3       UNC1        UNC2       UMFS      NEH1       NEH2       MEH2
                   "#f9476b", "#fb90a6","#fddae1", "#cf7fbc",  "#e2b2d6", "#fec155", "#e1bb94", "#fbd0a5", "#b58383")

#col_gradient = viridis_pal(option = "C")(17)  # n = number of colors seeked
#col_gradient = rep(c("#049DD9", "#F25C05"), c(9, 8))
# PC1-2 for individual populations
order1 = c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "DBX1", "DBX2", "DBX3",  "UNC1", "UNC2", "UMFS", "NEH1", "NEH2", "MEH2")
tab_pop$pop <-factor(tab_pop$pop, levels=order1)
# for 2 clusters
col_gradient = rep(c("#325A98", "#F62A00"), c(9, 8))
order2 = c("Wild", "Selected")
tab_pop$pop1 <-factor(tab_pop$pop1, levels=order2)

p1 <- ggplot(tab_pop, aes(x = EV1, y = EV2)) + 
  geom_mark_ellipse(aes(fill = pop, label = pop, color=pop), show.legend = FALSE, alpha = 0.2, label.fontface = c("plain"), con.cap=0, label.buffer = unit(0, 'mm'))+
  geom_point(size=2, aes(color=pop, shape=pop1), alpha = 0.8)+
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
  #stat_ellipse(aes(fill = pop), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = F) +
  #scale_fill_manual(values = c("#F25C05", "#049DD9"))
  scale_fill_manual(values = col_gradient) +
  guides(colour = guide_legend(order = 2), 
         shape = guide_legend(override.aes = list(shape = c(0,2)),
                              order = 1))
p2

graph2ppt(file="PCA_509_red_blue",width=8,height=6)

# Fst matrix 

library(ggplot2)
library(ComplexHeatmap)
library(circlize)
setwd("~/Dropbox/Mac/Documents/HG/Domestication/05_pairwise_fst")
dt = read.delim("pairwise.fst.txt", header = TRUE, sep='\t')
head(dt)

plot_dt <- as.matrix(dt)
rownames(plot_dt) <- c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "DBX1", "DBX2", "DBX3",  "UNC1", "UNC2", "UMFS", "NYH1", "NEH1", "NEH2", "MEH2")
colnames(plot_dt) <- c("MEW1", "MEW2", "LIW1", "LIW2", "DBW1", "DBW2", "NCW1", "NCW2", "DBX1", "DBX2", "DBX3",  "UNC1", "UNC2", "UMFS", "NYH1", "NEH1", "NEH2", "MEH2")
plot_dt <- as.matrix(plot_dt)
class(plot_dt)

fa = rep(c("WILD", "DOM"), times = c(9, 9))
fa_col = c("WILD" = "#aecdc2", "DOM" = "#f0b8b8")
dend1 = cluster_between_groups(plot_dt, fa)

Heatmap(plot_dt, name = "Pairwise Fst",
        col = colorRamp2(c(0, 0.1, 0.3), c("white", "skyblue", "orange")),
        rect_gp = gpar(col = "white", lwd = 2),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        width = unit(14, "cm"), height = unit(8, "cm"),
        row_names_gp = gpar(fontsize = 12),
        column_names_gp = gpar(fontsize = 12),
        top_annotation = HeatmapAnnotation(Group = fa, col = list(Group = fa_col)))

graph2ppt(file="Pairwise_fst",width=10,height=6)

