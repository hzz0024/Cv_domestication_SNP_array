library(dplyr)
library(stringr)

setwd("~/Dropbox/Mac/Documents/HG/Domestication/22_enrichment/enrichment")

vcftools  = "/Users/HG/Dropbox/Mac/Documents/HG/Github/BioinfoTools/vcftools_0.1.13/bin/vcftools";
plink  = "/Users/HG/Dropbox/Mac/Documents/HG/Domestication/14_ROH/plink";
system(paste(plink, " --vcf genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe.recode.vcf --allow-extra-chr --make-bed --out genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe", sep=""))


##############
# annotation #
##############
# it is better to run command below directly in the terminal
system("python3 extract_gene_V2.py -i Gene_annotation_all.txt -t genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe.bed -o genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe.bed.gene.txt")

###################
# Extract GO data #
###################
# for whole SNPs
LOC_all <- read.table("genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe.bed.gene.txt", header = FALSE, sep = "\t", quote="", col.names = paste0("V",seq_len(6)), fill = TRUE)
#LOC_all <- read.delim("All_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers.bed.gene.txt", header = FALSE, sep='\t')
dim(LOC_all)
LOC_all <- LOC_all[which(LOC_all$V6 != ""),]
dim(LOC_all)
# from NCBI protein database
#emapper = read.delim("out.emapper.annotation.withLOC.uniq.txt", header = TRUE, sep='\t')
emapper = read.table("out.emapper.annotation.withLOC.uniq.txt", header = FALSE, sep='\t', quote="", comment.char="#", col.names = paste0("V",seq_len(21)), fill = TRUE)
dim(emapper)
#emapper$X.Locus
# extract the genes mapped by whole SNPs
extract_df <-emapper[which(emapper$V1 %in% paste0("LOC",LOC_all$V6)),]
#extract_df <-emapper[which(emapper$X.Locus %in% paste0("LOC",LOC_all$V6)),]
dim(extract_df)
extract_df_filter <- extract_df[which(extract_df$V10 != "-"),]
dim(extract_df_filter)
write.table(extract_df_filter, file = "domestication.wholeSNPs.emapper.annotation.withLOC.uniq.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#######################
# Format outliers     #
#######################
library(stringr)
library(dplyr)
# format_outlier <- function(fname){
#   output_list <- read.delim(fname, header = TRUE, sep='\t')
#   df <- data.frame(output_list$chr, output_list$pos, output_list$pos)
#   df <- as.data.frame(df)
#   write.table(df, file = paste0(fname, ".bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
# }

#fname = 'permutation_outliers_n_45968.txt.bed.gene.txt'
format_LOC <- function(fname){
  gene_list <- read.delim(fname, header = FALSE, sep='\t')
  dim(gene_list)
  gene_list <- gene_list[which(gene_list$V6 != ""),]
  dim(gene_list)
  all_list <- str_split(gene_list$V6, ";")
  unique_list <- sort(paste0("LOC", unique(unlist(all_list))))
  write.table(unique_list, file = paste0(fname, ".unique.LOC.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

system(paste("python3 extract_gene_V2.py -i Gene_annotation_all.txt -t Predictor_loci_classification_best_10perc_unique_37.txt -o Predictor_loci_classification_best_10perc_unique_37.gene.txt"))
format_LOC("Predictor_loci_classification_best_10perc_unique_37.gene.txt")

system(paste("python3 extract_gene_V2.py -i Gene_annotation_all.txt -t Predictor_loci_classification_best_10perc_unique_37_10K.txt -o Predictor_loci_classification_best_10perc_unique_37_10K.gene.txt"))
format_LOC("Predictor_loci_classification_best_10perc_unique_37_10K.gene.txt")

###################
# Build R package #
###################
system(paste("Rscript 17_create_orgdb_from_emapper.R domestication.wholeSNPs.emapper.annotation.withLOC.uniq.txt"))

#######################
# enrichment analysis #
#######################

## load packages
library(tidyverse)
library(clusterProfiler)
library(stringr)
library(export)
# set up working location

dir.create('R_Library', recursive = T)

## prepare GO and KEGG lib
install.packages('./org.My.eg.db_1.0.tar.gz', 
                 repos = NULL,
                 lib = 'R_Library') 

# load OrgDB
library(org.My.eg.db, lib = 'R_Library')

# load gene ID
shared_SNP <- read.delim("Predictor_loci_classification_best_10perc_unique_37_10K.gene.txt.unique.LOC.txt", header = F, sep='\t')$V1
shared_SNP <- read.delim("Predictor_loci_classification_best_10perc_unique_37.gene.unique.LOC.txt", header = F, sep='\t')$V1

# GO enrichment
shared_SNP_go <- enrichGO(gene = shared_SNP,
                             OrgDb = org.My.eg.db,
                             keyType = 'GID',
                             ont = 'ALL',
                             qvalueCutoff = 0.1,
                             pvalueCutoff = 0.1)

shared_SNP_go_df <- as.data.frame(shared_SNP_go)

write.table(shared_SNP_go_df, file = "shared_SNP_go_df.txt", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

head(SGS_singleSNP_go_df)

require(DOSE)
require(enrichplot)
require(viridis)
## plot

barplot(shared_SNP_go, showCategory = 10, split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale="free")

# key results for GO
dotplot(shared_SNP_go, showCategory = 10, split="ONTOLOGY") + 
  scale_color_viridis(option="plasma") + 
  facet_grid(ONTOLOGY~., scale="free")+ 
  scale_y_discrete(labels=function(x) str_wrap(x, width=60)) + 
  labs(size="Counts",col="FDR")
graph2ppt(file="permutation_outliers_n_45968_go.pptx", width=10, height=10)

cnetplot(shared_SNP_go, 
         showCategory = 5,
         node_label = "all", # category | gene | all | none
         circular = TRUE, 
         colorEdge = TRUE)

###################################
#########  KEGG enrichment ########
###################################

emapper <- read_delim('./domestication.wholeSNPs.emapper.annotation.withLOC.uniq.txt', 
                      "\t", escape_double = FALSE, col_names = FALSE, 
                      comment = "#", trim_ws = TRUE) %>%
  dplyr::select(GID = X1, 
                KO = X12, 
                Pathway = X13)

# load gene ID
SGS_single_SNP <- read.delim("pop_n_477_pcadapt_outflank.shared.10k.outlier.igv_Dom_Wild.sliding.zfst.outlier.merged.igv.outlier.merged.igv.gene.txt.unique.LOC.txt", header = F, sep='\t')$V1

pathway2gene <- dplyr::select(emapper, Pathway, GID) %>%
  separate_rows(Pathway, sep = ',', convert = F) %>%
  filter(str_detect(Pathway, 'ko')) %>%
  mutate(Pathway = str_remove(Pathway, 'ko'))

library(magrittr)
get_path2name <- function(){
  keggpathid2name.df <- clusterProfiler:::kegg_list("pathway")
  keggpathid2name.df[,1] %<>% gsub("path:map", "", .)
  colnames(keggpathid2name.df) <- c("path_id","path_name")
  return(keggpathid2name.df)
}

pathway2name <- get_path2name()

SGS_singleSNP_ekp <- enricher(SGS_single_SNP,
                              TERM2GENE = pathway2gene,
                              TERM2NAME = pathway2name,
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.05)

SGS_singleSNP_ekp_df <- as.data.frame(SGS_singleSNP_ekp)

head(SGS_singleSNP_ekp_df)

## plot
barplot(SGS_singleSNP_ekp, showCategory = 10)

dotplot(SGS_singleSNP_ekp, showCategory = 10) + 
  scale_color_viridis(option="plasma") + 
  labs(size="Counts",col="FDR") +
  scale_size_area(max_size = 14)

graph2ppt(file="SGS_singleSNP_KEGG.pptx", width=8, height=8)

cnetplot(SGS_singleSNP_ekp, 
         #foldChange = geneList, 
         showCategory = 3,
         node_label = "all", # category | gene | all | none
         circular = FALSE, 
         colorEdge = TRUE)
## -----------------------------------------------------------
save(SGS_singleSNP_go, SGS_singleSNP_ekp, file = './SGS_singleSNP_enrich.rdata')

