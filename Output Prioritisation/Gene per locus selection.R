#################### Gene Per Locus Selection After Gene Priotiziation ####################
#*Restart R after running STRINGdb PPI.R if running this code straight after

#1. Identify direct and secondary PPIs with known BP genes for all genes
#2. Calculate statistics from model scores per locus (mean, SD, variance, SD+1 etc.)
#3. Apply filtering rules to select genes per locus
#4. Plotting distribution of model scores for selected genes
#5. Count number of genes filtered at each step


setwd("~/Documents/PhD Year 2/Ranked Gene Loci/Gene Selection Per Locus")
library(data.table)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)

########################################################################
#1. Identify direct and secondary PPIs with known BP genes for all input scored genes

#Read in PPI counts from interactions in stringdb
interactors <- fread("150_Experimental_interactions.txt")  #file is made from downloadable STRINGdb PPI data which has been processed to list interactors per gene

#Read in file of genes grouped by loci
bp_loci <- fread("bp_loci_ranked_Evangelou_XGB.txt")

#Merge files to give each BP loci gene its matching interactors listed in an additonal column
bp_loci2  <- merge(bp_loci, interactors, by = 'Gene', all.x = TRUE)

#Identify known BP genes - to later find their interactors (to be used in further gene filtering below)
bp_known <- filter(bp_loci2, Training_Score == 1)

#Identify genes with a direct interaction with known BP genes (creating 'direct_PPI_count' column)   
bp_PPI1 <- bp_loci2 %>%
  mutate(direct_PPI_count = str_count(interactors, str_c(bp_known$Gene, collapse = '|')))

#Identify genes that have interactors interacting with interactors of BP genes (creating 'secondary_PPI_count' column)

#Interactors need to be separated to be individually counted (initially all interactors are comma separated in one cell per gene)
sep_rows <- . %>% separate_rows(interactors, sep = ", ")

#Code for bp_PPI2 will produce an empty object if ran straight after STRINGdb PPI code - need to restart R
bp_PPI2 <- bp_PPI1 %>% 
  sep_rows() %>% 
  mutate(
    found = !is.na(match(interactors, sep_rows(bp_known)$interactors))) %>% 
  group_by(Gene) %>% 
  summarise(
    interactors = toString(interactors), 
    secondary_PPI_count = sum(found))

#Combining direct and secondary PPI counts into 1 dataset and dropping interactors column
bp_PPI_total <- merge(bp_PPI1, bp_PPI2, by = 'Gene', all.x = TRUE)
bp_PPI_total <- select(bp_PPI_total, Gene, direct_PPI_count, secondary_PPI_count)

#Getting genes at their loci with their PPI counts
bp_PPI_df <- merge(bp_loci, bp_PPI_total, by = 'Gene', all.x = TRUE)
bp_PPI_df <- unique(bp_PPI_df)

#Re-ordering columns to more easily view genes and their PPI counts
df <- select(bp_PPI_df, loci, Gene, XGB_Score, Training_Score, direct_PPI_count, secondary_PPI_count, rsID_LEAD_SNP,
             Median_Pvalue,  minTRAIT,  Gene_type)

#########################################################################################
#2. Calculate statistics from model scores per locus (to get +1SD - Upper_SD_Threshold' column for filtering)

#Creating statistical columns from model predictions
df <- df %>% group_by(loci) %>% mutate(AvgScore_Per_Locus=mean(XGB_Score,na.rm = T),
                                       diff = max(XGB_Score) - XGB_Score,
                                       Avg_TopGene_Diff_Per_Locus = sum(diff) / (n() - 1),
                                       Variance_Per_Locus = var(XGB_Score,na.rm = T),
                                       IQR_Per_Locus = IQR(XGB_Score,na.rm = T),
                                       SD_Per_Locus = sd(XGB_Score,na.rm = T),
                                       LocusGenes = toString(Gene))

df$Upper_SD_Threshold <- df$AvgScore_Per_Locus + df$SD_Per_Locus

###################################################
#3. Apply filtering rules to select genes per locus

#Filtering to select gene(s) per locus if they are:
#1. >Upper_SD Threshold (+1 SD) OR
#2. If multiple genes between average score - Upper SD then filter by PPI
#3. Gene with highest direct PPI chosen OR if matching direct_PPI_count:
#4. Gene with highest secondary PPI chosen OR if matching secondary_PPI_count:
#5. all genes left chosen

df$direct_PPI_count[is.na(df$direct_PPI_count)] <- 0

#Identifying training genes scored at 1 (known BP genes that automatically are selected per loci)
known_BP_genes <- filter(df, Training_Score==1)

#Identifying loci where only 1 gene is present
one_gene <-  df %>%
  group_by(loci) %>%
 filter(n() == 1)

#Identifying genes that meet the first condition of being > +1SD (Upper_SD_Threshold):
UpperSD_Genes <- df %>%
  group_by(loci) %>%
  #filtering within loci with multiple genes (loci with already singular genes not entering filter):
  filter(if(n() > 1) {(XGB_Score > Upper_SD_Threshold) } else TRUE) %>%
  filter(if(n() > 1) {XGB_Score > Upper_SD_Threshold & direct_PPI_count == max(direct_PPI_count)} else TRUE) %>%
  filter(if(n() > 1) {XGB_Score > Upper_SD_Threshold & secondary_PPI_count == max(secondary_PPI_count)} else TRUE)

#Identifying loci where only 1 gene is > average model score (and therefore no further PPI filtering is needed):
Avg_Genes <- df %>% 
  group_by(loci) %>% 
  subset(!(loci %in% UpperSD_Genes$loci)) %>%
  subset(!(Gene %in% known_BP_genes$Gene)) %>%
  filter(XGB_Score > AvgScore_Per_Locus) %>%
  filter(n() == 1)

#Selecting gene(s) per locus depending direct and secondary PPI filters:
PPI_Genes <- df %>% 
  group_by(loci) %>% 
  subset(!(loci %in% UpperSD_Genes$loci)) %>%
  subset(!(loci %in% Avg_Genes$loci)) %>%
  filter(XGB_Score >= AvgScore_Per_Locus) %>%
  filter(direct_PPI_count == max(direct_PPI_count)) %>%
  filter(secondary_PPI_count == max(secondary_PPI_count)) 

#Combining selected genes per loci:
new_df <- rbind(UpperSD_Genes, Avg_Genes, PPI_Genes, known_BP_genes, one_gene)
new_df <- new_df %>% arrange(loci)

loci_ranked_df <- select(new_df, loci, Gene, XGB_Score, Training_Score, LocusGenes, AvgScore_Per_Locus,
                SD_Per_Locus,Upper_SD_Threshold, direct_PPI_count, secondary_PPI_count,
                Avg_TopGene_Diff_Per_Locus,
                Variance_Per_Locus, IQR_Per_Locus, rsID_LEAD_SNP,
                minTRAIT, Median_Pvalue, Gene_type)

colnames(loci_ranked_df)[1] <- "Loci"
colnames(loci_ranked_df)[3] <- "XGB Score"
colnames(loci_ranked_df)[5] <- "LocusGene(s)"
loci_ranked_df <- loci_ranked_df[!duplicated(loci_ranked_df[2:17]),]
loci_ranked_df2 <- loci_ranked_df[gtools::mixedorder(loci_ranked_df$Loci), ]
write.csv(loci_ranked_df2, "./all_loci_ordered_and_selected.csv", row.names = FALSE)


loci_ranked_df2 <- loci_ranked_df[!duplicated(loci_ranked_df[,c('Gene')]),]
loci_ranked_df2 <- loci_ranked_df2[gtools::mixedorder(loci_ranked_df2$Loci), ]
loci_ranked_df2$Loci <- cumsum(!duplicated(loci_ranked_df2$Loci))

count_loci <- filter(loci_ranked_df, `XGB Score` <= 0.5)


write.csv(loci_ranked_df2, "./Gene_Per_Locus_Selection_Evangelou_XGB.csv", row.names = FALSE)

sexbias <- filter(loci_ranked_df2, GTEx_signif_sexbias==1)
sexbias <- sexbias %>% ungroup(Loci) %>% select(Gene, `XGB Score`)
fwrite(sexbias, "GTExsexbiased_selected_genes_evangelou.txt", sep = '\t',row.names = FALSE)

low <- filter(loci_ranked_df, `XGB Score` <0.5)
low <- low %>% ungroup(Loci) %>% select(Gene, `XGB Score`)
write.csv(low, "lowscored_selected_genes.txt", row.names = FALSE)

mismatchextract <- subset(bp_loci, !(Gene %in% loci_ranked_df2$Gene)) 
mismatchextract <- select(mismatchextract, Gene, XGB_Score)
mismatchextract <- mismatchextract[!duplicated(mismatchextract[,c('Gene')]),]
fwrite(mismatchextract, "unselected_genes.txt", sep = '\t', row.names = FALSE)

selected <- loci_ranked_df2 %>% ungroup(Loci) %>% select(Gene, `XGB Score`)
fwrite(selected, "selected_genes.txt", sep = '\t', row.names = FALSE)

genes05 <- df %>% ungroup(loci) %>% filter(XGB_Score >=0.5) %>% select(Gene, XGB_Score)
genes05  <- genes05[!duplicated(genes05[,c('Gene')]),]
fwrite(genes05, "genes_0.5_threshold.txt", sep = '\t', row.names = FALSE)

genes06 <- df %>% ungroup(loci) %>% filter(XGB_Score >=0.6) %>% select(Gene, XGB_Score)
genes06  <- genes06[!duplicated(genes06[,c('Gene')]),]
fwrite(genes06, "genes_0.6_threshold.txt", sep = '\t', row.names = FALSE)

genes07 <- df %>% ungroup(loci) %>% filter(XGB_Score >=0.7) %>% select(Gene, XGB_Score)
genes07  <- genes07[!duplicated(genes07[,c('Gene')]),]
fwrite(genes07, "genes_0.7_threshold.txt", sep = '\t', row.names = FALSE)
genes07_training <- merge(genes07, known_BP_genes, all=TRUE)
genes07_training  <- genes07_training[!duplicated(genes07_training[,c('Gene')]),]
genes07_training <- select(genes07_training, Gene, XGB_Score, Training_Score)
fwrite(genes07_training, "genes_0.7_and_training.txt", sep = '\t', row.names = FALSE)


genes08 <- df %>% ungroup(loci) %>% filter(XGB_Score >=0.8) %>% select(Gene, XGB_Score)
genes08  <- genes08[!duplicated(genes08[,c('Gene')]),]
fwrite(genes08, "genes_0.8_threshold.txt", sep = '\t', row.names = FALSE)
genes08_training <- merge(genes08, known_BP_genes, all=TRUE)
genes08_training  <- genes08_training[!duplicated(genes08_training[,c('Gene')]),]
genes08_training <- select(genes08_training, Gene, XGB_Score, Training_Score)
fwrite(genes08_training, "genes_0.8_and_training.txt", sep = '\t', row.names = FALSE)


#########################################################
#4. Plotting selected genes per locus score distribution:

barfill <- "#4271AE"
barlines <- "#1F3552"

p <- ggplot(loci_ranked_df, aes(x=`XGB Score`))+
  geom_histogram(color="#1F3552", fill="steelblue") + 
  scale_x_continuous(breaks = seq(0, 1.2, by = 0.1)) + #limits = c(0, 1)
  scale_y_continuous(breaks = seq(0, 300, by = 10)) +
  theme(axis.text.x = element_text(color = "grey20", size = 16, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain"))+
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title=element_text(size=14), legend.text=element_text(size=14)) +
  ggtitle("Extreme Gradient Boosting Score Distribution (Selected Genes Per Locus)")  + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Gene Count") 

p

ggsave("bp_loci_Distribution_selected_Evangelou_XGB.png",width = 10, height = 8, dpi=300,
       limitsize = FALSE)

###############################################
#5. Count number of genes filtered at each step

#Counting number of loci with >1 gene selected:
multigeneloci <- loci_ranked_df2[duplicated(loci_ranked_df2$Loci) | duplicated(loci_ranked_df2$Loci, fromLast = TRUE), ] 
length(unique(multigeneloci$Loci)) #2 loci with multiple genes still selected
#write.csv(multigeneloci , "./duploci.csv", row.names = FALSE) 
avgcount <- Avg_Genes[!duplicated(Avg_Genes[,c('Gene')]),]
sdcount <- UpperSD_Genes[!duplicated(UpperSD_Genes[,c('Gene')]),]
#Number of genes that were selected by 1st filter (> +1 SD model score) are in UpperSD_Genes object 
nrow(sdcount)
#Number of genes that were selected by being the only gene > the average score are in Avg_Genes object 
nrow(avgcount)

length(unique(PPI_Genes$loci)) #loci entering PPI filtering

#Count number of loci (seeing how many still have >1 gene past direct PPI filter - those entering loci_ranked_df filter)
count1stPPI <- df %>% 
  group_by(loci) %>% 
  subset(!(loci %in% UpperSD_Genes$loci)) %>%
  subset(!(loci %in% Avg_Genes$loci)) %>%
  filter(XGB_Score >= AvgScore_Per_Locus) %>%
  filter(direct_PPI_count == max(direct_PPI_count))

count1stPPI1 <- select(count1stPPI, Gene, loci)
PPI1_dup <- unique(count1stPPI1[duplicated(count1stPPI1$loci), "loci"])
count1stPPI2 <-count1stPPI1[!count1stPPI1$loci %in% PPI1_dup$loci,]

