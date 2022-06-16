######Ordering Gene into Loci by LD######
#1. Load all necessary files (KnownBP variants, Gene annotated file, and XGBoost predictions
#2. Idenitfy lead SNPs from KnownBP variant file (KnownBP_Apr2020_LDr2-8_500kb_gwasRes.txt) and sort into individual groups
#3. Grouping variants with Gene annotated if they are within 500kb+/- of lead SNPs - identifying Gene at loci
#4. Merge in extra columns of information to Gene now grouped at their loci

setwd("~/Documents/PhD Year 2/Ranked Gene Loci")
library(tidygraph)
library(tidyverse)
library(plyr)
library(data.table)
library(magrittr)
library(janitor)
options(scipen = 999)

###################################################################################################
#1. Load all necessary files (BP associated variants, Gene annotated file, and XGB predictions)

#Creating list of all lead snp ids and their CPs with their annotated genes:
bp_loci <- fread('KnownBP_Apr2020_LDr2-8_500kb_gwasRes.txt')
bp_loci <- select(bp_loci, CP, CP_LEADSNP, rsID_LEAD_SNP, 
                 Source, Type, LocusName,minP, minTRAIT)

#Need the gene annotations for all BP loci to be able to merge with XGB gene predictions by gene:
BP_genes <- fread("Genes_bp_loci.bed", fill = TRUE)
colnames(BP_genes)[1] <- 'Chr'
colnames(BP_genes)[2] <- 'Start'
colnames(BP_genes)[3] <- 'End'
colnames(BP_genes)[12] <- 'CP'
#p-value and minTRAIT for adding to all_loci_and_scores ordered loci dataset later
colnames(BP_genes)[9] <-'minP'
colnames(BP_genes)[10] <-'minTRAIT' 
colnames(BP_genes)[16] <-'Gene'

total_CP <- select(BP_genes, Gene, CP, Chr, Start, End)

total_sentinel_CP <- join(bp_loci, total_CP)
total_sentinel_CP  <- total_sentinel_CP[!duplicated(total_sentinel_CP), ]

#Loading in XGB predictions for all genes and combining into one dataset
XGB <- fread('BPgene_ranking_xgb_evangelou.txt')
training <- fread('xgb_training_evangelou.csv')

#Removing 'least likely' genes which were not associated in GWAS and so are not in BP loci
training <- filter(training, BPlabel_encoded != "0.1" | is.na(BPlabel_encoded))

colnames(XGB)[1] <- 'Gene'
colnames(training)[1] <- 'Gene'
colnames(training)[2] <- 'XGB_training_prediction'

BP_genes_scores1 <- merge(total_CP , XGB, by = "Gene", all.x = TRUE)
BP_genes_scores2 <- merge(BP_genes_scores1, training, by = "Gene", all.x = TRUE)

df1 <- filter(BP_genes_scores2, !is.na(BP_genes_scores2$XGB_Score))
df2<- filter(BP_genes_scores2, !is.na(BP_genes_scores2$XGB_training_prediction))
df <- rbind(df1, df2)
df <- unique(df)

###########################################################################################################################
#2. Idenitfy lead SNPs from Known BP variant file (KnownBP_Apr2020_LDr2-8_500kb_gwasRes.txt) and sort into individual groups

#Selecting lead SNPs only:
leads <- filter(total_sentinel_CP, !is.na(total_sentinel_CP$CP_LEADSNP))
leads <- filter(leads, !is.na(leads$Start)) 

#Separating multiple inputs in CP_LEADSNP column to individual CPs to use as lead SNP boundaries
#Individual lead CPs grouped by originally being in the same row in Known BP variant file
#If other rows share any same CPs they also join that group with highest and lowest CPs being the boundaries for that locus
#Groups assigned numbers for each locus
leads <- data.table(leads)

#Tidygraph package codes each lead SNP as a node within a graph
#Connected lead SNPs (those with matching lead SNP CPs or multiple lead CPs per row) then form a group
lead_group <- leads %>% 
  separate_rows(CP_LEADSNP, sep = ", ") %>% 
  as_tbl_graph() %>% 
  activate(nodes) %>% 
  #mutate produces error if dplyr package isn't named to be used:
  dplyr::mutate(loci = group_components()) %>% 
  as_tibble()

#lead_group head:
#name         loci
#1:10000476	  106
#1:100829685	344

colnames(lead_group)[1] <- 'CP'

#Combine loci groupings with all other columns for lead SNPs
lead_group <- merge(lead_group, leads, by = 'CP', all.x = T)

lead_LD <- lead_group %>% arrange(loci)
#Create Chr and Start columns from CP column (to use Start position as boundaries per loci)
lead_LD$CP <- gsub(":", ",", lead_LD$CP)
lead_df1 <- data.frame(do.call('rbind', strsplit(as.character(lead_LD$CP), ',', fixed = TRUE)))
#Append loci column to Chr and Start columns
lead_df <- cbind(lead_df1, lead_LD[2])
#lead_df:
#	X1  X2       loci
# 17	43572419	1
#	17	43572896	1

#Identify largest and smallest CPs per locus group
max_ld <- aggregate(.~loci, lead_df, FUN = max, na.rm = TRUE, na.action = NULL)
max_ld$max <- max_ld$X2
min_ld <- aggregate(.~loci, lead_df, FUN = min, na.rm = TRUE, na.action = NULL)
min_ld$min <- min_ld$X2

#Combine min and max columns together in the same dataset
min_max_ld <-cbind(min_ld, max_ld)
#Remove repeated columns
min_max_ld <- min_max_ld[ -c(5:7)]

#############################################################################################################
#3. Grouping variants with Gene annotated if they are within 500kb+/- of lead SNPs - identifying Gene at loci

#Rename Chr and Start to match BP loci dataset
colnames(min_max_ld)[2] <- 'chrom'
min_max_ld$chrom <- as.numeric(min_max_ld$chrom)
min_max_ld$max <- as.numeric(min_max_ld$max)
min_max_ld$min <- as.numeric(min_max_ld$min)

df$chrom <- as.numeric(df$Chr)
df$position <- as.numeric(df$Start)

#Using lead snp CPs to set 500kb +/- from min and max CPs per loci to create boundaries:
setDT(min_max_ld)
min_max_ld[, c("low", "high") := .(min - 500000, max  + 500000)]

#Find prioritised gene's variant matches on chromosome, with position between low&high

#Create ID column for a reference point when grouping variants into loci (ID number will match loci number)
min_max_ld$ID <- seq.int(nrow(min_max_ld))
df$ID <- seq.int(nrow(df))

#Order genes by loci - genes fitting into overlapping loci position will become one group
df[min_max_ld, loci := i.ID, on = .(chrom, position > low, position < high ) ]

#Order loci/groups 
dt <- unique(df)
dt <- dt %>% arrange(loci)
dt$loci <- cumsum(!duplicated(dt$loci))

#Check for NAs - Genes with CP_LEADSNP in one study and NA in another (creating NAs in loci column)
na_test <- filter(df, is.na(loci))
na_Gene <- filter(na_test, !duplicated(Gene))
na_snps <-  total_sentinel_CP[Gene %in% na_Gene$Gene]

#write.csv(df, "./AllLoci_ranked.csv", row.names = FALSE)

ordered_loci_dt <- select(dt, loci, Gene, XGB_Score, BPlabel_encoded,
                          XGB_training_prediction)
dt_ordered = unique(ordered_loci_dt)
colnames(dt_ordered)[2] <- 'Gene'

###########################################################################
#4. Merge in extra columns of information to Gene now grouped at their loci

#Merge in extra column and comparative method scores per gene
gene_length <- fread('GeneLength.txt')
colnames(gene_length)[1] <- 'Gene'
colnames(dt_ordered)[2] <- 'Gene'
gtex_sex_bias <- fread('GTEx_sexbias.txt')
CPs_Per_Gene <-  as.data.table(df)[, toString(CP), by = Gene]
colnames(CPs_Per_Gene)[1] <- 'Gene'
colnames(CPs_Per_Gene)[2] <- 'CPs'

df_rank <- Reduce(function(x, y) merge(x, y, by = 'Gene', all.x = TRUE), 
                  list(dt_ordered, CPs_Per_Gene, gene_length, gtex_sex_bias))
#Combine training model scores and training scores into one prediction score column:
df_rank$XGB_Score <- coalesce(df_rank$XGB_Score, df_rank$XGB_training_prediction)

data1 <- fread('GPrior3label_unknown.txt')
data2 <- fread('GPrior3label_training.txt')
data3 <- fread('toppgene_results.csv')
mantisml <- fread('mantisml.csv')

all_scores <- Reduce(function(x, y) merge(x, y, by = 'Gene', all = TRUE), 
                     list(data1,data2, data3,mantisml))
all_scores  <- all_scores[!duplicated(all_scores$Gene), ]

all_scores$GPrior3label <- coalesce(all_scores$Gprior3label_training, 
                                    all_scores$Gprior3label_unknown)

scores <- all_scores

x <- scores$GPrior3label
scores$GPrior3label <- x/100

scores_matches <- scores[Gene %in% df_rank$Gene]

all_loci_and_scores <- merge(df_rank, scores_matches, by = 'Gene', all.x = T)

#Get all lead snp information per locus and add to all_loci_and_scores dataset
summ_leads <- leads %>% group_by(Gene) %>% summarise_all(toString)
summ_leads <- summ_leads[, c(1, 3:7)]
colnames(summ_leads)[1] <-'Gene'

all_loci_and_scores <- merge(all_loci_and_scores, summ_leads, by = 'Gene', all.x = TRUE)

#Get minimum p-values per gene and minTRAIT per gene:

med_pvalue <- fread('median_pvalues.csv')
all_loci_and_scores2 <- merge(all_loci_and_scores, med_pvalue, by = 'Gene', all.x = TRUE)


#minTRAIT:

bp_loci_genes <- select(BP_genes, Gene, minTRAIT)
bp_loci_genes <- setDT(bp_loci_genes)[, lapply(.SD, paste, collapse = ", "), by = Gene]

bp_loci_trait <- select(bp_loci_genes, Gene, minTRAIT)
all_loci_and_scores2 <- merge(all_loci_and_scores2, bp_loci_trait, by ='Gene', all.x = TRUE)

#Clean minTRAIT column to remove dupliXGBes
all_loci_and_scores2$minTRAIT <- sapply(all_loci_and_scores2$minTRAIT , function(x) paste(unique(unlist(str_split(x,", "))), collapse = ", "))
all_loci_and_scores2$minTRAIT <-gsub(" , ", "", all_loci_and_scores2$minTRAIT)
all_loci_and_scores2$minTRAIT <-gsub("^, ", "", all_loci_and_scores2$minTRAIT)
all_loci_and_scores2$minTRAIT <-gsub(", $", "", all_loci_and_scores2$minTRAIT)

all_loci_and_scores3 <- select(all_loci_and_scores2, loci, Gene, BPlabel_encoded,
                XGB_Score, Median_Pvalue, GPrior3label,
                ToppGeneScore, mantis_ml_proba, GeneLength, GTEx_signif_sexbias, rsID_LEAD_SNP, Source, Type, LocusName, CP_LEADSNP, CPs, 
                minTRAIT)

colnames(all_loci_and_scores3)[3]  <- 'Training_Score'
colnames(all_loci_and_scores3)[6]  <- 'GPrior_Score'
colnames(all_loci_and_scores3)[7]  <- 'ToppGene_Score'
colnames(all_loci_and_scores3)[8]  <- 'Mantis_probability'

gene_type <- fread('hg19Rel92_AllgeneTypes_0kb.txt')
colnames(gene_type)[4] <- 'Gene'
colnames(gene_type)[5] <- 'Gene_type'
gene_type <- select(gene_type, Gene, Gene_type)
all_loci_and_scores3 <- merge(all_loci_and_scores3, gene_type, by ='Gene', all.x = TRUE)


all_loci_and_scores3 <- all_loci_and_scores3[order(as.numeric(sub(':.*', '', all_loci_and_scores3$loci))), ]
all_loci_and_scores3 <- all_loci_and_scores3[gtools::mixedorder(all_loci_and_scores3$loci), ]
all_loci_and_scores3$loci <- cumsum(!duplicated(all_loci_and_scores3$loci))
all_loci_and_scores3 <- dplyr::filter(all_loci_and_scores3, Gene != 'Y_RNA')
fwrite(all_loci_and_scores3, "./generanked_allLoci.txt", row.names = FALSE)


#Removing genes with unknown gene IDs for other methods:
all_loci_and_scores_df <- all_loci_and_scores3 %>% 
  filter_at(vars(GPrior_Score, ToppGene_Score, Mantis_probability), any_vars(!is.na(.)))
all_loci_and_scores_df <- all_loci_and_scores_df[order(as.numeric(sub(':.*', '', all_loci_and_scores_df$loci))), ]
all_loci_and_scores_df <- all_loci_and_scores_df[gtools::mixedorder(all_loci_and_scores_df$loci), ]
all_loci_and_scores_df$loci <- cumsum(!duplicated(all_loci_and_scores_df$loci))

all_loci_and_scores_pval <- filter(all_loci_and_scores_df, Median_Pvalue <= 0.00000005)
fwrite(all_loci_and_scores_pval, "./EvangelouGenes_pvalue_insignif.txt", row.names = FALSE)


all_loci_and_scores_df <- all_loci_and_scores_df[,-c(12:16)]
all_loci_and_scores_df$rsID_LEAD_SNP <-gsub("NA, ", "", all_loci_and_scores_df$rsID_LEAD_SNP)
all_loci_and_scores_df$rsID_LEAD_SNP <-gsub(", NA", "", all_loci_and_scores_df$rsID_LEAD_SNP)
all_loci_and_scores_df <- select(all_loci_and_scores_df, 'Gene',	'loci',	'Training_Score',	'XGB_Score',
                                 'Median_Pvalue',	'Gene_type',	'rsID_LEAD_SNP',	'minTRAIT')

all_loci_and_scores_df <- all_loci_and_scores_df[!duplicated(all_loci_and_scores_df), ]

fwrite(all_loci_and_scores_df, "./bp_loci_ranked_Evangelou_XGB.txt", row.names = FALSE)

dups <- all_loci_and_scores_df[!duplicated(all_loci_and_scores_df[,c('Gene')]),]
dups <- dups[gtools::mixedorder(dups$loci), ]
dups$loci <- cumsum(!duplicated(dups$loci))

dups <- all_loci_and_scores_df %>% group_by(Gene) %>% filter(n() > 1)


test <- all_loci_and_scores_df %>%
    group_by(loci) %>%
   filter(n_distinct(Gene) > 1) %>%
    ungroup()
dups2 <- test %>% group_by(Gene) %>% filter(n() > 1)


fwrite(dups, "./no_dupliXGBeLoci_bp_loci_ranked_Evangelou_XGB.txt", row.names = FALSE)


ggplot(all_loci_and_scores3, aes(x=XGB_Score))+
  geom_histogram(color="darkblue", fill="lightblue") +
  ggtitle("XGBoost Prioritisation Score Distribution")  +  
  ylab("Gene Count") +
  theme(text = element_text(size=16), plot.title = element_text(size = 16, hjust = 0.5))

ggsave("bp_loci_Distribution_Evangelou_XGB.tiff",width = 9, height = 4, dpi=300,
       limitsize = FALSE)