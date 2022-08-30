Data Preprocessing
================
Hannah L Nicholls
2021-07-18

## Contents

[1. Filtering insignificant
variants](#1.%20Filtering%20insignificant%20variants)<br /> [2. GWAS
variant p-value
filtering](#2.%20GWAS%20variant%20p-value%20filtering)<br /> [3. GWAS
variant data filtering](#2.%20GWAS%20variant%20data%20filtering)<br />
[4. Protein-protein interaction
filtering](#3.%20Protein-protein%20interaction%20filtering)<br /> [5.
Variant-level feature
processing](#4.%20Variant-level%20feature%20processing)<br /> [6.
Feature pre-procesing](#5.%20Feature%20pre-procesing)<br /> [7. Merging
databases and dividing training
data](#6.%20Merging%20databases%20and%20dividing%20training%20data)<br />
[8. Gene length correlation](#8.%20Gene%20length%20correlation)<br />
[9. SessionInfo](#9.%20SessionInfo)

## Necessary Libraries

``` r
library(stringr)
library(plyr)
library(dplyr)
library(tidyverse)
library(data.table)
library(janitor)
library(sqldf)
library(tidyr)
library(splitstackshape)
library(matrixStats) 
library(GenomicRanges)
options(scipen=999)
memory.limit(size=56000)
```

    ## [1] Inf

## <a id="1. Filtering insignificant variants"></a>1. Filtering insignificant variants

-   Insignificant variants defined as having a p-value \>0.15, no
    linkage disquilibrium r2 measures, and are not within 500kb +/-
    blood pressure loci
-   Insignificant variants once identified and annotated to genes later
    go further filtering by protein-protein interactions (point 3.)
-   All variants from the whole GWAS are also taken in this first
    section to be annotated in unix scripts (providing gene annotations
    per variants and ANNOVAR pathogenic variant annotations)

### BP associated variants across all GWAS as of April 2020:

-   Includes Evangelou et al. variants and enables identification of
    newer associated variants to avoid them being identified as
    insignificant variants

``` r
#Selecting non-BP loci variants only

#Read in signficant gwas variants file
cp_bp <- fread("KnownBP_Apr2020_LDr2-8_500kb_gwasRes.txt", fill = TRUE)
head(cp_bp)
```

    ##             CP             CP_LEADSNP               CP_LDSNP rsID_LEAD_SNP
    ## 1:  1:10000476 1:10000476, 1:10357806 1:10000476, 1:10000476   rs182770070
    ## 2:  1:10081736             1:10000476             1:10081736          <NA>
    ## 3: 1:100829685            1:100829685            1:100829685   rs114558965
    ## 4: 1:100829720            1:100829685            1:100829720          <NA>
    ## 5: 1:100833868            1:100829685            1:100833868          <NA>
    ## 6:  1:10089376             1:10000476             1:10089376          <NA>
    ##    dupCount LocusName                       Source      Type CHR       POS
    ## 1:        1       n/a BPICE_UKBfinemapCommon(1210) Secondary   1  10000476
    ## 2:       NA      <NA>                         <NA>      <NA>  NA        NA
    ## 3:        1    CDC14A GxL_Fuentes_Educ_T2nov_TRANS      Lead   1 100829685
    ## 4:       NA      <NA>                         <NA>      <NA>  NA        NA
    ## 5:       NA      <NA>                         <NA>      <NA>  NA        NA
    ## 6:       NA      <NA>                         <NA>      <NA>  NA        NA
    ##    snpid chr bpos   a1   a2 freq BETAsbp Psbp BETAdbp Pdbp BETApp Ppp minP
    ## 1:  <NA>  NA   NA <NA> <NA>   NA      NA   NA      NA   NA     NA  NA   NA
    ## 2:  <NA>  NA   NA <NA> <NA>   NA      NA   NA      NA   NA     NA  NA   NA
    ## 3:  <NA>  NA   NA <NA> <NA>   NA      NA   NA      NA   NA     NA  NA   NA
    ## 4:  <NA>  NA   NA <NA> <NA>   NA      NA   NA      NA   NA     NA  NA   NA
    ## 5:  <NA>  NA   NA <NA> <NA>   NA      NA   NA      NA   NA     NA  NA   NA
    ## 6:  <NA>  NA   NA <NA> <NA>   NA      NA   NA      NA   NA     NA  NA   NA
    ##    minTRAIT BETAmean
    ## 1:     <NA>       NA
    ## 2:     <NA>       NA
    ## 3:     <NA>       NA
    ## 4:     <NA>       NA
    ## 5:     <NA>       NA
    ## 6:     <NA>       NA

### Whole GWAS data:

``` r
head(full_gwas)
```

    ##              snpid chr     bpos a1 a2   freq BETAsbp    Psbp BETAdbp   Pdbp
    ## 1:   1:2556125_C_T   1  2556125  t  c 0.3255 -0.0262 0.41300 -0.0113 0.5388
    ## 2:   1:2556548_C_T   1  2556548  t  c 0.3261 -0.0274 0.39270 -0.0121 0.5096
    ## 3:   1:2556709_G_A   1  2556709  a  g 0.3257 -0.0263 0.41210 -0.0116 0.5266
    ## 4: 12:11366987_C_T  12 11366987  t  c 0.9443  0.0355 0.61460  0.0019 0.9631
    ## 5: 17:21949792_C_A  17 21949792  a  c 0.4570 -0.0384 0.20690 -0.0043 0.8065
    ## 6: 17:21955349_T_G  17 21955349  t  g 0.5253  0.0505 0.09562  0.0103 0.5574
    ##     BETApp    Ppp    minP minTRAIT BETAmean
    ## 1: -0.0157 0.4690 0.41300      SBP -0.01875
    ## 2: -0.0160 0.4615 0.39270      SBP -0.01975
    ## 3: -0.0155 0.4749 0.41210      SBP -0.01895
    ## 4:  0.0185 0.7007 0.61460      SBP  0.01870
    ## 5: -0.0212 0.3050 0.20690      SBP -0.02135
    ## 6:  0.0248 0.2303 0.09562      SBP  0.03040

### Identifying insignificant variants

Identifying variants not present in the BP associated variants file:

``` r
full_gwas$CP <- paste(full_gwas$chr, full_gwas$bpos, sep = ":")
fwrite(full_gwas, "whole-gwas-variants.txt")

CP_matches <- cp_bp[CP %in% full_gwas$CP]

#Find insignificant variants (variants/rows not present in significant variant file)
mismatch_extract <- subset(full_gwas, !(CP %in% cp_bp$CP))
fwrite(mismatch_extract, "mismatched-genes.csv")
```

Filtering by p-value \>0.15 and by variants not within 500kb+/- of blood
pressure loci:

``` r
#Filtering again by only variants not 500kb +/- sentinel SNPs:

#Idenitfy lead SNPs in significant variants file
sentinels <- filter(cp_bp, Type == "Lead")

#Select relevant columns and rename to match insignificant variants column names
sentinels <- select(sentinels, CP, CHR, POS)
colnames(sentinels)[2] <- "chrom"
colnames(sentinels)[3] <- "position"
colnames(mismatch_extract)[2] <- "chrom"
colnames(mismatch_extract)[3] <- "position"
sentinels$chrom <-  as.numeric(sentinels$chrom)
sentinels <- data.table(sentinels)
mismatch_extract$chrom <-  as.numeric(mismatch_extract$chrom)
mismatch_extract <- data.table(mismatch_extract)

#Create boundaries 500kb +/- from lead SNPs
sentinels[, c("low", "high") := .(position - 500000, position  + 500000)]

#Find matches on chrom, with position between low&high for insignificant variants
mismatch_extract[sentinels, match := i.CP,
     on = .(chrom, position > low, position < high)]

#Discard all found matches and then drop the match-column
df <- mismatch_extract[is.na(match)][, match := NULL][]

#Output is only the insignificant variants not within 500kb+/- lead snps
fwrite(df, "genes-filtered_featuredist5kb.csv")

#Filtering with Giri Variants:
#Repeats the same 500kb+/- filtering code as above for further sentinel SNP data

dt2 <- fread("GiriGenes.csv")

#Create boundaries 500kb +/- from lead SNPs
dt2[, c("low", "high") := .(position - 500000, position  + 500000)]

#Find matches on chromosome, with position between low&high for insignificant variants
df[dt2, match := i.CP, on = .(chrom, position > low, position < high)]
#Discard all found matches and then drop the match-column
df <- df[is.na(match)][, match := NULL][]

#Merge two filtering insignificant variant files by all columns, so merging other columns back in:

dt <- merge(mismatch_extract, df)

#Create End column for file to suit bedtools format
dt$End <- dt$position
fwrite(dt, "gwas_insignificant-variants_FD5kb.txt")
```

Filtering by removing any insignificant variants in linkage
disequilibrum with associated blood pressure variants:

``` r
#Check if any insignificant variants have R2 - and remove:
insign <- fread("gwas_insignificant-variants_FD5kb.txt")
ld_1 <- fread("BP_Apr2020_ld_r2-1_500KB.txt")#all SNPs in ld with BP lead SNPs in ld 0.1 or greater r2
ld_1$CP1 <- paste(ld_1$CHR_B, ld_1$BP_B, sep = ":")
ld_1$CP2 <- paste(ld_1$CHR_A, ld_1$BP_A, sep = ":")
ld_2 <- select(ld_1, CP1, R2)
colnames(ld_2)[1] <- "CP"
ld_3 <- select(ld_1, CP2, R2)
colnames(ld_3)[1] <- "CP"
ld <- rbind(ld_2, ld_3)

CP_matches <- insign[CP %in% ld$CP]

final_insignif <- subset(insign, !(CP %in% ld$CP))

#Filtered file that next enters PPI filtering after gene annotation
fwrite(final_insignif, "gwas_insignificant-variants-filtered.txt")
```

Writing file to enter unix annotation in bedtools and ANNOVAR:

``` r
#Merge new insignificant variants with BP loci file

#Renaming columns to be accepted by bedotools and annovar
dt <- select(dt, "chrom", "position", "End", "a1", "a2", "BETAsbp", "BETAdbp", "BETApp",  "minP", "minTRAIT", "BETAmean", "CP")
colnames(dt)[3] <- "End"
colnames(dt)[4] <- "Ref"
colnames(dt)[5] <- "Alt"

#Also renaming KnownBP_Apr2020_LDr2-8_500kb_gwasRes.txt (significant variants/BP loci file) to match
bploci <- cp_bp
colnames(bploci)[14] <- "Ref"
colnames(bploci)[15] <- "Alt"

#Creating new chr and start columns from CP
dat <- setDT(bploci)[, tstrsplit(CP, "[:]", type.convert = TRUE)]
colnames(dat)[1] <- "Chr"
colnames(dat)[2] <- "Start"
dat$End <- dat$Start
#Adding chr and start columns back into dataset
bploci <- cbind(bploci, dat)

#Selecting only the columns needed for downstream processing
bpdf <- select(bploci, "Chr", "Start", "End",  "Ref", "Alt", "BETAsbp", "BETAdbp", "BETApp", "minP", "minTRAIT", "BETAmean", "CP")

#Combining with filtering insignificant variants data
df$End <- df$position
df <- select(df, "chrom", "position", "End",  "a1", "a2", "BETAsbp", "BETAdbp", "BETApp", "minP", "minTRAIT", "BETAmean", "CP")
colnames(df)[1] <- "Chr"
colnames(df)[2] <- "Start"
colnames(df)[4] <- "Ref"
colnames(df)[5] <- "Alt"

#Combining the variants together into one dataset
totaldata <- rbind(bpdf, df)

#Arranging Chr to be ordered 1-22 to allow data to be accepted by bedtools:
totaldata <- totaldata %>% arrange(Chr, Start)

#Needs to be tab delimited format for bedtools:
write.table(totaldata, file = "associated-and-insignificant-variants.txt", append = TRUE, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
#totaldata file enters unix scripts for gene annotation

head(totaldata)
```

    ##    Chr  Start    End Ref Alt BETAsbp BETAdbp  BETApp      minP minTRAIT
    ## 1:   1 752566 752566   a   g  0.0687  0.0605  0.0177 0.0171400      DBP
    ## 2:   1 885689 885689   a   g  0.2350  0.1369  0.0807 0.0009269      DBP
    ## 3:   1 885699 885699   a   g -0.2403 -0.1341 -0.0863 0.0009015      SBP
    ## 4:   1 886006 886006   t   c -0.2408 -0.1346 -0.0871 0.0008683      SBP
    ## 5:   1 887801 887801   a   g -0.2382 -0.1306 -0.0874 0.0009413      SBP
    ## 6:   1 888639 888639   t   c -0.2472 -0.1330 -0.0933 0.0006506      SBP
    ##    BETAmean       CP
    ## 1:  0.06460 1:752566
    ## 2:  0.18595 1:885689
    ## 3: -0.18720 1:885699
    ## 4: -0.18770 1:886006
    ## 5: -0.18440 1:887801
    ## 6: -0.19010 1:888639

## Data pre-processing following bedtools and ANNOVAR annotation in bash scripts

## <a id="2. GWAS variant p-value filtering"></a>2. GWAS variant p-value filtering

-   ANNOVAR variant feature filtering and selecting insigifnicant genes
    for PPI filtering:

``` r
df <- fread("Total_Genes_Annotated.ANN.hg19_multianno.txt", fill=TRUE)
#Annovar puts its header into 2nd row, needing rows 1 and 2 combined
tmp <- row_to_names(df, row_number = 1)
tmp2 <- row_to_names(tmp[,200:213], row_number = 1)
tmp <- tmp[-1, ]
df <- cbind(tmp[,1:199], tmp2)#


#Select ANNOVAR features from meta-predictors and features that do not use overlapping data:
df <- select(df, Gene, CP, Feature.Source, Feature.Chr, Feature.Start,Feature.End, 
            Feature.Type, Feature.Source, Feature.Dist, minP, BETAdbp, BETAsbp,
            BETApp, BETAmean, minTRAIT, REVEL, MetaSVM_rankscore, MetaLR_rankscore, MCAP, wgEncodeBroadHmmHuvecHMM, EFIN_SP_score)

df <- filter(df, !is.na(df$Gene))

#Write to be used in feature processing of epigenetic features (needs both variants and genes)
fwrite(df, "Variants_GenesAnnotated.txt", sep = "\t", row.names = FALSE)
```

``` r
df <- fread("Variants_GenesAnnotated.txt")

insignif <- fread("gwas_insignificant-variants-filtered.txt")

CP_matches <- df[CP %in% insignif$CP]

cp_bp <- fread("KnownBP_Apr2020_LDr2-8_500kb_gwasRes.txt", fill = TRUE)
sentinels <- filter(cp_bp, Type == "Lead")

#Select relevant columns and rename to match insignificant variants column names
sentinels <- select(sentinels, CP, CHR, POS)
colnames(sentinels)[2] <- "chr"
colnames(CP_matches)[4] <- "chr"
sentinels$chrom <-  as.numeric(sentinels$chrom)
sentinels <- data.table(sentinels)
CP_matches$chrom <-  as.numeric(CP_matches$chrom)
CP_matches <- data.table(CP_matches)

#Create boundaries 500kb +/- from lead SNPs
sentinels[, c("low", "high") := .(POS - 500000, POS  + 500000)]

#Find matches on chrom, with position between low&high for insignificant variants

gr1 <- makeGRangesFromDataFrame(
  data.frame(
    chr=sentinels$chr,
    start=sentinels$low,
    end=sentinels$high),
  keep.extra.columns=TRUE)

gr2 <- makeGRangesFromDataFrame(
  data.frame(
    chr=CP_matches$chr,
    start=CP_matches$Feature.Start,
    end=CP_matches$Feature.End,
    Gene = CP_matches$Gene),
  keep.extra.columns=TRUE)

no_overlaps <- gr2[-queryHits(findOverlaps(gr2, gr1, type="any")),] 

no_overlap_genes <- unique(no_overlaps$Gene)

matches <- df[Gene %in% no_overlap_genes]

#Condense variants to genes:
data <- setDT(matches)[, lapply(.SD, paste, collapse = ", "), by = Gene]

data <- select(data, Gene, minP, CP)

#Get min p-value per gene:
get_range <- function(x) {
  x <- type.convert(str_split(x, ",\\s+", simplify = TRUE), na.strings = c(".", "NA"))
  x <- t(apply(x, 1L, function(i) {
    i <- i[!is.na(i)]
    if (length(i) < 1L) c(NA_real_, NA_real_) else range(i)
  }))
  dimnames(x)[[2L]] <- c("min", "max")
  x
}

dt <- data[, c(Gene = .(Gene), lapply(.SD, get_range)), .SDcols = -"Gene"]

dt <- select(dt, Gene, minP.min)

dt2 <- dt

pval = dt2[dt2$minP.min > 0.0000005]

pval01 = dt2[dt2$minP.min > 0.1] 
pval015 = dt2[dt2$minP.min > 0.15] 
pval02 = dt2[dt2$minP.min > 0.2] 
pval025 = dt2[dt2$minP.min > 0.25] 
pval03 = dt2[dt2$minP.min > 0.3] 
pval035 = dt2[dt2$minP.min > 0.35] 
pval04 = dt2[dt2$minP.min > 0.4] 
pval045 = dt2[dt2$minP.min > 0.45] 
pval05 = dt2[dt2$minP.min > 0.5] 
pval055 = dt2[dt2$minP.min > 0.55] 
pval06 = dt2[dt2$minP.min > 0.6] 
pval065 = dt2[dt2$minP.min > 0.65] 
pval07 = dt2[dt2$minP.min > 0.7] 
pval075 = dt2[dt2$minP.min > 0.75] 
pval08 = dt2[dt2$minP.min > 0.8] 
pval085 = dt2[dt2$minP.min > 0.85] 
pval09 = dt2[dt2$minP.min > 0.9] 

fwrite(pval, 'pval_forPPI.txt', row.names = FALSE)
fwrite(pval01, 'pval01_forPPI.txt', row.names = FALSE)
fwrite(pval015, 'pval015_forPPI.txt', row.names = FALSE)
fwrite(pval02, 'pval02_forPPI.txt', row.names = FALSE)
fwrite(pval025, 'pval025_forPPI.txt', row.names = FALSE)
fwrite(pval03, 'pval03_forPPI.txt', row.names = FALSE)
fwrite(pval035, 'pval035_forPPI.txt', row.names = FALSE)
fwrite(pval04, 'pval04_forPPI.txt', row.names = FALSE)
fwrite(pval045, 'pval045_forPPI.txt', row.names = FALSE)
fwrite(pval05, 'pval05_forPPI.txt', row.names = FALSE)
fwrite(pval055, 'pval055_forPPI.txt', row.names = FALSE)
fwrite(pval06, 'pval06_forPPI.txt', row.names = FALSE)
fwrite(pval065, 'pval065_forPPI.txt', row.names = FALSE)
fwrite(pval07, 'pval07_forPPI.txt', row.names = FALSE)
fwrite(pval075, 'pval075_forPPI.txt', row.names = FALSE)
fwrite(pval08, 'pval08_forPPI.txt', row.names = FALSE)
fwrite(pval085, 'pval085_forPPI.txt', row.names = FALSE)
fwrite(pval09, 'pval09_forPPI.txt', row.names = FALSE)
```

## <a id="3. GWAS variant data filtering"></a>3. GWAS variant data filtering

-   ANNOVAR variant feature filtering and selecting insigifnicant genes
    for PPI filtering:

``` r
df <- fread("Variants_GenesAnnotated.txt")

#Identify insignificant genes from variants
insignif <-  fread('pval015_forPPI.txt')
insignif$insignif <- 1
insign_df <-  merge(df, insignif, by = 'Gene', all.y = TRUE)
insign_df <- filter(insign_df, !is.na(insign_df$Gene))
insign_df_genes <- setDT(insign_df)[, lapply(.SD, paste, collapse = ", "), by = Gene]
insign_df_genes  <- insign_df_genes[!duplicated(insign_df_genes$Gene), ]
#insignif column then used in further filtering by PPI later

#Identify Evangelou et al. genes from whole GWAS variants
evangelou <- fread('WholeDataFDLD08.csv')
colnames(evangelou)[149] <-'Gene'
evangelou  <- select(evangelou, CP, Gene, minP, BETAsbp, BETAdbp, BETApp, minTRAIT, Feature.Type)

BPlocigenes <- df[Gene %in% evangelou$Gene]

#Compress per variant-level features to their genes (variant values become comma separated values within one cell/row per gene)

data <- setDT(BPlocigenes)[, lapply(.SD, paste, collapse = ", "), by = Gene]

#Combine with insignificant genes
data <- rbind(data, insign_df_genes, fill=TRUE)
fwrite(data, 'Genes_evangel_and_insignif.txt', row.names = FALSE)

#Identify gene list to undergo PPI filtering in later code
dataforPPI <- select(data, Gene, insignif)

get_range <- function(x) {
  x <- type.convert(str_split(x, ",\\s+", simplify = TRUE), na.strings = c(".", "NA"))
  x <- t(apply(x, 1L, function(i) {
    i <- i[!is.na(i)]
    if (length(i) < 1L) c(NA_real_, NA_real_) else range(i)
  }))
  dimnames(x)[[2L]] <- c("min", "max")
  x
}

dataforPPI <- dataforPPI[, c(Gene = .(Gene), lapply(.SD, get_range)), .SDcols = -"Gene"]
dataforPPI <- select(dataforPPI, Gene, insignif.max)
colnames(dataforPPI)[2] <- 'insignif'

fwrite(dataforPPI, 'Genelist_forPPI.txt', row.names = FALSE)
```

## <a id="4. Protein-protein interaction filtering"></a>4. Protein-protein interaction filtering

``` r
gc()
```

    ##             used   (Mb) gc trigger    (Mb) limit (Mb)   max used    (Mb)
    ## Ncells   3455240  184.6   37684940  2012.6         NA   47106174  2515.8
    ## Vcells 142920123 1090.4 3342379942 25500.4      65536 4177222844 31869.7

``` r
memory.limit(9999999999)
```

    ## [1] Inf

``` r
#Match genes to string identifiers

df <- fread("Variants_GenesAnnotated.txt")

#Condense variants to genes:
data <- setDT(df)[, lapply(.SD, paste, collapse = ", "), by = Gene]

Genelist <- select(data, Gene)

protein_alias <- fread('9606.protein.aliases.v11.0.txt')
colnames(protein_alias)[1] <- 'protein_id'
colnames(protein_alias)[2] <- 'Gene'
protein_alias <- select(protein_alias, Gene, protein_id)
protein_info <- fread('9606.protein.info.v11.0.txt')
colnames(protein_info)[1] <- 'protein_id'
colnames(protein_info)[2] <- 'Gene'
protein_info <- select(protein_info, Gene, protein_id)
PPIs <- rbind(protein_alias, protein_info)


#Total data downloaded from string including all scores per interactions loaded:

proteindf <- fread("9606.protein.links.full.v11.0.txt")
dim(proteindf)
```

    ## [1] 11759454       16

``` r
names(proteindf)
```

    ##  [1] "protein1"                 "protein2"                
    ##  [3] "neighborhood"             "neighborhood_transferred"
    ##  [5] "fusion"                   "cooccurence"             
    ##  [7] "homology"                 "coexpression"            
    ##  [9] "coexpression_transferred" "experiments"             
    ## [11] "experiments_transferred"  "database"                
    ## [13] "database_transferred"     "textmining"              
    ## [15] "textmining_transferred"   "combined_score"

``` r
colnames(proteindf)[1] <- "STRING_id1"
colnames(proteindf)[2] <- "STRING_id2"

#Interaction scores only from coexpression, experiments and database used:
protein_df <- select(proteindf, STRING_id1, STRING_id2, 
                     experiments)
protein_df <- filter(protein_df, experiments != 0)
protein_df$Max <- apply(protein_df[,3], 1, FUN = max)
protein_df <- filter(protein_df, Max >= 150) 

protein_df <- select(protein_df, STRING_id1, STRING_id2)
```

``` r
#Matching Gene and string identifiers to interactors 

#Duplicating STRING_id to have renamed matching column STRING_id1 
#(to match with named STRING_id1 interactors in protdf dataset, finding the gene names for those interactors)
PPIs$STRING_id1 <- PPIs$protein_id
PPIs$STRING_id2 <- PPIs$protein_id
#Having both STRING_id1 and STRING_id means they will correspond to the same columns in protdf respectively

mapped1 <- select(PPIs, STRING_id1, Gene)
colnames(mapped1)[2]<- "Gene1"
mapped1 <- unique(mapped1)
string1_df <- merge(protein_df, mapped1, by='STRING_id1', all.x = TRUE, allow.cartesian=TRUE)
string1_df  <- unique(string1_df)
colnames(Genelist)[1] <- 'Gene1'
string1_df_genes <- merge(string1_df, Genelist, by='Gene1', all.y = TRUE, allow.cartesian=TRUE)
string1_df_genes <- unique(string1_df_genes)

head(string1_df_genes)
```

    ##    Gene1           STRING_id1           STRING_id2
    ## 1:   7SK                 <NA>                 <NA>
    ## 2:  A1BG 9606.ENSP00000263100 9606.ENSP00000351697
    ## 3:  A1BG 9606.ENSP00000263100 9606.ENSP00000284116
    ## 4:  A1BG 9606.ENSP00000263100 9606.ENSP00000294119
    ## 5:  A1BG 9606.ENSP00000263100 9606.ENSP00000384792
    ## 6:  A1BG 9606.ENSP00000263100 9606.ENSP00000416515

``` r
#Finding HGNC gene symbols for other STRING id column:
mapped2 <- select(PPIs, STRING_id2, Gene)
mapped2 <- unique(mapped2)
string2_df <- merge(protein_df, mapped2, by='STRING_id2', all.x =TRUE, allow.cartesian=TRUE)
string2_df <- unique(string2_df)
colnames(Genelist)[1] <- 'Gene'
string2_df_genes <- merge(string2_df, Genelist, by='Gene', all.y = TRUE, allow.cartesian=TRUE)
string2_df_genes <- unique(string2_df_genes)

all_string_IDs <- inner_join(string1_df_genes, string2_df_genes)
all_string_IDs <- select(all_string_IDs, STRING_id1, Gene1, STRING_id2, Gene)

colnames(all_string_IDs)[4] <- 'Gene2'
all_string_IDs <- filter(all_string_IDs, !is.na(STRING_id1))

fwrite(all_string_IDs, 'STRING_IDs_Experimental_150.txt', row.names = FALSE)

#Compress gene's interactors for first column:
genes1 <- all_string_IDs[, c(2,4)]
genes1 <- genes1[!(genes1$Gene1=="NA"), ]
genes1 <- data.table(genes1)
genes1 <- genes1[,lapply(.SD, function(col) paste(col, collapse =", ")), 
                 by=.(Gene1)]
genes1$Gene2 <- gsub("NA,", "", genes1$Gene2)


head(genes1)
```

    ##     Gene1
    ## 1:   A1BG
    ## 2:   A1CF
    ## 3:    A2M
    ## 4:  A2ML1
    ## 5: A4GALT
    ## 6:  A4GNT
    ##                                                                                                                                                                                                                                                                                  Gene2
    ## 1:                                                                                                                                                                                         REV3L, GDPD1, UBXN1, WDR62, FDXR, CRISP3, KCMF1, PSME4, ZBTB40, RNF123, UBR4, GLIPR1, UBAC1
    ## 2:                                                                                                                                                                                                             SDPR, FBP2, KHSRP, PSMG2, EIF4G1, PSMG3, CELF2, APOBEC1, TNPO2, SYNCRIP
    ## 3: GDPD1, KLK13, AMBP, ADAMTS7, IL1B, PDGFB, IL10, CDK7, ADAM19, MMP2, ADAMTS12, LCAT, TGIF1, PZP, APOE, HAPLN3, ANXA6, NGF, LEP, ADAMTS1, HSPA5, MYOC, KLK2, PTX3, MAPK6, LRP1, IL4, TGFBI, PRADC1, METRN, CELA1, SPACA3, HIPK1, KLK3, LYZ, APP, PAEP, B2M, METRNL, CTSB, CPB2, LOXL1
    ## 4:                                                                                                                                                                                                                                          L2HGDH, ECH1, THRA, ZFAND4, FCRL5, UGT1A10
    ## 5:                                                                                                                                                                                                                                                ABHD17B, ABHD14A, ICK, PTGFRN, FCGRT
    ## 6:                                                                                                                                                                                                                     ICK, FAM3C, IMPA2, ARF5, ABHD14A, FAM213A, JMJD8, RNF167, FANK1

``` r
#Compress gene's interactors for second column:
genes2 <- all_string_IDs[, c(2,4)]
genes2 <- genes2[!(genes2$Gene2=="NA"),]
genes2 <- data.table(genes2)
genes2 <- genes2[,lapply(.SD, function(col) paste(col, collapse = ", ")), 
                 by=.(Gene2)]
genes2$Gene1 <- gsub("NA,", "", genes2$Gene1)

#Combine into one dataset:
colnames(genes1)[1] <- "Gene"
colnames(genes2)[1] <- "Gene"

interactors_df <- merge(genes1, genes2, by='Gene', all=T)
#Combine interactors from both gene1 and gene2 and remove duplicates:
PPI_df  <- transform(interactors_df , interactors = paste(Gene1, Gene2, sep = ", "))
PPI_df$interactors <- sapply(PPI_df$interactors , function(x) paste(unique(unlist(str_split(x,", "))), collapse = ", "))

interactors <- select(PPI_df, Gene, interactors)

fwrite(interactors, "150_Experimental_interactions.txt", sep="\t", row.names = F, quote = F)
```

### PPI filtering of insignificant genes by identifying genes with no close interactions with known blood pressure genes

-   Known blood pressure genes identified as such if they have
    interacting blood pressure drugs or mechanisms (provided in
    Genelist_drug_and_databases.txt file).
-   Genes with no direct PPIs with known blood pressures or no secondary
    PPIs (interactors interacting with known blood pressure gene
    interactors) selected as final geneset of least likely/insignificant
    blood pressure genes to be scored at 0.1 on machine learning
    regression classification.

``` r
rm(list = ls())


Gene_list <- fread("Genelist_forPPI.txt")
gene_types <- fread('hg19Rel92_AllgeneTypes_0kb.txt')
colnames(gene_types)[c(4,5)] <- c('Gene', 'Type')
gene_types <- select(gene_types, Gene, Type)
Genelist <- merge(Gene_list, gene_types, by='Gene', all.x=TRUE)
Genelist <- filter(Genelist, Type == 'protein_coding')
LL_count <- filter(Genelist, !is.na(insignif))


interactors <- fread("150_Experimental_interactions.txt")
drugdf <- fread("Genelist_drug_and_databases.txt")
textminingdf <- fread('genie_textminingBP.txt')
colnames(textminingdf)[3] <-'Gene'
genes_to_label <- Reduce(function(x, y) merge(x, y, by = 'Gene', all.x = TRUE),
                         list(Genelist, drugdf, textminingdf))
genes_to_label <- select(genes_to_label, Gene, BPmed, Mechanism,
                         SideeffectFreq, Rank, insignif)

genes_to_label[, BP_templabel := fcase( 
  Mechanism == 1 | BPmed >= 1 | SideeffectFreq > 1 | Rank >= 1, 'most likely', 
  insignif == 1 , 'least likely', 
  default = 'other')]

bp_labels <- select(genes_to_label, Gene, BP_templabel)
#Genes with BP drug or BP mechanism (known BP genes)
ml_df <- filter(bp_labels, BP_templabel == "most likely")
bp_genes <- fread('Genes_BPloci.bed')
colnames(bp_genes)[12] <-'CP'
colnames(bp_genes)[16] <-'Gene'
bp_genes$CP <- str_replace_all(bp_genes$CP, ',', ':')
bp_total_loci <- fread('KnownBP_Apr2020_LDr2-8_500kb_gwasRes.txt')
bp_genes2 <- merge(bp_genes, bp_total_loci, by='CP', all.x=TRUE)
bp_genes2 <- filter(bp_genes2, !is.na(Type))
bp_genes3 <- setDT(bp_genes2)[, lapply(.SD, paste, collapse = ", "), by = Gene]
bp_genes3 <- select(bp_genes3, Gene)
bp_genes3 <- select(bp_genes, Gene)

ml_df <- rbind(ml_df, bp_genes3, fill=TRUE)

ml_df[is.na(ml_df)] <- "most likely"
ml_df <- ml_df[!duplicated(ml_df$Gene), ] #3,786 genes

#Insignificant genes from GWAS (p-value, no LD, no within 500kb+/- bploci)
ll_df <- filter(bp_labels, BP_templabel == "least likely") #10,471 genes

bp_loci <- ml_df 
bp_genes_PPI <- merge(bp_loci, interactors, by = 'Gene', all.x = TRUE)

#Get secondary interactors:
interactors_bp_PPI <- bp_genes_PPI %>% 
  mutate(interactors = strsplit(as.character(interactors), ", ")) %>% 
  unnest(interactors)
interactors_bp_PPI <- select(interactors_bp_PPI, interactors)
colnames(interactors_bp_PPI)[1] <- 'Gene'
interactors_bp_PPI <- unique(interactors_bp_PPI)

interactors_bp_PPI <- interactors_bp_PPI %>%
  mutate_if(is.character, trimws)

secondary_bp_PPI <- merge(interactors_bp_PPI, interactors, by = 'Gene', all.x = TRUE)

#Directly connected PPIs with known BP genes (matched_direct_interactor column)
#Least likely genes that interact with a known BP gene
ll_df$matched_direct_interactor <-  ll_df$Gene %in% unlist(strsplit
                                                           (bp_genes_PPI$interactors,
                                                             ", "))

ll_df$matched_secondary_interactor <- ll_df$Gene %in% unlist(strsplit
                                                             (secondary_bp_PPI$interactors,
                                                               ", "))

string_ids <- fread('9606.protein.info.v11.0.txt')
colnames(string_ids)[2] <- 'Gene'
ll_matches <- ll_df[Gene %in% string_ids$Gene]

#Select insignificant genes with no direct or secondary BP gene interactors

exp150_false_interactions1 <- filter(ll_df, matched_direct_interactor == "FALSE")

exp150_false_interactions2 <- filter(ll_df, matched_direct_interactor == "FALSE"
                                 & matched_secondary_interactor == "FALSE") 

bp_true_interactions1 <- filter(ll_df, matched_direct_interactor == "TRUE") 
bp_true_interactions2 <- filter(ll_df, matched_secondary_interactor == "TRUE") 


#read in genes with >0.1 r2 with BP loci to ensure genes are filtered:
ld_genes <- fread('LD01_genes.txt')
exp150_no_ld_genes <- subset(exp150_false_interactions2, !(Gene %in% ld_genes$Gene))

write.table(exp150_no_ld_genes,"150_Experimental_pval015_filtered_unassociated_genes.txt",sep="\t", row.names = F, quote = F)
```

## <a id="5. Variant-level feature processing"></a>5. Variant-level feature processing

-   Min or max value per gene for each variant feature selected with
    min/max depending on feature significance direction

``` r
df <- fread('Genes_evangel_and_insignif.txt')
ppifiltered1 <- fread('150_Experimental_pval015_filtered_unassociated_genes.txt', select=c('Gene'))
ppifiltered1$insignif <- 1

df$insignif[df$insignif==""] <- NA
df1 <- filter(df, is.na(insignif))

data <- rbind(df1, ppifiltered1, fill=TRUE)
data <- data[!duplicated(data$Gene), ]

#Select min or max value for each feature per gene between compressed values

get_range <- function(x) {
  x <- type.convert(str_split(x, ",\\s+", simplify = TRUE), na.strings = c(".", "NA"))
  x <- t(apply(x, 1L, function(i) {
    i <- i[!is.na(i)]
    if (length(i) < 1L) c(NA_real_, NA_real_) else range(i)
  }))
  dimnames(x)[[2L]] <- c("min", "max")
  x
}

dt <- data[, c(Gene = .(Gene), lapply(.SD, get_range)), .SDcols = -"Gene"]


dt <- subset(dt, select=-c(wgEncodeBroadHmmHuvecHMM.min, wgEncodeBroadHmmHuvecHMM.max,
                             MetaSVM_rankscore.min, MetaLR_rankscore.min,
                             MCAP.min, REVEL.min, EFIN_SP_score.max, minP.max,
                            minTRAIT.min, minTRAIT.max, BETAmean.min,BETAmean.max,
                           BETAdbp.min, BETAdbp.max, BETApp.min, BETApp.max, BETAsbp.min, BETAsbp.max,
                           CP.min,CP.max, minP.max, minP.min, Feature.Dist.max,
                           Feature.Dist.min, Feature.Type.max, Feature.Type.min,
                           Feature.End.min, Feature.Start.max, Feature.Start.min,
                           Feature.Chr.max, Feature.Source.min, Feature.End.max,
                           Feature.Source.max, Feature.Chr.min, insignif.min, insignif.max,
                           CP.min, CP.max
                           ))

dt2 <- select(data, Gene, wgEncodeBroadHmmHuvecHMM)

dat <- select(data, Gene, BETAsbp, BETAdbp, BETApp, BETAmean)
dat <- dat[!duplicated(dat$Gene),]

#Columns most columns need maximum value per gene select.
#Columns in between and not selected either need text removal first or no selection performed (minTRAIT column)

#For columns with text, remove text or convert to numeric values to select between values per gene/row

dt2$wgEncodeBroadHmmHuvecHMM <- gsub("[A-Z]|[a-z]|[=]|\\s|_|/", "", dt2$wgEncodeBroadHmmHuvecHMM)

dt2 <- dt2 %>% 
  tidyr::separate_rows(wgEncodeBroadHmmHuvecHMM, sep = ",", convert = TRUE) %>%
  dplyr::group_by(Gene) %>% 
  dplyr::mutate(wgEncodeBroadHmmHuvecHMM.count = n()) %>% 
  dplyr::slice(1) %>% 
  ungroup(Gene)

dt2_wg <- select(dt2, wgEncodeBroadHmmHuvecHMM.count)
dt1 <- cbind(dt, dt2_wg)


#Select the maximum absolute value between all beta phenotypes per gene:

#Select maximum absolute beta value per gene (between sbp, dbp, and pp)
#Selecting one phenotype at a time then combining columns:
#betadbp:
dat2 <- select(dat, Gene, BETAdbp)%>%
  tidyr::separate_rows(BETAdbp, sep = ",", convert = TRUE)

dat2$BETAdbp <- as.numeric(dat2$BETAdbp)
dat2 <- dat2 %>% 
  group_by(Gene) %>%
  summarise(Negatives1 = +(min(BETAdbp, na.rm = TRUE) == -max(abs(BETAdbp), na.rm = TRUE)),
            BETAdbp = max(abs(BETAdbp), na.rm = TRUE))

#Replace infinite values with NA
dat2 <- do.call(data.frame,lapply(dat2, function(x) replace(x, is.infinite(x),NA)))

#Repeat for BETApp, BETAsbp, and BETAmean
#BETApp:
dat3 <- select(dat, Gene, BETApp)%>%
  tidyr::separate_rows(BETApp, sep = ",", convert = TRUE)

dat3$BETApp <- as.numeric(dat3$BETApp)
dat3 <- dat3 %>% 
  group_by(Gene) %>%
  summarise(Negatives2 = +(min(BETApp, na.rm = TRUE) == -max(abs(BETApp), na.rm = TRUE)),
            BETApp = max(abs(BETApp), na.rm = TRUE))

dat3 <- do.call(data.frame,lapply(dat3, function(x) replace(x, is.infinite(x),NA)))

#BETAsbp:
dat4 <- select(dat, Gene, BETAsbp)%>%
  tidyr::separate_rows(BETAsbp, sep = ",", convert = TRUE)

dat4$BETAsbp <- as.numeric(dat4$BETAsbp)
dat4 <- dat4 %>% 
  group_by(Gene) %>%
  summarise(Negatives3 = +(min(BETAsbp, na.rm = TRUE) == -max(abs(BETAsbp), na.rm = TRUE)),
            BETAsbp = max(abs(BETAsbp), na.rm = TRUE))

dat4 <- do.call(data.frame,lapply(dat4, function(x) replace(x, is.infinite(x),NA)))

#Add created beta columns into df datframe and then removing duplicate named columns:
df <- cbind(dat, dat2[,2:3], dat3[,2:3], dat4[,2:3])
df <- subset(df, select = -c(BETAdbp, BETAsbp,BETApp, BETAmean))

#Return negative beta values to those which were negative per phenotype
df <- df %>%
  mutate(BETAdbp = ifelse(Negatives1, -BETAdbp, BETAdbp))
df <- df %>%
  mutate(BETApp = ifelse(Negatives2, -BETApp, BETApp))
df <- df %>%
  mutate(BETAsbp = ifelse(Negatives3, -BETAsbp, BETAsbp))

#Selecting only the 3 phenotype beta value columns for selecting the maximum absolute value between them:
df1sel <- df[,c(3, 5, 7)]
dat_df <- tibble::as_tibble(df1sel)

#Function to select maximum absolute value rowwise out of the 3 phenotypes:
f <- function(data) {
  tmp <- Filter(is.numeric, data)
  if(inherits(data, "tbl_df")) {
    tmp <- as.matrix(tmp)
  }
  tmp[cbind(1:nrow(tmp),
            max.col(replace(x <- abs(tmp), is.na(x), -Inf)))]
}

#Performing the f function to select the maximum absolute value into a new 'betamax' column
df1_beta_selected2 <- f(dat_df)
df1_beta_selected2 <- as.list(df1_beta_selected2)
df_beta <- data.frame(matrix(unlist(df1_beta_selected2), nrow=length(df1_beta_selected2), byrow = T))
colnames(df_beta)[1] <- 'betamax'
#Adding betamax column into main dataset:
df_new <- cbind(dt1, df_beta)
df_new <- unique(df_new)
#Drop any columns missing all values
df_all_features <- df_new %>% select_if(~!all(is.na(.)))

insignif <- select(data, Gene, insignif)
df_all_features <- merge(df_all_features, insignif, all.x = TRUE)
#Write final file:
fwrite(df_all_features, "pval015_genelist_filtered_evangelou.txt", 
       sep = "\t", row.names = FALSE)

head(df_all_features)
```

    ##       Gene REVEL.max MetaSVM_rankscore.max MetaLR_rankscore.max MCAP.max
    ## 1:   AAMDC        NA                    NA                   NA       NA
    ## 2:   ABCB9        NA                    NA                   NA       NA
    ## 3:   ABCC8        NA                    NA                   NA       NA
    ## 4:   ABCC9        NA                    NA                   NA       NA
    ## 5: ABHD16A        NA                    NA                   NA       NA
    ## 6: ABHD17C        NA                    NA                   NA       NA
    ##    wgEncodeBroadHmmHuvecHMM.count betamax insignif
    ## 1:                            250  0.7121     <NA>
    ## 2:                              7 -0.4107     <NA>
    ## 3:                            280 -0.4145     <NA>
    ## 4:                             56  0.4561     <NA>
    ## 5:                            671 -0.3551     <NA>
    ## 6:                             52  0.4113     <NA>

## <a id="6. Feature pre-procesing"></a>6. Feature pre-procesing

### IPA blood pressure gene curation:

``` r
df1 <- fread('CVD-BP-IPA.txt')
colnames(df1)[1] <-'Gene'
colnames(df1)[3] <-'Activity_BP_IPA'
df1 <- select(df1, Gene, Activity_BP_IPA)
df1 <- unique(df1)
df2 <- fread('familial-HTN-IPA.txt')
colnames(df2)[1] <-'Gene'
colnames(df2)[3] <-'Activity_FamilialBP_IPA'
df2 <- select(df2, Gene, Activity_FamilialBP_IPA)
df2 <- unique(df2)
df3 <- fread('DBP-IPA.txt')
colnames(df3)[1] <-'Gene'
colnames(df3)[3] <-'Activity_DBP_IPA'
df3 <- select(df3, Gene, Activity_DBP_IPA)
df3 <- unique(df3)
df4 <- fread('SBP-IPA.txt')
colnames(df4)[1] <-'Gene'
colnames(df4)[3] <-'Activity_SBP_IPA'
df4 <- select(df4, Gene, Activity_SBP_IPA)
df4 <- unique(df4)
df5 <- fread('BP-disorder-IPA.txt')
colnames(df5)[3] <-'Activity_BPdisorder_IPA'
colnames(df5)[1] <-'Gene'
df5 <- select(df5, Gene, Activity_BPdisorder_IPA)
df5 <- unique(df5)

df6 <- fread('IPA-Pressure.txt')
colnames(df6)[4] <-'Gene'
df6 <- select(df6, Gene)
df6 <- separate_rows(df6, Gene, sep = ",")
df6 <- unique(df6)
df6$IPA_Pressure <- 1

df7 <- fread('maternal-HTN-IPA.txt')
colnames(df7)[4] <-'Gene'
df7 <- select(df7, Gene)
df7 <- separate_rows(df7, Gene, sep = ",")
df7 <- unique(df7)
df7$IPA_Matneral <- 1

df8 <- fread('MAP-IPA.txt')
colnames(df8)[4] <-'Gene'
df8 <- select(df8, Gene)
df8 <- separate_rows(df8, Gene, sep = ",")
df8 <- unique(df8)
df8$MAP_IPA <- 1

df9 <- fread('PP-IPA.txt')
colnames(df9)[4] <-'Gene'
df9 <- select(df9, Gene)
df9 <- separate_rows(df9, Gene, sep = ",")
df9 <- unique(df9)
df9$PP_IPA <- 1


merged <- Reduce(function(x, y) merge(x, y,by='Gene', all=TRUE), list(df1,df2, df3, df4, df5, df6, df7, df8, df9))

### Getting gene's presence in IPA under BP terms
merged <- merged[,1]
merged$IPA_BP <- 1 

genes <- fread('Variants_GenesAnnotated.txt')
genes <- select(genes, Gene)
genes <- unique(genes)
df <- merge(genes, merged, by='Gene', all.x=TRUE)


#### Getting gene's activity in BP (increase or decrease)
levels( factor(df1$Activity_BP_IPA) )
```

    ##  [1] "[decreased activity, increased activity]"                        
    ##  [2] "[decreased activity, unknown change in activity]"                
    ##  [3] "[increased activity, unknown change in activity]"                
    ##  [4] "decreased activity"                                              
    ##  [5] "decreased activity,increased activity"                           
    ##  [6] "decreased activity,increased activity,unknown change in activity"
    ##  [7] "decreased activity,unknown change in activity"                   
    ##  [8] "increased activity"                                              
    ##  [9] "increased activity,unknown change in activity"                   
    ## [10] "unknown change in activity"

``` r
#"[decreased activity, increased activity]"                        
#[2] "[decreased activity, unknown change in activity]"                
#[3] "[increased activity, unknown change in activity]"                
#[4] "decreased activity"                                              
#[5] "decreased activity,increased activity"                           
#[6] "decreased activity,increased activity,unknown change in activity"
#[7] "decreased activity,unknown change in activity"                   
#[8] "increased activity"                                              
#[9] "increased activity,unknown change in activity"                   
#[10] "unknown change in activity"  

filt1 <-filter(df1, Activity_BP_IPA=='decreased activity,increased activity')
filt2 <-filter(df1, Activity_BP_IPA=='decreased activity,increased activity,unknown change in activity')
#filt3 <-filter(df1, Activity_BP_IPA=='unknown change in activity')

df1_2 <-subset(df1, !(Gene %in% filt1$Gene)) 
df1_2 <-subset(df1_2, !(Gene %in% filt2$Gene)) 
#df1_2 <-subset(df1_2, !(Gene %in% filt3$Gene)) 

levels( factor(df1_2$Activity_BP_IPA) )
```

    ## [1] "[decreased activity, increased activity]"        
    ## [2] "[decreased activity, unknown change in activity]"
    ## [3] "[increased activity, unknown change in activity]"
    ## [4] "decreased activity"                              
    ## [5] "decreased activity,unknown change in activity"   
    ## [6] "increased activity"                              
    ## [7] "increased activity,unknown change in activity"   
    ## [8] "unknown change in activity"

``` r
filter1 <- filter(df1_2, Activity_BP_IPA=='decreased activity')
filter1$IPA_Activity <- 1
filter2 <- filter(df1_2, Activity_BP_IPA=='increased activity')
filter2$IPA_Activity <- 2
df1_final <- rbind(filter1, filter2)

final <- merge(df, df1_final, by='Gene', all.x=TRUE)
final$IPA_BP[is.na(final$IPA_BP)] <- 0
final$IPA_Activity[is.na(final$IPA_Activity)] <- 0
final <- select(final, Gene, IPA_BP, IPA_Activity)
fwrite(final, 'IPA_activity.txt', row.names=FALSE)
```

### PanglaoDB cell-type ubiquitous index:

``` r
df <- fread('PanglaoDB_markers.txt')

get_range <- function(x) {
  x <- type.convert(str_split(x, ",\\s+", simplify = TRUE), na.strings = c(".", "NA"))
  x <- t(apply(x, 1L, function(i) {
    i <- i[!is.na(i)]
    if (length(i) < 1L) c(NA_real_, NA_real_) else range(i)
  }))
  dimnames(x)[[2L]] <- c("min", "max")
  x
}

dt <- df[, c(Gene = .(Gene), lapply(.SD, get_range)), .SDcols = -"Gene"]

dt <- select(dt, Gene, ubiquitousness.index.max)
dt <- unique(dt)

fwrite(dt, 'panglaodb_cells_highestscores.txt', row.names=FALSE)
```

### CpG Island counts per gene from UCSC data:

``` r
#Read in Variant file including gene names from GWAS variant filtering code:
file1 <- fread("Variants_GenesAnnotated.txt")
file1 <- select(file1, Gene, Feature.Chr, Feature.Start, Feature.End)
colnames(file1)[2] <- "Chromosome"
colnames(file1)[3] <- "Start"
colnames(file1)[4] <- "End"
file1$Position <- file1$Start

#Read in gene length form hg19 ref file (used in bedtools step):
#Gene length as originally going to be used for correlation but was removed to be used a ML feature instead
#genelength <- fread("GeneLength.txt")

file2<-fread("cpgIslandExt.txt")
colnames(file2)[2] <- "Chromosome"
colnames(file2)[3] <- "Start"
colnames(file2)[4] <- "End"
file2$count <- 1
file2$Chromosome <- gsub('chr', '', file2$Chromosome)
file1$Chromosome <- as.integer(file1$Chromosome)
file2$Chromosome <- as.integer(file2$Chromosome)

#Match CpG Chr Start and End to those in file1:
df <- file1[,.(Chromosome,Position,Start, End, Gene)][file2,
                                                .(Chromosome,Position=x.Position, Start, End,Gene,count),
                                                on=.(Chromosome,Position>=Start,Position<=End),nomatch=0L]
#Compress columns from variant to gene level:

CpGcount <- df %>%
  group_by(Gene) %>%
  summarize(text = str_c(count, collapse = ", "))

colnames(CpGcount)[1] = "Gene"
colnames(CpGcount)[2] = "CpGcount"

#Sum rows to get count per gene:
CpGcount$CpGcount <- sapply(CpGcount$CpGcount, function(row) sum(as.numeric(strsplit(row, ", ")[[1]])), USE.NAMES = FALSE)

CpGcount <- CpGcount[complete.cases(CpGcount), ]
CpGcount <- CpGcount[!duplicated(CpGcount$Gene), ]
fwrite(CpGcount, "Genes_CpG.txt") 

head(CpGcount)
```

    ## # A tibble: 6 × 2
    ##   Gene     CpGcount
    ##   <chr>       <dbl>
    ## 1 A1BG-AS1       12
    ## 2 A3GALT2        36
    ## 3 AACS          209
    ## 4 AAMDC         149
    ## 5 AAR2           34
    ## 6 AATF          124

### DNase cluster counts per gene from UCSC data:

``` r
file1 <- fread("Variants_GenesAnnotated.txt")
file1 <- select(file1, Gene, Feature.Chr, Feature.Start, Feature.End)
colnames(file1)[2] <- "Chromosome"
colnames(file1)[3] <- "Start"
colnames(file1)[4] <- "End"
file1$Position <- file1$Start

file2<-fread("DNaseClusters.txt")
file2$chrom <- gsub('chr', '', file2$chrom)
colnames(file2)[2] <- "Chr"
colnames(file2)[3] <- "Start"
colnames(file2)[4] <- "End"

file2$Chr <- as.integer(file2$Chr)
file1$Chr <- as.integer(file1$Chr)

file2$count <- 1

df <- file1[,.(Chr,Position,Start, End, Gene)][file2,
                                      .(Chr,Position=x.Position, Start, End,Gene, name, score, count),
                                      on=.(Chr,Position>=Start,Position<=End),nomatch=0L]


score <- df %>%
  group_by(Gene) %>%
  summarize(text = str_c(score, collapse = ", "))

cluster <- df %>%
  group_by(Gene) %>%
  summarize(text = str_c(count, collapse = ", "))

merged <- Reduce(function(x, y) merge(x, y, by='Gene', all=TRUE), list(cluster, score))

colnames(merged)[2] = "DNaseCluster_count"
colnames(merged)[3] = "DNase_highestscore"

merged$DNaseCluster_count <- sapply(merged$DNaseCluster_count, function(row) sum(as.numeric(strsplit(row, ", ")[[1]])), USE.NAMES = FALSE)


DNaseClusters <- select(merged, Gene, DNaseCluster_count)
DNaseClusters <- DNaseClusters[complete.cases(DNaseClusters), ]
fwrite(DNaseClusters, "DNase_Gene.txt") 

head(DNaseClusters)
```

    ##        Gene DNaseCluster_count
    ## 1       7SK                  2
    ## 2  A1BG-AS1                 12
    ## 3       A2M                122
    ## 4   A2M-AS1                 10
    ## 5     A2ML1                197
    ## 6 A2ML1-AS2                  1

### Enhancer site counts per gene from UCSC data:

``` r
file1 <- fread("Variants_GenesAnnotated.txt")
file1 <- select(file1, Gene, Feature.Chr, Feature.Start, Feature.End)
colnames(file1)[2] <- "Chromosome"
colnames(file1)[3] <- "Start"
colnames(file1)[4] <- "End"
file1$Position <- file1$Start

file2<-fread("GeneHancer.txt")
file2$chrom <- gsub('chr', '', file2$chrom)
colnames(file2)[1] <- "Chr"
colnames(file2)[2] <- "Start"
colnames(file2)[3] <- "End"

file2$Chr <- as.integer(file2$Chr)
file1$Chr <- as.integer(file1$Chr)
file2$Start <- as.integer(file2$Start)
file1$Start <- as.integer(file1$Start)
file2$End <- as.integer(file2$End)
file1$End <- as.integer(file1$End)

test <- file1[,.(Chr,Position,Gene)][file2,
                                      .(Chr,Position=x.Position, Start, End,Gene, score,name,elementType),
                                      on=.(Chr,Position>=Start,Position<=End),nomatch=0L]


setDT(test)[, count := uniqueN(name), by =Gene]


score <- test %>%
  group_by(Gene) %>%
  summarize(text = str_c(score, collapse = ", "))

count <- test %>%
  group_by(Gene) %>%
  summarize(text = str_c(count, collapse = ", "))



count <- Reduce(function(x, y) merge(x, y, by='Gene', all=TRUE), list(score, count))

colnames(count)[2] = "GeneHancerHighestScore"
colnames(count)[3] = "EnhancerCount"

count$EnhancerCount <- sapply(count$EnhancerCount, function(row) sum(as.numeric(strsplit(row, ", ")[[1]])), USE.NAMES = FALSE)


count <- select(count, Gene,EnhancerCount)
count <- count[complete.cases(count), ]
fwrite(count, "GeneHancer_Gene.txt") 

head(count)
```

    ##       Gene EnhancerCount
    ## 1      7SK             2
    ## 2 A1BG-AS1            12
    ## 3  A2M-AS1            10
    ## 4   A4GALT           190
    ## 5     AACS           209
    ## 6    AADAC            78

### H3K4me1 site counts per gene from UCSC data (Huvec cell line):

``` r
file1 <- fread("Variants_GenesAnnotated.txt")
file1 <- select(file1, Gene, Feature.Chr, Feature.Start, Feature.End)
colnames(file1)[2] <- "Chromosome"
colnames(file1)[3] <- "Start"
colnames(file1)[4] <- "End"
file1$Position <- file1$Start


file2<-fread("wgEncodeBroadHistoneHuvecH3k4me1StdPk.txt")
colnames(file2)[2] <- "Chr"
file2$Chr <- gsub('chr', '', file2$Chr)
colnames(file2)[3] <- "Start"
colnames(file2)[4] <- "End"
colnames(file2)[8] <- "SignalValue_H3k4me1"

file2$Chr <- as.numeric(file2$Chr)
file1$Chr <- as.integer(file1$Chr)
file2$ID <- seq.int(nrow(file2))

file2$count <- 1
file2[, c("low", "high") := .(Start - 5000, End  + 5000)]

df <- file1[,.(Chr,Position,Start, End, Gene)][file2,
                                      .(Chr,Position=x.Position, Start, End,Gene,ID,SignalValue_H3k4me1, count),
                                      on=.(Chr,Position>=low,Position<=high),nomatch=0L]

df2 <- df %>%
  group_by(Gene) %>%
  summarize(text = str_c(SignalValue_H3k4me1, collapse = ", "))
colnames(df2)[2] = "SignalValue_H3k4me1"

df2$SignalValue_H3k4me1_median <- sapply(strsplit(df2$SignalValue_H3k4me1, ","), function(x) median(as.numeric(x)))

genelength <- fread("GeneLength.txt")

count <- df %>%
  group_by(Gene) %>%
  summarize(text = str_c(count, collapse = ", "))

count$text <- sapply(strsplit(count$text, ", "), function(x) sum(as.numeric(x)))

colnames(count)[2] = "H3k4me1_count"

merge_df <- merge(df2, count, by='Gene')
final <- select(merge_df, Gene, H3k4me1_count, SignalValue_H3k4me1_median)
final <- final[complete.cases(final), ]
fwrite(final, "H3k4me1_Gene.txt")

head(final)
```

    ##        Gene H3k4me1_count SignalValue_H3k4me1_median
    ## 1       7SK            10                  10.253500
    ## 2      A1BG            16                   5.427715
    ## 3  A1BG-AS1            24                   5.427715
    ## 4       A2M           122                  10.072100
    ## 5   A2M-AS1            10                  10.072100
    ## 6 A2ML1-AS1           120                   6.268290

### H3K4me3 site counts per gene from UCSC data (Huvec cell line):

``` r
file1 <- fread("Variants_GenesAnnotated.txt")
file1 <- select(file1, Gene, Feature.Chr, Feature.Start, Feature.End)
colnames(file1)[2] <- "Chromosome"
colnames(file1)[3] <- "Start"
colnames(file1)[4] <- "End"
file1$Position <- file1$Start

file2<-fread("wgEncodeBroadHistoneHuvecH3k4me3StdPk.txt")
colnames(file2)[2] <- "Chr"
file2$Chr <- gsub('chr', '', file2$Chr)
colnames(file2)[3] <- "Start"
colnames(file2)[4] <- "End"
colnames(file2)[8] <- "SignalValue_H3k4me3"

file2$Chr <- as.numeric(file2$Chr)
file1$Chr <- as.integer(file1$Chr)
file2$ID <- seq.int(nrow(file2))

file2$count <- 1
file2[, c("low", "high") := .(Start - 5000, End  + 5000)]

df <- file1[,.(Chr,Position,Start, End, Gene)][file2,
                                               .(Chr,Position=x.Position, Start, End,Gene,ID,SignalValue_H3k4me3, count),
                                               on=.(Chr,Position>=low,Position<=high),nomatch=0L]

df2 <- df %>%
  group_by(Gene) %>%
  summarize(text = str_c(SignalValue_H3k4me3, collapse = ", "))
colnames(df2)[2] = "SignalValue_H3k4me3"

df2$SignalValue_H3k4me3_median <- sapply(strsplit(df2$SignalValue_H3k4me3, ","), function(x) median(as.numeric(x)))

genelength <- fread("GeneLength.txt")

count <- df %>%
  group_by(Gene) %>%
  summarize(text = str_c(count, collapse = ", "))

count$text <- sapply(strsplit(count$text, ", "), function(x) sum(as.numeric(x)))

colnames(count)[2] = "H3k4me3_count"

merge_df <- merge(df2, count, by='Gene')
final <- select(merge_df, Gene, H3k4me3_count, SignalValue_H3k4me3_median)
final <- final[complete.cases(final), ]
fwrite(final, "H3k4me3_Gene.txt")

head(final)
```

    ##       Gene H3k4me3_count SignalValue_H3k4me3_median
    ## 1      7SK             2                   19.01710
    ## 2     A1BG             8                    3.54593
    ## 3 A1BG-AS1            12                    3.54593
    ## 4      A2M           122                   11.12310
    ## 5  A2M-AS1            10                   11.12310
    ## 6  A3GALT2            36                    4.36892

### H3k27Ac site counts per gene from UCSC data (Huvec cell line):

``` r
file1 <- fread("Variants_GenesAnnotated.txt")
file1 <- select(file1, Gene, Feature.Chr, Feature.Start, Feature.End)
colnames(file1)[2] <- "Chromosome"
colnames(file1)[3] <- "Start"
colnames(file1)[4] <- "End"
file1$Position <- file1$Start

file2<-fread("wgEncodeBroadHistoneHuvecH3k27AcStdPk.txt")
colnames(file2)[2] <- "Chr"
file2$Chr <- gsub('chr', '', file2$Chr)
colnames(file2)[3] <- "Start"
colnames(file2)[4] <- "End"
colnames(file2)[8] <- "SignalValue_H3k27Ac"

file2$Chr <- as.numeric(file2$Chr)
file1$Chr <- as.integer(file1$Chr)
file2$ID <- seq.int(nrow(file2))

file2$count <- 1
file2[, c("low", "high") := .(Start - 5000, End  + 5000)]

df <- file1[,.(Chr,Position,Start, End, Gene)][file2,
                                      .(Chr,Position=x.Position, Start, End,Gene,ID, count, SignalValue_H3k27Ac),
                                      on=.(Chr,Position>=low,Position<=high),nomatch=0L]

df2 <- df %>%
  group_by(Gene) %>%
  summarize(text = str_c(SignalValue_H3k27Ac, collapse = ", "))
colnames(df2)[2] = "SignalValue_H3k27Ac"

df2$SignalValue_H3k27Ac_median <- sapply(strsplit(df2$SignalValue_H3k27Ac, ","), function(x) median(as.numeric(x)))

genelength <- fread("GeneLength.txt")

count <- df %>%
  group_by(Gene) %>%
  summarize(text = str_c(count, collapse = ", "))

count$text <- sapply(strsplit(count$text, ", "), function(x) sum(as.numeric(x)))

colnames(count)[2] = "H3k27Ac_count"

merge_df <- merge(df2, count, by='Gene')
final <- select(merge_df, Gene, H3k27Ac_count, SignalValue_H3k27Ac_median)
final <- final[complete.cases(final), ]
fwrite(final, "H3k27Ac_Gene.txt")

head(final)
```

    ##      Gene H3k27Ac_count SignalValue_H3k27Ac_median
    ## 1     7SK             2                   16.89090
    ## 2     A2M           122                   10.02010
    ## 3 A2M-AS1            10                   10.02010
    ## 4    AACS           209                    6.06894
    ## 5   AAGAB           100                   17.44700
    ## 6   AAMDC           149                    7.38085

### DeepSEA variant annotation

-   Functional significance scores taken per variants associated to
    blood pressure (found by multiple querys in DeepSEA’s online
    analysis tool)
-   DeepSEA variant file then merged with gene list by chromosome
    position and minimum value per gene selected

``` r
file1 <- fread('infile.vcf.out_file1.funsig')
file2 <- fread('infile.vcf.out_file2.funsig')
file3 <- fread('infile.vcf.out_file3.funsig')
file4 <- fread('infile.vcf.out_file4.funsig')
file5 <- fread('infile.vcf.out_file5.funsig')
file6 <- fread('infile.vcf.out_file6.funsig')
file7 <- fread('infile.vcf.out_file7.funsig')
file8 <- fread('infile.vcf.out_file8.funsig')
file9 <- fread('infile.vcf.out_file9.funsig')
file10 <- fread('infile.vcf.out_file10.funsig')
file11<- fread('infile.vcf.out_file11.funsig')
file12<- fread('infile.vcf.out_file12.funsig')
file13<- fread('infile.vcf.out_file13.funsig')
file14<- fread('infile.vcf.out_file14.funsig')
file16<- fread('infile.vcf.out_file16.funsig')
file17 <- fread('infile.vcf.out_file17.funsig')
file18<- fread('infile.vcf.out_file18.funsig')
file19<- fread('infile.vcf.out_file19.funsig')
file20 <- fread('infile.vcf.out_file20.funsig')
file21<- fread('infile.vcf.out_file21.funsig')


df <- rbind(file1,file2,file3,file4,file5,file6,file7,file8,file9,file10,
            file11,file12,file13,file14,file16,file17,file18,file19,file20,file21)

#Create CP column:
df$chr <- gsub("^.{0,3}", "", df$chr)
df$CP <- paste(df$chr, df$pos,sep=":")
colnames(df)[7]<-'DeepSEA_Functional_Significance'
df <- dplyr::select(df, CP, DeepSEA_Functional_Significance)

loci <- fread('Variants_GenesAnnotated.txt') #Made in GWAS variant data filtering
#Get genes per variant and their DeepSEA score:
deepsea_df <- merge(loci, df, by='CP',  all.x=TRUE)

deepsea_df <- subset(deepsea_df, select=c(Gene, DeepSEA_Functional_Significance))

temp <- deepsea_df[, lapply(.SD, paste0, collapse = ", "), by = Gene]

#Get min values per gene

min2 = function(x) if(all(is.na(x))) NA else min(x,na.rm = T)
getmin = function(col) str_extract_all(col,"[0-9\\.-]+") %>%
  lapply(.,function(x)min2(as.numeric(x)) ) %>%
  unlist() 

deepsea <- temp %>%
  mutate_at(names(temp)[2],getmin)

write.table(deepsea,"deepsea_genes.txt",sep="\t",row.names=FALSE)
head(deepsea)
```

    ##            Gene DeepSEA_Functional_Significance
    ## 1:         UPF2                              NA
    ## 2:      R3HCC1L                       0.0386700
    ## 3:        LOXL4                       0.0062338
    ## 4: RP11-34A14.3                       0.0613090
    ## 5:      PYROXD2                              NA
    ## 6:         HPS1                              NA

### GTEx tissue fold change data merged to combine all individual tissue files

Downloaded GTEx eQTL data (GTEx_Analysis_v8_eQTL.tar) and unzipped to
get egene files (\*.egenes.txt.gz files contain data for all genes
tested) 1. Read all files into environment and select only Gene name and
FC columns 2. Rename columns 3. Merge datasets together with genelist
for final file

``` r
#Read in all files in GTEx directory
temp = list.files(pattern="*.txt")
for (i in 1:length(temp)) assign(temp[i], fread(temp[i]))

#Rename Gene column for every dataset in environment
my_func <- function(x) {
  x <- x %>%
    select(gene_name,log2_aFC) 
}

e <- .GlobalEnv
nms <- ls(pattern = ".txt$", envir = e)
for(nm in nms) e[[nm]] <- my_func(e[[nm]])

#Rename all FC columns to unqiue names
colnames(Adipose_Subcutaneous.v8.egenes.txt)[2] <- 'Adipose_Subcutaneous_GTExFC'
colnames(Adipose_Visceral_Omentum.v8.egenes.txt)[2] <- 'Adipose_Visceral_Omentum_GTExFC'
colnames(Adrenal_Gland.v8.egenes.txt)[2] <- 'Adrenal_Gland_GTExFC'
colnames(Artery_Aorta.v8.egenes.txt)[2] <- 'Artery_Aorta_GTExFC'
colnames(Artery_Coronary.v8.egenes.txt)[2] <- 'Artery_Coronary_GTExFC'
colnames(Artery_Tibial.v8.egenes.txt)[2] <- 'Artery_Tibial_GTExFC'
colnames(Brain_Amygdala.v8.egenes.txt)[2] <- 'Brain_Amygdala_GTExFC'
colnames(Brain_Anterior_cingulate_cortex_BA24.v8.egenes.txt)[2] <- 'Brain_Anterior_cingulate_cortex_BA24_GTExFC'
colnames(Brain_Caudate_basal_ganglia.v8.egenes.txt)[2] <- 'Brain_Caudate_basal_ganglia_GTExFC'
colnames(Brain_Cerebellar_Hemisphere.v8.egenes.txt)[2] <- 'Brain_Cerebellar_Hemisphere_GTExFC'
colnames(Brain_Cerebellum.v8.egenes.txt)[2] <- 'Brain_Cerebellum_GTExFC'
colnames(Brain_Cortex.v8.egenes.txt)[2] <- 'Brain_Cortex_GTExFC'
colnames(Brain_Frontal_Cortex_BA9.v8.egenes.txt)[2] <- 'Brain_Frontal_Cortex_BA9_GTExFC'
colnames(Brain_Hippocampus.v8.egenes.txt)[2] <- 'Brain_Hippocampus_GTExFC'
colnames(Brain_Hypothalamus.v8.egenes.txt)[2] <- 'Brain_Hypothalamus_GTExFC'
colnames(Brain_Nucleus_accumbens_basal_ganglia.v8.egenes.txt)[2] <- 'Brain_Nucleus_accumbens_basal_ganglia_GTExFC'
colnames(Brain_Putamen_basal_ganglia.v8.egenes.txt)[2] <- 'Brain_Putamen_basal_ganglia_GTExFC'
colnames(`Brain_Spinal_cord_cervical_c-1.v8.egenes.txt`)[2] <- '`Brain_Spinal_cord_cervical_c-1`_GTExFC'
colnames(Brain_Substantia_nigra.v8.egenes.txt)[2] <- 'Brain_Substantia_nigra_GTExFC'
colnames(Breast_Mammary_Tissue.v8.egenes.txt)[2] <- 'Breast_Mammary_Tissue_GTExFC'
colnames(Cells_Cultured_fibroblasts.v8.egenes.txt)[2] <- 'Cells_Cultured_fibroblasts_GTExFC'
colnames(`Cells_EBV-transformed_lymphocytes.v8.egenes.txt`)[2] <- '`Cells_EBV-transformed_lymphocytes`_GTExFC'
colnames(Colon_Sigmoid.v8.egenes.txt)[2] <- 'Colon_Sigmoid_GTExFC'
colnames(Colon_Transverse.v8.egenes.txt)[2] <- 'Colon_Transverse_GTExFC'
colnames(Esophagus_Gastroesophageal_Junction.v8.egenes.txt)[2] <- 'Esophagus_Gastroesophageal_Junction_GTExFC'
colnames(Esophagus_Mucosa.v8.egenes.txt)[2] <- 'Esophagus_Mucosa_GTExFC'
colnames(Esophagus_Muscularis.v8.egenes.txt)[2] <- 'Esophagus_Muscularis_GTExFC'
colnames(Heart_Atrial_Appendage.v8.egenes.txt)[2] <- 'Heart_Atrial_Appendage_GTExFC'
colnames(Heart_Left_Ventricle.v8.egenes.txt)[2] <- 'Heart_Left_Ventricle_GTExFC'
colnames(Kidney_Cortex.v8.egenes.txt)[2] <- 'Kidney_Cortex_GTExFC'
colnames(Liver.v8.egenes.txt)[2] <- 'Liver_GTExFC'
colnames(Lung.v8.egenes.txt)[2] <- 'Lung_GTExFC'
colnames(Minor_Salivary_Gland.v8.egenes.txt)[2] <- 'Minor_Salivary_Gland_GTExFC'
colnames(Muscle_Skeletal.v8.egenes.txt)[2] <- 'Muscle_Skeletal_GTExFC'
colnames(Nerve_Tibial.v8.egenes.txt)[2] <- 'Nerve_Tibial_GTExFC'
colnames(Ovary.v8.egenes.txt)[2] <- 'Ovary_GTExFC'
colnames(Pancreas.v8.egenes.txt)[2] <- 'Pancreas_GTExFC'
colnames(Pituitary.v8.egenes.txt)[2] <- 'Pituitary_GTExFC'
colnames(Prostate.v8.egenes.txt)[2] <- 'Prostate_GTExFC'
colnames(Skin_Not_Sun_Exposed_Suprapubic.v8.egenes.txt)[2] <- 'Skin_Not_Sun_Exposed_Suprapubic_GTExFC'
colnames(Skin_Sun_Exposed_Lower_leg.v8.egenes.txt)[2] <- 'Skin_Sun_Exposed_Lower_leg_GTExFC'
colnames(Small_Intestine_Terminal_Ileum.v8.egenes.txt)[2] <- 'Small_Intestine_Terminal_Ileum_GTExFC'
colnames(Spleen.v8.egenes.txt)[2] <- 'Spleen_GTExFC'
colnames(Stomach.v8.egenes.txt)[2] <- 'Stomach_GTExFC'
colnames(Testis.v8.egenes.txt)[2] <- 'Testis_GTExFC'
colnames(Thyroid.v8.egenes.txt)[2] <- 'Thyroid_GTExFC'
colnames(Uterus.v8.egenes.txt)[2] <- 'Uterus_GTExFC'
colnames(Vagina.v8.egenes.txt)[2] <- 'Vagina_GTExFC'
colnames(Whole_Blood.v8.egenes.txt)[2] <- 'Whole_Blood_GTExFC'

rm(temp, nm, nms, my_func, i, e)
setwd("~/Documents/PhD Year 2/BP-GWAS-Predict/All GWAS data sorting")

genes <-  fread("pval015_genelist_filtered_evangelou.txt")
genes <- select(genes, Gene)
colnames(genes)[1] <- 'gene_name'

#Remove duplicate genes in any datasets (that were causing the following merge to be too large)
df_varnames <- ls()[ grep(".txt$", ls()) ]
for (dfname in df_varnames) {
  df <- get(dfname)
  assign(dfname, df[! duplicated(df$gene_name), ])
}

GTEx_df <- Reduce(function(x, y) merge(x, y, by='gene_name', all.x=TRUE), list(genes, Adipose_Subcutaneous.v8.egenes.txt,
                                                                       Adipose_Visceral_Omentum.v8.egenes.txt,
                                                                       Adrenal_Gland.v8.egenes.txt,
                                                                       Artery_Aorta.v8.egenes.txt,
                                                                       Artery_Coronary.v8.egenes.txt,
                                                                       Artery_Tibial.v8.egenes.txt,
                                                                       Brain_Amygdala.v8.egenes.txt,
                                                                       Brain_Anterior_cingulate_cortex_BA24.v8.egenes.txt,
                                                                       Brain_Caudate_basal_ganglia.v8.egenes.txt,
                                                                       Brain_Cerebellar_Hemisphere.v8.egenes.txt,
                                                                       Brain_Cerebellum.v8.egenes.txt,
                                                                       Brain_Cortex.v8.egenes.txt,
                                                                       Brain_Frontal_Cortex_BA9.v8.egenes.txt,
                                                                       Brain_Hippocampus.v8.egenes.txt,
                                                                       Brain_Hypothalamus.v8.egenes.txt,
                                                                       Brain_Nucleus_accumbens_basal_ganglia.v8.egenes.txt,
                                                                       Brain_Putamen_basal_ganglia.v8.egenes.txt,
                                                                       `Brain_Spinal_cord_cervical_c-1.v8.egenes.txt`,
                                                                       Brain_Substantia_nigra.v8.egenes.txt,
                                                                       Breast_Mammary_Tissue.v8.egenes.txt,
                                                                       Cells_Cultured_fibroblasts.v8.egenes.txt,
                                                                       `Cells_EBV-transformed_lymphocytes.v8.egenes.txt`,
                                                                       Colon_Sigmoid.v8.egenes.txt,
                                                                       Colon_Transverse.v8.egenes.txt,
                                                                       Esophagus_Gastroesophageal_Junction.v8.egenes.txt,
                                                                       Esophagus_Mucosa.v8.egenes.txt,
                                                                       Esophagus_Muscularis.v8.egenes.txt,
                                                                       Heart_Atrial_Appendage.v8.egenes.txt,
                                                                       Heart_Left_Ventricle.v8.egenes.txt,
                                                                       Kidney_Cortex.v8.egenes.txt,
                                                                       Liver.v8.egenes.txt,
                                                                       Lung.v8.egenes.txt,
                                                                       Minor_Salivary_Gland.v8.egenes.txt,
                                                                       Muscle_Skeletal.v8.egenes.txt,
                                                                       Nerve_Tibial.v8.egenes.txt,
                                                                       Ovary.v8.egenes.txt,
                                                                       Pancreas.v8.egenes.txt,
                                                                       Pituitary.v8.egenes.txt,
                                                                       Prostate.v8.egenes.txt,
                                                                       Skin_Not_Sun_Exposed_Suprapubic.v8.egenes.txt,
                                                                       Skin_Sun_Exposed_Lower_leg.v8.egenes.txt,
                                                                       Small_Intestine_Terminal_Ileum.v8.egenes.txt,
                                                                       Spleen.v8.egenes.txt,
                                                                       Stomach.v8.egenes.txt,
                                                                       Testis.v8.egenes.txt,
                                                                       Thyroid.v8.egenes.txt,
                                                                       Uterus.v8.egenes.txt,
                                                                       Vagina.v8.egenes.txt,
                                                                       Whole_Blood.v8.egenes.txt
))


colnames(GTEx_df)[1] <- 'Gene'
#write.table(GTEx_df,"./GTEx_FC.txt", row.names=FALSE)

dim(GTEx_df)
```

    ## [1] 2097   50

``` r
names(GTEx_df)
```

    ##  [1] "Gene"                                        
    ##  [2] "Adipose_Subcutaneous_GTExFC"                 
    ##  [3] "Adipose_Visceral_Omentum_GTExFC"             
    ##  [4] "Adrenal_Gland_GTExFC"                        
    ##  [5] "Artery_Aorta_GTExFC"                         
    ##  [6] "Artery_Coronary_GTExFC"                      
    ##  [7] "Artery_Tibial_GTExFC"                        
    ##  [8] "Brain_Amygdala_GTExFC"                       
    ##  [9] "Brain_Anterior_cingulate_cortex_BA24_GTExFC" 
    ## [10] "Brain_Caudate_basal_ganglia_GTExFC"          
    ## [11] "Brain_Cerebellar_Hemisphere_GTExFC"          
    ## [12] "Brain_Cerebellum_GTExFC"                     
    ## [13] "Brain_Cortex_GTExFC"                         
    ## [14] "Brain_Frontal_Cortex_BA9_GTExFC"             
    ## [15] "Brain_Hippocampus_GTExFC"                    
    ## [16] "Brain_Hypothalamus_GTExFC"                   
    ## [17] "Brain_Nucleus_accumbens_basal_ganglia_GTExFC"
    ## [18] "Brain_Putamen_basal_ganglia_GTExFC"          
    ## [19] "`Brain_Spinal_cord_cervical_c-1`_GTExFC"     
    ## [20] "Brain_Substantia_nigra_GTExFC"               
    ## [21] "Breast_Mammary_Tissue_GTExFC"                
    ## [22] "Cells_Cultured_fibroblasts_GTExFC"           
    ## [23] "`Cells_EBV-transformed_lymphocytes`_GTExFC"  
    ## [24] "Colon_Sigmoid_GTExFC"                        
    ## [25] "Colon_Transverse_GTExFC"                     
    ## [26] "Esophagus_Gastroesophageal_Junction_GTExFC"  
    ## [27] "Esophagus_Mucosa_GTExFC"                     
    ## [28] "Esophagus_Muscularis_GTExFC"                 
    ## [29] "Heart_Atrial_Appendage_GTExFC"               
    ## [30] "Heart_Left_Ventricle_GTExFC"                 
    ## [31] "Kidney_Cortex_GTExFC"                        
    ## [32] "Liver_GTExFC"                                
    ## [33] "Lung_GTExFC"                                 
    ## [34] "Minor_Salivary_Gland_GTExFC"                 
    ## [35] "Muscle_Skeletal_GTExFC"                      
    ## [36] "Nerve_Tibial_GTExFC"                         
    ## [37] "Ovary_GTExFC"                                
    ## [38] "Pancreas_GTExFC"                             
    ## [39] "Pituitary_GTExFC"                            
    ## [40] "Prostate_GTExFC"                             
    ## [41] "Skin_Not_Sun_Exposed_Suprapubic_GTExFC"      
    ## [42] "Skin_Sun_Exposed_Lower_leg_GTExFC"           
    ## [43] "Small_Intestine_Terminal_Ileum_GTExFC"       
    ## [44] "Spleen_GTExFC"                               
    ## [45] "Stomach_GTExFC"                              
    ## [46] "Testis_GTExFC"                               
    ## [47] "Thyroid_GTExFC"                              
    ## [48] "Uterus_GTExFC"                               
    ## [49] "Vagina_GTExFC"                               
    ## [50] "Whole_Blood_GTExFC"

## <a id="7. Merging databases and dividing training data"></a>7. Merging databases and dividing training data

``` r
#Load all collected feature data
df1 <- fread('pval015_genelist_filtered_evangelou.txt')
df2 <- fread('GTEx_TPM.txt')
df2 <- df2[,2:56] #For GTEx median TPM data
df3 <- fread('GTEx_FC.txt')
df4<-fread("gwascata_genes.txt", fill=TRUE)
df4$logpval_gwascatalog[is.na(df4$logpval_gwascatalog)] <- 0
#Getting partial gene matches for longer strings for genes in GWAScatalog:
df1$Genes_Extended <- 
  df4$Gene[sapply(df1$Gene, 
                         function(x) match(x, substr(df4$Gene, 1, nchar(x))))]
colnames(df4)[1] <- 'Genes_Extended'
df1 <- merge(df1, df4, by='Genes_Extended', all.x=T)
df1 <- subset(df1, select = -c(Genes_Extended, gwastrait) )
df5 <- fread('deepsea_genes.txt')
df6 <- fread('exomiser_scores.txt') 
df7 <- fread('Genelist_drug_and_databases.txt')
df8 <- fread("genic_intolerance_v3_12Mar16.txt")
df8 <- select(df8, Gene, `ALL_0.1%`)
colnames(df8)[2] <-'RVIS_Score'
df9 <- fread('genie_textminingBP.txt')
colnames(df9)[3] <-'Gene'
df9 <- select(df9, Gene, Rank)
df10 <- fread('panglaodb_cells_highestscores.txt')
df11 <- fread('Genes_CpG.txt')
df12 <- fread('DNase_Gene.txt')
df13 <- fread('GeneHancer_Gene.txt')
df14<- fread('SDI_genes.txt')
df15<- fread('H3K4me1_Gene.txt')
df16 <- fread('H3K4me3_Gene.txt')
df17 <- fread('H3K27Ac_Gene.txt')
df18 <- fread('exac_march16_pLI.txt')
colnames(df18)[2] <- 'Gene'
colnames(df18)[20] <- 'pLI_ExAC'
colnames(df18)[21] <- 'CNVs_ExAC'
df18<-df18[,c(2, 20, 21)]
df19 <- fread('GeneLength.txt')
colnames(df19)[1] <- 'Gene'
df20 <- fread('MGImarkerQuery_hypertension.txt', fill=TRUE)
colnames(df20)[8] <- 'Gene'
df20 <- select(df20, Gene)
df20$MGI_Gene <- 1
df21 <- fread('Genes_ratmodels.txt')
df21 <- select(df21, Gene)
df21$Ratmodel_Gene <- 1
df22 <- fread('GDI.txt')
df22 <- df22[,1:2]
df23 <- fread('GTEx_sexbias.txt')
df24 <- fread('HIPred.tsv')
df25 <- fread('IPA_Activity.txt')
df26 <- fread('O_GLcNAc_Score.txt')
df27 <- fread('EMS.txt')
df28 <- fread('PharmGKB_Score.txt')
df29 <- fread('ExAc_constraint.txt')

#Merge all loaded data
all_features <- Reduce(function(x, y) merge(x, y, by = 'Gene', all.x = TRUE),
                       list(df1, df2, df3, df5, df6, df7, df8, df9, df10, df11,
                            df12, df13, df14, df15, df16, df17, df18, df19,df20,
                            df21, df22, df23, df24, df25, df26, df27, df28, df29))
all_features <- all_features[!duplicated(all_features$Gene), ]

#Label training genes and write final file
all_features$BPlabel <- NA
all_features$BPmed[all_features$BPmed==0] <- NA

all_features[, BPlabel := fcase( 
  Mechanism == 1 | BPmed >= 3 , 'most likely', 
  BPmed < 3 | SideeffectFreq > 1 | Rank >= 1, 'probable',
  insignif == 1, 'least likely',
  default = 'unknown')]

all_features = subset(all_features, select = -c(BPmed, Mechanism, Others, insignif, SideeffectFreq, Sideeffectseverity, Rank) )


testbp1 <- filter(all_features, BPlabel == 'least likely') 
testbp2 <- filter(all_features, BPlabel == 'probable')  
testbp3 <- filter(all_features, BPlabel == 'most likely') 

training <- filter(all_features, BPlabel != 'unknown')
unknown <- filter(all_features, BPlabel == 'unknown')

write.table(training,"pval015_BP_training.txt", sep = "\t", row.names = F, quote = F)
write.table(unknown,"pval015_BP_unknown.txt", sep = "\t", row.names = F, quote = F)

write.table(all_features,"pval015_all_genes_all_features.txt", sep = "\t", row.names = F, quote = F)


dim(training)
```

    ## [1] 293 140

``` r
names(training)
```

    ##   [1] "Gene"                                             
    ##   [2] "REVEL.max"                                        
    ##   [3] "MetaSVM_rankscore.max"                            
    ##   [4] "MetaLR_rankscore.max"                             
    ##   [5] "MCAP.max"                                         
    ##   [6] "wgEncodeBroadHmmHuvecHMM.count"                   
    ##   [7] "betamax"                                          
    ##   [8] "logpval_gwascatalog"                              
    ##   [9] "Adipose - Subcutaneous_GTExTPM"                   
    ##  [10] "Adipose - Visceral (Omentum)_GTExTPM"             
    ##  [11] "Adrenal Gland_GTExTPM"                            
    ##  [12] "Artery - Aorta_GTExTPM"                           
    ##  [13] "Artery - Coronary_GTExTPM"                        
    ##  [14] "Artery - Tibial_GTExTPM"                          
    ##  [15] "Bladder_GTExTPM"                                  
    ##  [16] "Brain - Amygdala_GTExTPM"                         
    ##  [17] "Brain - Anterior cingulate cortex (BA24)_GTExTPM" 
    ##  [18] "Brain - Caudate (basal ganglia)_GTExTPM"          
    ##  [19] "Brain - Cerebellar Hemisphere_GTExTPM"            
    ##  [20] "Brain - Cerebellum_GTExTPM"                       
    ##  [21] "Brain - Cortex_GTExTPM"                           
    ##  [22] "Brain - Frontal Cortex (BA9)_GTExTPM"             
    ##  [23] "Brain - Hippocampus_GTExTPM"                      
    ##  [24] "Brain - Hypothalamus_GTExTPM"                     
    ##  [25] "Brain - Nucleus accumbens (basal ganglia)_GTExTPM"
    ##  [26] "Brain - Putamen (basal ganglia)_GTExTPM"          
    ##  [27] "Brain - Spinal cord (cervical c-1)_GTExTPM"       
    ##  [28] "Brain - Substantia nigra_GTExTPM"                 
    ##  [29] "Breast - Mammary Tissue_GTExTPM"                  
    ##  [30] "Cells - Cultured fibroblasts_GTExTPM"             
    ##  [31] "Cells - EBV-transformed lymphocytes_GTExTPM"      
    ##  [32] "Cervix - Ectocervix_GTExTPM"                      
    ##  [33] "Cervix - Endocervix_GTExTPM"                      
    ##  [34] "Colon - Sigmoid_GTExTPM"                          
    ##  [35] "Colon - Transverse_GTExTPM"                       
    ##  [36] "Esophagus - Gastroesophageal Junction_GTExTPM"    
    ##  [37] "Esophagus - Mucosa_GTExTPM"                       
    ##  [38] "Esophagus - Muscularis_GTExTPM"                   
    ##  [39] "Fallopian Tube_GTExTPM"                           
    ##  [40] "Heart - Atrial Appendage_GTExTPM"                 
    ##  [41] "Heart - Left Ventricle_GTExTPM"                   
    ##  [42] "Kidney - Cortex_GTExTPM"                          
    ##  [43] "Kidney - Medulla_GTExTPM"                         
    ##  [44] "Liver_GTExTPM"                                    
    ##  [45] "Lung_GTExTPM"                                     
    ##  [46] "Minor Salivary Gland_GTExTPM"                     
    ##  [47] "Muscle - Skeletal_GTExTPM"                        
    ##  [48] "Nerve - Tibial_GTExTPM"                           
    ##  [49] "Ovary_GTExTPM"                                    
    ##  [50] "Pancreas_GTExTPM"                                 
    ##  [51] "Pituitary_GTExTPM"                                
    ##  [52] "Prostate_GTExTPM"                                 
    ##  [53] "Skin - Not Sun Exposed (Suprapubic)_GTExTPM"      
    ##  [54] "Skin - Sun Exposed (Lower leg)_GTExTPM"           
    ##  [55] "Small Intestine - Terminal Ileum_GTExTPM"         
    ##  [56] "Spleen_GTExTPM"                                   
    ##  [57] "Stomach_GTExTPM"                                  
    ##  [58] "Testis_GTExTPM"                                   
    ##  [59] "Thyroid_GTExTPM"                                  
    ##  [60] "Uterus_GTExTPM"                                   
    ##  [61] "Vagina_GTExTPM"                                   
    ##  [62] "Whole Blood_GTExTPM"                              
    ##  [63] "Adipose_Subcutaneous_GTExFC"                      
    ##  [64] "Adipose_Visceral_Omentum_GTExFC"                  
    ##  [65] "Adrenal_Gland_GTExFC"                             
    ##  [66] "Artery_Aorta_GTExFC"                              
    ##  [67] "Artery_Coronary_GTExFC"                           
    ##  [68] "Artery_Tibial_GTExFC"                             
    ##  [69] "Brain_Amygdala_GTExFC"                            
    ##  [70] "Brain_Anterior_cingulate_cortex_BA24_GTExFC"      
    ##  [71] "Brain_Caudate_basal_ganglia_GTExFC"               
    ##  [72] "Brain_Cerebellar_Hemisphere_GTExFC"               
    ##  [73] "Brain_Cerebellum_GTExFC"                          
    ##  [74] "Brain_Cortex_GTExFC"                              
    ##  [75] "Brain_Frontal_Cortex_BA9_GTExFC"                  
    ##  [76] "Brain_Hippocampus_GTExFC"                         
    ##  [77] "Brain_Hypothalamus_GTExFC"                        
    ##  [78] "Brain_Nucleus_accumbens_basal_ganglia_GTExFC"     
    ##  [79] "Brain_Putamen_basal_ganglia_GTExFC"               
    ##  [80] "`Brain_Spinal_cord_cervical_c-1`_GTExFC"          
    ##  [81] "Brain_Substantia_nigra_GTExFC"                    
    ##  [82] "Breast_Mammary_Tissue_GTExFC"                     
    ##  [83] "Cells_Cultured_fibroblasts_GTExFC"                
    ##  [84] "`Cells_EBV-transformed_lymphocytes`_GTExFC"       
    ##  [85] "Colon_Sigmoid_GTExFC"                             
    ##  [86] "Colon_Transverse_GTExFC"                          
    ##  [87] "Esophagus_Gastroesophageal_Junction_GTExFC"       
    ##  [88] "Esophagus_Mucosa_GTExFC"                          
    ##  [89] "Esophagus_Muscularis_GTExFC"                      
    ##  [90] "Heart_Atrial_Appendage_GTExFC"                    
    ##  [91] "Heart_Left_Ventricle_GTExFC"                      
    ##  [92] "Kidney_Cortex_GTExFC"                             
    ##  [93] "Liver_GTExFC"                                     
    ##  [94] "Lung_GTExFC"                                      
    ##  [95] "Minor_Salivary_Gland_GTExFC"                      
    ##  [96] "Muscle_Skeletal_GTExFC"                           
    ##  [97] "Nerve_Tibial_GTExFC"                              
    ##  [98] "Ovary_GTExFC"                                     
    ##  [99] "Pancreas_GTExFC"                                  
    ## [100] "Pituitary_GTExFC"                                 
    ## [101] "Prostate_GTExFC"                                  
    ## [102] "Skin_Not_Sun_Exposed_Suprapubic_GTExFC"           
    ## [103] "Skin_Sun_Exposed_Lower_leg_GTExFC"                
    ## [104] "Small_Intestine_Terminal_Ileum_GTExFC"            
    ## [105] "Spleen_GTExFC"                                    
    ## [106] "Stomach_GTExFC"                                   
    ## [107] "Testis_GTExFC"                                    
    ## [108] "Thyroid_GTExFC"                                   
    ## [109] "Uterus_GTExFC"                                    
    ## [110] "Vagina_GTExFC"                                    
    ## [111] "Whole_Blood_GTExFC"                               
    ## [112] "DeepSEA_Functional_Significance"                  
    ## [113] "ExomiserScore"                                    
    ## [114] "RVIS_Score"                                       
    ## [115] "ubiquitousness.index.max"                         
    ## [116] "CpGcount"                                         
    ## [117] "DNaseCluster_count"                               
    ## [118] "EnhancerCount"                                    
    ## [119] "SDI"                                              
    ## [120] "H3k4me1_count"                                    
    ## [121] "SignalValue_H3k4me1_median"                       
    ## [122] "H3k4me3_count"                                    
    ## [123] "SignalValue_H3k4me3_median"                       
    ## [124] "H3k27Ac_count"                                    
    ## [125] "SignalValue_H3k27Ac_median"                       
    ## [126] "pLI_ExAC"                                         
    ## [127] "CNVs_ExAC"                                        
    ## [128] "Gene_length"                                      
    ## [129] "MGI_Gene"                                         
    ## [130] "Ratmodel_Gene"                                    
    ## [131] "GDI_Score"                                        
    ## [132] "GTEx_signif_sexbias"                              
    ## [133] "HIPred"                                           
    ## [134] "IPA_BP"                                           
    ## [135] "IPA_Activity"                                     
    ## [136] "O_GlcNAc_Score"                                   
    ## [137] "EMS_Max"                                          
    ## [138] "PharmGKB_Score"                                   
    ## [139] "ExAc_constraint_obs_exp"                          
    ## [140] "BPlabel"

## <a id="8. Gene length correlation"></a>8. Gene length correlation

``` r
train_variants <- select(training, Gene_length, CpGcount, EnhancerCount, DNaseCluster_count, H3k4me1_count,
                         H3k4me3_count, H3k27Ac_count, betamax, DeepSEA_Functional_Significance,
                         wgEncodeBroadHmmHuvecHMM.count, logpval_gwascatalog)
all_variants <- select(all_features, Gene_length, CpGcount, EnhancerCount, DNaseCluster_count, H3k4me1_count,
                         H3k4me3_count, H3k27Ac_count, betamax, DeepSEA_Functional_Significance,
                         wgEncodeBroadHmmHuvecHMM.count, logpval_gwascatalog)

training_table <- t(sapply(train_variants[, -1], function(x) {
  c(Pearsons = cor(train_variants$Gene_length, x, method = "pearson", use = "complete.obs"), 
    Spearman = cor(train_variants$Gene_length, x, method = "spearman", use = "complete.obs"))
}))

all_variants_table <- t(sapply(all_variants[, -1], function(x) {
  c(Pearsons = cor(all_variants$Gene_length, x, method = "pearson", use = "complete.obs"), 
    Spearman = cor(all_variants$Gene_length, x, method = "spearman", use = "complete.obs"))
}))

print('Training data gene length correlation for variant level features:')
```

    ## [1] "Training data gene length correlation for variant level features:"

``` r
print(training_table)
```

    ##                                      Pearsons    Spearman
    ## CpGcount                         0.9826259497  0.89924606
    ## EnhancerCount                    0.9780680462  0.87029991
    ## DNaseCluster_count               0.9795220050  0.90017240
    ## H3k4me1_count                    0.8322488066  0.88116252
    ## H3k4me3_count                    0.9609045404  0.87476255
    ## H3k27Ac_count                    0.8863011342  0.87776844
    ## betamax                          0.0176626649 -0.03526183
    ## DeepSEA_Functional_Significance -0.2026521477 -0.30311480
    ## wgEncodeBroadHmmHuvecHMM.count   0.0242076148 -0.03672767
    ## logpval_gwascatalog             -0.0006769054  0.17709528

``` r
print('Total data gene length correlation for variant level features:')
```

    ## [1] "Total data gene length correlation for variant level features:"

``` r
print(all_variants_table)
```

    ##                                     Pearsons    Spearman
    ## CpGcount                         0.002047835  0.93671627
    ## EnhancerCount                   -0.014141909  0.91657671
    ## DNaseCluster_count              -0.006560648  0.92443591
    ## H3k4me1_count                   -0.006066534  0.88000832
    ## H3k4me3_count                   -0.008147886  0.89163857
    ## H3k27Ac_count                   -0.008658441  0.88019883
    ## betamax                         -0.029448609  0.02278911
    ## DeepSEA_Functional_Significance -0.020286745 -0.38983447
    ## wgEncodeBroadHmmHuvecHMM.count  -0.006386533 -0.03509751
    ## logpval_gwascatalog              0.013937852  0.30100994

## <a id="9. SessionInfo"></a>9. SessionInfo

``` r
sessionInfo()
```

    ## R version 4.2.1 (2022-06-23)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur ... 10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] GenomicRanges_1.48.0  GenomeInfoDb_1.32.3   IRanges_2.30.1       
    ##  [4] S4Vectors_0.34.0      BiocGenerics_0.42.0   matrixStats_0.62.0   
    ##  [7] splitstackshape_1.4.8 sqldf_0.4-11          RSQLite_2.2.16       
    ## [10] gsubfn_0.7            proto_1.0.0           janitor_2.1.0        
    ## [13] data.table_1.14.2     forcats_0.5.2         purrr_0.3.4          
    ## [16] readr_2.1.2           tidyr_1.2.0           tibble_3.1.8         
    ## [19] ggplot2_3.3.6         tidyverse_1.3.2       dplyr_1.0.9          
    ## [22] plyr_1.8.7            stringr_1.4.1        
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] httr_1.4.4             bit64_4.0.5            jsonlite_1.8.0        
    ##  [4] modelr_0.1.9           assertthat_0.2.1       blob_1.2.3            
    ##  [7] GenomeInfoDbData_1.2.8 googlesheets4_1.0.1    cellranger_1.1.0      
    ## [10] yaml_2.3.5             pillar_1.8.1           backports_1.4.1       
    ## [13] glue_1.6.2             chron_2.3-57           digest_0.6.29         
    ## [16] XVector_0.36.0         rvest_1.0.3            snakecase_0.11.0      
    ## [19] colorspace_2.0-3       htmltools_0.5.3        pkgconfig_2.0.3       
    ## [22] broom_1.0.1            haven_2.5.1            zlibbioc_1.42.0       
    ## [25] scales_1.2.1           tzdb_0.3.0             googledrive_2.0.0     
    ## [28] generics_0.1.3         ellipsis_0.3.2         cachem_1.0.6          
    ## [31] withr_2.5.0            cli_3.3.0              magrittr_2.0.3        
    ## [34] crayon_1.5.1           readxl_1.4.1           memoise_2.0.1         
    ## [37] evaluate_0.16          fs_1.5.2               fansi_1.0.3           
    ## [40] xml2_1.3.3             tools_4.2.1            hms_1.1.2             
    ## [43] gargle_1.2.0           lifecycle_1.0.1        munsell_0.5.0         
    ## [46] reprex_2.0.2           compiler_4.2.1         rlang_1.0.4           
    ## [49] RCurl_1.98-1.8         grid_4.2.1             rstudioapi_0.14       
    ## [52] bitops_1.0-7           rmarkdown_2.16         gtable_0.3.0          
    ## [55] DBI_1.1.3              R6_2.5.1               lubridate_1.8.0       
    ## [58] knitr_1.40             fastmap_1.1.0          bit_4.0.4             
    ## [61] utf8_1.2.2             stringi_1.7.8          Rcpp_1.0.9            
    ## [64] vctrs_0.4.1            dbplyr_2.2.1           tidyselect_1.1.2      
    ## [67] xfun_0.32
