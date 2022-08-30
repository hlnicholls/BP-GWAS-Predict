#!/bin/sh
#$ -cwd           # Set the working directory for the job to the current directory
#$ -pe smp 1      # Request 1 core
#$ -l h_rt=120:0:0 # Request 24 hour runtime
#$ -l h_vmem=15G   # Request 1GB RAM

module load bedtools

	# ANNOTATION USING BED TOOLS ALL SNPS AND PROXIES  #########################################################
	# 1. First annotate the SNPs within the genes boundries (start and end of transciption)
	# 2. Annotate remaining SNPs to within Feature
	# 3. Annotate remaining SNPs to closest gene/feature
	# 4. Merge lists
	# 5. Annotate using ANNOVAR

#ALLGWAS_variants.txt columns:
#1.  Chr  
#2.  Start    
#3.  End 
#4.  Ref 
#5.  Alt 
#6.  BETAsbp 
#7.  BETAdbp  
#8.  BETApp      
#9.  minP 
#10. minTRAIT 
#11. BETAmean       
#12. CP

#ALLGWAS_variants.txt is just GWASresults750k_minP_meanBETA.csv with CP column created
head ALLGWAS_variants.txt
#Chr     Start   End     Ref     Alt     BETAsbp BETAdbp BETApp  minP    minTRAIT        BETAmean        CP
#1       752566  752566  a       g       0.0687  0.0605  0.0177  0.01714 DBP     0.0646  1:752566
#1       885689  885689  a       g       0.235   0.1369  0.0807  0.0009269       DBP     0.18595 1:885689
#1       885699  885699  a       g       -0.2403 -0.1341 -0.0863 0.0009015       SBP     -0.1872 1:885699
#1       886006  886006  t       c       -0.2408 -0.1346 -0.0871 0.0008683       SBP     -0.1877 1:886006
#1       887801  887801  a       g       -0.2382 -0.1306 -0.0874 0.0009413       SBP     -0.1844 1:887801
#1       888639  888639  t       c       -0.2472 -0.133  -0.0933 0.0006506       SBP     -0.1901 1:888639
#1       888659  888659  t       c       -0.2323 -0.1277 -0.0845 0.001183        SBP     -0.18   1:888659
#1       889238  889238  a       g       -0.2282 -0.1275 -0.0812 0.00167 SBP     -0.17785        1:889238
#1       892745  892745  a       g       0.2401  0.1317  0.0868  0.0008727       SBP     0.1859  1:892745



#1) Filter hg19Rel92_AllgeneTypes_EnsemblHavana_0kb.txt for antisense, processed_transcript, protein_coding, and pseudogenes

#For some reason filtering in R and writing a tab delim file later gives an error from bedtools that file is not in tab format - so not filtering before bedtools
#hg19Rel92_AllgeneTypes_0kb.txt and AllgeneTypes_filtered.txt are both listed as same file type in unix and running cat -t does not work either

    #Ensure chromosomes are in 1-22 order
	sort -n -k1,1 -k2,2 hg19Rel92_AllgeneTypes_0kb.txt  > AllgeneTypes_sorted.txt


	# 1) annotate within genes (bedtools v2.17.0)
	sed 's/ /\t/g' ALLGWAS_variants.txt | tail -n +2 | intersectBed -a - -b AllgeneTypes_sorted.txt -wa -wb -loj > genes0kb.bed

	
	# get those not in genes
		awk 'BEGIN{OFS="\t"}{if($13 ==".")print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' genes0kb.bed > toFeatures.avinput
		wc -l genes0kb.bed   #7083535 
		wc -l toFeatures.avinput #7083535 


	# 2) Annotate remaining SNPs to closest gene/feature (within 5KB)

	sort -n -k1,1 -k2,2 toFeatures.avinput > toFeaturesSorted.avinput

	closestBed -d -t all -a toFeaturesSorted.avinput -b AllgeneTypes_sorted.txt  > toFeaturesSorted.bed


###########################################################################

head toFeaturesSorted.bed  #-d


#1       918384  918384  t       g       0.0533  0.0365  0.0155  0.05253 DBP     0.0449  1:918384        1       910579  917497  C1orf170        protein_coding  ensembl_havana       #887
#1       918573  918573  a       g       -0.0518 -0.0381 -0.0132 0.04347 DBP     -0.04495        1:918573        1       910579  917497  C1orf170        protein_coding       ensembl_havana  1076
#1       924898  924898  a       c       0.3172  0.1081  0.2183  0.0003453       PP      0.21265 1:924898        1       931346  933431  RP11-54O7.17    lincRNA havana       6448
#1       926351  926351  t       c       0.3092  0.1043  0.2149  0.0004249       PP      0.20675 1:926351        1       931346  933431  RP11-54O7.17    lincRNA havana       4995
#1       926431  926431  a       t       -0.3046 -0.1033 -0.2119 0.0005051       PP      -0.20395        1:926431        1       931346  933431  RP11-54O7.17    lincRNA      havana  #4915
#1       926621  926621  a       c       -0.3099 -0.1046 -0.2153 0.0004156       PP      -0.20725        1:926621        1       931346  933431  RP11-54O7.17    lincRNA      havana  4725
#1       928578  928578  a       g       0.3071  0.1027  0.2149  0.0004225       PP      0.2049  1:928578        1       931346  933431  RP11-54O7.17    lincRNA havana       2768
#1       928836  928836  t       c       0.305   0.1049  0.2106  0.0005334       PP      0.20495 1:928836        1       931346  933431  RP11-54O7.17    lincRNA havana       2510
#1       929327  929327  a       g       -0.3015 -0.1    -0.2127 0.0004839       PP      -0.20075        1:929327        1       931346  933431  RP11-54O7.17    lincRNA      havana  2019
#1       940005  940005  a       g       0.0344  0.0351  -0.0012 0.06487 DBP     0.03475 1:940005        1       943456  943609  RP11-54O7.10    pseudogene      havana       3451

#-D *****ERROR: -D option must be followed with "ref", "a", or "b"  - not sure how to apply yet
 

		##use ALL DISTANCE ..FILTER FOR DISTANCE WHEN PERFORMING ANALYSIS
		
		# LIST WITHIN 50KB
		awk '{if($19 <= 5000) print $0}' toFeaturesSorted.bed > ClosestGeneFeature5kb.bed
		
		#not actually used:
		# # # # LIST TO ANNOTATE TO NAME as CHR-POS
		# # # awk '{if($31 > 50000) print $0}' ClosestGeneFeature.bed |cut -f 1-26 - | uniq > toReNameLocusCHRPOS.txt
		# # # wc -l ClosestGeneFeature50kb.bed
		# # # 169791 ClosestGeneFeature50kb.bed
		# # # wc -l toReNameLocusCHRPOS.txt
		# # # 30393 toReNameLocusCHRPOS.txt
		
	# 4) CONCATENATE LISTS
	awk 'BEGIN{OFS="\t"}{if($13 != ".")print $0,0}' genes0kb.bed > T1

	# ClosestGeneFeature.bed (as in its current format)
	cat T1  ClosestGeneFeature5kb.bed | sort -k1,1 -k2,2 -n - > Total_Genes.bed
	
head Total_Genes.bed


#1       752566  752566  a       g       0.0687  0.0605  0.0177  0.01714 DBP     0.0646  1:752566        1       745489  753092  RP11-206L10.10  processed_transcripthavana   0
#1       885689  885689  a       g       0.235   0.1369  0.0807  0.0009269       DBP     0.18595 1:885689        1       879584  894689  NOC2L   protein_coding  ensembl_havana       0
#1       885699  885699  a       g       -0.2403 -0.1341 -0.0863 0.0009015       SBP     -0.1872 1:885699        1       879584  894689  NOC2L   protein_coding  ensembl_havana       0
#1       886006  886006  t       c       -0.2408 -0.1346 -0.0871 0.0008683       SBP     -0.1877 1:886006        1       879584  894689  NOC2L   protein_coding  ensembl_havana       0
#1       887801  887801  a       g       -0.2382 -0.1306 -0.0874 0.0009413       SBP     -0.1844 1:887801        1       879584  894689  NOC2L   protein_coding  ensembl_havana       0
#1       888639  888639  t       c       -0.2472 -0.133  -0.0933 0.0006506       SBP     -0.1901 1:888639        1       879584  894689  NOC2L   protein_coding  ensembl_havana       0
#1       888659  888659  t       c       -0.2323 -0.1277 -0.0845 0.001183        SBP     -0.18   1:888659        1       879584  894689  NOC2L   protein_coding  ensembl_havana       0
#1       889238  889238  a       g       -0.2282 -0.1275 -0.0812 0.00167 SBP     -0.17785        1:889238        1       879584  894689  NOC2L   protein_coding  ensembl_havana       0
#1       892745  892745  a       g       0.2401  0.1317  0.0868  0.0008727       SBP     0.1859  1:892745        1       879584  894689  NOC2L   protein_coding  ensembl_havana       0
#1       893280  893280  a       g       0.2378  0.1332  0.086   0.001053        SBP     0.1855  1:893280        1       879584  894689  NOC2L   protein_coding  ensembl_havana       0

module load R 
R 
library('data.table')
library('dplyr')
df <- fread('Total_Genes.bed')

colnames(df)[1] <-'Chr'
colnames(df)[2] <-'Start'
colnames(df)[3] <-'End'
colnames(df)[4] <-'Ref'
colnames(df)[5] <-'Alt'
colnames(df)[6] <-'BETAsbp'
colnames(df)[7] <-'BETAdbp'
colnames(df)[8] <-'BETApp'
colnames(df)[9] <-'minP'
colnames(df)[10] <-'minTRAIT'
colnames(df)[11] <-'BETAmean'
colnames(df)[12] <-'CP'
colnames(df)[13] <-'Feature.Chr'
colnames(df)[14] <-'Feature.Start'
colnames(df)[15] <-'Feature.End'
colnames(df)[16] <-'Gene'
colnames(df)[17] <-'Feature.Type'
colnames(df)[18] <-'Feature.Source'
colnames(df)[19] <-'Feature.Dist'

fwrite(df, 'Total_genes.txt', sep='\t', row.names=FALSE, col.names=TRUE)

head(df)

df1 <- filter(df, Feature.Type=='protein_coding')
df2 <- filter(df, Feature.Type=='processed_transcript')
df3 <- filter(df, Feature.Type=='pseudogene')
df4 <- filter(df, Feature.Type=='antisense')
final <- rbind(df1, df2, df3, df4)

nrow(final)
#[1] 3980785

tmp <- final[-which(duplicated(final$Gene)), ]
nrow(tmp)
#nrow: 33847

fwrite(final, 'Total_genes_final_selectTypes.txt', sep='\t', row.names=FALSE, col.names=TRUE)
fwrite(df, 'Total_genes_final_allTypes.txt', sep='\t', row.names=FALSE, col.names=TRUE)
