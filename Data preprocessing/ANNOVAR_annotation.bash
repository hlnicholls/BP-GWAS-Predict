#!/bin/sh
#$ -cwd           # Set the working directory for the job to the current directory
#$ -pe smp 1      # Request 1 core
#$ -l h_rt=120:0:0 # Request 24 hour runtime
#$ -l h_vmem=15G   # Request 1GB RAM


module load annovar/2018Apr16

table_annovar.pl Total_genes.txt /data/WHRI-Bioinformatics/RefData/Human/annovar/humandb --buildver hg19 -out Total_Genes_SelectTypes_Annotated.ANN -remove -protocol refGene,genomicSuperDups,wgRna,gwasCatalog,wgEncodeBroadHmmHuvecHMM,dbnsfp33a,clinvar_20200316,intervar_20180118,gnomad_genome,mcap,revel,mitimpact24 -operation g,r,r,r,r,f,f,f,f,f,f,f -nastring . -otherinfo
