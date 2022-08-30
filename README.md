# BP-GWAS-Predict
Jupyter+R: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/hlnicholls/BP-GWAS-Predict/HEAD)<br />
RStudio: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/hlnicholls/BP-GWAS-Predict/HEAD?urlpath=urlpath%3Drstudio)
(binder will only work when the repo is made public)

BP-GWAS-Predict is a machine learning pipeline for predicting disease-causing genes post-genome-wide association study in blood pressure.


# Contents
[1. Introduction](#introduction)<br />
[2. Installation](#installation)<br />
[3. Step-by-step pipeline guide](#step-by-step-pipeline-guide)<br />
[4. References](#references)<br />
[5. Contact](#contact)<br />

# Introduction

This machine learning method aims to apply an optimized model for identifying real genome-wide association study (GWAS) signals from an individual blood pressure (BP) GWAS by Evangelou et al. in 2018 [1].

# Installation

Data pre-processing was applied in R and machine learning was applied in Python via the Anaconda environment. The setup of the Python Anaconda environment can be installed by running in the Anaconda command prompt:

```conda env create -f environment.yml```

Another option available is to install required R and Python packages by running ```install.R``` and ```pip install -r requirements.txt```

Gene annotation was applied in unix using ```bedtools``` and ```ANNOVAR```. Installation for bedtools is described here (https://bedtools.readthedocs.io/en/latest/content/installation.html) and installation for ANNOVAR is described here (https://annovar.openbioinformatics.org/en/latest/user-guide/download/)

## Step-by-step Pipeline Guide

![overview](https://user-images.githubusercontent.com/53306752/174135239-da7c48e5-67d4-4f3a-8cc1-a227465353e2.png)

**Overview of the Machine Learning Framework** 

GWAS variants from Evangelou et al. were annotated to genes and evaluated by machine learning. Genes were scored according to the their likely associations to BP based on their BP-drug relationship (genes interacting with BP drugs scored at 1), publication significance to BP (genes significantly named in BP publications scored at 0.75) and protein-protein interactions with known BP genes (genes least likely to affect BP with no interactions with BP genes scored at 0.1). Associated genes alongside insignificant least likely BP genes (p-value <0.15, not in 500kb+/- loci, no linkage disequilibrium) were annotated. Biological/functional data was collected from 20 databases. Further gene filtering identified genes that could be use as training data in the machine learning stage. Machine learning was applied benchmarking 14 models using 8 selected features. The top performing trained model was then used to score genes not in the training data, with the genes and their corresponding scores being then sorted into their loci to select the best genes per locus for gene enrichment analysis

### 1. Data preprocessing and integration:

An R markdown can be found [here](https://github.com/hlnicholls/BP-GWAS-Predict/blob/main/Data%20preprocessing/Data-Preprocessing-in-R.md) detailing the complete data-processing code steps outlined below.

#### Whole GWAS data preprocessing

Whole GWAS data was taken in order to annotate variant-level data and then take the most significant value per gene. The whole GWAS data was also used to identify variants with high p-values (>0.15) and no linkage disequilibrium measured with BP loci, provide an initial variants that could be identified as least likely to affect BP for further filtering. The further filtering involved selecting only variants not within 500kb+/- BP loci and not closely connected with known BP genes in protein-protein interaction (PPI) networks.

#### Variant-to-Gene annotation

Variants are then annotated to genes using bedtools within unix using a bash script ```Bedtools_Gene_Annotation.bash``` found [here](https://github.com/hlnicholls/BP-GWAS-Predict/blob/main/Data%20preprocessing/Bedtools_Gene_Annotation.bash)

Bedtools gene annotation involved using the hg19/GRCh37 reference genome from Ensembl (release 92, Homo_sapiens.GRCh37.87) to annotate variants to their closest genes with a cut off of 50kb distance accepted to annotate a variant to that gene.

The annotated genes were further fitered by gene type, selecting only: protein-coding, processed transcripts, pseudogenes, and antisense genes - to ensure genes would no be heavily missing feature data on machine learning.

#### Variant feature annotation and filtering
Variants were annotated in unix using ANNVOAR (2018Apr16). 

The ANNOVAR annotation was completed using ```ANNOVAR_annotation.bash``` found [here](https://github.com/hlnicholls/BP-GWAS-Predict/blob/main/Data%20preprocessing/ANNOVAR_annotation.bash)

Once ANNOVAR variant annotation was complete, the whole GWAS genes were filtered to only selected associated genes and genes found to be insignificant/least likely to affect BP (p-value >0.15, no linkage disequilibrium, not within 500kb+/- BP loci) - filtering out other genes in the whole GWAS data.

#### Least likely BP gene filtering

Genes identified as least likely to affect BP by p-value and linkage disequilibrium were further filtered by PPI distance to known BP genes, creating a final list of least likely genes to be scored at 0.1 on training.

#### Variant-level features
All variant-level annotations per gene had the gene's most significant variant-level value selected to represent that gene, excluding for beta values in the GWAS summary statistics, which had the maximum absolute value taken across all three BP phenotypes (systolic, diastolic and pulse pressure) for each gene.

#### Gene-level features
GTEx foldchange data was the only gene-level annotation that required further processing before merging into a final file, requiring all GTEx tissues to be merged into one file. All other collected feaures required no further processing and could be directly merged into a final file.

#### Training data curation
The final integrated file was then divided into training data and genes to be predicted by the trained model. This was done by identifying three gene groups that would make up the training data:
1. known BP genes (genes with BP drugs and BP mechanisms)
2. Genes that are probable to affect BP (genes with significance in BP publications in text-mining or had drug interactions known to have a BP side effect)
3. Least likely genes (genes fitlered by distance to BP loci, high GWAS p-value, no linkage disequilibrum r2 measured, and no close PPIs with known BP genes)



### 2. Machine learning:

#### Feature pre-processing

Features in the training data were removed if there were >25% missing or had >0.9 correlation in the ```EDA-training-data.ipynb``` script found [here](https://github.com/hlnicholls/BP-GWAS-Predict/blob/main/Machine%20learning/EDA-training-data.ipynb)

Feature imputation was performed using random forest imputation followed by feature selection using BorutaShap in the ```All-Features-ML.ipynb``` script found [here](https://github.com/hlnicholls/BP-GWAS-Predict/blob/main/Machine%20learning/All-Features-ML.ipynb)

#### Machine learning benchmarking
 
Regression classification was used with all models benchmarked using repeated nested 5-fold nested cross validation, with each model undergoing Bayesian hyperparamter tuning (code for benchmarking can be found [here](https://github.com/hlnicholls/BP-GWAS-Predict/blob/main/Machine%20learning/Selected-Features-ML.ipynb)).  Median model performances across the repeated folds were then taken to assess model performance. The top performaning model using all features had its best model parameters applied in BorutaShap feature selection to select features.

After feature selection, the models were then benchmarked against each other in repeated nested cross-validation again, with the top performing model used for further analysis. Final model benchmarking that selected the top performing model can be found [here](https://github.com/hlnicholls/BP-GWAS-Predict/blob/main/Machine%20learning/Model%20Benchmarking/Selected-Features-ML.ipynb).

The models benchmarked (extreme gradient boosting, catboost, lightgbm gradient boosting, random forest, extratrees, decision tree, k-nearest neighbors, LASSO and elasticnet logistic regrssion) were further compared with voting and stacking regressors using all other models and a bagging regressor using extreme gradient boosting within the ```Meta-estimators benchmarking.ipynb``` script found [here](https://github.com/hlnicholls/BP-GWAS-Predict/blob/main/Machine%20learning/Meta-estimators%20benchmarking.ipynb).

The trained top performing model was given non-training data genes to predict their likelihood of affecting BP within the ```XGB-Predictions.ipynb``` script found [here](https://github.com/hlnicholls/BP-GWAS-Predict/blob/main/Machine%20learning/XGB-Predictions.ipynb).

The top performing model was further anlysed using SHAP in the ```XGB-SHAP-interpretation.ipynb``` script found [here](https://github.com/hlnicholls/BP-GWAS-Predict/blob/main/Machine%20learning/XGB-SHAP-interpretation.ipynb).
 
### 3. Prioritization analysis:

#### Gene selection per locus 

After machine learning prioritization the genes were grouped to their loci using the ```Loci-ordering.R``` script found [here](https://github.com/hlnicholls/BP-GWAS-Predict/blob/main/Output%20Prioritisation/Loci%20ordering%20by%20LD.R)

A gene per locus selection strategy was then employed to select the best genes per each locus, which was carried out in the ```Gene-per-locus-selection.R``` script found [here](https://github.com/hlnicholls/BP-GWAS-Predict/blob/main/Output%20Prioritisation/Gene%20per%20locus%20selection.R).

The strategy to select the top gene per locus involved testing each loci's genes under filtering conditions:

1. Machine learning prediction (XGB_Score) is > +1 standard deviation (SD) OR the gene is already a known BP gene (scored 1 in training).<br />
If the conditions in 1. are not met by any of the genes in a locus then:<br />
2. Genes > the average score are selected for further filtering by protein-protein interactions (PPI).<br />
3. Gene with highest direct PPI (direct_PPI_count) with known disease-genes is chosen OR if direct PPI is matching between genes they enter 4.<br />
4. Gene with highest secondary PPI (secondary_PPI_count) with known disease-genes is chosen OR if there are still matching PPI counts genes enter step 5.<br />
5. All genes matching PPIs are chosen


# References
1. Evangelou, E., Warren, H.R., Mosen-Ansorena, D. et al. Genetic analysis of over 1 million people identifies 535 new loci associated with blood pressure traits. Nat Genet 50, 1412–1425 (2018). https://doi.org/10.1038/s41588-018-0205-x

# Contact

To report any issues or for further information, please contact: 

- Hannah Nicholls: hannahnicholls@qmul.ac.uk
