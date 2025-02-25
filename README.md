# Sputum Neutrophil Elastase associates with microbiota and P.aeruginosa in bronchiectasis 
*Manuscript Analysis*

## Prerequisites
* MacOS/Windows/Linux
* R (version 3.5.1)

## Installing R [typical time = 5-10 minutes]
1. Open an internet browser and go to www.r-project.org.
2. Click the "download R" link in the middle of the page under "Getting Started."
3. Select a CRAN location (a mirror site) and click the corresponding link.
4. Click on the "Download R for [your operating system]" link at the top of the page.
5. Click on the file containing the latest version of R under "Files."
6. Save the .pkg file, double-click it to open, and follow the installation instructions.

## Running the Analysis
* [Human 16S rRNA gene Sequencing](https://github.com/segalmicrobiomelab/bronchiectasis.neutrophil.elastase/tree/master/Microbiome.Analysis)

## For the Human 16S rRNA Gene Sequencing
The Required Input Data
* Microbe Abundance (no-miss-table-dada2.qza):An abundance table with microbial species in rows and samples in columns.
* Taxonomy Table (taxonomy.qza): Table of the taxonomic annotation
* MetaData (BX_Sputum_convert.txt): A spreadsheet with samples in rows and metadata in columns
* Tree File (rooted-tree_quality.qza): 

1. Open R
1. Follow Code Bronchiectasis.Neutrophil.Elastase.r 
  1. At the start of the code you will need to load packages required for analysis, for example:
```
library(vegan)
```
* If you don't have this package you will need to install it, for example:
```
install.pakcages("vegan")
```
