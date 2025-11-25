---
title: 'Introduction'
teaching: 10
exercises: 2
---

:::::::::::::::::::::::::::::::::::::: questions 

- What are the main types of functional enrichment analysis approaches, and how do they differ?
- When should you choose one enrichment strategy over another for RNA-seq data?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Understand the conceptual differences between over-representation analysis (ORA) and functional class scoring (FCS)
- Learn how enrichment tools (e.g. `clusterProfiler`, `fgsea`, `RegEnrich` and `STRINGdb`) implement these approaches using pathway and gene-set databases 

::::::::::::::::::::::::::::::::::::::::::::::::

## Introduction

Sometimes, there is an extensive list of genes to interpret after differential gene expres-sion analysis, and it is not feasible to go through the biological function of each gene one at a time. A common downstream procedure is functional enrichment analysis (or gene set testing), which aims to determine which pathways or gene networks the differ-entially expressed genes are implicated in. There are many gene set testing methods available, and it is useful to try several of them.

The purpose of this tutorial is to demonstrate how to perform functional enrichment analysis/gene set testing using various tools/packages in R. We will use data from the Nature Cell Biology paper, EGF-mediated induction of Mcl-1 at the switch to lactation is essential for alveolar cell survival (<https://www.ncbi.nlm.nih.gov/pubmed/25730472>). This study examined the expression profiles of basal and luminal cells in the mammary gland of virgin, pregnant and lactating mice.


## Load and read required libraries

We begin by loading the required packages. Please read the following libraries:


``` r
library(edgeR)
library(goseq)
library(fgsea)
library(EGSEA)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(enrichplot)
library(pathview)
library(edgeR)
library(impute)
library(preprocessCore)
library(RegEnrich)
```
## Inspect Datasets

We will use several files for this workshop:

- Results from differential expression analysis `debasal` and `deluminal` with genes in rows and logFC/p-values in columns 
- Sample information file `factordata` – gives details of sample ID and groups
-	Gene lengths file `seqdata`
-	Filtered counts file `filteredcounts` – genes in rows and counts for each sample in columns, lowly expressed genes removed
-	Hallmarks gene set file for mouse from MSigDB loaded in .RData format – `Mm.H`

Let's inspect the files:


``` r
debasal <- read.csv("data/limma-voom_basalpregnant-basallactate", header = TRUE, sep = "\t")
deluminal <- read.csv("data/limma-voom_luminalpregnant-luminallactate", header = TRUE, sep = "\t")
factordata <- read.table("data/factordata", header = TRUE, sep = "\t")

#To view the first 5 rows of the dataset
head(debasal)
```

``` output
  ENTREZID        SYMBOL                                     GENENAME     logFC
1    24117          Wif1                      Wnt inhibitory factor 1  1.819943
2   381290        Atp2b4 ATPase, Ca++ transporting, plasma membrane 4 -2.143885
3   226101          Myof                                    myoferlin -2.329744
4    16012        Igfbp6 insulin-like growth factor binding protein 6 -2.896115
5   231830       Micall2                                 MICAL-like 2  2.253400
6    78896 1500015O10Rik                   RIKEN cDNA 1500015O10 gene  2.807548
   AveExpr         t      P.Value    adj.P.Val        B
1 2.975545  19.85403 5.722034e-11 5.366685e-07 15.55490
2 3.944066 -19.07173 9.406224e-11 5.366685e-07 15.09463
3 6.223525 -18.30281 1.562524e-10 5.366685e-07 14.55585
4 1.978449 -18.21558 1.657202e-10 5.366685e-07 14.13954
5 4.760597  18.00994 1.905713e-10 5.366685e-07 14.33472
6 3.036519  18.60321 2.037466e-10 5.366685e-07 14.35640
```

You can also view the entire file in a different tab using `View()`:


``` r
View(debasal)
```


::::: challenge

-	How many columns are there in `debasal` and `deluminal` objects?

-	What are the different types of samples in this analysis? Hint: Look at `factordata` file.

:::::

## Summary

::::::::::::::::::::::::::::::::::::: keypoints 

## Commonly used analyses following differenital gene expression (DGE)

-	Over-representation analysis (ORA): Tests whether DGE list contains more genes from a specific pathway or gene set

-	Functional class scoring (FCS): Evaluates coordinated shifts in expression across all gene sets 

-	Protein-protein interactions (PPI): Maps the functional connections between proteins to reveal network structure or pathways involved

::::::::::::::::::::::::::::::::::::::::::::::::


