---
title: Setup
---

This workshop provides a practical introduction to functional enrichment analysis following differential expression in RNA-seq studies. 
We will compare two major enrichment strategies, **over-representation analysis (ORA)** and **functional class scoring (FCS)**, and discuss when each approach is most appropriate. Participants will learn how to implement these methods in R using packages including `clusterProfiler`, `fgsea`, `Reg-Enrich` and `STRINGdb`, drawing on pathway and gene-set resources such as **Gene Ontology**, **KEGG Pathway Database** and **Molecular Signatures Database**. By the end of this workshop, you will have a clear understanding of how to interpret enriched pathways in RNA-seq data. 

This workshop is largely based on the Galaxy-based workshop [RNA-seq genes to pathways](https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/rna-seq-genes-to-pathways/tutorial.html) on the Galaxy Training Network. 

:::: prereq

-	Installed R and RStudio

-	Have basic R knowledge

-	Completed ‘Intro to R for Biologists’ and ‘RNA-seq: From reads to counts to genes’, or equivalent

::::

## R Packages & Datasets

In this workshop, we will learn how to use `clusterProfiler`, `fgsea`, `RegEnrich` and `STRINGdb` tools,
along with related dependencies `org.Mm.eg.db`, `impute` and `preprocessCore`. 

Please install the following packages:

```r

## Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

## List of Bioconductor packages
bioc_packages <- c(
    "edgeR",
    "goseq",
    "fgsea",
    "EGSEA",
    "clusterProfiler",
    "org.Mm.eg.db",
    "enrichplot",
    "pathview",
    "preprocessCore",
    "RegEnrich",
    "STRINGdb"
)

## Install Bioconductor packages
BiocManager::install(bioc_packages, ask = FALSE, update = TRUE)

## Install CRAN packages
cran_packages <- c(
    "ggplot2",
    "impute"
)

install.packages(cran_packages)


```



The original dataset is from [Fu et al., 2015](<https://www.ncbi.nlm.nih.gov/pubmed/25730472>) and described in the
[this tutorial](https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/rna-seq-genes-to-pathways/tutorial.html).

We have provided the data we will be using within the GitHub repo. To access them:

- [limma-voom_basalpregnant-basallactate](episodes/data/limma-voom_basalpregnant-basallactate)
- [limma-voom_luminalpregnant-luminallactate](episodes/data/limma-voom_luminalpregnant-luminallactate)
- [seqdata](episodes/data/seqdata)
- [mouse_hallmark_sets](episodes/data/mouse_hallmark_sets.RData)
- [factordata](episodes/data/factordata)
- [filteredcounts](episodes/data/limma-voom_filtered_counts)
- [mouseTFs](episodes/data/BrowseTF_TcoF-DB.csv)

The code assumes that these files are in a folder called "data" but you can adjust the code to the correct download path as needed.


You can load some of the data directly from Zenodo into your RStudio enviornment with the following code:

```r

# To download files from Zenodo

dataurl <- "https://zenodo.org/record/2596382/files/"

debasal <- read.csv(paste0(dataurl,"limma-voom_basalpregnant-basallactate"), header = TRUE, sep = "\t")
deluminal <- read.csv(paste0(dataurl,"limma-voom_luminalpregnant-luminallactate"), header = TRUE, sep = "\t")
seqdata <- read.csv(paste0(dataurl,"seqdata"), header = TRUE, sep = "\t")
load(paste0(dataurl,"mouse_hallmark_sets")) #loads as Mm.H
factordata <- read.table(paste0(dataurl,"factordata"), header = TRUE, sep = "\t")
filteredcounts <- read.csv(paste0(dataurl,"limma-voom_filtered_counts"), header = TRUE, sep = "\t")

```

## Summary

:::: checklist

Attendees are required to bring their own laptop computers. **Please ensure you have installed:**

- [Chrome](https://www.google.com/chrome/) or [FireFox](https://www.mozilla.org/en-US/)
- [R](https://cran.ms.unimelb.edu.au/) (Download and install the latest version of R using the UniMelb mirror)
- [RStudio](https://posit.co/download/rstudio-desktop/#download)
- R packages required for this workshop (see above)
- Datasets required for this workshop (see above)

::::


## Citations

- Maria Doyle, Belinda Phipson, 3: RNA-seq genes to pathways (Galaxy Training Materials). https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/rna-seq-genes-to-pathways/tutorial.html Online; accessed Tue Nov 25 2025

- Hiltemann, Saskia, Rasche, Helena et al., 2023 Galaxy Training: A Powerful Framework for Teaching! PLOS Computational Biology 10.1371/journal.pcbi.1010752

- Batut et al., 2018 Community-Driven Data Analysis Training for Biology Cell Systems 10.1016/j.cels.2018.05.012
