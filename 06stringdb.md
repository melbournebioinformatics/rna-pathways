---
title: 'Interaction networks with `StringDB`'
teaching: 10
exercises: 2
---

:::::::::::::::::::::::::::::::::::::: questions 

- How can we use `STRINGdb` to visualise protein–protein interaction networks for our DE genes?  
- How do we map our gene identifiers to the IDs used by STRING?  
- What information does STRING functional enrichment add beyond standard GO/KEGG analysis?  

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Load and initialise the `STRINGdb` object for mouse.  
- Map a set of differentially expressed genes to STRING identifiers.  
- Visualise a protein–protein interaction network for top DE genes.  
- Retrieve and inspect STRING functional enrichment results.  

::::::::::::::::::::::::::::::::::::::::::::::::


``` error
Error in `library()`:
! there is no package called 'edgeR'
```

``` error
Error in `library()`:
! there is no package called 'goseq'
```

``` error
Error in `library()`:
! there is no package called 'fgsea'
```

``` error
Error in `library()`:
! there is no package called 'EGSEA'
```

``` error
Error in `library()`:
! there is no package called 'clusterProfiler'
```

``` error
Error in `library()`:
! there is no package called 'org.Mm.eg.db'
```

``` error
Error in `library()`:
! there is no package called 'ggplot2'
```

``` error
Error in `library()`:
! there is no package called 'enrichplot'
```

``` error
Error in `library()`:
! there is no package called 'pathview'
```

``` error
Error in `library()`:
! there is no package called 'edgeR'
```

``` error
Error in `library()`:
! there is no package called 'impute'
```

``` error
Error in `library()`:
! there is no package called 'preprocessCore'
```

``` error
Error in `library()`:
! there is no package called 'RegEnrich'
```

``` error
Error in `library()`:
! there is no package called 'STRINGdb'
```


## Interaction networks with `StringDB`
So far, we have focused on pathway-level enrichment. Another useful way to interpret RNA-seq results is to look at ***protein–protein interaction (PPI) networks***:
Are our differentially expressed genes part of the same complexes or signalling modules?

The STRINGdb package provides an interface to the STRING database, which aggregates known and predicted PPIs from multiple sources (experiments, databases, text-mining, etc.).

In this lesson we will:

- Initialise a STRINGdb object for mouse.
- Map our top differentially expressed genes to STRING IDs.
- Plot an interaction network.
- Retrieve functional enrichment results from STRING.


``` r
# Initialize STRINGdb for mouse (taxonomy ID: 10090)
string_db <- STRINGdb$new(version = "12", species = 10090, score_threshold = 400, input_directory = "")
```

``` error
Error:
! object 'STRINGdb' not found
```

``` r
# Prepare data: select top 200 DE genes (by adjusted P value)
top200 <- debasal[order(debasal$adj.P.Val), ][1:200, ]
top200_mapped <- string_db$map(top200, "ENTREZID", removeUnmappedRows = TRUE)
```

``` error
Error:
! object 'string_db' not found
```

``` r
# Plot the protein interaction network
string_db$plot_network(top200_mapped$STRING_id)
```

``` error
Error:
! object 'string_db' not found
```

``` r
# Get functional enrichment (GO, KEGG, Reactome)
enrichment <- string_db$get_enrichment(top200_mapped$STRING_id)
```

``` error
Error:
! object 'string_db' not found
```

``` r
head(enrichment)
```

``` error
Error:
! object 'enrichment' not found
```

There are many available ways of exploring your data using the STRING database that can't be covered in one tutorial but you can learn more by reading the [vignette](https://www.bioconductor.org/packages/release/bioc/vignettes/STRINGdb/inst/doc/STRINGdb.pdf) and inspect available functions within the `STRINGdb` package by running:


``` r
STRINGdb$methods()
```

``` error
Error:
! object 'STRINGdb' not found
```

Read more about STRING:

- Szklarczyk, Damian et al. “The STRING database in 2025: protein networks with directionality of regulation.” Nucleic acids research vol. 53,D1 (2025): D730-D737. doi:10.1093/nar/gkae1113



::::::::::::::::::::::::::::::::::::: keypoints 

- `STRINGdb` links your genes to protein–protein interaction networks from the STRING database.

- Mapping from gene IDs (e.g. ENTREZ) to STRING IDs is a crucial first step.

- Network visualisation can reveal modules of interconnected DE genes that may not be obvious from lists or tables.

- STRING provides its own functional enrichment, which can complement results from `clusterProfiler` and `fgsea`.

::::::::::::::::::::::::::::::::::::::::::::::::

