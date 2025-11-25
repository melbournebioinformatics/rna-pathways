---
title: 'Conclusion'
teaching: 10
exercises: 2
---

:::::::::::::::::::::::::::::::::::::: questions 

- What have we learned about functional enrichment and pathway analysis?
- How do different methods complement one another when interpreting RNA-seq results?


::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Summarise the key concepts introduced across the lesson series.
- Understand how different gene set and network tools fit together in a typical analysis workflow.
- Recognise when and why to choose each enrichment method.

::::::::::::::::::::::::::::::::::::::::::::::::

## Conclusion

In this tutorial, we have explored several complementary approaches for interpreting RNA-seq results beyond differential expression alone. Through using these various R packages, we are able to get insights biological processes and pathways involved in the differential expression of genes observed.

Specifcally, we worked through:

- **Over-representation analysis (ORA)** with `clusterProfiler`  
- **Gene set enrichment analysis (GSEA)** using `fgsea`  
- **Regulatory network analysis** with `RegEnrich`  
- **Protein–protein interaction networks** via `STRINGdb`

Although each tool uses different assumptions and statistical frameworks, they all aim to answer a similar biological question:

> *Which biological processes, pathways, or regulators help explain the gene expression changes we observe?*

By applying multiple methods, you can cross-validate findings and gain a more complete picture of the molecular biology underlying your condition of interest.

You should now feel comfortable:

- preparing gene lists or ranked gene sets  
- running several types of enrichment analyses  
- visualising pathway-level patterns  
- integrating results from complementary tools  
- exploring interaction networks and regulatory drivers  

These approaches form a core part of transcriptomic interpretation and are widely used in modern functional genomics.

::::::::::::::::::::::::::::::::::::: keypoints 

- Enrichment methods help translate gene-level changes into biological meaning.  
- Different tools (ORA, GSEA, network-based methods) answer different but complementary questions.  
- Combining methods provides stronger and more interpretable biological insights.  
- Functional enrichment is an essential component of any RNA-seq analysis workflow.

::::::::::::::::::::::::::::::::::::::::::::::::

