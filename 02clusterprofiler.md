---
title: 'Gene Ontology testing with `clusterProfiler`'
teaching: 10
exercises: 2
---

:::::::::::::::::::::::::::::::::::::: questions 

- What are the different types of GO terms (BP, MF, CC)?
- How do we perform ORA using `enrichGO()` function?
-	How can we run GSEA-style functional class scoring with `gseGO()` function?


::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Apply GO-based enrichment methods using `clusterProfiler`
- Perform both ORA and GSEA using the GO terms database
- Build confidence in navigating GO resources and interpreting enriched terms 

::::::::::::::::::::::::::::::::::::::::::::::::

## Introduction

The Gene Ontology (GO) project is a major bioinformatics initiative that standardises how we describe gene functions across species, organising them into three categories: Biological Process, Molecular Function and Cellular Component. `clusterProfiler` is an R package that allows us to test whether these GO terms are associated with our RNA-seq results and gain insight into the pathways or functions represented in our data. This section demonstrates how to perform both **over-representation analysis (ORA)** and **functional class scoring (FCS)** with GO database, depending on whether you are working with a list of significant genes or full ranked expression data. 

## Over-Representation Analysis (ORA)

**ORA** tests whether a list of significant genes are linked to specific GO terms. The input is a vector of gene IDs (or list of genes) that passes your differential expression cut-off. ORA can be run separately for downregulated and upregulated genes to reveal which GO terms are enriched in each direction.

We first subset the `debasal` dataset to extract genes with adjusted p-value below 0.01 and store this set of significant genes in an object called `genes`. We then run enrichGO function using this gene list, specifying the organism database `org.Mm.eg.db`, the identifier type `ENTREZID` and the GO category of interest `CC` (for cellular component). The function is configured with standard p-value and q-value, using Benjamini-Hochberg correction. We use the function `head()` to check the first few lines of output. 






``` r
debasal$Status <- debasal$adj.P.Val < 0.01
gene <- debasal$ENTREZID[debasal$Status]

ego <- enrichGO(gene = gene,
                OrgDb = org.Mm.eg.db,
                keyType = 'ENTREZID',
                ont = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.05,
                readable = TRUE)
head(ego)
```

``` output
                   ID                       Description GeneRatio   BgRatio
GO:0030684 GO:0030684                       preribosome   63/2804 105/25804
GO:0022626 GO:0022626                cytosolic ribosome   66/2804 124/25804
GO:0032040 GO:0032040          small-subunit processome   44/2804  73/25804
GO:0005819 GO:0005819                           spindle  119/2804 450/25804
GO:0044391 GO:0044391                 ribosomal subunit   68/2804 185/25804
GO:0022627 GO:0022627 cytosolic small ribosomal subunit   28/2804  38/25804
           RichFactor FoldEnrichment   zScore       pvalue     p.adjust
GO:0030684  0.6000000       5.521541 16.21001 3.722343e-34 2.810369e-31
GO:0022626  0.5322581       4.898141 15.19256 2.316529e-31 8.744898e-29
GO:0032040  0.6027397       5.546753 13.58297 2.200476e-24 5.537865e-22
GO:0005819  0.2644444       2.433568 10.71178 1.051630e-20 1.984952e-18
GO:0044391  0.3675676       3.382565 11.35559 1.354371e-20 2.045101e-18
GO:0022627  0.7368421       6.780839 12.45143 1.432799e-19 1.802938e-17
                 qvalue
GO:0030684 1.754661e-31
GO:0022626 5.459899e-29
GO:0032040 3.457580e-22
GO:0005819 1.239309e-18
GO:0044391 1.276864e-18
GO:0022627 1.125669e-17
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          geneID
GO:0030684                                                                                                                                                                                                                                                                                                                                                                                   Wdr36/Rps15a/Riok2/Rrs1/Nop56/Rcl1/Utp14b/Utp18/Wdr75/Rps16/Riok1/Pno1/Nol6/Nat10/Utp25/Rps11/Rps5/Ltv1/Rrp1b/Rps17/Rps4x/Rps7/Mak16/Srfbp1/Wdr74/Pwp2/Ppan/Rps19/Bysl/Pes1/Wdr46/Rps8/Nip7/Rps9/Rrp15/Nop14/Utp15/Krr1/Mphosph10/Tbl3/Noc4l/Noc2l/Rps27a/Rrp9/Wdr43/Rps19bp1/Rps28/Tsr1/Rps12/Utp4/Rps3a1/Fbl/Heatr1/Ftsj3/Mrto4/Rps23/Prkdc/Dimt1/Nob1/Rps24/Ebna1bp2/Dhx37/Rps14
GO:0022626                                                                                                                                                                                                                                                                                                                                                             Rpl27/Rpl32/Rpl23a/Rplp0/Rps15a/Rpl10a/Rpl3/Rpl12/Rpl19/Rps16/Rps2/Rpl22/Rpl4/Rps20/Rplp2/Rps11/Rps5/Rps10/Rps15/Rpl18a/Rps26/Rps17/Rps4x/Rpl14/Rps7/Rpl7a/Rpl5/Rps29/Rps19/Rpl18/Rpl6/Rps8/Rpl11/Rpl36/Rps9/Rpl10/Rps18/Rpl15/Gspt1/Rpl24/Rpl31/Rps25/Rps27a/Abce1/Rpl21/Etf1/Rpl13a/Rps28/Rps12/Rps21/Rps3a1/Rpl34/Rpl36a/Rpl23/Rpl8/Usp10/Rps3/Rpl35/Rplp1/Rps23/Ppargc1a/Rpl41/Rps24/Rpl26/Rps14/Rpsa
GO:0032040                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 Wdr36/Rps15a/Nop56/Rcl1/Utp14b/Utp18/Wdr75/Rps16/Pno1/Nol6/Nat10/Utp25/Rps11/Rps5/Rps17/Rps4x/Rps7/Pwp2/Rps19/Wdr46/Rps8/Rps9/Nop14/Utp15/Krr1/Mphosph10/Tbl3/Noc4l/Rps27a/Rrp9/Wdr43/Rps19bp1/Rps28/Rps12/Utp4/Rps3a1/Fbl/Heatr1/Rps23/Prkdc/Dimt1/Rps24/Dhx37/Rps14
GO:0005819 Wdr5/Kntc1/Nedd9/Mical3/Espl1/Fam110a/Gsk3b/Ino80/Kif11/Sgo1/Kif22/Cep170/Bub1b/Lzts2/Adrb2/Kif18b/Sirt2/Eml4/Afg2a/Prc1/Knstrn/Champ1/Cdk1/Zzz3/Kat2b/Ska1/Kif20a/Diaph3/Birc5/Unc119/Cdc20/Gpsm2/Luzp1/Arhgef2/Shcbp1/Ccnb1/Pard3/Cspp1/Invs/Hmmr/Gem/Ccsap/Plk1/Cenpf/Cenpe/Dctn1/Spice1/Kifc1/Topors/Nek6/Clasp1/Cep350/Nup62/Anxa11/Cdc27/Ercc2/Rassf10/Tacc3/Pmf1/Dzip1l/Cdk5rap2/Spag5/Cdca8/Dpysl2/Kif2a/Ncor1/Rps6ka2/Hecw2/Kif14/Npm1/Ctdp1/Mtcl1/Ddx11/Racgap1/Kmt5b/Mapk14/Clasp2/Git1/Spdl1/Tppp/Ccdc66/Plekhg6/Ckap2/Arl8a/Parp4/Prpf19/Mad2l1/Ckap5/Mapre3/Rps3/Kif2c/Ikbkg/Kif23/Chmp3/Kif15/Nusap1/Hnrnpu/Ect2/Dlgap5/Aurkb/Cltc/Hspa2/Rangap1/Csnk1d/Slc25a5/Rcc2/Nsun2/Ralbp1/Map4/Tbl1x/Tubb2a/Ppp2cb/Cep63/Nek7/Haus4/Tpx2/Tbccd1/Ckap2l/Tmem201
GO:0044391                                                                                                                                                                                                                                                                                                                                              Rpl27/Rpl32/Rpl23a/Rplp0/Rps15a/Rpl10a/Rpl3/Rpl12/Rpl19/Rps16/Rps2/Rpl22/Rpl4/Rps20/Rplp2/Rps11/Rps5/Rps10/Mrpl17/Rps15/Rpl18a/Rps26/Rps17/Rps4x/Rpl14/Rps7/Rpl7a/Rpl5/Rps29/Rps19/Rpl18/Rpl6/Rack1/Rps8/Rpl11/Rpl36/Rps9/Rpl10/Rps18/Npm1/Rpl15/mt-Rnr2/Rpl24/Rpl31/Rps25/Rps27a/Mrpl12/Rpl21/Rpl13a/Rps28/Rps12/Rps21/Rps3a1/Rpl34/Rpl36a/Rpl23/Rpl8/Rps3/Rpl35/Mrpl52/Rplp1/Rps23/Rpl41/Rps24/Rpl26/Mrps30/Rps14/Rpsa
GO:0022627                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   Rps15a/Rps16/Rps2/Rps20/Rps11/Rps5/Rps10/Rps15/Rps26/Rps17/Rps4x/Rps7/Rps29/Rps19/Rps8/Rps9/Rps18/Rps25/Rps27a/Rps28/Rps12/Rps21/Rps3a1/Rps3/Rps23/Rps24/Rps14/Rpsa
           Count
GO:0030684    63
GO:0022626    66
GO:0032040    44
GO:0005819   119
GO:0044391    68
GO:0022627    28
```

We can then use `dotplot()` function to visualise the results in the form of a dot plot. From the plot below, we can see that GO term cellular component spindle, membrane microdomain and ribosome are top enriched terms.


``` r
dotplot(ego)
```

<img src="fig/02clusterprofiler-rendered-unnamed-chunk-3-1.png" alt="" style="display: block; margin: auto;" />

::::: challenge

Challenge! Can you identify enriched GO term biological process in `deluminal` dataset? Are the enriched pathways similar?

:::::
## Gene Set Enrichment Analysis (GSEA)

We can also perform **GSEA** using GO database. GSEA is a type of functional class scoring method that evaluates whether genes belonging to a GO term tend to appear at the top or bottom of a ranked gene list, rather than relying on a cut-off (i.e. adj.P.Val < 0.01). The input is a continuous ranking metric (e.g. log2FC) for all genes. This allows the detection of subtle but coordinated shifts in GO terms for both downregulated and upregulated pathways. 

We begin by creating a ranked gene list for GSEA by extracting the logFC values from `debasal` dataset and its corresponding `ENTREZID`. We then sort this vector in a decreasing order so that the upregulated genes appear at the top of the list and the downregulated genes at the bottom. Using this ranked gene list, we run `gseGO()` to perform GSEA on GO terms `CC`, by specifying the organism database, gene ID type, gene set limits and p-value cut-off for enrichment. 


``` r
debasal_genelist <- debasal$logFC
names(debasal_genelist) <- debasal$ENTREZID
debasal_genelist <- sort(debasal_genelist, decreasing = TRUE)

ego3 <- gseGO(gene          = debasal_genelist,
                OrgDb         = org.Mm.eg.db,
                keyType       = 'ENTREZID',
                ont           = "CC",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
head(ego3)
```

``` output
                   ID                              Description setSize
GO:0030684 GO:0030684                              preribosome     104
GO:0022626 GO:0022626                       cytosolic ribosome     110
GO:0000776 GO:0000776                              kinetochore     167
GO:0000779 GO:0000779 condensed chromosome, centromeric region     177
GO:0044391 GO:0044391                        ribosomal subunit     169
GO:0000775 GO:0000775           chromosome, centromeric region     242
           enrichmentScore      NES       pvalue     p.adjust       qvalue rank
GO:0030684       0.6665752 2.374098 1.000000e-10 3.700000e-09 3.700000e-09 3377
GO:0022626       0.6280927 2.282534 1.000000e-10 3.700000e-09 3.700000e-09 4038
GO:0000776       0.5753751 2.198080 1.000000e-10 3.700000e-09 3.700000e-09 1254
GO:0000779       0.5636384 2.146015 1.000000e-10 3.700000e-09 3.700000e-09 1254
GO:0044391       0.5552514 2.097483 2.925402e-10 9.019991e-09 9.019991e-09 4724
GO:0000775       0.5286434 2.059075 1.000000e-10 3.700000e-09 3.700000e-09 1417
                             leading_edge
GO:0030684 tags=64%, list=21%, signal=51%
GO:0022626 tags=72%, list=26%, signal=54%
GO:0000776  tags=25%, list=8%, signal=23%
GO:0000779  tags=24%, list=8%, signal=22%
GO:0044391 tags=62%, list=30%, signal=44%
GO:0000775  tags=21%, list=9%, signal=20%
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                core_enrichment
GO:0030684                                                                                                                                                                                                                                    72515/59028/14113/235036/67223/72462/104662/215193/59014/67973/67619/98956/27966/217995/56095/57741/67222/69902/66538/67134/213773/21771/105372/66164/53414/67920/78294/71340/208144/57750/107071/20085/110816/57315/52705/230082/20055/20115/20116/217109/20103/20088/64934/66249/68052/66254/100608/54127/20091/67045/20042/73674/353258/20068/76846/72554/267019/20102/20044/27207/69072/195434/225348/14791/57294/66475/27993
GO:0022626                                                                                                                                                              16785/56040/269261/22224/20084/19982/68436/22186/78294/67186/67097/20085/67671/20055/19951/11837/100503670/20115/27367/20116/54217/27370/76808/20103/270106/268449/20088/19896/67025/68052/20090/75617/24015/20054/27050/54127/26961/67115/67891/67945/114641/22121/19946/20091/19899/20042/66489/67427/66480/66481/65019/19921/20068/19988/19933/76846/267019/20102/20044/27207/19981/19942/14852/19941/57294/66475/19944/225363/27176/57808/16898/19934/110954/68193/11815/67281/207214/105083/319195
GO:0000776                                                                                                                                                                                                                                                                                                                                                                                                       20877/12235/66468/66977/66570/108000/67629/12236/76464/208628/26886/12615/268697/18817/102920/54141/67052/18005/60411/107995/72415/68549/70385/22137/11799/73804/51944/72155/229841/381318/71876/68014/67177/56150/69928/66934/66442/67037/19387/101994/236930
GO:0000779                                                                                                                                                                                                                                                                                                                                                                                                 20877/12235/66468/66977/54392/66570/108000/67629/12236/76464/208628/26886/12615/268697/18817/102920/54141/67052/18005/60411/107995/72415/68549/70385/22137/11799/73804/51944/72155/229841/381318/71876/68014/67177/56150/69928/66934/66442/67037/19387/101994/236930
GO:0044391 16785/56040/269261/20084/56282/19982/66973/68436/22186/78294/67186/67097/20085/67671/20055/18148/19951/11837/100503670/20115/14694/68836/27367/20116/54217/27370/76808/20103/270106/268449/20088/19896/67025/68052/20090/75617/69163/20054/27050/54127/26961/67115/67891/67945/114641/22121/19946/20091/19899/20042/66489/59054/67427/60441/66480/66481/65019/19921/27397/20068/118451/19988/19933/76846/267019/79044/20102/20044/27207/78523/19981/19942/66230/19941/57294/66475/19944/94063/27176/57808/16898/102060/66258/19934/110954/28028/68193/75398/67281/207214/319195/50529/26451/14109/19989/20104/64657/64655/66407/20005/94065/216767/67840/67308/19943
GO:0000775                                                                                                                                                                                                                                                                                                                                 20877/12235/66468/66977/54392/66570/108000/52276/67629/72107/12236/70645/76464/208628/26886/12615/71988/268697/18817/102920/54141/67052/18005/60411/217653/107995/72415/68549/70385/22137/11799/73804/51944/21973/72155/229841/381318/217578/71876/68014/67177/56150/17345/69928/66934/66442/67037/19387/101994/236930/218973/219114
             log2err
GO:0030684       NaN
GO:0022626       NaN
GO:0000776       NaN
GO:0000779       NaN
GO:0044391 0.8140358
GO:0000775       NaN
```

``` r
dotplot(ego3)
```

<img src="fig/02clusterprofiler-rendered-unnamed-chunk-4-1.png" alt="" style="display: block; margin: auto;" />

We can also use the `gseaplot()` function to visualise GSEA result for a specific gene set. In this example, we select the top-ranked enriched GO term (geneSetID = 1). The result-ing plot displays how genes contributing to the enrichment of this GO term are distributed in the ranked gene list. 


``` r
gseaplot(ego3, by = "all", title = ego3$Description[1], geneSetID = 1)
```

<img src="fig/02clusterprofiler-rendered-unnamed-chunk-5-1.png" alt="" style="display: block; margin: auto;" />


::::::::::::::::::::::::::::::::::::: keypoints 

- GO terms are divided into Biological Process (BP), Molecular Function (MF) and Cellular Component (CC), which can be analysed separately or together depending on the biological question.
- The `enrichGO()` and `gseGO()` functions in `clusterProfiler` allow users to perform **ORA** and **GSEA** using the GO database directly.
- GO testing results highlight gene sets or pathways that are overrepresented in your dataset, allowing interpretation of downregulated or upregulated genes.

::::::::::::::::::::::::::::::::::::::::::::::::

