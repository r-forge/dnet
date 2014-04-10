<a href="index.html"><IMG src="dnet_logo1.png" height="100px" id="logo"></a>

## Introduction

`dnet` is an open-source R package for integrative analysis of high-throughput digitised data (e.g. expression, replication timing and mutational profiles) in terms of `network`, `ontology` and `evolution`. 

`dnet` intends to analyse the biological network whose nodes/genes are associated with dynamic information such as expression levels across samples, for identifying gene-active dynamic subnetworks.

`dnet` also supports enrichment analysis using a variety of ontologies (their annotations) and using gene phylostratific age information. 

`dnet` is able to conduct analyses in human but also common model organisms: mouse, rat, chicken, c.elegans, fruitfly, zebrafish and arabidopsis.

`dnet` aims to deliver a tool with rich visuals but less inputs.

## Features

* Identification of gene-active dynamic networks
* Network-based sample classifications and visualisations on 2D sample landscape
* Random walk with restart for network affinity calculation
* Enrichment analysis (Fisher's exact test, Hypergeometric test, Binomial test or KS-like test) that can account for the hierarchy of the ontology
* A wide variety of data supported: ontologies (including `Gene Ontology`, `Disease Ontology`, `Human Phenotype` and `Mammalian Phenotype`), phylostratific age information and gene association networks in many organisms. For details, please refer to [RData](http://dnet.r-forge.r-project.org/rdata.html).
* This package can run on `Windows`, `Mac` and `Linux`

## User I/O

<a href="javascript:newWin('dnet_userIO.png', 'dnet_workflow.png', '800', '800')" title="Click to enlarge"><img style="max-width:95%;border:1px solid #0000FF;box-shadow:5px 5px 2px #C0C0FF;" src='dnet_userIO.png', width="400px" /></a>
