This directory (http://dnet.r-forge.r-project.org/data) contains curated/compiled data as part of an open-source R package 'dnet'. 

The data is stored as R data files (ended with '.RData' extension), and can be broadly grouped into the following four categories:

1. Ontologies
They locate in the subdirectory 'Obo', and contain information on ontology terms. These terms are organised as a direct acyclic graph (DAG), which is further stored as an object of the class 'igraph' (see http://igraph.sourceforge.net/doc/R/aaa-igraph-package.html). Ontologies include Gene Ontology (GO) and its three subontologies (BP: Biological Process; MF: Molecular Function; CC: Cellular Component), Human Phenotype (HP) and its three subontologies (PA: Phenotypic Abnormality; ON: ONset and clinical course; MI: Mode of Inheritance), Disease Ontology (DO), and Mammalian Phenotype (MP).

2. Databases in model organisms
Model organisms include human (Hs; tax_id=10090), mouse (Mm; tax_id=10090), arabidopsis (At; tax_id=3702), c.elegans (Ce; tax_id=6239), fruitfly (Dm; tax_id=7227), zebrafish (Da; tax_id=7955), rat (Rn; tax_id=10116) and chicken (Gg; tax_id=9031). They can be found in the subdirectory entitled with 2-letter species identifiers. For example, databases in human are available at http://dnet.r-forge.r-project.org/data/Hs, containing information on Entrez Genes, their annotations by ontologies, their phylostratific age (PS), and their interacting network derived from STRING.

3. Genesets in human
These genesets are derived from the molecular signatures database (Msigdb; http://www.broadinstitute.org/gsea/msigdb/index.jsp). They locate in the subdirectory 'Msigdb'.

4. Datasets used as demos
They locate in the subdirectory 'Datasets'.


For details, please refer to documentations at http://dnet.r-forge.r-project.org/docs.html


Fang H, Gough J. (2014) DNET: dynamic networks via an integrative analysis of network, expression, evolution and ontology data. R package version 1.0.0. http://dnet.r-forge.r-project.org