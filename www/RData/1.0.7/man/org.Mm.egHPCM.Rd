\name{org.Mm.egHPCM}
\alias{org.Mm.egHPCM}

\title{Annotations of Mouse Entrez Genes (EG) by Human Phenotype Clinical Modifier (HPCM).}


\usage{
org.Mm.egHPCM <- dRDataLoader(RData='org.Mm.egHPCM')
}

\description{
An R object that contains associations between HPCM terms and Mouse Entrez Genes. This data is first prepared based on \url{http://purl.obolibrary.org/obo/hp.obo} and  \url{http://compbio.charite.de/hudson/job/hpo.annotations.monthly/lastStableBuild/artifact/annotation/ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt}, which results in annotations of Human Entrez Genes. Then, these annotations are transferred to Mouse Entrez Genes based on \url{ftp://anonymous@ftp.ncbi.nih.gov/pub/HomoloGene/build67/homologene.data}.
}

\value{
    an object of class "GS", a list with following components:
    \itemize{
        \item{\code{set_info}: a matrix of nSet X 4 containing gene set information, where nSet is the number of gene sets (i.e. HPCM terms), and the 4 columns are "setID" (i.e. "Term ID"), "name" (i.e. "Term Name"), "namespace" and "distance"}
        \item{\code{gs}: a list of gene sets, each storing gene members thereof. Always, gene sets are identified by "setID" and gene members identified by "Entrez ID"}
    }
}

\references{
Robinson et al. (2012) The Human Phenotype Ontology: a tool for annotating and analyzing human hereditary disease. \emph{Am J Hum Genet}, 83:610-615.
\cr
Sayers et al. (2011) Database resources of the National Center for Biotechnology Information. \emph{Nucleic Acids Res}, 39:D38-51.
}

\keyword{datasets}

\examples{
org.Mm.egHPCM <- dRDataLoader(RData='org.Mm.egHPCM')
names(org.Mm.egHPCM)
}
