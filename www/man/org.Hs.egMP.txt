\name{org.Hs.egMP}
\alias{org.Hs.egMP}

\title{Annotations of Human Entrez Genes (EG) by Mammalian Phenotype (MP).}


\usage{
org.Hs.egMP <- dRDataLoader(RData='org.Hs.egMP')
}

\description{
An R object that contains associations between Mammalian Phenotype terms and Human Entrez Genes. This data is prepared based on \url{ftp://ftp.informatics.jax.org/pub/reports/MPheno_OBO.ontology} and \url{ftp://ftp.informatics.jax.org/pub/reports/MGI_PhenoGenoMP.rpt}, which results in annotations of Mouse Entrez Genes. Then, these annotations are transferred to Human Entrez Genes based on \url{ftp://anonymous@ftp.ncbi.nih.gov/pub/HomoloGene/build67/homologene.data}.
}

\value{
    an object of class "GS", a list with following components:
    \itemize{
        \item{\code{set_info}: a matrix of nSet X 4 containing gene set information, where nSet is the number of gene sets (i.e. MP terms), and the 4 columns are "setID" (i.e. "Term ID"), "name" (i.e. "Term Name"), "namespace" and "distance"}
        \item{\code{gs}: a list of gene sets, each storing gene members thereof. Always, gene sets are identified by "setID" and gene members identified by "Entrez ID"}
    }
}

\references{
Smith et al. (2009) The Mammalian Phenotype Ontology: enabling robust annotation and comparative analysis. \emph{Wiley Interdiscip Rev Syst Biol Med}, 1:390-399.
\cr
Sayers et al. (2011) Database resources of the National Center for Biotechnology Information. \emph{Nucleic Acids Res}, 39:D38-51.
}

\keyword{datasets}

\examples{
org.Hs.egMP <- dRDataLoader(RData='org.Hs.egMP')
names(org.Hs.egMP)
}
