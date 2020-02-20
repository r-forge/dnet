\name{org.Mm.egMP}
\alias{org.Mm.egMP}

\title{Annotations of Mouse Entrez Genes (EG) by Mammalian Phenotype (MP).}


\usage{
org.Mm.egMP <- dRDataLoader(RData='org.Mm.egMP')
}

\description{
An R object that contains associations between Mammalian Phenotype terms and Mouse Entrez Genes. This data is prepared based on \url{ftp://ftp.informatics.jax.org/pub/reports/MPheno_OBO.ontology} and \url{ftp://ftp.informatics.jax.org/pub/reports/MGI_PhenoGenoMP.rpt}.
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
}

\keyword{datasets}

\examples{
org.Mm.egMP <- dRDataLoader(RData='org.Mm.egMP')
names(org.Mm.egMP)
}
