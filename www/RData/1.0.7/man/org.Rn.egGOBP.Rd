\name{org.Rn.egGOBP}
\alias{org.Rn.egGOBP}

\title{Annotations of Rat Entrez Genes (EG) by Gene Ontology Biological Process (GOBP).}


\usage{
org.Rn.egGOBP <- dRDataLoader(RData='org.Rn.egGOBP')
}

\description{
An R object that contains associations between Gene Ontology Biological Process terms and Rat Entrez Genes. This data is prepared based on \url{http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology.1_2.obo} and \url{ftp://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz}.
}

\value{
    an object of class "GS", a list with following components:
    \itemize{
        \item{\code{set_info}: a matrix of nSet X 4 containing gene set information, where nSet is the number of gene sets (i.e. GOBP terms), and the 4 columns are "setID" (i.e. "Term ID"), "name" (i.e. "Term Name"), "namespace" and "distance"}
        \item{\code{gs}: a list of gene sets, each storing gene members thereof. Always, gene sets are identified by "setID" and gene members identified by "Entrez ID"}
    }
}

\references{
Maglott et al. (2011) Entrez Gene: gene-centered information at NCBI. \emph{Nucleic Acids Res}, 39:D52-57.
\cr
Ashburner et al. (2000) Gene ontology: tool for the unification of biology. \emph{Nat Genet}, 25:25-29.
}


\keyword{datasets}

\examples{
org.Rn.egGOBP <- dRDataLoader(RData='org.Rn.egGOBP')
names(org.Rn.egGOBP)
}
