\name{org.Ce.eg}
\alias{org.Ce.eg}

\title{C.elegans Entrez Genes (EG).}


\usage{
org.Ce.eg <- dRDataLoader(RData='org.Ce.eg')
}

\description{
An R object that contains Entrez Gene information for the c.elegans. This data is prepared based on \url{ftp://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz}.
}

\value{
    an object of class "EG", a list with following components:
    \itemize{
        \item{\code{gene_info}: a matrix of nGene X 7 containing gene information, where nGene is the number of Entrez Genes, and the 7 columns are "GeneID", "Symbol", "description", "chromosome", "map_location", "Synonyms" and "dbXrefs"}
    }
}

\references{
Maglott et al. (2011) Entrez Gene: gene-centered information at NCBI. \emph{Nucleic Acids Res}, 39:D52-57.
}

\keyword{datasets}

\examples{
org.Ce.eg <- dRDataLoader(RData='org.Ce.eg')
names(org.Ce.eg)
org.Ce.eg$gene_info[1:5,]
}
