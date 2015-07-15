\name{org.Dm.eg}
\alias{org.Dm.eg}

\title{Fruitfly Entrez Genes (EG).}


\usage{
org.Dm.eg <- dRDataLoader(RData='org.Dm.eg')
}

\description{
An R object that contains Entrez Gene information for the fruitfly. This data is prepared based on \url{ftp://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz}.
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
org.Dm.eg <- dRDataLoader(RData='org.Dm.eg')
names(org.Dm.eg)
org.Dm.eg$gene_info[1:5,]
}
