\name{org.Gg.eg}
\alias{org.Gg.eg}

\title{Chicken Entrez Genes (EG).}


\usage{
org.Gg.eg <- dRDataLoader(RData='org.Gg.eg')
}

\description{
An R object that contains Entrez Gene information for the chicken. This data is prepared based on \url{ftp://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz}.
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
org.Gg.eg <- dRDataLoader(RData='org.Gg.eg')
names(org.Gg.eg)
org.Gg.eg$gene_info[1:5,]
}
