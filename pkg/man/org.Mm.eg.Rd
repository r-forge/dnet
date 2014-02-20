\name{org.Mm.eg}
\alias{org.Mm.eg}

\title{Mouse Entrez Genes (EG).}


\usage{
data(org.Mm.eg)
}

\description{
An R object that contains Entrez Gene information for the mouse. This data is prepared based on \url{ftp://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz}.
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