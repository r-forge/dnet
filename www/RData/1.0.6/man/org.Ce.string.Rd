\name{org.Ce.string}
\alias{org.Ce.string}

\title{C.elegans functional protein association network from STRING (version 10).}

\usage{
org.Ce.string <- dRDataLoader(RData='org.Ce.string')
}

\description{
An igraph object that contains a functional protein association network in c.elegans. The network is extracted from the STRING database (version 9.1). Only those associations with medium confidence (score>=0.4) are retained. 
}

\value{
    an object of class "igraph" (see \url{http://igraph.org/r/doc/aaa-igraph-package.html}). It has attributes for both vertices and edges. 
    Below are attributes for the vertices:
    \itemize{
        \item{\code{name}: unique id for the vertices}
        \item{\code{seqid}: protein seqid for the vertices}
        \item{\code{geneid}: Entrez geneid (if any) for the vertices}
        \item{\code{symbol}: gene symbol (if any) for the vertices}
        \item{\code{description}: gene description (if any) for the vertices}        
    }
    Below are attributes for the edges:
    \itemize{
        \item{\code{neighborhood_score}: predictive score based on neighborhood data}
        \item{\code{fusion_score}: predictive score based on fusion data}
        \item{\code{cooccurence_score}: predictive score based on cooccurence data}
        \item{\code{coexpression_score}: predictive score based on coexpression}
        \item{\code{experimental_score}: predictive score based on experimental data}
        \item{\code{database_score}: predictive score based on database}
        \item{\code{textmining_score}: predictive score based on text mining}
        \item{\code{combined_score}: combined score from all above predictive scores}
    }
}

\references{
Szklarczyk et al. (2015) STRING v10: protein-protein interaction networks, integrated over the tree of life. \emph{Nucleic Acids Res}, 43:D447-52.
}

\keyword{datasets}

\examples{
org.Ce.string <- dRDataLoader(RData='org.Ce.string')
org.Ce.string
}
