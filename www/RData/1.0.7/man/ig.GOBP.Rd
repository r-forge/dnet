\name{ig.GOBP}
\alias{ig.GOBP}

\title{Gene Ontology Biological Process (GOBP).}


\usage{
ig.GOBP <- dRDataLoader(RData='ig.GOBP')
}

\description{
An R object that contains information on Gene Ontology Biological Process terms. These terms are organised as a direct acyclic graph (DAG), which is further stored as an object of the class 'igraph' (see \url{http://igraph.org/r/doc/aaa-igraph-package.html}). This data is prepared based on \url{http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology.1_2.obo}.
}

\value{
    an object of class "igraph". As a direct graph, it has attributes to vertices/nodes and edges:
    \itemize{
        \item{\code{vertex attributes}: "name" (i.e. "Term ID"), "term_id" (i.e. "Term ID"), "term_name" (i.e "Term Name") and "term_distance" (i.e. Term Distance: the distance to the root; always 0 for the root itself)}
        \item{\code{edge attributes}: "relation" (either 'is_a' or 'part_of')}
    }
}

\references{
Ashburner et al. (2000) Gene ontology: tool for the unification of biology. \emph{Nat Genet}, 25:25-29.
}


\keyword{datasets}

\examples{
ig.GOBP <- dRDataLoader(RData='ig.GOBP')
ig.GOBP
}
