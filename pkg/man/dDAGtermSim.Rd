\name{dDAGtermSim}
\alias{dDAGtermSim}
\title{Function to calculate pair-wise semantic similarity between input terms based on a direct acyclic graph (DAG) with annotated data}
\usage{
dDAGtermSim(g, terms = NULL, method = c("Resnik", "Lin", "Schlicker",
"Jiang", "Pesquita"), precomputed = TRUE, verbose = TRUE)
}
\arguments{
\item{g}{an object of class "igraph" or "graphNEL"}

\item{terms}{the terms/nodes between which pair-wise
semantic similarity is calculated. If NULL, all terms in
the input dag will be used for calcluation, which is very
prohibitively expensive!}

\item{method}{the method used to measure semantic
similarity between input terms. It can be "Resnik" for
information content (IC) of most informative information
ancestor (MICA) (see
\url{http://arxiv.org/pdf/cmp-lg/9511007.pdf}), "Lin" for
2*IC at MICA divided by the sum of IC at pairs of terms
(see
\url{http://webdocs.cs.ualberta.ca/~lindek/papers/sim.pdf}),
"Schlicker" for weighted version of 'Lin' by the
1-prob(MICA) (see
\url{http://www.ncbi.nlm.nih.gov/pubmed/16776819}),
"Jiang" for 1 - diference between the sum of IC at pairs
of terms and 2*IC at MICA (see
\url{http://arxiv.org/pdf/cmp-lg/9709008.pdf}),
"Pesquita" for graph information content similarity
related to Tanimoto-Jacard index (ie. summed information
content of common ancestors divided by summed information
content of all ancestors of term1 and term2 (see
\url{http://www.ncbi.nlm.nih.gov/pubmed/18460186}))}

\item{precomputed}{logical to indicate whether a sparse
matrix of children x ancestors is pre-computed for the
subsequent identification of common ancestors. By
default, it sets to true. It is advisable to have enabled
it when there are many input terms}

\item{verbose}{logical to indicate whether the messages
will be displayed in the screen. By default, it sets to
true for display}
}
\value{
It returns a sparse matrix containing pair-wise semantic
similarity between input terms
}
\description{
\code{dDAGtermSim} is supposed to calculate pair-wise
semantic similarity between input terms based on a direct
acyclic graph (DAG) with annotated data.
}
\note{
When there are many input terms, 'precomputed' shoulb be
true; otherwise keep it as FALSE for no more than a handful
of input terms
}
\examples{
\dontrun{
# 1) load HPPA as igraph object
data(ig.HPPA)
g <- ig.HPPA

# 2) load human genes annotated by HPPA
data(org.Hs.egHPPA)

# 3) prepare for ontology and its annotation information
dag <- dDAGannotate(g, annotations=org.Hs.egHPPA,
path.mode="all_paths", verbose=TRUE)

# 4) calculate pair-wise semantic similarity between 5 randomly chosen terms
terms <- sample(V(dag)$name, 5)
sim <- dDAGtermSim(g=dag, terms=terms,
method=c("Resnik","Lin","Schlicker","Jiang","Pesquita"),
precomputed=FALSE, verbose=TRUE)
sim
}
}
\seealso{
\code{\link{dDAGinduce}}, \code{\link{dDAGancestor}}
}

