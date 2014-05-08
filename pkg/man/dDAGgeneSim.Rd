\name{dDAGgeneSim}
\alias{dDAGgeneSim}
\title{Function to calculate pair-wise semantic similarity between input terms based on a direct acyclic graph (DAG) with annotated data}
\usage{
dDAGgeneSim(g, genes = NULL, method.gene = c("average", "max",
"BM.average",
"BM.max", "BM.complete")[3], method.term = c("Resnik", "Lin",
"Schlicker",
"Jiang", "Pesquita")[3], force = TRUE, verbose = TRUE)
}
\arguments{
\item{g}{an object of class "igraph" or "graphNEL"}

\item{genes}{the genes between which pair-wise semantic
similarity is calculated. If NULL, all genes annotatable
in the input dag will be used for calcluation, which is
very prohibitively expensive!}

\item{method.gene}{the method used for how to derive
semantic similarity between genes from semantic
similarity between terms. It can be "average" for average
similarity between any two terms (one from gene 1, the
other from gene 2), "max" for the maximum similarity
between any two terms, "BM.average" for best-matching
(BM) based average similarity (i.e. for each gene in a
gene pair, first calculate maximum similarity for each
term, then their average of maximum similarity; the final
BM-based average similiary is the pre-calculated average
between two genes in pair), "BM.max" for BM based maximum
similarity (i.e. for each gene in a gene pair, first
calculate maximum similarity for each term, then their
average of maximum similarity; the final BM-based maximum
similiary is the maximum of the pre-calculated average
between two genes in pair), "BM.complete" for BM based
complete-linkage similarity (inspired by complete-linkage
procedure: for each gene in a gene pair, first calculate
maximum similarity for each term; then take the mimumum
of pre-calculated maximum similarity between two genes in
pair). By default, it sets "BM.average"}

\item{method.term}{the method used to measure semantic
similarity between terms. It can be "Resnik" for
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

\item{force}{logical to indicate whether the only most
specific terms (for each gene) will be used. By default,
it sets to true. It is always advisable to use this since
it is computationally fast but without compromising
accuracy (considering the fact that true-path-rule has
been applied when running \code{\link{dDAGannotate}}).}

\item{verbose}{logical to indicate whether the messages
will be displayed in the screen. By default, it sets to
true for display}
}
\value{
It returns a sparse matrix containing pair-wise semantic
similarity between input terms
}
\description{
\code{dDAGgeneSim} is supposed to calculate pair-wise
semantic similarity between input terms based on a direct
acyclic graph (DAG) with annotated data.
}
\note{
For the mode "shortest_paths", the induced subgraph is the
most concise, and thus informative for visualisation when
there are many nodes in query, while the mode "all_paths"
results in the complete subgraph.
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
allgenes <- unique(unlist(V(dag)$annotations))
genes <- sample(allgenes,5)
sim <- dDAGgeneSim(g=dag, genes=genes,
method.gene=c("average","max","BM.average","BM.max","BM.complete")[3],
method.term=c("Resnik","Lin","Schlicker","Jiang","Pesquita")[3],
verbose=TRUE)
sim
}
}
\seealso{
\code{\link{dDAGtermSim}}, \code{\link{dDAGinduce}},
\code{\link{dDAGtip}}
}

