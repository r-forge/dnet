\name{dRWRdating}
\alias{dRWRdating}
\title{Function to estimate RWR dating based sample relationships from the input gene-sample matrix and graph}
\usage{
dRWRdating(data, g, normalise = c("laplacian", "row", "column",
"none"),
restart = 0.5, normalise.affinity.matrix = c("none", "quantile"),
num.permutation = 10, p.adjust.method = c("BH", "BY", "bonferroni",
"holm", "hochberg", "hommel"), adjp.cutoff = 0.05, verbose = T)
}
\arguments{
\item{data}{an input gene-sample data matrix used for
seeds}

\item{g}{an object of class "igraph" or "graphNEL"}

\item{normalise}{the way to normalise the adjacency
matrix of the input graph. It can be 'laplacian' for
laplacian normalisation, 'row' for row-wise
normalisation, 'column' for column-wise normalisation, or
'none'}

\item{restart}{the restart probability used for RWR. The
restart probability takes the value from 0 to 1,
controlling the range from the starting nodes/seeds that
the walker will explore. The higher the value, the more
likely the walker is to visit the nodes centered on the
starting nodes. At the extreme when the restart
probability is zero, the walker moves freely to the
neighbors at each step without restarting from seeds,
i.e., following a random walk (RW)}

\item{normalise.affinity.matrix}{the way to normalise the
output affinity matrix. It can be 'none' for no
normalisation, 'quantile' for quantile normalisation to
ensure that columns (if multiple) of the output affinity
matrix have the same quantiles}

\item{num.permutation}{the number of permutations used to
for generating the distribution of contact strength under
randomalisation}

\item{p.adjust.method}{the method used to adjust
p-values. It can be one of "BH", "BY", "bonferroni",
"holm", "hochberg" and "hommel". The first two methods
"BH" (widely used) and "BY" control the false discovery
rate (FDR: the expected proportion of false discoveries
amongst the rejected hypotheses); the last four methods
"bonferroni", "holm", "hochberg" and "hommel" are
designed to give strong control of the family-wise error
rate (FWER). Notes: FDR is a less stringent condition
than FWER}

\item{adjp.cutoff}{the cutoff of adjusted pvalue to
construct the contact graph}

\item{verbose}{logical to indicate whether the messages
will be displayed in the screen. By default, it sets to
true for display}
}
\value{
an object of class "dContact", a list with following
components: \itemize{ \item{\code{zscore}: a symmetric
matrix storing zscore between pairwise samples}
\item{\code{pval}: a symmetric matrix storing pvalue
between pairwise samples} \item{\code{adjpval}: a symmetric
matrix storing adjusted pvalue between pairwise samples}
\item{\code{cgraph}: the constructed contact graph (as a
'igraph' object) under the cutoff of adjusted value}
\item{\code{call}: the call that produced this result} }
}
\description{
\code{dRWRdating} is supposed to estmate sample
relationships from the given gene-sample matrix and graph
according to random walk restart (RWR) dating based method.
The method includes: 1) RWR in the graph using the input
matrix as seeds; 2) calculation of dating based contact
strength (inner products of RWR-smooth columns of intput
matrix); 3) estimation of the contact signficance by a
randomalisation procedure.
}
\note{
none
}
\examples{
# 1) generate a random graph according to the ER model
g <- erdos.renyi.game(100, 1/100)

# 2) produce the induced subgraph only based on the nodes in query
subg <- dNetInduce(g, V(g), knn=0)
V(subg)$name <- 1:vcount(subg)

# 3) obtain the pre-computated affinity matrix
PTmatrix <- dRWR(subg, normalise="laplacian", restart=0.75)
# visualise affinity matrix
visHeatmapAdv(PTmatrix, Rowv=FALSE, Colv=FALSE, colormap="wyr",
KeyValueName="Affinity")

# 3) obtain affinity matrix given sets of seeds
# define sets of seeds as data
# each seed with equal weight (i.e. all non-zero entries are '1')
aSeeds <- c(1,0,1,0,1)
bSeeds <- c(0,0,1,0,1)
data <- data.frame(aSeeds,bSeeds)
rownames(data) <- 1:5
# calcualte their two contact graph
dContact <- dRWRdating(data=data, g=subg)
dContact
}
\seealso{
\code{\link{dRWR}}
}

