\name{dRWRpipeline}
\alias{dRWRpipeline}
\title{Function to setup a pipeine to estimate RWR-based contact strength between samples from an input gene-sample data matrix and an input graph}
\usage{
dRWRpipeline(data, g, method = c("direct", "indirect"),
normalise = c("laplacian", "row", "column", "none"), restart = 0.5,
normalise.affinity.matrix = c("none", "quantile"),
permutation = c("degree", "none"), num.permutation = 10,
p.adjust.method = c("BH", "BY", "bonferroni", "holm", "hochberg",
"hommel"),
adjp.cutoff = 0.05, verbose = T)
}
\arguments{
\item{data}{an input gene-sample data matrix used for
seeds}

\item{g}{an object of class "igraph" or "graphNEL"}

\item{method}{the method used to calculate RWR. It can be
'direct' for directly applying RWR, 'indirect' for
indirectly applying RWR (first pre-compute affinity
matrix and then derive the affinity score)}

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

\item{permutation}{how to do permutation. It can be
'degree' for degree-preserving permutation, 'none' for
permutation in random}

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
components: \itemize{ \item{\code{ratio}: a symmetric
matrix storing ratio (the observed against the expected)
between pairwise samples} \item{\code{zscore}: a symmetric
matrix storing zscore between pairwise samples}
\item{\code{pval}: a symmetric matrix storing pvalue
between pairwise samples} \item{\code{adjpval}: a symmetric
matrix storing adjusted pvalue between pairwise samples}
\item{\code{cgraph}: the constructed contact graph (as a
'igraph' object) under the cutoff of adjusted value}
\item{\code{call}: the call that produced this result} }
}
\description{
\code{dRWRpipeline} is supposed to estimate sample
relationships (ie. contact strength between samples) from
an input gene-sample matrix and an input graph. The
pipeline includes: 1) random walk restart (RWR) of the
input graph using the input matrix as seeds; 2) calculation
of contact strength (inner products of RWR-smoothed columns
of input matrix); 3) estimation of the contact signficance
by a randomalisation procedure. It supports two methods how
to use RWR: 'direct' for directly applying RWR in the given
seeds; 'indirectly' for first pre-computing affinity matrix
of the input graph, and then deriving the affinity score.
}
\note{
The choice of which method to use RWR depends on the number
of seed sets and the number of permutations for statistical
test. If the total product of both numbers are huge, it is
better to use 'indrect' method (for a single run). However,
if the user wants to re-use pre-computed affinity matrix
(ie. re-use the input graph a lot), then it is highly
recommended to sequentially use \code{\link{dRWR}} and
\code{\link{dRWRcontact}} instead.
}
\examples{
# 1) generate a random graph according to the ER model
g <- erdos.renyi.game(100, 1/100)

# 2) produce the induced subgraph only based on the nodes in query
subg <- dNetInduce(g, V(g), knn=0)
V(subg)$name <- 1:vcount(subg)

# 3) estimate RWR dating based sample relationships
# define sets of seeds as data
# each seed with equal weight (i.e. all non-zero entries are '1')
aSeeds <- c(1,0,1,0,1)
bSeeds <- c(0,0,1,0,1)
data <- data.frame(aSeeds,bSeeds)
rownames(data) <- 1:5
# calcualte their two contact graph
dContact <- dRWRpipeline(data=data, g=subg)
dContact
}
\seealso{
\code{\link{dRWR}}, \code{\link{dRWRcontact}}
}

