\name{dRWR}
\alias{dRWR}
\title{Function to implement Random Walk with Restart (RWR) to pre-compute affinity matrix for the input graph}
\usage{
  dRWR(g,
    normalise = c("laplacian", "row", "column", "none"),
    restart = 0.75)
}
\arguments{
  \item{g}{an object of class "igraph" or "graphNEL"}

  \item{normalise}{the way to normalise the adjacency
  matrix of the input graph. It can be 'laplacian' for
  laplacian normalisation, 'row' for row-wise
  normalisation, 'column' for column-wise normalisation, or
  'none'}

  \item{restart}{the restart probability used for RWR}
}
\value{
  \itemize{ \item{\code{PTmatrix}: affinity matrix with the
  dimension of n X n, where n is the number of nodes in the
  input graph. Columns stand for starting nodes walking
  from, and rows for ending nodes walking to. Therefore, a
  column for a starting node represents a steady-state
  affinity vector that the starting node will visit all the
  ending nodes in the graph} }
}
\description{
  \code{dRWR} is supposed to implement Random Walk with
  Restart (RWR) to pre-compute affinity matrix for for
  nodes in the input graph with respect to the starting
  node (loop over every node in the graph)
}
\note{
  The input graph will treat as an unweighted graph if
  there is no 'weight' edge attribute assocaited
}
\examples{
# 1) generate a random graph according to the ER model
g <- erdos.renyi.game(100, 1/100)

# 2) produce the induced subgraph only based on the nodes in query
subg <- dNetInduce(g, V(g), knn=0)

# 3) calculate the affinity matrix
PTmatrix <- dRWR(subg, normalise="laplacian", restart=0.75)

# 4) visualise affinity matrix
graphics::image(PTmatrix, col=visColormap("wyr")(64), zlim=c(0,1))
}
\seealso{
  \code{\link{dNetInduce}}
}

