\name{afCEC}
\alias{afCEC}

\title{
afCEC
}

\description{
Performs the Active Function Cross-Entropy Clustering on the data set.
}

\usage{
afCEC (
points,
maxClusters,
initialLabels="k-means++",
cardMin=0.05,
costThreshold=-0.000001,
minIterations=0,
maxIterations=100,
numberOfStarts=1,
method="Hartigan",
values="quadratic",
interactive=FALSE
)
}

\arguments{

\item{points}{
A (#points x n) matrix of data containing n-dimensional points stored row-by-row.
}

\item{maxClusters}{
Indicates the maximal number of clusters, that afCEC can partition the data into.
}

\item{initialLabels}{
Initial labelling determining the data membership to the particular clusters. There are 3 options allowed:
\itemize{
\item{"random" - causes the labelling to be randomly generated.}
\item{"k-means++" - causes the labelling to be generated using the k-means++ heuristics.}
\item{A (#points x numberOfStarts) matrix of type integer containing values in the range: 0...maxClusters - 1. The value x in the row i and column j indicates that in the j-th start, the i-th point will be initially assigned to the (x + 1)-th cluster.}
}
The defalt value is "k-means++".
}

\item{cardMin}{
Value so that if the number of points in particular cluster relatively to the size of the whole data drops bellow, then that cluster gets removed and it's data is assigned to the other clusters. The default value is 0.05.
}

\item{costTreshold}{
Negative value so that if the difference between the calculated cost functions in the subsequent iterations is greater than it then the afCEC terminates returning the solution from the terminated iteration as the final solution (in the present start), provided that at least minIterations already passed. The default value is -0.000001;
}

\item{minIterations}{
Value indicating the minimal number of iterations needed before afCEC can terminate, provided that no error occurred. The default value is 0.
}

\item{maxIterations}{
Value indicating the maximal number of iterations in one start, that afCEC cannot exceed. The default value is 100.
}

\item{numberOfStarts}{
Value indicating the number of runs of the algorithm to be performed. The best solution is chosen out of the solutions given in the particular steps, so that it has the lowest value of the cost function among them. The default value is 1.
}

\item{method}{
Heuristics used to perform the clustering. There are 2 options allowed:
\itemize{
\item{"Lloyd" - indicates that the Lloyd heuristics will be used to perform the clustering.}
\item{"Hartigan" - indicates that the Hartigan heuristics will be used to perform the clustering.}
}
The default value is "Hartigan".
}

\item{values}{
Definition of the active function used to perform the clustering. There are 3 options allowed:
\itemize{
\item{"quadratic" - indicates that the function of the form:
\cr
\cr
f(x_1,...,x_(n-1))=x_1^2+...+x_(n-1)^2+x_1+...+x_(n-1)+1,
\cr
\cr
where n is the dimensionality of the data, will be used to perform the clustering.
}
\item{"String containing the user-defined formula. See more in \code{User-defined formulas}.}
\item{A ((m*n) x #points) matrix containing the values of the particular components of the active function placed according to the following layout (row by row):
\cr
\cr
First m rows:
\cr
\cr
[f_1(x_2_1,...,x_n_1), f_1(x_2_2,...,x_n_2), ..., f_1(x_2_#points,...,x_n_#points)]\cr
[f_2(x_2_1,...,x_n_1), f_2(x_2_2,...,x_n_2), ..., f_2(x_2_#points,...,x_n_#points)]\cr
...\cr
[f_m(x_2_1,...,x_n_1), f_m(x_2_2,...,x_n_2), ..., f_m(x_2_#points,...,x_n_#points)]\cr
\cr
Second m rows:
\cr
\cr
[f_1(x_1_1,x_3_1,...,x_n_1), f_1(x_1_2,x_3_2,...,x_n_2), ...]\cr
[f_2(x_1_1,x_3_1,...,x_n_1), f_2(x_1_2,x_3_2,...,x_n_2), ...]\cr
...\cr
[f_m(x_1_1,x_3_1,...,x_n_1), f_m(x_1_2,x_3_2,...,x_n_2), ...]\cr
\cr
Last m rows:
\cr
\cr
[f_1(x_1_1,...,x_(n-1)_1), f_1(x_1_2,...,x_(n-1)_2), ...]\cr
[f_2(x_1_1,...,x_(n-1)_1), f_2(x_1_2,...,x_(n-1)_2), ...]\cr
...\cr
[f_m(x_1_1,...,x_(n-1)_1), f_m(x_1_2,...,x_(n-1)_2), ...],\cr
\cr
where: n - dimensionality of the data, x_i_j - j-th coordinate of the i-th point of data, m - number of components of the active function.
}
}
The default value is "quadratic".
\cr
\cr
\strong{Remarks:}
\cr
\cr
In the case of the third way of defining the function, the returned object doesn't contain the information about the means coordinates corresponding to the active direction. Moreover, the plotting capabilities of the afCEC package are in that case severely impaired. See \code{plot}.
}

\item{interactive}{
Indicates if the algorithm runs in the "interactive" mode. The "interactive" mode allows to track the intermediate steps of the afCEC. Instead of one object of the afCEC class, the whole list of such objects is returned, where each item of the list corresponds to the intermediate step of the algorithm. See \code{value} section for more details. The default value is FALSE.
}

}
\value{
\itemize{
\item{Empty list - if clustering failed.}
\item{Object of class afCEC containing the parameters of the best fitted mixed generalised multivariate normal distribution - if argument interactive=FALSE and clustering succeded.}
\item{List of k objects of class afCEC containing the parameters of the best fitted mixed generalised multivariate normal distribution in the k intermediate steps of running of the algorithm - if argument interactive=TRUE and clustering succeeded.}
}
}

\seealso{
\code{\link[=UserDefinedFormulas]{User-defined formulas}}\cr
\code{\link{plot}}
}

\keyword{ afCEC }