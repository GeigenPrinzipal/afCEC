\name{afCEC-package}
\alias{afCEC-package}
\docType{package}
\title{
Active Function Cross-Entropy Clustering
}
\description{
Active Function Cross-Entropy Clustering (afCEC) partitions the n-dimensional data into the clusters by finding the parameters of the mixed generalised multivariate normal distribution, that optimally approximates the scattering of the data in the n-dimensional space, whose density function is of the form:

p_1*N(mi_1,^sigma_1,sigma_1,f_1)+...+p_k*N(mi_k,^sigma_k,sigma_k,f_k)

Additionally, the Active Function Cross-Entropy Clustering performs the automatic reduction of the unnecessary clusters. It's implemented in the form of a customizable function \code{afCEC}.
}
\details{
\tabular{ll}{
Package: \tab afCEC\cr
Type: \tab Package\cr
Version: \tab 0.9.1\cr
Date: \tab 2014-10-06\cr
License: \tab GPL-3 \cr
}
}
\author{
Krzysztof Byrski
}

\seealso{
\code{\link{afCEC}.}
}
\keyword{ clustering, entropy, gaussian, kmeans }
