\name{SampleSphereUniform}
\alias{SampleSphereUniform}

\title{
SampleSphereUniform
}

\usage{
SampleSphereUniform(n, dim)
}

\description{
Generates \code{n} uniformly distributed points on a \code{dim}-dimensional unit sphere embeded in the (\code{dim}+1)-dimensional space. The chart below depicts the result of executing the code placed in the "examples" section:

\tabular{c}{
\if{html}{\figure{SampleSphereUniform.svg}{options: width=300pt}}
\if{latex}{\figure{SampleSphereUniform.pdf}{options: width=300pt}}
}
}

\value{
A (\code{n} x \code{dim}+1) matrix of (\code{dim}+1)-dimensional points stored row-by-row
}

\examples{
library(afCEC);
points <- SampleSphereUniform(1000, 2);
plot3d(points, xlab="X", ylab="Y", zlab="Z");
}