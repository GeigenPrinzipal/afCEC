\name{SampleMeshUniform}
\alias{SampleMeshUniform}

\title{
SampleMeshUniform
}

\usage{
SampleMeshUniform(path, n)
}

\description{
Generates \code{n} uniformly distributed points on the surface of the mesh of the 3D object in the *.obj format with the path given by the parameter \code{path}. The chart below depicts the result of executing the code placed in the "examples" section:

\tabular{c}{
\if{html}{\figure{SampleMeshUniform.svg}{options: width=300pt}}
\if{latex}{\figure{SampleMeshUniform.pdf}{options: width=300pt}}
}
}

\value{
A (\code{n} x 3) matrix of three-dimensional points stored row-by-row
}

\examples{
library(afCEC);
points <- SampleMeshUniform("neptune.obj", 10000);
plot3d(points, aspect=FALSE, axes=FALSE, xlab="", ylab="", zlab="");
axes3d(box=TRUE, xlen=0, ylen=0, zlen=0);
}