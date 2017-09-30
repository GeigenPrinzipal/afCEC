\name{SampleTorusUniform}
\alias{SampleTorusUniform}

\title{
SampleTorusUniform
}

\usage{
SampleTorusUniform(R, r, n)
}

\description{
Generates \code{n} approximately uniformly distributed points on the surface of the three-dimensional torus obtained by revolution of the circle with the radius \code{r} around the Z axis. The parameter \code{R} corresponds to the distance from the torus center to the center of any circular cross section. The chart below depicts the result of executing the code placed in the "examples" section:

\tabular{c}{
\if{html}{\figure{SampleTorusUniform.svg}{options: width=300pt}}
\if{latex}{\figure{SampleTorusUniform.pdf}{options: width=300pt}}
}
}

\section{Remarks}{
The function uses the inverse cumulative distribution function method to obtain the approximately uniformly distributed samples. As one of the marginal cdfs corresponding to one of the parameters of the torus parametric equation doesn't have the inverse function that can be expressed in terms of the elementary functions, the routine \code{SampleTorusUniform} uses its reasonable approximation. Therefore, the resulting sampling is exact only for approximately close values of \code{R} and \code{r}.
}

\value{
A (\code{n} x 3) matrix of three-dimensional points stored row-by-row
}

\examples{
library(afCEC);
points <- SampleTorusUniform(0.75, 0.25, 1000);
plot3d(points, aspect=FALSE, xlab="X", ylab="Y", zlab="Z");
}