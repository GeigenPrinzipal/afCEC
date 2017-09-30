\name{SampleBallUniform}
\alias{SampleBallUniform}

\title{
SampleBallUniform
}

\usage{
SampleBallUniform(n, dim)
}

\description{
Generates \code{n} uniformly distributed points inside a \code{dim}-dimensional unit ball. The chart below depicts the result of executing the code placed in the "examples" section:

\tabular{c}{
\if{html}{\figure{SampleBallUniform.svg}{options: width=300pt}}
\if{latex}{\figure{SampleBallUniform.pdf}{options: width=300pt}}
}
}

\value{
A (\code{n} x \code{dim}) matrix of \code{dim}-dimensional points stored row-by-row
}

\examples{
library(afCEC);
points <- SampleBallUniform(1000, 3);
plot3d(points, xlab="X", ylab="Y", zlab="Z");
}