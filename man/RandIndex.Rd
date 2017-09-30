\name{RandIndex}
\alias{RandIndex}

\title{
RandIndex
}

\usage{
RandIndex(labels1, labels2)
}

\description{
Computes the Rand similarity index between two clusterings determined by the indices of the clusters the particular data points belong to, given by the parameters: \code{labels1} and \code{labels2} respectively. Both parameters are the equal-length vectors of the type integer with values greater than or equal to 0. If one of the above-mentioned conditions isn't met, the exception is thrown and the function returns with errors.
}

\value{
A value in the range [0,1]
}

\examples{
library(afCEC);
points1 <- SampleBallUniform(1000,2);
points2 <- (SampleBallUniform(2000,2)*2) + matrix(rep(c(2,0),2000),2000,2,byrow=TRUE);
points <- rbind(points1,points2);
labels <- c(rep(0,1000),rep(1,2000));
result <- afCEC(points,2);
ri <- RandIndex(labels,result$labels);
print(
    paste(
         "The Rand index between the model clustering and the clustering obtained by the ",
        "afCEC function is: ",
        ri,
        ".",
        sep=""
    )
);
}