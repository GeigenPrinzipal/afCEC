\name{plot}
\alias{plot}

\title{
plot
}

\description{
Plots the chart showing the clustering of the data. Depending on the dimensionality of the data passed to the afCEC function, it draws either 2D or 3D chart. The function doesn't work with the other dimensions.
}

\usage{
plot (
x,
draw_points=TRUE,
draw_means=TRUE,
draw_ellipsoids=TRUE,
draw_surfaces=TRUE,
confidence=0.95,
grid_resolution=32
)
}

\arguments{

\item{x}{
Object of class afCEC returned from the afCEC function.
}

\item{draw_points}{
If the value is TRUE, function draws the points belonging to the particular clusters, so that the points belonging to the different clusters have the different colors. The default value is TRUE.
}

\item{draw_means}{
If the value is TRUE, function draws the means of the clusters as the big black dots / spheres. The default value is TRUE.\cr
\cr
\strong{Remarks:}\cr
\cr
This functionality doesn't work if object x was given after clustering the data with the active function defined by the matrix of it's components values.
}

\item{draw_ellipsoids}{
If the value is TRUE, function draws the curved ellipses / ellipsoids of confidence of the particular clusters. The default value is TRUE.\cr
\cr
\strong{Remarks:}\cr
\cr
This functionality doesn't work if object x was given after clustering the data with the active function defined by the matrix of it's components values.
}

\item{draw_surfaces}{
If the value is TRUE, function draws the line / curved surface spanned on the inactive principal axes of the curved ellipses / ellipsoids of confidence of the particular clusters. The default value is TRUE.\cr
\cr
\strong{Remarks:}\cr
\cr
This functionality doesn't work if object x was given after clustering the data with the active function defined by the matrix of it's components values.
}

\item{confidence}{
Determines the percentile of data belonging to the particular clusters, the corresponding ellipses / ellipsoids of confidence should contain. For example, if the value 0.5 is specified, then, for each cluster, it's ellipse / ellipsoid of confidence will contain 0.5 * 100 [\%] of assigned points. The default value is 0.95.
}

\item{grid_resolution}{
Determines the grid resolution using to approximate the curves / surfaces drawn in the plot function. The default value is 32.
}

}
\value{
None.
}

\seealso{
\code{\link{afCEC}}
}