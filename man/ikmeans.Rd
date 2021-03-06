% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ikmeans.R
\name{ikmeans}
\alias{ikmeans}
\title{Interval kmeans clustering.}
\usage{
ikmeans(x, centers, nstart = 10, distance = "euclid", trace = FALSE,
  iter.max = 20)
}
\arguments{
\item{x}{An 3D interval array.}

\item{centers}{A number or interval, number of cluster for clustering or pre init centers.}

\item{nstart}{A number, number of execution to find the best result.}

\item{distance}{A string ('euclid': Euclidian distance, 'hausdorff': Hausdorff distance).}

\item{trace}{A boolean, tracing information on the progress of the algorithm is produced.}

\item{iter.max}{the maximum number of iterations allowed.}
}
\description{
Culster interval data with kmeans algorithm.
}
\examples{
ikmeans(iaggregate(iris, col=5), 2)
ikmeans(iaggregate(iris, col=5), iaggregate(iris, col=5))
}
