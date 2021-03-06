% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interval.R
\name{plot.interval}
\alias{plot.interval}
\title{Interval plotting.}
\usage{
\method{plot}{interval}(x, ...)
}
\arguments{
\item{x}{An interval.}

\item{...}{:
- 'col' Colors of rectangles.
- 'add' Add to exsiting plot or not.
- Other options from plot.}
}
\description{
plot for interval.
}
\examples{
plot(iaggregate(iris, 5))
plot(iaggregate(iris, 5), col=4)
plot(iaggregate(iris, 5), add=TRUE)
}
