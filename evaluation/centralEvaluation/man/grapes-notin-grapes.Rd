% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/centralEvaluation.r
\name{\%notin\%}
\alias{\%notin\%}
\title{Helper function}
\usage{
lhs \%notin\% rhs
}
\arguments{
\item{lhs}{The left hand side argument}

\item{rhs}{The right hand side argument}
}
\value{
Elements of the \code{lhs} argument that are not in the \code{rhs} argument
}
\description{
This little helper function is from
stackoverflow https://stackoverflow.com/questions/38351820/negation-of-in-in-r.
It was easier to define a not-in function than redo all those lists

Function definition is based on this stackoverflow answer:
https://stackoverflow.com/questions/5831794/opposite-of-in-exclude-rows-with-values-specified-in-a-vector
}
