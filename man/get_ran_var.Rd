% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/r2_helpers.R
\name{get_ran_var}
\alias{get_ran_var}
\title{Extracts random effect variances}
\usage{
get_ran_var(mod, overdisp_name = NULL)
}
\arguments{
\item{mod}{An lme4 model object.}

\item{overdisp_name}{name of overdispersion term}
}
\description{
This function computes the sum of random effect variances where one
or more of the random effects are random slopes. It uses Paul Johnson' method
to compute the average group variance across the levels of a covariate.
This function extracts only grouping factors, no residual or overdispersion.
}
\keyword{internal}
