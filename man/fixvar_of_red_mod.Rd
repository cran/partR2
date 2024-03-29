% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/r2_helpers.R
\name{fixvar_of_red_mod}
\alias{fixvar_of_red_mod}
\title{Calculate fixed effect variance from a reduced model}
\usage{
fixvar_of_red_mod(partvar, mod, dat, expct, overdisp_name, R2_type)
}
\arguments{
\item{partvar}{One or more fixed effect variables which are taken out
of the model.}

\item{mod}{merMod object.}

\item{dat}{Data.frame to fit the model}

\item{expct}{Expectation}

\item{overdisp_name}{Name of overdispersion term}

\item{R2_type}{"marginal" or "conditional"}
}
\value{
R2 of reduced model.
}
\description{
Calculate fixed effect variance from a reduced model
}
\keyword{internal}
