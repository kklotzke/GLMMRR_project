% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.R
\name{summary.RRglmerMod}
\alias{summary.RRglmerMod}
\title{Summarizing GLMMRR fits for fixed-effect models}
\usage{
\method{summary}{RRglmerMod}(object, p1p2.digits = 2, ...)
}
\arguments{
\item{object}{an object of class RRglm.}

\item{p1p2.digits}{number of digits for aggregating data based on the level of protection (default: 2).}

\item{...}{further arguments passed to or from other methods.}
}
\value{
An object of class summary.RRglmerMod. Extends the class \code{summary.glmerMod} with Randomize Response data.
}
\description{
Summarizing GLMMRR fits for fixed-effect models
}
