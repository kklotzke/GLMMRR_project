% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RRglmer.R
\name{RRglmer}
\alias{RRglmer}
\title{Fitting Generalized Linear Mixed-Effects Models with binary Randomized Response data}
\usage{
RRglmer(
  formula,
  item,
  link,
  RRmodel,
  p1,
  p2,
  data,
  control = glmerControl(),
  na.action = "na.omit",
  ...
)
}
\arguments{
\item{formula}{a two-sided linear formula object describing both the fixed-effects and fixed-effects part of the model,
with the response on the left of a ~ operator and the terms, separated by + operators, on the right.
Random-effects terms are distinguished by vertical bars ("|") separating expressions for design matrices from grouping factors.}

\item{item}{optional item identifier for long-format data.}

\item{link}{a glm link function for binary outcomes. Must be a function name.
Available options: "RRlink.logit", "RRlink.probit", "RRlink.cloglog" and "RRlink.cauchit"}

\item{RRmodel}{the Randomized Response model, defined per case.
Available options: "DQ", "Warner", "Forced", "UQM", "Crosswise", "Triangular" and "Kuk"}

\item{p1}{the Randomized Response parameter p1, defined per case. Must be 0 <= p1 <= 1.}

\item{p2}{the Randomized Response parameter p2, defined per case. Must be 0 <= p2 <= 1.}

\item{data}{a data frame containing the variables named in \code{\link{formula}} as well as the Randomized Response model and parameters.
If the required information cannot be found in the data frame, or if no data frame is given, then the variables are taken
from the environment from which RRglmer is called.}

\item{control}{a list (of correct class, resulting from \code{\link{lmerControl}()} or \code{\link{glmerControl}()} respectively) containing control parameters,
including the nonlinear optimizer to be used and parameters to be passed through to the nonlinear optimizer,
see the \code{*lmerControl} documentation for details.}

\item{na.action}{a function that indicates what should happen when the data contain NAs.
The default action (\code{\link{na.omit}}, as given by \code{getOption("na.action"))})
strips any observations with any missing values in any variables.}

\item{...}{other potential arguments to be passed to \code{\link{glmer}}.}
}
\value{
An object of class RRglmerMod. Extends the class \code{glmerMod} with Randomize Response data,
for which many methods are available (e.g. \code{methods(class="glmerMod")}).
}
\description{
Fit a generalized linear mixed-effects model (GLMM) with binary Randomized Response data.
Both fixed effects and random effects are specified via the model formula.
Randomize response parameters can be entered either as single values or as vectors.
Implemented as a wrapper for \code{\link{glmer}}. Reference: Fox, J-P, Veen, D. and Klotzke, K. (2018).
Generalized Linear Mixed Models for Randomized Responses. \emph{Methodology.} https://doi.org/10.1027/1614-2241/a000153
}
\examples{
# Fit the model with fixed effects for gender, RR and pp
# and a random effect for age using the logit link function.
# The Randomized Response parameters p1, p2 and model
# are specified for each observation in the dataset.
out <- RRglmer(response ~ Gender + RR + pp + (1|age), link="RRlink.logit", RRmodel=RRmodel,
         p1=RRp1, p2=RRp2, data=Plagiarism, na.action = "na.omit",
         etastart = rep(0.01, nrow(Plagiarism)),
         control = glmerControl(optimizer = "Nelder_Mead", tolPwrss = 1e-03), nAGQ = 1)
summary(out)
}
\seealso{
\code{\link{lme4}}
}
