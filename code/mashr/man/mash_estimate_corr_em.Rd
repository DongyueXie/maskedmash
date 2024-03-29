% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/est_cor.R
\name{mash_estimate_corr_em}
\alias{mash_estimate_corr_em}
\title{Fit mash model and estimate residual correlations using EM algorithm}
\usage{
mash_estimate_corr_em(
  data,
  Ulist,
  init,
  max_iter = 30,
  tol = 1,
  est_cor = TRUE,
  track_fit = FALSE,
  prior = c("nullbiased", "uniform"),
  details = TRUE,
  ...
)
}
\arguments{
\item{data}{a mash data object, eg as created by \code{mash_set_data}}

\item{Ulist}{a list of covariance matrices to use}

\item{init}{the initial value for the residual correlation. If it is
not given, we use result from
\code{estimate_null_correlation_simple}}

\item{max_iter}{maximum number of iterations to perform}

\item{tol}{convergence tolerance}

\item{est_cor}{whether to estimate correlation matrix (TRUE) or the
covariance matrix (FALSE)}

\item{track_fit}{add an attribute \code{trace} to output that saves
current values of all iterations}

\item{prior}{indicates what penalty to use on the likelihood, if any}

\item{details}{whether to return details of the model, if it is
TRUE, the mash model, the number of iterations and the value of
objective functions will be returned}

\item{...}{other parameters pass to \code{mash}}
}
\value{
the estimated correlation matrix and the
fitted mash model \cr

\item{V}{estimated residual correlation matrix}

\item{mash.model}{fitted mash model}
}
\description{
Estimates a residual correlation matrix from data using an ad hoc EM
algorithm.
}
\details{
Returns the estimated residual correlation matrix among conditions.
We estimate the residual correlation matrix using an ad hoc em algorithm.
The update in the ad hoc M step is not guaranteed to increase the likelihood,
therefore, the EM algorithm is stopped before the likelihood drops.
The residual correlation matrix V is estimated using the posterior
second moment of the noise.

Warning: This method could take some time.  The
\code{\link{estimate_null_correlation_simple}} gives a quick
approximation for the null correlation matrix.
}
\examples{
simdata = simple_sims(100,5,1)
m.1by1 = mash_1by1(mash_set_data(simdata$Bhat,simdata$Shat))
strong.subset = get_significant_results(m.1by1,0.05)
random.subset = sample(1:nrow(simdata$Bhat),20)
data.strong = mash_set_data(simdata$Bhat[strong.subset,], simdata$Shat[strong.subset,])
data.tmp = mash_set_data(simdata$Bhat[random.subset,], simdata$Shat[random.subset,])
U_pca = cov_pca(data.strong, 3)
U_ed = cov_ed(data.strong, U_pca)
Vhat = mash_estimate_corr_em(data.tmp, U_ed)
}
