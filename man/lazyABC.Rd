\name{lazyABC}
\alias{lazyABC}
\title{Perform lazy ABC algorithm}
\usage{
lazyABC(yobs, n.its, eps, stopstep, alpha = NULL, parallel = TRUE,
  S0 = 1e+06, I0 = 1, R0 = 0, n.subsample = 100)
}
\arguments{
  \item{yobs}{Observed number recovered}

  \item{n.its}{Number of iterations to perform}

  \item{eps}{ABC threshold}

  \item{stopstep}{When to consider early stopping (set to
  Inf for ordinary ABC)}

  \item{alpha}{Continuation probability function: a
  function of one non-negative integer (number infectious)
  which returns a scalar in [0,1]. Leave NULL for ordinary
  ABC.}

  \item{parallel}{If TRUE iterations are performed in
  parallel}

  \item{S0}{Initial number susceptible}

  \item{I0}{Initial number infectious}

  \item{R0}{Initial number recovered}

  \item{n.subsample}{Size of subsample}
}
\value{
A list comprising: ABCsample - dataframe of R0 and weight;
time - sum of elapsed time in each core.
}
\description{
Perform rejection sampling lazy ABC on the SIR epidemic
example or, as a special case, perform ordinary ABC
}
\details{
Prior is gamma=1 and beta~Gamma(5,1) (equivalently
R0~Gamma(5,1))
}
