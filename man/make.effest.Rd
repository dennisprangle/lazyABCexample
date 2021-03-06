\name{make.effest}
\alias{make.effest}
\title{Create efficiency estimator function}
\usage{
make.effest(phi, gamma, T2, T1bar)
}
\arguments{
  \item{phi}{Decision statistics in training set (assumed
  scalar)}

  \item{gamma}{Estimated gamma values (probability of ABC
  acceptance) corresponding to phi}

  \item{T2}{Estimated T2 values (time to complete
  simulation) corresponding to phi}

  \item{T1bar}{Mean T1 value (time for initial simulation
  stage)}
}
\value{
A function \code{effest(lambda)} which estimates relative
efficiency of lazy ABC compared to standard ABC.
}
\description{
Outputs a function of lambda (scalar) which estimates
relative efficiency
}

