\name{SIRsample}
\alias{SIRsample}
\title{Sample from a SIR output}
\usage{
SIRsample(N, R, n)
}
\arguments{
  \item{N}{Total population size}

  \item{R}{Number recovered in population}

  \item{n}{Subsample size}
}
\value{
Number recovered in subsample
}
\description{
Number recovered in a subsample of SIR output
}
\details{
Takes a sample of n from a population of N of whom R are
recovered, and returns the number of the sample who are
found to be recovered.
}

