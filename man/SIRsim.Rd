\name{SIRsim}
\alias{SIRsim}
\title{Markovian SIR epidemic model}
\usage{
SIRsim(beta, gamma, S0, I0, R0, t0, norm_beta = TRUE, t_end = Inf,
  step_end = Inf)
}
\arguments{
  \item{beta}{The infection parameter}

  \item{gamma}{Recovery rate}

  \item{S0}{Initial number susceptible}

  \item{I0}{Initial number infectious}

  \item{R0}{Initial number recovered}

  \item{t0}{Initial time}

  \item{norm_beta}{Whether beta should be normalised (i.e.
  divided) by population size}

  \item{t_end}{Time at which the simulation will terminate}

  \item{step_end}{Number of steps after which the
  simulation will terminate}
}
\value{
A matrix with S, I, R, t (these 4 are the model output) and
elapsed (computer time taken) columns is returned.  This
includes the initial state.
}
\description{
Functions to simulate from a Markovian SIR epidemic model.
}
\details{
\code{SIRsim} performs a simulation.  This terminates when
(a) the specified time is reached or (b) the specified
number of steps have been performed or (c) zero infectious
remain.  A simulation terminated for (a) or (b) can be
continued by starting another simulation from the final
state. \code{SIRsample} supposes a sample of n are taken
from a population of N of whom R are recovered, and returns
the number of the sample who are found to be recovered.
}

