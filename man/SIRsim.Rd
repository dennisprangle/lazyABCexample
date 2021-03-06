\name{SIRsim}
\alias{SIRsim}
\title{Markovian SIR epidemic model}
\usage{
SIRsim(beta, gamma, S0, I0, R0, t0, norm_beta = TRUE, t_end = Inf,
  step_end = Inf, thinning = 0)
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

  \item{thinning}{An integer determining what type of
  output to return. Zero gives the final state. A positive
  value thins the output by this factor i.e.
  \code{output=1} gives all output, \code{output=2} gives
  results for step 0,2,4,... The final state is always
  returned.}
}
\value{
A dataframe with columns step (number of steps performed),
S, I, R (state at this step) and elapsed (computer time
taken). The steps reported depend on the thinning variable.
}
\description{
Simulate from a Markovian SIR epidemic model.
}
\details{
This function performs a simulation.  This terminates when
(a) the specified time is reached or (b) the specified
number of steps have been performed or (c) zero infectious
remain.  A simulation terminated for (a) or (b) can be
continued by starting another simulation from the final
state.
}

