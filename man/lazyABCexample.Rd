\docType{package}
\name{lazyABCexample}
\alias{lazyABCexample}
\alias{lazyABCexample-package}
\title{Example lazy ABC implementation}
\description{
Implements lazy ABC on a SIR epidemic application.
}
\details{
This R package implements an epidemics example from "Lazy
ABC" Dennis Prangle (2014), available at
http://arxiv.org/abs/1405.7867 **nb this example appears in
v2 only (still in preparation)**

The package contains commands to simulate from the model of
interest (SIRsim) and to perform lazy ABC, or standard ABC
as a special case (lazyABC). Two scripts are also supplied
as demos. Both perform the analysis in the paper.
\code{SIRproduction} is the exact code used in the paper.
\code{SIRexample} is a more user-friendly version (less
code, more explanatory plots). The scripts can be run via
\code{demo(SIRexample)} or \code{demo(SIRproduction)}. They
use parallel processing via \code{mclapply}, so the user
may wish to run \code{options(mc.cores=X)} for some
appropriate value of X before execution.
}
\author{
Dennis Prangle \email{d.b.prangle@reading.ac.uk}
}

