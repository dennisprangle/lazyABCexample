lazyABCexample
==============

This R package implements an epidemics example from

"Lazy ABC" Dennis Prangle (2014)

Available at http://arxiv.org/abs/1405.7867

**nb this example appears from v2 of this paper**


The package contains commands to simulate from the model of interest (SIRsim) and to perform lazy ABC, or standard ABC as a special case (lazyABC).
Two scripts are also supplied as demos.
Both perform the analysis in the paper.
"SIRproduction" is the exact code used in the paper.
"SIRexample" is a more user-friendly version (less code, more explanatory plots).
The scripts can be run by "demo(SIRexample)" or "demo(SIRproduction)".
They use parallel processing via "mclapply", so the user may wish to run "options(mc.cores=X)" for some appropriate value of X before execution.

