#' Markovian SIR epidemic model
#'
#' Functions to simulate from a Markovian SIR epidemic model.
#'
#' \code{SIRsim} performs a simulation.  This terminates when (a) the specified time is reached or (b) the specified number of steps have been performed or (c) zero infectious remain.  A simulation terminated for (a) or (b) can be continued by starting another simulation from the final state.
#' \code{SIRsample} supposes a sample of n are taken from a population of N of whom R are recovered, and returns the number of the sample who are found to be recovered.
#'
#' @param beta The infection parameter
#' @param gamma Recovery rate
#' @param S0 Initial number susceptible
#' @param I0 Initial number infectious
#' @param R0 Initial number recovered
#' @param t0 Initial time
#' @param norm_beta Whether beta should be normalised (i.e. divided) by population size
#' @param t_end Time at which the simulation will terminate
#' @param step_end Number of steps after which the simulation will terminate
#' @param thinning An integer determining what type of output to return. Zero gives the final state. A positive value thins the output by this factor i.e. \code{output=1} gives all output, \code{output=2} gives results for step 0,2,4,... The final state is always returned.
#' @return If \code{thinning>0} a matrix is returned with columns step (number of steps performed), S, I, R (state at this step) and elapsed (computer time taken). The steps reported depend on the thinning variable. If \code{thinning==0} then only the final step is returned, as a vector.

#' @export
SIRsim <- function(beta, gamma, S0, I0, R0, t0, norm_beta=TRUE, t_end=Inf, step_end=Inf, thinning=0) {
    start_time <- proc.time()[3]
    if (norm_beta) beta <- beta / (S0+I0+R0)
    ##Initialise output
    if (thinning>0) {
        output <- matrix(nrow=ceiling((2*S0+I0+1)/thinning), ncol=6) ##n.b. 2*S0+I0+1 is maximum number of events, and 1 is added for initial state
        colnames(output) <- c("step","S","I","R","t","elapsed")
        output[1,] <- c(0,S0,I0,R0,t0,0)
        thincounter <- 0 ##How many steps since output last recorded
        recordedsofar <- 1 ##How many outputs recorded so far
    }
    ##Initialise current state
    S <- S0; I <- I0; R <- R0; t <- t0; step <- 0
    ##Main loop
    while(TRUE) {
        SI_rate <- S*I*beta
        IR_rate <- I*gamma
        total_rate <- SI_rate + IR_rate
        t <- t+rexp(1, total_rate)
        if (t > t_end) break;
        if (runif(1)*total_rate < SI_rate) {
          ##SI event
            S <- S-1; I <- I+1
        } else {
          ##IR event
            I <- I-1; R <- R+1
        }       
        step <- step+1
        if (thinning > 0) {
            thincounter <- thincounter+1
            if (thincounter == thinning) {
                output[recordedsofar+1,] <- c(step,S,I,R,t,proc.time()[3]-start_time)
                recordedsofar <- recordedsofar+1
                thincounter <- 0
            }
        }
        if (step==step_end || I==0) break;
    }
    if (thinning==0) {
        output <- c(step,S,I,R,t,proc.time()[3]-start_time)
        names(output) <- c("step","S","I","R","t","elapsed")
    } else {
        if (thincounter>1) {
            output[recordedsofar+1,] <- c(step,S,I,R,t,proc.time()[3]-start_time)
            recordedsofar <- recordedsofar+1
        }        
        output <- output[seq(1,recordedsofar),,drop=FALSE]
    }
    return(output)
}

#' @export
SIRsample <- function(N, R, n) {
    rhyper(1, R, N-R, n)
}
