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
#' @return A matrix with S, I, R, t (these 4 are the model output) and elapsed (computer time taken) columns is returned.  This includes the initial state.

#' @export
SIRsim <- function(beta, gamma, S0, I0, R0, t0, norm_beta=TRUE, t_end=Inf, step_end=Inf) {
    start_time <- proc.time()[3]
    if (norm_beta) beta <- beta / (S0+I0+R0)
    ##Initialise output
    output <- matrix(nrow=2*S0+I0+1, ncol=5) ##Memory allocated for max number of events plus first state
    colnames(output) <- c("S","I","R","t","elapsed")
    output[1,] <- c(S0,I0,R0,t0,0)
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
        output[step+1,] <- c(S,I,R,t,proc.time()[3]-start_time)
        if (step==step_end || I==0) break;
    }
    return(output[seq(1,step+1),])
}

#' @export
SIRsample <- function(N, R, n) {
    rhyper(1, R, N-R, n)
}
