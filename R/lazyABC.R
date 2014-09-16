#' Perform lazy ABC algorithm
#'
#' Perform rejection sampling lazy ABC on the SIR epidemic example or, as a special case, perform ordinary ABC
#'
#' @param yobs Observed number recovered
#' @param n.its Number of iterations to perform
#' @param eps ABC threshold
#' @param stopstep When to consider early stopping (set to Inf for ordinary ABC)
#' @param alpha Continuation probability function: a function of one non-negative integer (number infectious) which returns a scalar in [0,1]. Leave NULL for ordinary ABC.
#' @param parallel If TRUE iterations are performed in parallel
#' @param S0 Initial number susceptible
#' @param I0 Initial number infectious
#' @param R0 Initial number recovered
#' @param n.subsample Size of subsample
#'
#' @details Prior is gamma=1 and beta~Gamma(5,1) (equivalently R0~Gamma(5,1))
#' 
#' @return A list comprising: ABCsample - dataframe of R0 and weight; time - sum of elapsed time in each core.
#'
#' @export
lazyABC <- function(yobs, n.its, eps, stopstep, alpha=NULL, parallel=TRUE,
                    S0=1E6, I0=1, R0=0, n.subsample=100) {
    doiteration <- function() { ##A single lazy ABC iteration
        t0 <- proc.time()[3]
        betastar <- rgamma(1, shape=3, scale=1)
        gammastar <- 1
        sim1 <- SIRsim(betastar, gammastar, S0, I0, R0, t0=0, step_end=stopstep, thinning=0)
        if (sim1$I == 0) {
            ##Simulation finished before condition to consider early stopping met
            early_stopping <- FALSE
            astar <- 1
            ystar <- SIRsample(N=S0+I0+R0, R=sim1$R, n=n.subsample)
        } else {
            ##Consider early stopping
            astar <- alpha(sim1$I)
            if (runif(1)<=astar) { ##Continuation
                early_stopping <- FALSE
                sim2 <- SIRsim(betastar, gammastar, S0=sim1$S, I0=sim1$I, R0=sim1$R, t0=0, step_end=Inf, thinning=0)
                ystar <- SIRsample(N=S0+I0+R0, R=sim2$R, n=n.subsample)
            } else { ##Early stopping
                early_stopping <- TRUE
            }
        }
        if (early_stopping) {
            t1 <- proc.time()[3]
            out <- c(R0=betastar, weight=0, time=t1-t0)
        } else {
            lstarABC <- abs(ystar-yobs) <= eps
            w <- lstarABC / astar ##Assume importance weight equals prior
            t1 <- proc.time()[3]
            out <- c(R0=betastar, weight=w, time=t1-t0)
        }
        return(out)
    }
    if (parallel) {
        samp <- mclapply(1:n.its, function(i){ doiteration() },
                         mc.preschedule=FALSE)
        samp <- do.call(rbind, samp)
    } else {
        samp <- replicate(n=n.its, doiteration())
        samp <- t(samp)
    }
    samp <- data.frame(samp)
    ttotal <- sum(samp$time)
    samp <- samp[samp$weight>0,-3]
    return(list(ABCsample=samp, time=ttotal))
}
