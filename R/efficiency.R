#' Create efficiency estimator function
#'
#' Outputs a function of lambda (scalar) which estimates relative efficiency
#'
#' @param phi Decision statistics in training set (assumed scalar)
#' @param gamma Estimated gamma values (probability of ABC acceptance) corresponding to phi
#' @param T2 Estimated T2 values (time to complete simulation) corresponding to phi
#' @param T1bar Mean T1 value (time for initial simulation stage)
#'
#' @return
#' A function \code{effest(lambda)} which estimates relative efficiency of lazy ABC compared to standard ABC.
#' @export
make.effest <- function(phi, gamma, T2, T1bar) {
    n.train <- length(phi)
    effest <- function(lambda, plot=FALSE) {
        alpha <- pmin(1, lambda*sqrt(gamma/T2)) #Continuation probabilities
        if (plot) plot(phi, alpha)
        sumw2.base <- mean(gamma) ##Mean of squared weights under standard ABC
        sumw2 <- mean(gamma/alpha, na.rm=TRUE) ##Expected mean of weights^2.  na.rm takes care of cases where gamma=alpha=0.
        T2bar <- mean(T2)
        t.base <- (T1bar+T2bar)*n.train ##Expected time under standard ABC
        t.lazy <- T1bar*n.train + T2bar*sum(alpha) ##E(1st stage time)*number of iterations + E(2nd stage time)*expected number of continuations
        t.base*sumw2.base / (t.lazy*sumw2) ##Relative efficiency
    }
    return(effest)
}
