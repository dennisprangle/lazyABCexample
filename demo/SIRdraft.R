require(magrittr)
require(mgcv)
require(parallel)
##Simulate true data
R0.true <- 2
set.seed(1)
Rtrue <- SIRsim(R0.true,1,1E5-1E3,1E3,0,0)$R
Robs <- SIRsample(1E5, Rtrue, 100) ##84

##Ordinary ABC
myeps <- 10
res.ord <- lazyABC(Robs, 1E3, eps=myeps, stopstep=Inf,
                   S0=1E5-1E3, I0=1E3, R0=0)

##Ad-hoc lazy ABC
res.lazy1 <- lazyABC(Robs, 1E3, eps=myeps, stopstep=1000,
                     alpha=function(I1000){ if(I1000<1000) { 0.1 } else { 1 } },
                     S0=1E5-1E3, I0=1E3, R0=0)

##Sample training data
n.train <- 1E3
R0.train <- rgamma(n.train, shape=3, scale=1)
trainingsim <- function(i) {
    out <- SIRsim(R0.train[i],1,1E5-1E3,1E3,0,0,thinning=1000)
    out <- cbind(out, train.it=i)
    return(out)
}
temp <- mclapply(1:n.train, trainingsim, mc.preschedule=FALSE)
train.final <- sapply(1:n.train,
                      function(i) {
                          x <- temp[[i]]
                          myR <- tail(x, n=1)$R
                          SIRsample(1E5+1, myR, 100)
                      })
elapsed.final <- sapply(1:n.train,
                        function(i) {
                            x <- temp[[i]]
                            tail(x, n=1)$elapsed
                      })
train.all <- do.call(rbind, temp)

##Compare training data at a particular time step to observations
train.subset <- train.all[train.all$step==1E3,]
subset.I <- train.subset$I
subset.Rfinal <- train.final[train.subset$train.it]
subset.diff <- abs(subset.Rfinal - Robs)
plot(subset.I, subset.Rfinal)

##Fit a non-linear regression of observable data on decision statistic
resp <- cbind(success=subset.Rfinal, fail=100-subset.Rfinal) ##Correct form for glm-like function
icov <- subset.I ##Simple variable name for glm-like function
fit <- gam(resp ~ s(icov), family=binomial)
zz <- predict(fit, type="response")
points(icov, 100*zz, col="red")
points(subset.I, qbinom(p=0.025, size=100, prob=zz), col="green")
points(subset.I, qbinom(p=0.975, size=100, prob=zz), col="green")
##Estimate gamma (or "hit probability") given decision statistic
hitPr <- function(icov, Robs, eps) {
    p <- predict(fit, newdata=data.frame(icov=icov), type="response")
    sapply(p, function(q) dbinom(seq(Robs-eps,Robs+eps), size=100, prob=q) %>% sum)
}
xx <- seq(min(icov), max(icov), length.out=100)
plot(xx, hitPr(xx, Robs, myeps), type='l')
##Estimate gamma just for training data
temp <- predict(fit, type="response")
gamma.train <- sapply(temp, function(q) dbinom(seq(Robs-myeps,Robs+myeps), size=100, prob=q) %>% sum)

##Fit a non-linear regression of time remaining given decision statistic
respT <- elapsed.final - train.subset$elapsed
plot(subset.I, respT)
fit.tbar <- gam(respT ~ s(icov))
lines(xx, predict(fit.tbar, data.frame(icov=xx)), col="red")

##Function to estimate efficiency
eff.est <- make.effest(phi=icov, gamma=gamma.train, T2=predict(fit.tbar), T1bar=mean(train.subset$elapsed))


eff.est <- function(lambda) {
    alpha.train <- pmin(1, lambda*sqrt(gamma.train/T2bar)) #Continuation probability
    plot(icov, alpha.train)
    sumw2.base <- sum(gamma.train)
    sumw2 <- sum(gamma.train/alpha.train, na.rm=TRUE) ##Expected sum of weights^2.  na.rm takes care of cases where gamma.train=alpha.train=0.
    t.base <- Tstan*n.train
    t.lazy <- T1bar*n.train + mean(T2bar)*sum(alpha.train) ##E(1st stage time)*number of iterations + E(2nd stage time)*expected number of continuations
    ##cat(lambda, t.base, t.lazy, sumw, sumw2, sumw/sumw2, t.base*sumw / (t.lazy*sumw2), "\n")
    t.base*sumw2.base / (t.lazy*sumw2)
}

tomin <- function(lambda) -eff.est(lambda)
optimise(tomin, lower=0, upper=1E4)

##Do lazy ABC with optimal alpha
alpha.opt <- function(I1000) {
    T2bar <- predict(fit.tbar, data.frame(icov=I1000))
    gamma <- hitPr(I1000, Robs, 10)
    min(1, 1.539693 * sqrt(gamma/T2bar))
}

res.lazy2 <- lazyABC(Robs, 1E3, eps=10, stopstep=1000,
                     alpha=alpha.opt,
                     S0=1E5-1E3, I0=1E3, R0=0)
