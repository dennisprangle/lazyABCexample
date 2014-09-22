require(magrittr)
require(mgcv)
require(parallel)
require(ggplot2)
##Simulate true data
R0.true <- 2
set.seed(1)
Rtrue <- SIRsim(R0.true,1,1E5-1E3,1E3,0,0)$R
Robs <- SIRsample(1E5, Rtrue, 100) ##73

##ORDINARY ABC
myeps <- 1
res.ord <- lazyABC(Robs, 1E4, eps=myeps, stopstep=Inf,
                   S0=1E5-1E3, I0=1E3, R0=0)

##AD-HOC LAZY ABC
res.lazy1 <- lazyABC(Robs, 1E4, eps=myeps, stopstep=1000,
                     alpha=function(I1000){ if(I1000<=1000) { 0.1 } else { 1 } },
                     S0=1E5-1E3, I0=1E3, R0=0)

##Sample training data
n.train <- 1E3
R0.train <- rgamma(n.train, shape=3, scale=1)
trainingsim <- function(i) {
    set.seed(i+n.train) ##(a) Reproducable (b) Different seeds to lazyABC code
    out <- SIRsim(R0.train[i],1,1E5-1E3,1E3,0,0,thinning=1000)    
    out <- cbind(out, train.it=i)
    return(out)
}
train.samp <- mclapply(1:n.train, trainingsim, mc.preschedule=FALSE)

##STANDARD TUNING
##Process training data
train.final <- sapply(1:n.train,
                      function(i) {
                          x <- train.samp[[i]]
                          myR <- tail(x, n=1)$R
                          SIRsample(1E5+1, myR, 100)
                      }) ##Observation from each training run
elapsed.final <- sapply(1:n.train,
                        function(i) {
                            x <- train.samp[[i]]
                            tail(x, n=1)$elapsed
                      }) ##Total elapsed time of each training run
train.1000 <- lapply(1:n.train,
                      function(i) {
                          x <- train.samp[[i]]
                          x[2,] ##n.b. process runs at least 1000 steps so it's guaranteed that this row correspond to step 1000
                      }) ##State after 1000 steps from each training run
train.1000 <- do.call(rbind, train.1000)

##Compare training data at a particular time step to observations
I1000 <- train.1000$I ##I1000 from training data
train.obs <- train.final[train.1000$train.it] ##Observations corresponding to I1000
train.diff <- abs(train.obs - Robs)
plot(I1000, train.obs)

##Fit a non-linear regression of observable data on decision statistic
resp <- cbind(success=train.obs, fail=100-train.obs) ##Form needed for glm-like function
fit.resp <- gam(resp ~ s(I1000), family=binomial)

##Check non-linear regression by plotting predictions and CIs
zz <- predict(fit.resp, type="response")
points(I1000, 100*zz, col="red")
points(I1000, qbinom(p=0.025, size=100, prob=zz), col="green")
points(I1000, qbinom(p=0.975, size=100, prob=zz), col="green")

##Estimate gamma (aka "hit probability") given decision statistic
hitPr <- function(I1000, Robs, eps) {
    p <- predict(fit.resp, newdata=data.frame(I1000=I1000), type="response")
    sapply(p, function(q) dbinom(seq(Robs-eps,Robs+eps), size=100, prob=q) %>% sum)
}

##Plot hit probability
xx <- seq(min(I1000), max(I1000), length.out=100)
plot(xx, hitPr(xx, Robs, myeps), type='l')

##Estimate gamma on training data
temp <- predict(fit.resp, type="response")
gamma.train <- sapply(temp, function(q) dbinom(seq(Robs-myeps,Robs+myeps), size=100, prob=q) %>% sum)

##Fit a non-linear regression of time remaining given decision statistic
respT <- elapsed.final - train.1000$elapsed
fit.tbar <- gam(respT ~ s(I1000), gaussian(link = "log"))
T2bar.train <- predict(fit.tbar, type="response")

##Plot time regression against data as a check
plot(I1000, respT)
lines(xx, predict(fit.tbar, data.frame(I1000=xx), type="response"), col="red")

##Create function to estimate efficiency
eff.est <- make.effest(phi=I1000, gamma=gamma.train, T2=T2bar.train, T1bar=mean(train.1000$elapsed))

##Choose lambda to optimise efficiency
tomin <- function(lambda) -eff.est(lambda)
temp <- optimise(tomin, lower=0, upper=1E4)
temp
lambda.opt <- temp$minimum
##Plot optimal choice of alpha
eff.est(lambda.opt, plot=TRUE)

##Do lazy ABC with optimal alpha
alpha.opt <- function(I1000) {
    T2bar <- predict(fit.tbar, data.frame(I1000=I1000), type="response")
    gamma <- hitPr(I1000, Robs, 10)
    min(1, lambda.opt * sqrt(gamma/T2bar))
}
res.lazy2 <- lazyABC(Robs, 1E4, eps=myeps, stopstep=1000,
                     alpha=alpha.opt,
                     S0=1E5-1E3, I0=1E3, R0=0)

##CONSERVATIVE TUNING
##Estimate gamma
train.z <- (train.diff <= 3)
fit.cons <- gam(train.z ~ s(I1000), family=binomial)
gammacons.train <- predict(fit.cons, type="response")
##Plot gamma estimate
plot(I1000, predict(fit.cons, type="response"))
rug(I1000[train.z==0], side=1)
rug(I1000[train.z==1], side=3)

##Create function to estimate efficiency
eff.est.cons <- make.effest(phi=I1000, gamma=gammacons.train, T2=T2bar.train, T1bar=mean(train.1000$elapsed))

##Choose lambda to optimise efficiency
tomin <- function(lambda) -eff.est.cons(lambda)
temp <- optimise(tomin, lower=0, upper=1E4)
temp
lambda.opt.cons <- temp$minimum
##Plot optimal alpha
eff.est.cons(lambda.opt.cons, plot=TRUE)

##Do lazy ABC with optimal alpha
alpha.opt.cons <- function(I1000) {
    T2bar <- predict(fit.tbar, data.frame(I1000=I1000), type="response")
    gamma <- predict(fit.cons, data.frame(I1000=I1000), type="response")
    min(1, lambda.opt.cons * sqrt(gamma/T2bar))
}
res.lazy3 <- lazyABC(Robs, 1E4, eps=myeps, stopstep=1000,
                     alpha=alpha.opt.cons,
                     S0=1E5-1E3, I0=1E3, R0=0)

##EXAMINE OUTPUT
##Compare number of acceptances
nrow(res.ord$ABCsample)
nrow(res.lazy1$ABCsample)
nrow(res.lazy2$ABCsample)
nrow(res.lazy3$ABCsample)

##Doublecheck conservative tuning acceptances are a subset of ordinary ABC ones
table(rownames(res.lazy3$ABCsample) %in% rownames(res.ord$ABCsample))

##Compare ESS values
ess(res.ord$ABCsample$weight)
ess(res.lazy1$ABCsample$weight)
ess(res.lazy2$ABCsample$weight)
ess(res.lazy3$ABCsample$weight)

##Compare times
res.ord$time
res.lazy1$time
res.lazy2$time
res.lazy3$time

##Relative efficiencies
eff.ord <- ess(res.ord$ABCsample$weight) / res.ord$time
eff.lazy1 <- ess(res.lazy1$ABCsample$weight) / res.lazy1$time
eff.lazy2 <- ess(res.lazy2$ABCsample$weight) / res.lazy2$time
eff.lazy3 <- ess(res.lazy3$ABCsample$weight) / res.lazy3$time
eff.lazy1 / eff.ord
eff.lazy2 / eff.ord
eff.lazy3 / eff.ord

##What weights did lazy ABC have?
table(res.lazy1$ABCsample$weight)
table(res.lazy2$ABCsample$weight)
table(res.lazy3$ABCsample$weight)

##Estimates of posterior mean and variance
post <- function(res) {
    x <- res$ABCsample$R0
    w <- res$ABCsample$weight
    Ex <- sum(x*w) / sum(w)
    Ex2 <- sum(x^2*w) / sum(w)
    c(mean=Ex, sd=sqrt(Ex2-Ex^2))
}
post(res.ord)
post(res.lazy1)
post(res.lazy2)
post(res.lazy3)

##Estimated posterior histograms
qplot(x=R0, weight=weight, data=res.ord$ABCsample, binwidth=0.1)
qplot(x=R0, weight=weight, data=res.lazy1$ABCsample, binwidth=0.1)
qplot(x=R0, weight=weight, data=res.lazy2$ABCsample, binwidth=0.1)
qplot(x=R0, weight=weight, data=res.lazy3$ABCsample, binwidth=0.1)
