require(magrittr)
require(mgcv)
require(parallel)
##Simulate true data
R0.true <- 2
set.seed(1)
Rtrue <- SIRsim(R0.true,1,1E5-1E3,1E3,0,0)$R
Robs <- SIRsample(1E5, Rtrue, 100) ##73

##Ordinary ABC
myeps <- 1
res.ord <- lazyABC(Robs, 1E4, eps=myeps, stopstep=Inf,
                   S0=1E5-1E3, I0=1E3, R0=0)

##Ad-hoc lazy ABC
res.lazy1 <- lazyABC(Robs, 1E4, eps=myeps, stopstep=1000,
                     alpha=function(I1000){ if(I1000<1000) { 0.1 } else { 1 } },
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
temp <- mclapply(1:n.train, trainingsim, mc.preschedule=FALSE)

##Process training data
train.final <- sapply(1:n.train,
                      function(i) {
                          x <- temp[[i]]
                          myR <- tail(x, n=1)$R
                          SIRsample(1E5+1, myR, 100)
                      }) ##Observation from each training run
elapsed.final <- sapply(1:n.train,
                        function(i) {
                            x <- temp[[i]]
                            tail(x, n=1)$elapsed
                      }) ##Total elapsed time of each training run
train.1000 <- lapply(1:n.train,
                      function(i) {
                          x <- temp[[i]]
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
fit.tbar <- gam(respT ~ s(I1000))

##Plot time regression against data as a check
plot(I1000, respT)
lines(xx, predict(fit.tbar, data.frame(I1000=xx)), col="red")

##Function to estimate efficiency
eff.est <- make.effest(phi=I1000, gamma=gamma.train, T2=predict(fit.tbar), T1bar=mean(train.1000$elapsed))

##Choose lambda to optimise efficiency
tomin <- function(lambda) -eff.est(lambda)
temp <- optimise(tomin, lower=0, upper=1E4)
temp
lambda.opt <- temp$minimum
##Plot optimal choice of alpha
eff.est(lambda.opt, plot=TRUE)

##Do lazy ABC with optimal alpha
alpha.opt <- function(I1000) {
    T2bar <- predict(fit.tbar, data.frame(I1000=I1000))
    gamma <- hitPr(I1000, Robs, 10)
    min(1, lambda.opt * sqrt(gamma/T2bar))
}
res.lazy2 <- lazyABC(Robs, 1E4, eps=myeps, stopstep=1000,
                     alpha=alpha.opt,
                     S0=1E5-1E3, I0=1E3, R0=0)
