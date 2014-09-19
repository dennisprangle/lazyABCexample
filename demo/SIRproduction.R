require(magrittr)
require(mgcv)
require(parallel)
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
train.starttime <- proc.time()[3] ##nb simulation time will be added later
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

##Fit a non-linear regression of observable data on decision statistic
resp <- cbind(success=train.obs, fail=100-train.obs) ##Form needed for glm-like function
fit.resp <- gam(resp ~ s(I1000), family=binomial)

##Estimate gamma (aka "hit probability") given decision statistic
hitPr <- function(I1000, Robs, eps) {
    p <- predict(fit.resp, newdata=data.frame(I1000=I1000), type="response")
    sapply(p, function(q) dbinom(seq(Robs-eps,Robs+eps), size=100, prob=q) %>% sum)
}

##Estimate gamma on training data
temp <- predict(fit.resp, type="response")
gamma.train <- sapply(temp, function(q) dbinom(seq(Robs-myeps,Robs+myeps), size=100, prob=q) %>% sum)

##Fit a non-linear regression of time remaining given decision statistic
respT <- elapsed.final - train.1000$elapsed
fit.tbar <- gam(respT ~ s(I1000), gaussian(link = "log"))
T2bar.train <- predict(fit.tbar, type="response")

##Create function to estimate efficiency
eff.est <- make.effest(phi=I1000, gamma=gamma.train, T2=T2bar.train, T1bar=mean(train.1000$elapsed))

##Choose lambda to optimise efficiency
tomin <- function(lambda) -eff.est(lambda)
temp <- optimise(tomin, lower=0, upper=1E4)
temp
lambda.opt <- temp$minimum

##How long did tuning take?
proc.time()[3] - train.starttime + sum(elapsed.final)

##Do lazy ABC with optimal alpha
alpha.opt <- function(I1000) {
    T2bar <- predict(fit.tbar, data.frame(I1000=I1000), type="response")
    gamma <- hitPr(I1000, Robs, 10)
    min(1, lambda.opt * sqrt(gamma/T2bar))
}
res.lazy2 <- lazyABC(Robs, 1E4, eps=myeps, stopstep=1000,
                     alpha=alpha.opt,
                     S0=1E5-1E3, I0=1E3, R0=0)

##CONSERVATIVE TUNING (nb repeat everything necessary so timing correct)
train.starttime <- proc.time()[3] ##nb simulation time will be added later
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

##Estimate gamma
train.z <- (train.diff <= 2)
fit.cons <- gam(z ~ s(I1000), family=binomial)
gammacons.train <- predict(fit.cons, type="response")

##Fit a non-linear regression of time remaining given decision statistic
respT <- elapsed.final - train.1000$elapsed
fit.tbar <- gam(respT ~ s(I1000), gaussian(link = "log"))
T2bar.train <- predict(fit.tbar, type="response")

##Create function to estimate efficiency
eff.est.cons <- make.effest(phi=I1000, gamma=gammacons.train, T2=T2bar.train, T1bar=mean(train.1000$elapsed))

##Choose lambda to optimise efficiency
tomin <- function(lambda) -eff.est.cons(lambda)
temp <- optimise(tomin, lower=0, upper=1E4)
temp
lambda.opt.cons <- temp$minimum

##How long did tuning take?
proc.time()[3] - train.starttime + sum(elapsed.final)

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

##Plot results
xx <- seq(min(I1000), max(I1000), length.out=200)
aa <- sapply(xx, alpha.opt) ##tuned alpha
bb <- ifelse(xx<1000, 0.1, 1) ##ad-hoc alpha
zz <- predict(fit.resp, newdata=data.frame(I1000=xx), type="response") ##y predictions
uu <- predict(fit.tbar, data.frame(I1000=xx), type="response") ##T2 predictions
cairo_pdf(file="SIRout.pdf", width=6, height=5, pointsize=10)
par(mfrow=c(2,2), mar=c(4, 4.5, 2, 1.5) + 0.1 )
plot(I1000, train.obs, pch=".", xlab="I(1000)", ylab="y", col=gray(0.2), main="A")
lines(xx, 100*zz, lwd=2)
lines(xx, qbinom(p=0.025, size=100, prob=zz), lty=3, lwd=2)
lines(xx, qbinom(p=0.975, size=100, prob=zz), lty=3, lwd=2)
abline(h=Robs, lty=3)
plot(I1000, respT, pch=".", col=gray(0.2), xlab="I(1000)", ylab=expression("T"[2]), main="B")
lines(xx, uu, lwd=2)
plot(xx, hitPr(xx, Robs, myeps), type='l', xlab="I(1000)", ylab=expression(hat(gamma)), main="C")
plot(xx, aa, type='l', xlab="I(1000)", ylab=expression(alpha), yaxt="n", main="D")
axis(2, at=0:4/4, labels=c("0","0.25","0.5","0.75","1"))
lines(xx, bb, lty=2)
rug(res.lazy1$ABCsample$phi)
dev.off()
