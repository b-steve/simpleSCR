## Code in this file is similar to scr-ll.r, but includes
## model-fitting by RTMB.
library(RTMB)

## A function to calculate distances.
calc.dists <- function(points1, points2){
    apply(points1, 1, function(x) sqrt((x[1] - points2[, 1])^2 + (x[2] - points2[, 2])^2))
}

## Objects that need to be hanging around in the workspace:
## capt : A matrix of capture histories.
## dists: A matrix of distances between mask points and traps.
## a:     Area of a mask cell.
scr.nll <- function(pars){
    ## Unlinking parameters
    D <- exp(pars[1])
    g0 <- plogis(pars[2])
    sigma <- exp(pars[3])
    ## Number of animals detected.
    n <- nrow(capt)
    ## Detection probabilities.
    det.probs <- g0*exp(-dists^2/(2*sigma^2))
    ## Probability of non-detection at each mask cell.
    p.avoid <- apply(1 - det.probs, 1, prod)
    ## Probability of detection at each mask cell.
    p.det <- 1 - p.avoid
    ## Effective sampling area.
    esa <- a*sum(p.det)
    ## A vector to save likelihood contributions.
    f.capt <- numeric(n)
    ## Looping over capture histories.
    out <- 0
    for (i in 1:n){
        out <- out + log(sum(exp(apply(det.probs, 1, function(x) sum(dbinom(capt[i, ], 1, x, log = TRUE))))*a) + .Machine$double.xmin)
    }
    ## Contribution from the number of capture hisories, after some cancellation.
    out <- out + n*log(D) - D*esa - log(factorial(n))
    ##f.n <- dpois(n, D*esa, log = TRUE)
    -out
}

## Loading in data and getting the workspace ready.
load("test-data.RData")
capt <- test.data$bin.capt
dists <- calc.dists(test.data$traps, test.data$mask)
a <- attr(test.data$mask, "area")

## Start values.
par.start <- c(log(0.1), qlogis(0.5), log(50))
## Evaluating the log-likelhood at the start values.
scr.nll(par.start)

## Now trying it with RTMB.
obj <- MakeADFun(scr.nll, parameters = par.start)
## Evaluating the log-likelihood at the start values with the RTMB
## function.
obj$fn(par.start)
## Fitting model. Wow, so fast!
fit <- optim(par.start, obj$fn, gr = obj$gr)

## Unlinking estimates.
## D:
exp(fit$par[1])
## g0:
plogis(fit$par[2])
## sigma:
exp(fit$par[3])
## Negative log-likelihood at the estimates.
fit$value

## Let's try a one-at-a-time approach. Note that here "pars" includes
## both the actual parameters, followed by the capture history(!!).
n <- nrow(capt)
n.traps <- ncol(capt)
dcapt <- function(pars){
    ## Unlinking parameters
    D <- exp(pars[1])
    g0 <- plogis(pars[2])
    sigma <- exp(pars[3])
    capti <- pars[4:(n.traps + 3)]
    ## Number of animals detected.
    n <- nrow(capti)
    ## Detection probabilities.
    det.probs <- g0*exp(-dists^2/(2*sigma^2))
    ## Negative log-likelihood contribution.
    -log(sum(exp(apply(det.probs, 1, function(x) sum(dbinom(capti, 1, x, log = TRUE))))*a) + .Machine$double.xmin)
}

dc.tape <- MakeTape(dcapt, numeric(3 + n.traps))
data.tape <- MakeTape(function(i) DataEval(function(i) capt[i, ] , i), n.traps)
dcapts <- function(pars){
    out <- 0
    for (i in 1:n){
        capti <- data.tape(i)
        out <- out + dc.tape(c(pars, capti))
    }
    out
}
dcapts.tape <- MakeTape(dcapts, numeric(3))
dcapts.tape(par.start)

## Here's the contribution from the number of detected animals. It's
## not quite a dpois() because I've cancelled a constant with a term
## that comes from the contributions from capture histories.
dn <- function(pars){
    ## Unlinking parameters
    D <- exp(pars[1])
    g0 <- plogis(pars[2])
    sigma <- exp(pars[3])
    ## Number of animals detected.
    n <- nrow(capt)
    ## Detection probabilities.
    det.probs <- g0*exp(-dists^2/(2*sigma^2))
    ## Probability of non-detection at each mask cell.
    p.avoid <- apply(1 - det.probs, 1, prod)
    ## Probability of detection at each mask cell.
    p.det <- 1 - p.avoid
    ## Effective sampling area.
    esa <- a*sum(p.det)
    -n*log(D) + D*esa + log(factorial(n))
}
dn.tape <- MakeTape(dn, numeric(3))
nll.fun <- function(pars){
    dn.tape(pars) + dcapts.tape(pars)
}
nll.tape <- MakeTape(nll.fun, numeric(3))

optim(par.start, nll.tape, nll.tape$jacobian)


