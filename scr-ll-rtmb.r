## Code in this file is similar to scr-ll.r, but includes
## model-fitting by RTMB
library(TMB)
compile("scr_nll.cpp")
library(RTMB)

## A function to calculate distances.
calc.dists <- function(points1, points2){
    apply(points1, 1, function(x) sqrt((x[1] - points2[, 1])^2 + (x[2] - points2[, 2])^2))
}

###############################################################
## Beginning of functions for a regular RTMB implementation. ##
###############################################################
scr.nll <- function(pars, capt, dists, a){
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
## A closure so that we can specify data explicitly.
scr.nll.closure <- function(capt, dists, a){
    function(pars) scr.nll(pars, capt, dists, a)
}

## A function to fit the model with regular RTMB.
fit.scr.rtmb <- function(capt, traps, mask, par.start){
    ## Computing distances.
    dists <- calc.dists(traps, mask)
    ## Area of a mask cell.
    a <- attr(test.data$mask, "area")
    message("Running MakeADFun...")
    mad.time <- system.time({obj <- MakeADFun(scr.nll.closure(capt, dists, a), par.start)})
    message("Completed in ", round(mad.time[3], 2), " seconds.\n")
    message("Fitting model...")
    fit.time <- system.time({fit <- optim(par.start, obj$fn, gr = obj$gr)})
    message("Completed in ", round(fit.time[3], 2), " seconds.\n")
    c(D = exp(fit$par[1]), g0 =  plogis(fit$par[2]), sigma = exp(fit$par[3]))
}

#############################################################
## Beginning of functions for a taped RTMB implementation. ##
#############################################################
dprob <- function(pars, dists){
    ## Unlinking parameters
    D <- exp(pars[1])
    g0 <- plogis(pars[2])
    sigma <- exp(pars[3])
    ## Computing probabilities and turning into a matrix.
    c(g0*exp(-dists^2/(2*sigma^2)))
}
## A closure so that we specify distances explicitly.
dprob.closure <- function(dists){
    function(pars) dprob(pars, dists)
}

## A function from which to create a tape for the likelihood
## contribution of a capture history.
dcapt <- function(pars, dists, a, n.traps, n.mask){
    ## Obtaining the detection probabilities.
    det.probs <- matrix(pars[1:(n.traps*n.mask)], nrow = n.mask)
    ## The capture history.
    capti <- pars[(n.traps*n.mask + 1):(n.traps*n.mask + n.traps)]
    ## Negative log-likelihood contribution.
    -log(sum(exp(apply(det.probs, 1, function(x) sum(dbinom(capti, 1, x, log = TRUE))))*a) + .Machine$double.xmin)
}
## A closure so that we can specify data explicitly.
dcapt.closure <- function(dists, a, n.traps, n.mask){
    function(pars) dcapt(pars, dists, a, n.traps, n.mask)
}

## A function from which to create a tape for the likelihood
## contributions of all capture histories, based on data and
## individual capture history tapes.
dcapts <- function(pars, n, data.tape, dcapt.tape, dprob.tape){
    out <- 0
    for (i in 1:n){
        out <- out + dcapt.tape(c(pars, data.tape(i)))
    }
    out
}
## A closure so that we can specify n and the tapes explicitly.
dcapts.closure <- function(n, data.tape, dcapt.tape, dprob.tape){
    function(pars) dcapts(pars, n, data.tape, dcapt.tape, dprob.tape)
}

## A function from which to create a tape for the contribution of n.
dn <- function(pars, n, dists, a, n.traps, n.mask){
    ## Unlinking parameters.
    D <- exp(pars[1])
    ## Detection probabilities.
    det.probs <- matrix(pars[2:(n.traps*n.mask + 1)], nrow = n.mask)
    ## Probability of non-detection at each mask cell.
    p.avoid <- apply(1 - det.probs, 1, prod)
    ## Probability of detection at each mask cell.
    p.det <- 1 - p.avoid
    ## Effective sampling area.
    esa <- a*sum(p.det)
    -n*log(D) + D*esa + log(factorial(n))
}
## A closure so that we can specify the data explicitly.
dn.closure <- function(n, dists, a, n.traps, n.mask){
    function(pars) dn(pars, n, dists, a, n.traps, n.mask)
}

## A function from which to compute the likelihood.
nll.pretape <- function(pars, dn.tape, dcapts.tape, dprob.tape){
    det.probs <- dprob.tape(pars)
    dn.tape(c(pars[1], det.probs)) + dcapts.tape(det.probs)
}
## A closure so that we can specify the tapes explicitly.
nll.closure <- function(dn.tape, dcapts.tape, dprob.tape){
    function(pars) nll.pretape(pars, dn.tape, dcapts.tape, dprob.tape)
}

## A function to fit the model with a "regular" tape setup.
fit.scr.tape <- function(capt, traps, mask, par.start, optimize.tape = TRUE){
    ## Number of capture histories.
    n <- nrow(capt)
    ## Computing distances.
    dists <- calc.dists(traps, mask)
    ## Area of a mask cell.
    a <- attr(test.data$mask, "area")
    ## Number of detectors.
    n.traps <- ncol(capt)
    ## Number of mask points.
    n.mask <- nrow(mask)
    message("Making tapes...")
    tape.time <- system.time({
        ## Making a tape for the detection probabilities.
        dprob.tape <- MakeTape(dprob.closure(dists), numeric(3))
        ## Making tape for a single capture history's likelihood
        ## contribution, provided the detection probabilities.
        dcapt.tape <- MakeTape(dcapt.closure(dists, a, n.traps, n.mask), numeric(n.traps*n.mask + n.traps))
        ## Making a data tape.
        data.tape <- MakeTape(function(i) DataEval(function(i) capt[i, ] , i), 1)
        ## Making a tape for the sum of likelihood contributions for all
        ## capture histories.
        dcapts.tape <- MakeTape(dcapts.closure(n, data.tape, dcapt.tape, dprob.tape), numeric(n.traps*n.mask))
        ## Making a tape for the likelihood contribution of n.
        dn.tape <- MakeTape(dn.closure(n, dists, a, n.traps, n.mask), numeric(n.traps*n.mask + 1))
        ## Making a tape for the full likelihood.
        nll.tape <- MakeTape(nll.closure(dn.tape, dcapts.tape, dprob.tape), numeric(3))
        if (optimize.tape){
            nll.tape$simplify("optimize")
        }
    })
    message("Completed in ", round(tape.time[3], 2), " seconds.\n")
    message("Fitting model...")
    fit.time <- system.time({fit <- optim(par.start, nll.tape, gr = nll.tape$jacobian)})
    message("Completed in ", round(fit.time[3], 2), " seconds.\n")
    c(D = exp(fit$par[1]), g0 =  plogis(fit$par[2]), sigma = exp(fit$par[3]))
}

########################################################################
## Beginning of functions for a one-by-one taped RTMB implementation. ##
########################################################################

## One-by-one likelihood function.
nll.obo <- function(pars, capt, n.mask, dprob.tape, dcapt.tape, dn.tape){
    ## Number of capture histories.
    n <- nrow(capt)
    ## Detection probabilities.
    det.probs <- dprob.tape(pars)
    ## Looping over capture histories.
    out <- 0
    for (i in 1:n){
        out <- out + dcapt.tape(c(det.probs, capt[i, ]))
    }
    out <- out + dn.tape(c(pars[1], det.probs))
    out
}

## Model-fitting. Note that we're not actually using gradients here
## because we don't have an all-in-one tape.
fit.scr.tape.obo <- function(capt, traps, mask, par.start, optimize.tape = TRUE){
    ## Number of capture histories.
    n <- nrow(capt)
    ## Computing distances.
    dists <- calc.dists(traps, mask)
    ## Area of a mask cell.
    a <- attr(test.data$mask, "area")
    ## Number of detectors.
    n.traps <- ncol(capt)
    ## Number of mask points.
    n.mask <- nrow(mask)
    message("Making tapes...")
    tape.time <- system.time({
        ## Making a tape for the detection probabilities.
        dprob.tape <- MakeTape(dprob.closure(dists), numeric(3))
        ## Making tape for a single capture history's likelihood
        ## contribution, provided the detection probabilities.
        dcapt.tape <- MakeTape(dcapt.closure(dists, a, n.traps, n.mask), numeric(n.traps*n.mask + n.traps))
        ## Making a tape for the likelihood contribution of n.
        dn.tape <- MakeTape(dn.closure(n, dists, a, n.traps, n.mask), numeric(n.traps*n.mask + 1))
        if (optimize.tape){
            dprob.tape$simplify("optimize")
            dcapt.tape$simplify("optimize")
            dn.tape$simplify("optimize")
        }
    })
    message("Completed in ", round(tape.time[3], 2), " seconds.\n")
    message("Fitting model...")
    fit.time <- system.time({fit <- optim(par.start, nll.obo, capt = capt, n.mask = n.mask,
                                          dprob.tape = dprob.tape, dcapt.tape = dcapt.tape,
                                          dn.tape = dn.tape)})
    message("Completed in ", round(fit.time[3], 2), " seconds.\n")
    c(D = exp(fit$par[1]), g0 =  plogis(fit$par[2]), sigma = exp(fit$par[3]))
}

########################################################################
## Beginning of regular TMB implementation. ##
########################################################################

fit.scr.tmb <- function(capt, traps, mask, par.start){
    ## Number of capture histories.
    n <- nrow(capt)
    ## Computing distances.
    dists <- calc.dists(traps, mask)
    ## Area of a mask cell.
    a <- attr(test.data$mask, "area")
    ## Number of detectors.
    n.traps <- ncol(capt)
    ## Number of mask points.
    n.mask <- nrow(mask)
    ## Loading .so object.
    message("Running MakeADFun...")
    mad.time <- system.time({
        dyn.load("scr_nll.so")
        scr.obj <- TMB::MakeADFun(data = list(capt = capt,
                                              mask_dists = dists,
                                              n = n,
                                              n_traps = n.traps,
                                              n_mask = n.mask,
                                              mask_area = a),
                                  parameters = list(pars = par.start),
                                  DLL = "scr_nll")
    })
    message("Completed in ", round(mad.time[3], 2), " seconds.\n")
    message("Fitting model...")
    fit.time <- system.time({fit <- optim(par.start, scr.obj$fn, gr = scr.obj$gr)})
    message("Completed in ", round(fit.time[3], 2), " seconds.\n")
    c(D = exp(fit$par[1]), g0 =  plogis(fit$par[2]), sigma = exp(fit$par[3]))
}

## Loading in data.
load("test-data.RData")
capt <- test.data$bin.capt
traps <- test.data$traps
mask <- test.data$mask
## Start values for optimisation.
par.start <- c(log(0.1), qlogis(0.5), log(50))
## Fitting model using regular TMB.
fit.scr.rtmb(capt, traps, mask, par.start)
## Doubling the amount of data doubles MakeADFun() time.
fit.scr.rtmb(rbind(capt, capt), traps, mask, par.start)

## Fitting model with the taped version.
fit.scr.tape(capt, traps, mask, par.start)
## Same again: doubling data doubles tape-making time.
fit.scr.tape(rbind(capt, capt), traps, mask, par.start)
## Is it worth tape optimization?
fit.scr.tape(capt, traps, mask, par.start, optimize.tape = FALSE)
fit.scr.tape(rbind(capt, capt), traps, mask, par.start, optimize.tape = FALSE)

## Fitting model with the one-by-one taped version.
fit.scr.tape.obo(capt, traps, mask, par.start)
fit.scr.tape.obo(capt, traps, mask, par.start, optimize.tape = FALSE)

## Fitting the model with the TMB version.
fit.scr.tmb(capt, traps, mask, par.start)
## Doubling the amount of data.
fit.scr.tmb(rbind(capt, capt), traps, mask, par.start)


## Direct comparison of RTMB and TMB.
system.time({fit.scr.rtmb(capt, traps, mask, par.start)})
system.time({fit.scr.tmb(capt, traps, mask, par.start)})

## Same again with three copies of the data (42 capture histories),
system.time({fit.scr.rtmb(rbind(capt, capt, capt), traps, mask, par.start)})
system.time({fit.scr.tmb(rbind(capt, capt, capt), traps, mask, par.start)})

## Same with ten copies of the data (140 capture histories).
system.time({fit.scr.rtmb(rbind(capt, capt, capt, capt, capt, capt, capt, capt, capt, capt),
                          traps, mask, par.start)})
system.time({fit.scr.tmb(rbind(capt, capt, capt, capt, capt, capt, capt, capt, capt, capt),
                         traps, mask, par.start)})
