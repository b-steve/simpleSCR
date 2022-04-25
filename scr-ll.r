## A required package.
library(spatstat)

## Computes the negative log-likelihood function for a simple SCR
## model.
##
## Arguments are as follows:
##
## pars: A vector of parameters at which to compute the
##       log-likelihood. They should appear in the following order:
##       (1) log(D), (2) logit(g0), (3) log(sigma). Note that density,
##       D, is animals per hectare.
##
## capt: Matrix of capture histories.
## 
## traps: Matrix containing x- and y-coordinates of trap locations.
##
## mask: Matrix containing x- and y-coordinates of mask point
##       locations. Needs to have attribute 'area', providing the area
##       of a single pixel.
##
## Note that trap and mask coordinates are given in metres.
scr.nll <- function(pars, capt, traps, mask){
    ## Unlinking parameters
    D <- exp(pars[1])
    g0 <- plogis(pars[2])
    sigma <- exp(pars[3])
    ## Number of animals detected.
    n <- nrow(capt)
    ## Number of traps.
    n.traps <- nrow(traps)
    ## Number of mask points.
    n.mask <- nrow(mask)
    ## Area of a single mask pixel.
    a <- attr(mask, "area")
    ## Constructing a distance matrix. The element (i, j) gives the
    ## distance between the ith mask point and the jth trap.
    mask.dists <- crossdist(mask[, 1], mask [, 2],
                            traps[, 1], traps[, 2])
    ## Constructing a detection probability matrix. The element (i, j)
    ## gives the probability of an animal located at the ith mask
    ## point being detected at the jth trap.
    mask.probs <- g0*exp(-mask.dists^2/(2*sigma^2))
    ## Constructing a detection probability vector. The ith element
    ## gives the probability of an animal located at the ith mask
    ## point being detected by *at least one* trap.
    p.avoid <- apply(1 - mask.probs, 1, prod)
    p.det <- 1 - p.avoid
    ## Calculating the effective sampling area.
    esa <- a*sum(p.det)
    ##Calculating likelihood contribution due to each
    ## detected animal's capture history.
    f.capt <- numeric(n)
    for (i in 1:n){
        ## Calculating the log-probability of the animal's capture
        ## history, conditional on being at the jth mask point.
        log.f.capt.given.s <- numeric(n.mask)
        for (j in 1:n.mask){
            log.f.capt.given.s[j] <- sum(dbinom(capt[i, ], 1, mask.probs[j, ], log = TRUE))
        }
        ## Summing the integrand over all mask points.
        f.capt[i] <- sum(exp(log.f.capt.given.s))/esa
    }
    ## Log-likelihood contribution from all capture histories
    ## calculated by the log of the sum of the individual likelihood
    ## contributions.
    log.f.capt <- sum(log(f.capt + .Machine$double.xmin))
    ## Log-likelihood contribution from the number of animals
    ## detected.
    log.f.n <- dpois(n, D*esa, log = TRUE)
    ## Overall log-likelihood. The last part accounts for the fact
    ## that we cannot observe an all-zero capture history.
    ll <- log.f.n + log.f.capt
    ## Returning negative log-likelihood.
    -ll
}

## Loading some data to test things out.
load("test-data.RData")

## Some start values: D = 0.1, g0 = 0.5, sigma = 50.
par.start <- c(log(0.1), qlogis(0.5), log(50))

## Fitting the model.
fit <- optim(par.start, scr.nll, capt = test.data$bin.capt, traps = test.data$traps, mask = test.data$mask)

## Unlinking estimates.
## D:
exp(fit$par[1])
## g0:
plogis(fit$par[2])
## sigma:
exp(fit$par[3])

fit$value
