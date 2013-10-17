
# simulation set up

## model predictions
Fm <- matrix(exp(Xf %*% bf), nages, nyears)
Qm <- matrix(exp(Xq %*% bq), nages, nyears)

# build n at age
n <- matrix(NA, nages, nyears)
n[1,] <- c(Xr %*% br)
n[-1,1] <- c(Xn1 %*% bn1)
Z <- Fm + M
for (i in 2:nages) {
  n[i,-1] <- n[i-1, -nyears] - Z[i-1, -nyears]
}
# do plus group
#for (j in 2:nyears) n[nages,j] <- log( exp(n[nages,j]) + exp(n[nages,j-1] - Z[nages, j-1]) )


# construct catches and indices
sim.catch <- log(Fm) - log(Z) + log(1 - exp(-Z)) + n
if (biomass) {
  sim.index <- log( colSums(Wt * Qm * exp(n) * exp(-Z*pz)) )
} else {
  sim.index <- log(Qm) + n - Z*pz
}

simout <- vector("list", nsim)
for (s. in 1:nsim) {

## make observations
catchobs <- sim.catch + rnorm(length(sim.catch), 0, exp(lsigc))
indexobs <- sim.index + rnorm(length(sim.index), 0, exp(lsigi))

## build FLR objects
require(FLCore)
stock <- FLStock(
           catch.n      = FLQuant(exp(catchobs), dimnames=list(age=ages, year=years)),
           catch.wt     = FLQuant(Wt, dimnames=list(age=ages, year=years)),
           mat          = FLQuant(Mat, dimnames=list(age=ages, year=years)),
           m            = FLQuant(M, dimnames=list(age=ages, year=years)),
           harvest      = FLQuant(NA, dimnames=list(age=ages, year=years), units = "f"),
           harvest.spwn = FLQuant(pz, dimnames=list(age=ages, year=years)),
           m.spwn       = FLQuant(pz, dimnames=list(age=ages, year=years)),
           stock.wt     = FLQuant(Wt, dimnames=list(age=ages, year=years))
         )
catch(stock) <- computeCatch(stock)
range(stock)[c("minfbar","maxfbar")] <- fbar
range(stock)["plusgroup"] <- NA 

indices <- FLIndices(simsurv = 
             FLIndex(
               index = FLQuant(exp(indexobs), dimnames=list(age= if (biomass) "all" else ages, year=years))
             ))
range(indices[[1]])[c("startf","endf")] <- pz

## fit model
vmodel = list(~1, ~1)

out <- a4aInternal(fmodel, qmodel, srmodel, n1model, 
                   stock = stock, indices = indices, center = FALSE)

## assess the fit
est <-
  rbind(pars(out)@qmodel[[1]]@params@.Data,
        pars(out)@vmodel[[1]]@params@.Data,
        pars(out)@vmodel[[2]]@params@.Data,
        pars(out)@stkmodel@params@.Data)

summ <- as.data.frame(est)
names(summ) <- "est"

summ $ sim <- c(bq, lsigc, lsigi, bf, bn1, br)
summ $ diff <- with(summ, (est - sim))
summ $ sdiff <- 100 * with(summ, diff / abs(sim))
summ $ parname <- rownames(summ)
summ $ type <- sapply(strsplit(summ $ parname, ":"), "[", 1)
rownames(summ) <- NULL

simout[[s.]] <- summ
}

# summarise sims
sims <- do.call(rbind, simout)



