

doSimulation <- function(fmodel, qmodel, srmodel, n1model, stock, indices, nsim = 10) {

  simstock <- stock
  simindices <- indices

  simout <- vector("list", nsim)
  for (s. in 1:nsim) {

## make observations
    catch.n(simstock) <- catch.n(stock) * exp( rnorm(length(catch.n(stock)), 0, exp(0.1)) )
    catch(simstock) <- computeCatch(simstock)
    index(simindices[[1]]) <- index(indices[[1]]) + rnorm(length(index(indices[[1]])), 0, exp(0.1))


## fit model
    out <- a4aInternal(fmodel, qmodel, srmodel, n1model, 
                       stock = simstock, indices = simindices, center = FALSE)

## assess the fit
    est <-
      rbind(pars(out)@qmodel[[1]]@params@.Data,
            pars(out)@vmodel[[1]]@params@.Data,
            pars(out)@vmodel[[2]]@params@.Data,
            pars(out)@stkmodel@params@.Data)

    summ <- as.data.frame(est)
    names(summ) <- "est"
    summ $ parname <- rownames(summ)
    summ $ type <- sapply(strsplit(summ $ parname, ":"), "[", 1)
    rownames(summ) <- NULL

    simout[[s.]] <- summ
  }

  do.call(rbind, simout)
}

#summ $ sim <- c(bq, lsigc, lsigi, bf, bn1, br)
#summ $ diff <- with(summ, (est - sim))
#summ $ sdiff <- 100 * with(summ, diff / abs(sim))




