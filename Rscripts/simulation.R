
# simulation set up

nsim <- 100

years <- 1990:2010
ages <- 1:10

nages <- length(ages)
nyears <- length(years)

# model data
data <- expand.grid(age = ages, year = years)

# model formulas
fmodel <- ~ 1
qmodel <- list(~ 1)
#qmodel <- list(~ factor(as.numeric(age > 2)))
srmodel <- ~ 1
n1model <- ~ 1

# design matrices
opts <- options(contrasts = c(unordered = "contr.sum", ordered = "contr.poly")) 
Xf <- model.matrix(fmodel, data)
Xq <- model.matrix(qmodel[[1]], data)
Xr <- model.matrix(srmodel, subset(data, age == ages[1]))
Xn1 <- model.matrix(n1model, subset(data, age > ages[1] & year == years[1]))
options(opts)

# model parameters
bf <- log(0.4)
bq <- log(0.001)
#bq <- log(c(0.0000001, 0.001))
br <- log(10000)
bn1 <- log(1000)
lsigc <- log(0.1)
lsigi <- log(0.1)

# m, wts and maturity
M <- matrix(0.2, nages, nyears)
Wt <- matrix(1, nages, nyears)
Mat <- matrix(rep(0:1, c(2, nages-2)), nages, nyears)
pz <- 0.5

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
n[nages,-1] <- log( exp(n[nages,-1]) + exp(n[nages,-nyears] - Z[i-1, -nyears]) )

# construct catches and indices
catch <- log(Fm) - log(Z) + log(1 - exp(-Z)) + n
index <- log( colSums(Wt * Qm * exp(n) * exp(-Z*pz)) )


simout <- vector("list", nsim)
for (s. in 1:nsim) {

## make observations
catchobs <- catch + rnorm(length(catch), 0, exp(lsigc))
indexobs <- index + rnorm(length(index), 0, exp(lsigi))

## build FLR objects
require(FLCore)
stock <- FLStock(
           catch.n      = FLQuant(exp(catchobs), dimnames=list(age=ages, year=years)),
           catch.wt     = FLQuant(Wt, dimnames=list(age=ages, year=years)),
           mat          = FLQuant(Mat, dimnames=list(age=ages, year=years)),
           m            = FLQuant(M, dimnames=list(age=ages, year=years)),
           harvest      = FLQuant(NA, dimnames=list(age=ages, year=years), units = "f"),
           harvest.spwn = FLQuant(0.6, dimnames=list(age=ages, year=years)),
           m.spwn       = FLQuant(0.6, dimnames=list(age=ages, year=years)),
           stock.wt     = FLQuant(Wt, dimnames=list(age=ages, year=years))
         )
catch(stock) <- computeCatch(stock)
range(stock)[c("minfbar","maxfbar")] <- c(4,9)

indices <- FLIndices(simsurv = 
             FLIndex(
               index = FLQuant(exp(indexobs), dimnames=list(age="all", year=years))
             ))
range(indices[[1]])[c("startf","endf")] <- 0.5

## fit model
covar = NULL
wkdir = NULL
verbose = FALSE
fit = "assessment"
vmodel = list(~1, ~1)

source("fitting_code_snippet.R")

## assess the fit

est <- out $ par.est[1:5]

sim <- list(fpar = bf, qpar = bq, vpar = c(lsigc, lsigi), ny1par = bn1, rpar = br)

summ <- cbind.data.frame(sim = unlist(sim), est = unlist(est))
summ $ diff <- with(summ, (est - sim))
summ $ sdiff <- with(summ, diff / abs(sim) )
summ $ parname <- rownames(summ)

simout[[s.]] <- summ
}

boxplot(sdiff ~ parname, data = do.call(rbind, simout), type = "g")



