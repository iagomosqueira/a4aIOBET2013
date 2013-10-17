


#
# a function to setup a simulated fishery conditioned on a current assessment
#

setupStock <- function(fmodel, srmodel, n1model, stock)
{

# set up dimension of problem
  dms <- dims(stock)

  years <- dms$minyear:dms$maxyear 
  ages <- dms$min:dms$max

  nages <- length(ages)
  nyears <- length(years)


# model data
  data <- expand.grid(age = ages, year = years)

# design matrices
  require(FLa4a)
  opts <- options(contrasts = c(unordered = "contr.sum", ordered = "contr.poly")) 
  Xf <- getX(fmodel, data)
  Xr <- getX(srmodel, subset(data, age == ages[1]))
  Xn1 <- getX(n1model, subset(data, age > ages[1] & year == years[1]))
  options(opts)

# model parameters
  bf <- c(solve(t(Xf) %*% Xf) %*% t(Xf) %*% log(c(harvest(bet))))
  br <- c(solve(t(Xr) %*% Xr) %*% t(Xr) %*% log(c(rec(bet))))
  bn1 <- c(solve(t(Xn1) %*% Xn1) %*% t(Xn1) %*% log(c(stock.n(bet)[-1,1])))

# m, wts and maturity
  M <- matrix(c(m(bet)[,1]), nages, nyears)
  Wt <- matrix(c(stock.wt(bet)[,1]), nages, nyears)
  Mat <- matrix(c(mat(bet)[,1]), nages, nyears)
  pz <- 0.6
  fbar <- c(4,6)

## model predictions
  Fm <- matrix(exp(Xf %*% bf), nages, nyears)

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


# construct catches
  sim.catch <- log(Fm) - log(Z) + log(1 - exp(-Z)) + n

  require(FLCore)
  stock <- FLStock(
             catch.n      = FLQuant(exp(sim.catch), dimnames=list(age=ages, year=years)),
             catch.wt     = FLQuant(Wt, dimnames=list(age=ages, year=years)),
             landings.n   = FLQuant(exp(sim.catch), dimnames=list(age=ages, year=years)),
             landings.wt  = FLQuant(Wt, dimnames=list(age=ages, year=years)),
             mat          = FLQuant(Mat, dimnames=list(age=ages, year=years)),
             m            = FLQuant(M, dimnames=list(age=ages, year=years)),
             harvest      = FLQuant(Fm, dimnames=list(age=ages, year=years), units = "f"),
             harvest.spwn = FLQuant(pz, dimnames=list(age=ages, year=years)),
             m.spwn       = FLQuant(pz, dimnames=list(age=ages, year=years)),
             stock.wt     = FLQuant(Wt, dimnames=list(age=ages, year=years)),
             stock.n      = FLQuant(exp(n), dimnames=list(age=ages, year=years)),

           )
  catch(stock) <- computeCatch(stock)
  landings(stock) <- computeLandings(stock)
  range(stock)[c("minfbar","maxfbar")] <- fbar
  range(stock)["plusgroup"] <- NA 


# return stuff
  list(fmodel = fmodel, srmodel = srmodel, n1model = n1model, stock = stock)

}





setupIndices <- function(stock, biomass = TRUE, qtype = "flat", pz = 0.6) {

# set up dimension of problem
  dms <- dims(stock)

  years <- dms$minyear:dms$maxyear 
  ages <- dms$min:dms$max

  nages <- length(ages)
  nyears <- length(years)

# define a set of catchability ogives
  qtype <- match.arg(qtype, c("dome","logistic","flat"))
  if (qtype == "dome") {
    qest <- rep(1e-5 * c(.1,.1,.9,1,1,.9,.5, 0.3, .3, .2), nyears)
  } else if (qtype == "logistic") {
    qest <- rep(1e-5 * c(.1,.1,.4,.5,.6,.9,1, 1, 1, 1), nyears)
  } else {
    qest <- rep(1e-5, nyears * nages)  
  }

# model data
  data <- expand.grid(age = ages, year = years)

# design matrices
  require(FLa4a)
  opts <- options(contrasts = c(unordered = "contr.sum", ordered = "contr.poly")) 
  Xq <- getX(qmodel[[1]], data)
  options(opts)

# model parameters
  bq <- c(solve(t(Xq) %*% Xq) %*% t(Xq) %*% log(qest))

## model predictions
  Qm <- matrix(exp(Xq %*% bq), nages, nyears)

# construct indices
  Z <- harvest(stock) + m(stock)
  if (biomass) {
    sim.index <- log( colSums(Qm * drop(stock.wt(stock) * stock.n(stock) * exp(-Z*pz))) )
  } else {
    sim.index <- log(Qm) + drop(stock.n(stock) - Z*pz)
  }

  indices <- FLIndices(simsurv = 
               FLIndex(
                 index = FLQuant(exp(sim.index), dimnames=list(age= if (biomass) "all" else ages, year=years))
               ))
  range(indices[[1]])[c("startf","endf")] <- pz

# return stuff
  list(qmodel = qmodel, indices = indices)

}

