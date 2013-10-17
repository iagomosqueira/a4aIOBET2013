
require(FLCore)
load('../data/bet_iotc_2011.RData')
range(bet)[c("minfbar","maxfbar")] <- c(4,9)
units(harvest(bet)) <- "f"


load('../data/betSS3.RData')


bet <- window(bet, 1952, 2011)

dms <- dims(bet)

years <- dms$minyear:dms$maxyear 
ages <- dms$min:dms$max

nages <- length(ages)
nyears <- length(years)

# for now make up some F values and N values
harvest[harvest == 0] <- min(harvest[harvest > 0])
harvest(bet) <- c(harvest)

stock.n(bet) <- c(stock.n[,-ncol(stock.n)])

# define a set of catchability ogives for the biomass survey
if (qtype == "dome") {
  qest <- rep(1e-5 * c(.1,.1,.9,1,1,.9,.5, 0.3, .3, .2), nyears)
} else if (qtype == "logistic") {
  qest <- rep(1e-5 * c(.1,.1,.4,.5,.6,.9,1, 1, 1, 1), nyears)
} else {
  qest <- rep(1e-5, nyears * nages)  
}


# model data
data <- expand.grid(age = ages, year = years)

# model formulas
fmodel <- ~ s(age, k = 4) + s(year, k = 20)
srmodel <- ~ factor(year)
n1model <- ~ factor(age)
qmodel <- list(~ s(age, k=4))
if (biomass) {
  a50 <- 4
  a90 <- 5

  a <- log(1/9) * a50 / (a90 - a50)
  b <- - a / a50
  fixedq <- 1/(1 + exp(-a - b*ages))
  
  qmodel <- list( ~ 1 )
}

# design matrices
require(FLa4a)
opts <- options(contrasts = c(unordered = "contr.sum", ordered = "contr.poly")) 
Xf <- getX(fmodel, data)
Xq <- getX(qmodel[[1]], data)
Xr <- getX(srmodel, subset(data, age == ages[1]))
Xn1 <- getX(n1model, subset(data, age > ages[1] & year == years[1]))
options(opts)

# model parameters
bf <- c(solve(t(Xf) %*% Xf) %*% t(Xf) %*% log(c(harvest(bet))))
bq <- c(solve(t(Xq) %*% Xq) %*% t(Xq) %*% log(qest))
br <- c(solve(t(Xr) %*% Xr) %*% t(Xr) %*% log(c(rec(bet))))
bn1 <- c(solve(t(Xn1) %*% Xn1) %*% t(Xn1) %*% log(c(stock.n(bet)[-1,1])))
lsigc <- log(0.1)
lsigi <- log(0.1)

# m, wts and maturity
M <- matrix(c(m(bet)[,1]), nages, nyears)
Wt <- matrix(c(stock.wt(bet)[,1]), nages, nyears)
Mat <- matrix(c(mat(bet)[,1]), nages, nyears)
pz <- 0.6
fbar <- c(4,6)


