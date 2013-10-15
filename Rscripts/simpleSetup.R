require(FLa4a)

years <- 1990:2010
ages <- 1:5

nages <- length(ages)
nyears <- length(years)

# model data
data <- expand.grid(age = ages, year = years)

# model formulas
fmodel <- ~ I(year-2000)
qmodel <- list(~ 1)
#qmodel <- list(~ factor(as.numeric(age > 2)))
srmodel <- ~ 1
n1model <- ~ 1

# design matrices
opts <- options(contrasts = c(unordered = "contr.sum", ordered = "contr.poly")) 
Xf <- getX(fmodel, data)
Xq <- getX(qmodel[[1]], data)
Xr <- getX(srmodel, subset(data, age == ages[1]))
Xn1 <- getX(n1model, subset(data, age > ages[1] & year == years[1]))
options(opts)

# model parameters
bf <- c(log(0.4), 0.1)
bq <- log(0.001)
#bq <- log(c(0.0000001, 0.001))
br <- log(10000)
bn1 <- log(1000)
lsigc <- log(0.05)
lsigi <- log(0.1)

# m, wts and maturity
M <- matrix(0.2, nages, nyears)
Wt <- matrix(1, nages, nyears)
Mat <- matrix(rep(0:1, c(2, nages-2)), nages, nyears)
pz <- 0.5
fbar <- c(2,4)


