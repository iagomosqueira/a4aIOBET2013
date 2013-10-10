# load.R - DESC
# load.R

# Copyright 2003-2013 FLR Team. Distributed under the GPL 2 or later
# Maintainer: Iago Mosqueira, JRC
# Soundtrack:
# Notes:

library(FLCore)
library(reshape2)
library(plyr)
library(ggplotFL)

# LOAD .dat

caa <- read.table("caa.dat", sep="\t", header=TRUE)
nc <- read.table("nc.dat", sep="\t", header=TRUE)

yrs <- sort(unique(caa$Year))

# CREATE catch

catch <- FLQuant(nc$NCmtFish, dimnames=list(year=yrs), quant='age', units='t')

# CREATE catch.n

caam <- melt(caa, id.vars = c("Method","ID","Species","Fishery","Year","Quarter","Tno"))
caam <- ddply(caam, c("Year","variable"), summarise, sum = sum(value))
caam <- acast(caam, Year ~ variable)

catch.n <- FLQuant(t(caam), dimnames=list(age=0:9, year=yrs), units='num')
catch.n[catch.n == 0] <- 1e-8

# CREATE catch.wt

# From IOTC–2011–WPTT13–41
Linf <- 169.06
K <- 0.32
t0 <- -0.34
length.n <- Linf * (1 - exp(-K*((0:9) - t0)))

catch.wt <- FLQuant(3.661 * 1e-5 * length.n ^ 2.901, dimnames=list(age=0:9,
	year=yrs), units='kg')

# CREATE m

m <- FLQuant(c(0.8, 0.8, rep(0.4, 8)), dimnames=list(age=0:9, year=yrs))

# CREATE mat

mat <- FLQuant(1 / (1 + exp(-0.25 * (length.n -110.888))), dimnames=list(age=0:9,
	year=yrs))

# CREATE FLStock

bet <- FLStock(name='BET', desc='Indian Ocean BET 2011',
	catch=catch,
	catch.n=catch.n,
	catch.wt=catch.wt,
	landings=catch,
	landings.n=catch.n,
	landings.wt=catch.wt,
	stock.wt=catch.wt,
	m=m,
	mat=mat)

# ASSIGN m.spwn & harvest.spwn

m.spwn(bet) <- 0.6
harvest.spwn(bet) <- 0.6

# LOAD cpue data

cpue <- read.table("cpue.dat", header=T, sep="\t")

index <- FLQuant(cpue$trop, dimnames=list(year=yrs), quant='age')

japll <- FLIndex(name='JAPLL', desc='STD. CPUE JAL LL BET IOTC–2011–WPTT13–52',
	index=index)

range(japll)[c("startf","endf")] <- 0.5

# SAVE file

save(bet, japll, file='bet_iotc_2011.RData')
