# a4aBET.R - DESC
# a4aBET.R

# Copyright 2003-2013 FLR Team. Distributed under the GPL 2 or later
# Maintainer: Iago Mosqueira, JRC
# Soundtrack:
# Notes:

# This is what I managed to get started with today

# - Changes in srmodel (geomean vs. beholt, CV value) appears to have a large effect

# - I am not clear what I am doing with regards to selectivity of CPUE. It is a LL, so
# should fully select ages 4-9, or maybe peak at 4-5 and then lower. I would like to
# test those two alternatives.

# - As you can see the CPUE is crazy, so no wonder any SA struggles


library(FLa4a)
library(ggplotFL)

# BET FLStock and JAP LL FLIndex
load('data/bet_iotc_2011.RData')

# CHANGE Fully selected ages
range(bet)[c("minfbar","maxfbar")] <- c(4,9)

# CREATE FLindices using JAPLL
idx <- FLIndices(japll=japll)

# run1 - Constant q
run1 <- a4a(fmodel = ~ s(age, k = 3) + s(year, k = 10), qmodel = list( ~1), stock = bet,
	indices = idx)

c(log(index(run1)[[1]] / index(idx[[1]])))
plot((bet + run1))

# run2 - q by year
run2 <- a4a(fmodel = ~ s(age, k = 3) + s(year, k = 10), qmodel = list(~factor(year)), stock = bet, indices = idx)

c(log(index(run2)[[1]] / index(idx[[1]])))
plot((bet + run2))

# run3 - geomean SR
run3 <- a4a(fmodel = ~ s(age, k = 3) + s(year, k = 10), qmodel = list(~factor(year)), stock = bet, indices = idx, srmodel = ~ geomean(CV=5))

log(index(run3)[[1]] / index(idx[[1]]))
plot((bet + run3))

# run4 - 
run4 <- a4a(fmodel = ~te(age, year, k=c(4, 10)), qmodel = list(~factor(year)), stock = bet, indices = idx, srmodel = ~ bevholt(CV=0.5))

log(index(run4)[[1]] / index(idx[[1]]))
plot((bet + run4))

# run5 - 
run4 <- a4a(fmodel = ~ s(age, k=6) + factor(year), qmodel = list(japll=~factor(year)), stock = bet, indices = idx, srmodel = ~ bevholt(CV=0.5))

log(index(run4)[[1]] / index(idx[[1]]))
plot((bet + run4))
