## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----install, eval=FALSE, message=FALSE, warning=FALSE------------------------
#  install.packages("fastAFT")

## ----faft, eval=TRUE, message=FALSE, warning=FALSE----------------------------
## Mayo PBC data
library(survival)
pbc_analy <- as.matrix(na.omit(pbc[,c("time","status","age","edema","bili","albumin","protime")]))
# log transformation for time, bili, albumin, and protime
pbc_analy[,c(1,5:7)] <- log(pbc_analy[,c(1,5:7)])
colnames(pbc_analy)[c(1,5:7)] <- paste("log",colnames(pbc_analy)[c(1,5:7)])
# convert status to censoring indicator
pbc_analy[,2] <- pbc_analy[,2]>1

## Fast censored linear regression
# Gehan weight
library(fastAFT)
fit.g <- faft(pbc_analy[,1],pbc_analy[,2],pbc_analy[,-c(1,2)],weight="Gehan")
fit.g
# logrank weight
fit.l <- faft(pbc_analy[,1],pbc_analy[,2],pbc_analy[,-c(1,2)],weight="logrank")
fit.l

