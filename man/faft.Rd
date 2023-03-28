\name{faft}
\alias{faft}
\title{Fast censored linear regression for the accelerated failure time (AFT)
  model}
\description{An implementation of the fast censored linear regression in
  Huang (2013).}
\usage{faft(x,dlt,z,weight="logrank",ynci=0,epl=0.95,epu=0.05)}
\arguments{
  \item{x}{follow-up time.}
  \item{dlt}{censoring indicator: 1 - event, 0 - censored.}
  \item{z}{matrix of covariates: each column corresponds to a covariate.}
  \item{weight}{either "logrank" or "Gehan" estimating function.}
  \item{ynci}{compute test inversion-based 95\% CI's? 1 - yes, 0 - no.}
  \item{epl}{parameter in (0,1) for determining the lower quantile from
  censored quantile regression (Huang 2010) for
  the preparatory estimation: sum of squared covariates for at-risk uncensored
  individuals is about $epl^(dim(z)[2]+1)$ in determinant.}
  \item{epu}{parameter in (0,1) for determining the upper quantile from
  censored quantile regression (Huang 2010) for
  the preparatory estimation: sum of squared covariates for at-risk uncensored
  individuals is about $epu^(dim(z)[2]+1)$ in determinant.}
}
\value{
  \item{weight}{either "logrank" or "Gehan" estimating function.}
  \item{beta}{estimated regression coefficient (the proposed).}
  \item{va}{sandwich variance estimate for beta.}
  \item{qif}{quadratic score statistic at beta.}
  \item{ci95}{test inversion-based 95\% CI's, only available if requested and
  successful.}
  \item{message}{point estimation: "success", "error - algorithm fails",
  or "warning - singular hessian".}
  \item{imsg}{numerical code for point and test inversion-based interval
  estimation: 0 - success in point and interval, 1 - error in point where
  algorithm fails, 2 - warning in point with singular hessian, 3 - success
  in point but failure in interval.}
  \item{beta1stp}{the one-step estimator.}
  \item{qif1stp}{quadratic score statistic at beta1stp.}
  \item{betainit}{the initial estimator.}
  \item{qifinit}{quadratic score statistic at betainit.}
}
\examples{
## simulate a dataset of size 100 with 2 covariates
size <- 100
npred <- 2
beta <- rep(1,npred)

cvt <- matrix(rnorm(size*npred),ncol=npred)
resid <- log(rexp(size))
event.t <- resid + cvt \%*\% beta
censr.t <- log(runif(size, 0, 6))
x <- pmin(event.t, censr.t)
dlt <- as.numeric(event.t<=censr.t)

## run censored linear regression
fit.g <- faft(x,dlt,cvt,weight="Gehan")
fit.l <- faft(x,dlt,cvt,weight="logrank")
}

\references{Huang, Y. (2010) Quantile calculus and censored
  regression, \emph{The Annals of Statistics} 38, 1607--1637.

  Huang, Y. (2013) Fast censored linear regression. \emph{Scandinavian Journal
  of Statistics} 40, 789--806.}
\author{Yijian Huang}
