faft <- function(x,dlt,z,weight="logrank",ynci=0,epl=0.95,epu=0.05)
{ size <- length(x)
  npred <- dim(as.matrix(z))[2]
  lrgh <- 1
  if(weight == "logrank") lrgh <- 0

  fit <- .Fortran("faft", as.double(x), as.integer(dlt), as.double(z),
         double(size), as.integer(size), as.integer(npred),

         as.double(epl), as.double(epu),

         double(npred+1), integer(npred+1), integer(size), double(size),
         double(npred+1), double(npred+1), double((npred+1)^2),
         double((npred+1)^2), double(npred+1), double(size),
         double((npred+1)^2),

         as.integer(lrgh), iew=integer(1), bt=double(npred*3),
         va=double(npred^2), chis=double(3), ci=double(npred*2),
         as.integer(ynci),

         double(npred), double(npred^2), double(npred*size*2),
         double(npred^2*size*2), double(npred), double(npred),
         double(npred^2), double(npred^2))

  if(ynci == 1 & fit$iew == 0) obj <-
         list(weight=switch(lrgh+1,"logrank","Gehan"),
              beta=matrix(fit$bt,nrow=npred)[,1],
              va=matrix(fit$va,nrow=npred),qif=fit$chis[1],
              ci95=t(matrix(fit$ci,nrow=npred,
                        dimnames=list(NULL,c("lower","upper")))),
              message="success",
              imsg=fit$iew,
              beta1stp=matrix(fit$bt,nrow=npred)[,2],qif1stp=fit$chis[2],
              betainit=matrix(fit$bt,nrow=npred)[,3],qifinit=fit$chis[3])
  else obj <-
         list(weight=switch(lrgh+1,"logrank","Gehan"),
              beta=matrix(fit$bt,nrow=npred)[,1],
              va=matrix(fit$va,nrow=npred),qif=fit$chis[1],
              message=switch(fit$iew+1,"success","error - algorithm fails",
                                       "warning - singular hessian",
                                       "success"),
              imsg=fit$iew,
              beta1stp=matrix(fit$bt,nrow=npred)[,2],qif1stp=fit$chis[2],
              betainit=matrix(fit$bt,nrow=npred)[,3],qifinit=fit$chis[3])
  obj
}
