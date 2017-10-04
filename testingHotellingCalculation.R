dat=as.matrix(non.necrotic.groups[which(benoit.gene.names.abund==unique(benoit.gene.names.abund[2])),])[1:4,]

  y1=dat[,1:3]### split data set into different treatment groups
  y2=dat[,4:6]
  y3=dat[,7:9]
  y4=dat[,10:12]
  
  n1 = ncol(y1)### sample size for each group
  n2 = ncol(y2) 
  n3 = ncol(y3) 
  n4 = ncol(y4) 
  
  ybar1= apply(y1, 1, mean)## get a mean for each variable for each group
  ybar2= apply(y2, 1, mean)
  ybar3= apply(y3, 1, mean)
  ybar4= apply(y4, 1, mean)
  
  totbar=apply(dat,1,mean)
  
  SS1=var(t(y1))
  SS2=var(t(y2))
  SS3=var(t(y3))
  SS4=var(t(y4))
  
  # calculate the variance of the variables with each other for each group, multiply by sample size
  vw = ((n1-1)*SS1+(n2-1)*SS2+(n3-1)*SS3+(n4-1)*SS4)/(n1+n2+n3+n4-4)
  vb= cbind(ybar1-totbar,ybar2-totbar,ybar3-totbar,ybar4-totbar)
  SSCP=3*(vb%*%t(vb))/(4-1)
  Evar=vw%*%t(vw)
  Hvar=vb%*%t(vb)
  A=Hvar%*%solve(Evar)
  

fstat=c()
hotelling.trace=c()
hotelling.pval.test=c()
for(ii in 1:length(unique(benoit.gene.names.abund))){
xx=t(as.matrix(non.necrotic.groups[which(benoit.gene.names.abund==unique(benoit.gene.names.abund[ii])),]))
n=ncol(xx)## Number of observations
p=nrow(xx)
if(p<3){
  next
}else{
xx=manova(xx~benoit.groups[1:12])
object=xx
ss.resid =crossprod(object$residuals)
ss.effect=crossprod(object$effects[2:4,])
D <- diag(1/sqrt(sss <- diag(ss.resid)))
sss <- sss + diag(ss.effect)
rss.qr <- qr(D %*% ss.resid %*% D, tol = 1e-07)
if(p>8){
  A1=ss.effect%*%condreg(ss.resid,8)$invS
}else{
  A1=ss.effect%*%solve(ss.resid)
}
hotelling.trace[ii]=tr(A)
fstat[ii]=(n-p-1)*tr(A1)/((n-4)*p)
hotelling.pval.test[ii]=1-pf(fstat, p, n-p)
}
}

mrows=c("benoit.groups[1:12]","Residuals" )
for (i in 1) {
  eigs <- array(NA, c(1, 9), dimnames = list(nmrows[-2], NULL))
  stats <- matrix(NA, 2, 5, dimnames = list(nmrows, c(test, "approx F", "num Df", "den Df", "Pr(>F)")))
  A1 <- qr.coef(rss.qr, D %*% ss[[i]] %*% D)
  eigs[i, ] <- Re(eigen(A1, symmetric = FALSE,only.values = TRUE)$values)
  stats[i, 1L:4L] <- switch(test, Pillai = Pillai(eigs[i,], df[i], df[nt]), Wilks = Wilks(eigs[i, ], df[i], df[nt]), `Hotelling-Lawley` = HL(eigs[i,], df[i], df[nt]), Roy = Roy(eigs[i, ], df[i],df[nt]))
  ok <- stats[, 2L] >= 0 & stats[, 3L] > 0 & stats[,4L] > 0
  ok <- !is.na(ok) & ok
  stats[ok, 5L] <- pf(stats[ok, 2L], stats[ok,3L], stats[ok, 4L], lower.tail = FALSE)
}
x <- list(row.names = nmrows, SS = ss, Eigenvalues = eigs,stats = cbind(Df = df, stats = stats))
}



############ ANOVA with Holling's T^2 test #############################

library(DescTools)
library(mvtnorm)
muH0 <- c(-1, 2)
HotellingsT2Test(Y1, mu=muH0)

library(CondReg)
library(PDSCE)

hotelling.pval=matrix(nrow=length(unique(benoit.gene.names.abund)),ncol=4)
hotelling.Rsq=matrix(nrow=length(unique(benoit.gene.names.abund)),ncol=4)
for (ii in 1:length(unique(benoit.gene.names.abund))){
  xx=as.matrix(non.necrotic.groups[which(benoit.gene.names.abund==unique(benoit.gene.names.abund[2])),])
  if(ncol(E.var)<3){
    xx=manova(xx~benoit.groups[1:12])
    xx=summary(xx,test="Hotelling-Lawley")
    hotelling.pval[ii,]=xx$`Response Y1`$coefficients[,4]
    hotelling.Rsq[ii,]=c(xx$`Response Y1`$r.squared,xx$`Response Y1`$fstatistic)
    next
  }
  if(ncol(E.var)>8){
    tot.means=apply(xx,1,mean)
    group.means=matrix(nrow=nrow(xx),ncol=4)
    mean.dif=matrix(ncol=12,nrow=nrow(xx))
    for (iii in 1:(nrow(xx))){
      group.means[iii,]=tapply(xx[iii,],benoit.groups[1:12],mean) 
      mean.dif[iii,]=c((xx[iii,1:3]-group.means[iii,1]),(xx[iii,4:6]-group.means[iii,2]),(xx[iii,7:9]-group.means[iii,3]),(xx[iii,10:12]-group.means[iii,4]))
      H.var=group.means[iii,]-tot.means[iii]
    }
    H.var
    
    
  }
  else {
    xx=manova(t(xx)~benoit.groups[1:12])
  }
  xx=summary(xx,test="Hotelling-Lawley")
  hotelling.pval[ii,]=xx$`Response Y1`$coefficients[,4]
  hotelling.Rsq[ii,]=c(xx$`Response Y1`$r.squared,xx$`Response Y1`$fstatistic)
}

xx$`Response Y1`$r.squared
xx$`Response Y2`$coefficients[,4]

benoit.groups

#################### MANOVA "by hand" ####################
dat=as.matrix(non.necrotic.groups[which(benoit.gene.names.abund==unique(benoit.gene.names.abund[2])),])
Ew1=condreg(dat[,1:3],3)$S
Ew2=condreg(dat[,4:6],3)$S
Ew3=condreg(dat[,7:9],3)$S
Ew4=condreg(dat[,10:12],3)$S
Ewtot=(Ew1+Ew2+Ew3+Ew4)/4

Eb=cbind(apply(dat[,1:3],1,mean),apply(dat[,4:6],1,mean),apply(dat[,7:9],1,mean),apply(dat[,10:12],1,mean))
Eb=var(Eb)/3
Eb/Ew



for (ii in 1:length(unique(benoit.gene.names.abund))){
  xx=as.data.frame(t(non.necrotic.groups[which(benoit.gene.names.abund==unique(benoit.gene.names.abund[2])),]))
  xx=manova(as.matrix(xx)~benoit.groups[1:12])
  xx=summary.manova(xx,test="Hotelling-Lawley")
  hotelling.pval[ii,]=xx$`Response Y1`$coefficients[,4]
  hotelling.Rsq[ii,]=c(xx$`Response Y1`$r.squared,xx$`Response Y1`$fstatistic)
}  

# calculate the test statistic and associated quantities
q=qr(tr(A1))$rank
b=(p+2*n)*(q+2*n)/(2*(2*n+1)*(n-1))
c=(2+(p*q+2)/(b-1))/(2*n)
fstat=(tr(A1)/c)*((4+(p*q+2)/b-1)/(p*q))

##### Automatic Way
xx=as.matrix(non.necrotic.groups[which(benoit.gene.names.abund==unique(benoit.gene.names.abund[2])),])[1:4,]
  xx=manova(t(xx)~benoit.groups[1:12])
  xx=summary(xx,test="Hotelling-Lawley")
  hotelling.pval.test=xx$stats[1,6]
  hotelling.Fstat.test=xx$stats[1,3]

