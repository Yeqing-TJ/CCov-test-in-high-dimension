rm(list=ls(all=TRUE))
#library(xtable)
source("outputHCCov.R")
source("outputHMDD.R")
source("Chen2010JASA1.R")

DD<-c(0.0,0.25,0.5,0.75,1.0)


Mpower=matrix(0,5,3)

for (k in 1:5){

delta=DD[k]
n<-120
p<-floor(exp(n^{0.4}))+230

rhop<-8
rep_time<-1000



eta<-0.04
alpha=0.05


q=floor(0.5*p^{0.7})
non0<-sqrt(eta/q)
beta<-c(rep(non0,q),rep(0,(p-q)))



set.seed(1000)
rho1<-runif(rhop,0,1)
#rho<-c(rho1,rep(0,p-rhop))
#Gamma<-cbind(diag(rep(rho[1],p)),matrix(0,p,m-p))
#for (i in 2:(p-1))
#{Gamma<-Gamma+cbind(matrix(0,p,i-1),diag(rep(rho[i],p)),matrix(0,p,m-p-i+1))}
#Gamma<-Gamma+cbind(matrix(0,p,m-p),diag(rep(rho[p],p)))

CCovstat=seq(0,0,length=rep_time) 
siindCCov=seq(0,0,length=rep_time) 
MDDstat2018=seq(0,0,length=rep_time) 
siindMDD2018=seq(0,0,length=rep_time) 
Chenstat2010=seq(0,0,length=rep_time) 
siindChen2010=seq(0,0,length=rep_time) 

for (j in 1:rep_time){
  print(j)
  z=matrix(rnorm((p+rhop-1)*n,0,1),n,(p+rhop-1))
  X=matrix(0,n,p)
  for (i in 1:p){
    X[,i]=z[,i:(i+rhop-1)]%*%rho1*(i^{delta/2})
  }
  Y<-X%*%beta+rnorm(n,0,1)
  
  
  Chenstat2010[j]=Chen2010JASA(X,Y,alpha)$NewStatTC
  siindChen2010[j]=Chen2010JASA(X,Y,alpha)$NewTC
  MDDstat2018[j]=outputHMDD(X,Y)
  siindMDD2018[j]=(MDDstat2018[j]>qnorm(1-alpha,0,1))
  
  CCovstat[j]=outputHCCov(X,Y)
  siindCCov[j]=(CCovstat[j]>qnorm(1-alpha,0,1))
# print(c(siindCCov[j],siindChen2010[j],siindMDD2018[j]))
  print(c(mean(siindCCov[1:j]),mean(siindMDD2018[1:j]),mean(siindChen2010[1:j])))
}
Mpower[k,]=c(mean(siindCCov),mean(siindMDD2018),mean(siindChen2010))
}
save(list=ls(all=TRUE),file="NS_Linear_H1_n120_normal.Rdata")

