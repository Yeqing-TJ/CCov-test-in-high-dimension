Chen2010JASA<-function(X,Y,alpha)
{
n<-dim(X)[1]
p<-dim(X)[2]
n<-length(Y)
XStar=X
Z=Y
Xnew=X-matrix(1,n,1)%*%t(colMeans(X))
S=t(Xnew)%*%Xnew/(n-1)
Q=sum(diag(diag(diag(Xnew%*%t(Xnew)))^2))/(n-1)
trSig2<-((n-1)/(n*(n-2)*(n-3)))*((n-1)*(n-2)*sum(diag(S%*%S))+sum(diag(S))^2-n*Q)
Y2n=(sum((Y%*%t(Y))*(X%*%t(X)))-sum(diag((Y%*%t(Y))*(X%*%t(X)))))/(n*(n-1))
XXn2=(X%*%t(X))%*%(Y%*%t(Y))
XXn<-X%*%t(X)
YYn=Y%*%t(Y)
diagXX<-rep(diag(XXn),each=n-1)
diagYY<-rep(diag(YYn),each=n-1)
offdiagXX<-as.logical(lower.tri(XXn)+upper.tri(XXn))
offdiagYY<-as.logical(lower.tri(YYn)+upper.tri(YYn))
VecOffdiagXX<-matrix(Y%*%t(Y),n^2,1)[offdiagXX]
VecOffdiagYY<-matrix(X%*%t(X),n^2,1)[offdiagYY]
Y4n<-(sum(XXn2)-sum(diag(XXn2))-sum(diagXX*VecOffdiagXX)-sum(diagYY*VecOffdiagYY))/(n*(n-1)*(n-2))
Y3n=sum(Y)*apply(t(X),1,sum)-t(X)%*%Y
Y5n<-(t(Y3n)%*%Y3n-2*n*(n-1)*(n-2)*Y4n-sum(diag(X%*%t(X)))*(t(Y)%*%Y)+
sum(diag((X%*%t(X))*(Y%*%t(Y))))-n*(n-1)*Y2n-sum(diag(X%*%t(X)))*(sum(Y)^2-sum(Y^2))-
sum(Y^2)*(sum(X%*%t(X))-sum(diag(X%*%t(X))))+2*sum(diagXX*VecOffdiagXX)+2*sum(diagYY*VecOffdiagYY))/(n*(n-1)*(n-2)*(n-3))
Tnp1<-Y2n-2*Y4n+Y5n
#sigma2<-mean(Y*Y)-(mean(Y))^2
sigma2<-var(Y)
#TC=as.vector((1/sigma2)*(n/sqrt(2*trSig2))*Tnp1)                     #####S.X. Chen
TC=as.vector((1/sigma2)*(sqrt(n*(n-1))/sqrt(2*trSig2))*Tnp1)                     #####S.X. Chen
rejTC<-0
if (TC>qnorm(1-alpha,0,1))
 {rejTC<-1}
return(list(NewStatTC=TC,NewTC=rejTC))
}