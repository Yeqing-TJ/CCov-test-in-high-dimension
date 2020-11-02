outputHCCov<- function(x,y){
##estimate Snp^2 using  the plug-in method##
	n=dim(x)[1]
      p=dim(x)[2]
      eachcomp= matrix(0,n,n)
    	A<- matrix(0,n,n);
    	B<- matrix(0,n,n);
                c_n=(1-2/n+2/(n*n))^4

     
      B=as.matrix(dist(y, method = "euclidean", diag = TRUE, upper = TRUE))^2/2;
      rowsB=apply(B,1,sum)
	Btilde=B-(rowsB/(n-0))%*%matrix(1,1,n)-matrix(1,n,1)%*%t(rowsB/(n-0))+sum(rowsB)/((n-0)*(n-0))


      Results=rep(0,p)

xx=x
for(j in 1:p){
xx[,j]=apply(x[,j]%*%matrix(1,1,n)<=t(x[,j]%*%matrix(1,1,n)),2,mean)
}

      for(j in 1:p){
	Results[j] <-CCov_C(x[,j],y);   ## Compute the test-statistics for each component
## On computing the variance for each component
         A=as.matrix(dist(xx[,j], method = "euclidean", diag = TRUE, upper = TRUE));
        rowsA=apply(A,1,sum)
	Atilde=A-(1/2-xx[,j]+xx[,j]^2)%*%matrix(1,1,n)-t((1/2-xx[,j]+xx[,j]^2)%*%matrix(1,1,n))+1/3
    eachcomp<-eachcomp+Atilde; 
       }
AB=(eachcomp^2)*(Btilde^2)
Shat2=(sum(AB)-sum(diag(AB)))/(n*(n-1)*c_n)


Tn=2*sum(Results)/sqrt(Shat2)
	return(Tn);	
}


CCov_C<- function(x,y){
	
n <- length(y)

y=y-mean(y)
posit=order(x,decreasing=FALSE)
y=y[posit]

CovtestStat=((n*n-5*n+6)*sum((cumsum(y)-y)^2)-sum((n*n-3*n-2*n*(1:n-1)+4*(1:n-1))*(cumsum(y^2)-y^2))+
2*sum((n*(1:n-2)-2*(1:n)+2)*y*(cumsum(y)-y))+2*sum(y*y*(1:n-1)*(1:n-1))+
(n*(n-5)/2-(n-1)*n*(2*n-1)/6)*sum(y^2))/(n*(n-1)*(n-2)*(n-3)*(n-4))

CovtestStat=sqrt(n*(n-1)/2)*CovtestStat
     return(as.vector(CovtestStat));
}


