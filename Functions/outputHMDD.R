outputHMDD<- function(x,y){
	n=dim(x)[1]
      p=dim(x)[2]
      eachcomp= matrix(0,n,n)
    	A<- matrix(0,n,n);
    	B<- matrix(0,n,n);
    	c_n=(n-3)^4/((n-1)^4)+2*(n-3)^4/((n-1)^4*(n-2)^3)+2*(n-3)/((n-1)^4*(n-2)^3)
     
      B=as.matrix(dist(y, method = "euclidean", diag = TRUE, upper = TRUE))^2/2;
      rowsB=apply(B,1,sum)
	Btilde=B-(rowsB/(n-2))%*%matrix(1,1,n)-matrix(1,n,1)%*%t(rowsB/(n-2))+sum(rowsB)/((n-1)*(n-2))


      Results=rep(0,p)



      for(j in 1:p){
	Results[j] <-UMDC_C(x[,j],y);   ## Compute the test-statistics for each component
## On computing the variance for each component
         A=as.matrix(dist(x[,j], method = "euclidean", diag = TRUE, upper = TRUE));
        rowsA=apply(A,1,sum)
	Atilde=A-(rowsA/(n-2))%*%matrix(1,1,n)-matrix(1,n,1)%*%t(rowsA/(n-2))+sum(rowsA)/((n-1)*(n-2))
    eachcomp<-eachcomp+Atilde; 
       }
AB=(eachcomp^2)*(Btilde^2)
Shat2=(sum(AB)-sum(diag(AB)))/(n*(n-1)*c_n)

Tn=sum(Results)/sqrt(Shat2)
	return(Tn);	
}


UMDC_C<- function(X,Y){
	
	n <- length(Y); 
      v1=matrix(1,1,n)
    	K <- matrix(0,n,n);
    	L <- matrix(0,n,n);


      L=-as.matrix(dist(Y, method = "euclidean", diag = TRUE, upper = TRUE))^2/2;
	K =as.matrix(dist(X, method = "euclidean", diag = TRUE, upper = TRUE));


        UMDCtestStat=-(sum(diag(K%*%L))+(v1%*%K%*%t(v1))*(v1%*%L%*%t(v1))/((n-1)*(n-2))-2*(v1%*%(K%*%L)%*%t(v1))/(n-2))/(n*(n-3));
	  UMDCtestStat=sqrt(n*(n-1)/2)*UMDCtestStat
     return(as.vector(UMDCtestStat));
}


