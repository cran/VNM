Deff<-function(weight,P,dose,LB,UB,r=10,grid=0.01,epsilon=.001,epsilon_w=10^-6)
{
	np=length(P)
	if (np==2)
	{t1=-P[2]
	 t2=P[2]*log(P[1])
	 T=c(1,t1,t2,0)} else if(np==3)
	 {t1=P[1]
	  t2=-P[3]
	  t3=P[3]*log(P[2])
	  T=c(t1,t2,t3,0)} else
      {t1=P[2]
	   t2=-P[4]
	   t3=P[4]*log(P[3])
	   t4=P[1]
	   T=c(t1,t2,t3,t4)}

	e1=epsilon_w
	e2=epsilon
	nit=r
	gr=grid

      ginv<-function(X, tol = sqrt(.Machine$double.eps)) {
           dnx <- dimnames(X)
           if(is.null(dnx)) dnx <- vector("list", 2)
           s <- svd(X)
           nz <- s$d > tol * s$d[1]
           structure(
                if(any(nz)) s$v[, nz] %*% (t(s$u[, nz])/s$d[nz]) else X,
                dimnames = dnx[2:1])
      }

      LB=round(LB,2)
      UB=round(UB,2)
      x=seq(LB,UB,gr)

      k=length(T)
      nit=nit

      infor2 <- function(T, X) {
          f = matrix(c(exp(1/2*(T[2] * X + T[3]))/(1 + exp(T[2] * X + T[3])),X*exp(1/2*(T[2] * X + T[3]))/
          (1 + exp(T[2] *X + T[3]))), nrow = 2, ncol = 1, byrow = F)
          f %*% t(f)
      }     
      infor3 <- function(T, X) {
          f = matrix(c(1/(1 + exp(T[2] * X + T[3])), (-T[1] * X * 
              exp(T[2] * X + T[3]))/(1 + exp(T[2] * X + T[3]))^2, 
              (-T[1] * exp(T[2] * X + T[3]))/(1 + exp(T[2] * X + 
                T[3]))^2), nrow = 3, ncol = 1, byrow = F)
          f %*% t(f)
      }
      infor4 <- function(T, X) {
          f = matrix(c(1/(1 + exp(T[2] * X + T[3])), (-T[1] * X * 
              exp(T[2] * X + T[3]))/(1 + exp(T[2] * X + T[3]))^2, 
              (-T[1] * exp(T[2] * X + T[3]))/(1 + exp(T[2] * X + 
                T[3]))^2, 1), nrow = 4, ncol = 1, byrow = F)
          f %*% t(f)
      }
 
      infor <- function(T,X){
      	if(np==2) {infor2(T,X)} else if (np==3)
      	{infor3(T,X)} else {infor4(T,X)}}
     


      upinfor<-function(W,T,X) {
           k=length(X)
           last_infor=infor(T,X[k])
           infor=(1-sum(W))*last_infor
           for (i in 1:(k-1))
               infor=infor+W[i]*infor(T,X[i])
           infor
      }

      d1<-function(T,X,XL,inv)
           sum(diag(inv%*%(infor(T,X)-infor(T,XL))))

      dd1<-function(T,X1,X2,XL,inv)
           sum(diag(-inv%*%(infor(T,X2)-infor(T,XL))%*%inv%*%(infor(T,X1)-infor(T,XL))))

      D_weight<-function(W,T,X,d) {
           p=length(W)
           k=length(X)
           inv=ginv(upinfor(W,T,X))
           M=upinfor(W,T,X)
           f1=rep(0,p)
           f2=matrix(c(rep(f1,p)),nrow=p,ncol=p,byrow=F)
           for (i in 1:p) 
                f1[i]=d1(T,X[i],X[k],inv)
           for(i in 1:p) 
                for(j in 1:p) 
                     f2[i,j]=dd1(T,X[i],X[j],X[k],inv)
           newweight=W-d*(f1%*%ginv(f2))
           newweight
      }

      search_weight<-function(X,T) {
           diff=10
           k=length(X)
           W=rep(1/k,k-1)
           while(diff>e1) {
                d=1
                NW=D_weight(W,T,X,d)
                minW=min(min(NW),1-sum(NW))
                while(minW<0 & d>.0001) {
                      d=d/2
                      NW=D_weight(W,T,X,d)
                      minW=min(min(NW),1-sum(NW))
                }
                NW=c(NW,1-sum(NW))
                n=length(NW)
                minW=min(NW)
                diff=max(abs(W-NW[1:n-1]))
                if (abs(minW)<.0000001||minW<0) {
                      for(i in 1:n) {
                             if (NW[i]==minW)NW[i]=0
                      }
           
                }
                D=rbind(X,NW)
                for (i in 1:n) 
                      if (D[2,i]==0) D[,i]=NA
                X=D[1,]
                W=D[2,]
                X=na.omit(X)
                W=na.omit(W)
                k=length(X)
                W=W[1:k-1]
           }
           W=c(W,1-sum(W))
           D=rbind(X,W)
           D
      }

      X1=runif(2,LB,UB)
      X2=runif(3,LB,UB)
      X3=c(LB,LB+(UB-LB)/3,LB+2*(UB-LB)/3,UB)
      if (np==2) {X=X1} else if (np==3) {X=X2} else {X=X3}
      wX=length(X)
      W=rep(1/wX,wX-1)
      n=1
      p=1
      
      f2 <- function(T, X) {
          matrix(c(exp(1/2*(T[2] * X + T[3]))/(1 + exp(T[2] * X + T[3])),X*exp(1/2*(T[2] * X + T[3]))/
        (1 + exp(T[2] *X + T[3]))), nrow = 2, ncol = 1, byrow = F)
          }     
      f3 <- function(T, X) {
          matrix(c(1/(1 + exp(T[2] * X + T[3])), (-T[1] * X * 
              exp(T[2] * X + T[3]))/(1 + exp(T[2] * X + T[3]))^2, 
              (-T[1] * exp(T[2] * X + T[3]))/(1 + exp(T[2] * X + 
                T[3]))^2), nrow = 3, ncol = 1, byrow = F)
          }      
      f4=function(T,X) {
            matrix(c(1/(1+exp(T[2]*X+T[3])),(-T[1]*X*exp(T[2]*X+T[3]))/(1+exp(T[2]*X+T[3]))^2,(-T[1]*exp(T[2]*X+T[3]))/(1+exp(T[2]*X+T[3]))			^2,1),nrow=4,ncol=1,byrow=F)
      	  }
 
      f <- function(T,X){
      	if(np==2) {f2(T,X)} else if (np==3)
      	{f3(T,X)} else {f4(T,X)}}

      M=upinfor(W,T,X)
      while(n<nit){
           x=seq(LB,UB,gr)
           n1=length(x)
           ds=rep(0,n1)
           inv=ginv(M)
           for (i in 1:n1) 
                 ds[i]=t(f(T,x[i]))%*%inv%*%f(T,x[i])
           for (i in 1:n1) 
                 if(max(ds)==ds[i])x[i]=x[i] else x[i]=NA
           newX=na.omit(x)
           newX=round(newX[1],2)
           newds=max(ds)
           an=1/(n+1)
           p=abs(newds-np)
           newM=(1-an)*M+an*f(T,newX)%*%t(f(T,newX))
           M=newM
           X=c(X,newX)
           n=n+1
      }
      r=length(X)
      X=unique(X[(r-np):r])
      X=sort(X,decreasing=F)
      if (length(X)==1) X=c(X,runif(1,LB,UB))
      n=1
      p=1
      it=1
      while(p>e2) {
           x=seq(LB,UB,gr)
           n1=length(x)
           ds=rep(0,n1)
           D=search_weight(X,T)
           X=D[1,]
           k=length(X)
           W=D[2,1:k-1]
           inv=ginv(upinfor(W,T,X))
           for (i in 1:n1) 
                ds[i]=t(f(T,x[i]))%*%inv%*%f(T,x[i])
           for (i in 1:n1) 
                if(max(ds)==ds[i])x[i]=x[i] else x[i]=NA
           newX=na.omit(x)
           newX=round(newX[1],2)
           newds=max(ds)
           X=c(X,newX)
           X=sort(X,decreasing=F)
           X=unique(X)
           newp=abs(newds-np)
           if(abs(newp-p)<.0000000001) newp=10^-20
           if(it>20)newp=10^-20
           p=newp
           it=it+1
      }


      X=D[1,]
      n=length(X)
      W=D[2,1:n-1]
      M=upinfor(W,T,X)
      x1=seq(LB,UB,gr)
      n1=length(x1)
      ds=rep(0,n1)
      for (i in 1:n1) 
           ds[i]=t(f(T,x1[i]))%*%ginv(M)%*%f(T,x1[i])-k

      # plot(x1,ds,cex=.3,main="Verify the D-optimal design",ylab="Sensitivity function",xlab="Dose levels")

      kk=length(dose)
      weight=weight[1:kk-1]
      eff=(det(upinfor(weight,T,dose))/det(upinfor(W,T,X)))^(1/np)
    
      R=new("PAR",fid="Deff",LB=LB,UB=UB,grid=grid,ds=ds)
      return(new("OPT",Par=R,Opt=D,Eff=eff))
}



