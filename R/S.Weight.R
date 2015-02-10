S.Weight<-function(X,P,lambda,delta,epsilon_w=10^-6)
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
      e=epsilon_w
	dt=delta
	q1=lambda[1]
	q2=lambda[2]
      q3=1-q1-q2

      ginv<-function(X, tol = sqrt(.Machine$double.eps)){
           dnx <- dimnames(X)
           if(is.null(dnx)) dnx <- vector("list", 2)
           s <- svd(X)
           nz <- s$d > tol * s$d[1]
           structure(
                if(any(nz)) s$v[, nz] %*% (t(s$u[, nz])/s$d[nz]) else X,
                dimnames = dnx[2:1])
      }

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

      g21<-function(T,dt) {
          matrix(c(1/((T[1]-dt)*T[2]),(-log((T[1]-dt)/dt)+T[3])/T[2]^2,-1/T[2],0),nrow=4,ncol=1,byrow=F)
      }

      g22<-function(T,dt) {
          matrix(c(-1/((T[1]+dt)*T[2]),(-log(-dt/(T[1]+dt))+T[3])/T[2]^2,-1/T[2],0),nrow=4,ncol=1,byrow=F)
      }

      tg1<-function(T) {
          matrix(c(0,-T[3]/(T[2])^2,-1/T[2],0),nrow=4,ncol=1,byrow=F)
      } 
      g1<-function(T){tempT=tg1(T)
      	    if (np==2) {matrix(c(tempT[2],tempT[3]),nrow=2,ncol=1,byrow=F)}else if
      	    (np==3){matrix(c(0,tempT[2],tempT[3]),nrow=3,ncol=1,byrow=F)}else
      	    {tempT}}

      tg2<-function(T, dt) {
          if(T[2]<0) g21(T,dt) else g22(T,dt)
      }    
	g2<-function(T,dt){tempT=tg2(T,dt)
      	    if (np==2) {matrix(c(tempT[2],tempT[3]),nrow=2,ncol=1,byrow=F)}else if
      	    (np==3){matrix(c(tempT[1],tempT[2],tempT[3]),nrow=3,ncol=1,byrow=F)}else
      	    {tempT}}
     
      d1<-function(T,X,XL,inv)
           sum(diag(inv%*%(infor(T,X)-infor(T,XL))))

      d2<-function(T,X,XL,inv)
           -t(g1(T))%*%inv%*%(infor(T,X)-infor(T,XL))%*%inv%*%g1(T)%*%(t(g1(T))%*%inv%*%g1(T))^-1

      d3<-function(T,X,XL,inv,dt) 
           -t(g2(T,dt))%*%inv%*%(infor(T,X)-infor(T,XL))%*%inv%*%g2(T,dt)%*%(t(g2(T,dt))%*%inv%*%g2(T,dt))^-1


      dd1<-function(T,X1,X2,XL,inv)
           sum(diag(-inv%*%(infor(T,X2)-infor(T,XL))%*%inv%*%(infor(T,X1)-infor(T,XL))))

      dd2<-function(T,X1,X2,XL,inv) {
           ((t(g1(T))%*%inv%*%(infor(T,X2)-infor(T,XL))%*%inv%*%(infor(T,X1)-infor(T,XL))%*%
           inv%*%g1(T)+t(g1(T))%*%inv%*%(infor(T,X1)-infor(T,XL))%*%inv%*%(infor(T,X2)-infor(T,XL))%*%
           inv%*%g1(T)%*%t(g1(T))%*%inv%*%g1(T))-(t(g1(T))%*%inv%*%(infor(T,X1)-infor(T,XL))%*%inv%*%
           g1(T)%*%t(g1(T))%*%inv%*%(infor(T,X2)-infor(T,XL))%*%inv%*%g1(T)))%*%
           (t(g1(T))%*%inv%*%g1(T))^-2
      }
 
      dd3<-function(T,X1,X2,XL,inv,dt) {
           ((t(g2(T,dt))%*%inv%*%(infor(T,X2)-infor(T,XL))%*%inv%*%(infor(T,X1)-infor(T,XL))%*%
           inv%*%g2(T,dt)+t(g2(T,dt))%*%inv%*%(infor(T,X1)-infor(T,XL))%*%inv%*%(infor(T,X2)-infor(T,XL))%*%
           inv%*%g2(T,dt)%*%t(g2(T,dt))%*%inv%*%g2(T,dt))-(t(g2(T,dt))%*%inv%*%(infor(T,X1)-infor(T,XL))%*%
           inv%*%g2(T,dt)%*%t(g2(T,dt))%*%inv%*%(infor(T,X2)-infor(T,XL))%*%inv%*%g2(T,dt)))%*%
           (t(g2(T,dt))%*%inv%*%g2(T,dt))^-2
      }

      D_weight<-function(W,T,X,d) {
           p=length(W)
           k=length(X)
           inv=ginv(upinfor(W,T,X))
           M=upinfor(W,T,X)
           f1=rep(0,p)
           f2=matrix(c(rep(f1,p)),nrow=p,ncol=p,byrow=F)
           for (i in 1:p) 
                f1[i]=q1*d1(T,X[i],X[k],inv)/np-q2*d2(T,X[i],X[k],inv)-q3*d3(T,X[i],X[k],inv,dt)

           for(i in 1:p) 
                for(j in 1:p) 
                     f2[i,j]=q1*dd1(T,X[i],X[j],X[k],inv)/np-q2*dd2(T,X[i],X[j],X[k],inv)-q3*dd3(T,X[i],X[j],X[k],inv,dt)

           newweight=W-d*(f1%*%ginv(f2))
           newweight
      }
    
      diff=10
      k=length(X)
      W=rep(1/k,k-1)

      while(diff>e) {
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
                 for(i in 1:n) 
                       if (NW[i]==minW) NW[i]=0          
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
      #cat("Optimal weights", "\n")
      #print(D)
      k=ncol(D)
      W=D[2,1:k-1]
      X=D[1,]
      inv=ginv(upinfor(W,T,X))
      p=length(W)
      f1=rep(0,p)
      f2=matrix(c(rep(f1,p)),nrow=p,ncol=p,byrow=F)

      for (i in 1:p) 
           f1[i]=q1*d1(T,X[i],X[k],inv)/4-q2*d2(T,X[i],X[k],inv)-q3*d3(T,X[i],X[k],inv,dt)
      
      for(i in 1:p) 
           for(j in 1:p)
                 f2[i,j]=q1*dd1(T,X[i],X[j],X[k],inv)/4-q2*dd2(T,X[i],X[j],X[k],inv)-q3*dd3(T,X[i],X[j],X[k],inv,dt)
     
      return(new("SW",Opt.W=D,First.C=f1,Second.C=diag(f2)))
}
      

