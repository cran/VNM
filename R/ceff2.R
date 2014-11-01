ceff2 <-
function(weight,P,dose,LB,UB,delta,r=10,grid=0.01,epsilon=.001,epsilon_w=10^-6)
{ 
      t1=P[2]
      t2=-P[4]
      t3=P[4]*log(P[3])
      t4=P[1]
      T=c(t1,t2,t3,t4)
	e1=epsilon_w
	e2=epsilon
	dt=delta
	nit=r
	gr=grid    
     
      ginv <- function(X, tol = sqrt(.Machine$double.eps)) {
           dnx <- dimnames(X)
           if (is.null(dnx)) 
                dnx <- vector("list", 2)
           s <- svd(X)
           nz <- s$d > tol * s$d[1]
           structure(if (any(nz)) 
                s$v[, nz] %*% (t(s$u[, nz])/s$d[nz])
                else X, dimnames = dnx[2:1])
      }

      LB = round(LB, 2)
      UB = round(UB, 2)
      x = seq(LB, UB, gr)
      k = length(T)
      nit = nit
     
      infor <- function(T, X) {
           f = matrix(c(1/(1 + exp(T[2] * X + T[3])), (-T[1] * X * exp(T[2] * X + T[3]))/(1 + exp(T[2] * X + T[3]))^2, (-T[1] * exp(T[2] * X + T[3]))/(1 + exp(T[2] * X + T[3]))^2, 1), nrow = 4, ncol = 1, byrow = F)
           f %*% t(f)
      }

      upinfor <- function(W, T, X) {
           k = length(X)
           last_infor = infor(T, X[k])
           infor = (1 - sum(W)) * last_infor
           for (i in 1:(k - 1)) 
                infor = infor + W[i] * infor(T, X[i])
           infor
      }

      g1<-function(T,dt) 
           matrix(c(1/((T[1]-dt)*T[2]),(-log((T[1]-dt)/dt)+T[3])/T[2]^2,-1/T[2],0),nrow=4,ncol=1,byrow=F)

      g2<-function(T,dt) 
           matrix(c(-1/((T[1]+dt)*T[2]),(-log(-dt/(T[1]+dt))+T[3])/T[2]^2,-1/T[2],0),nrow=4,ncol=1,byrow=F)

      g<-function(T,dt) 
           if(T[2]<0) g1(T,dt) else g2(T,dt)

      d3<-function(T,X,XL,inv,dt) 
           -t(g(T,dt))%*%inv%*%(infor(T,X)-infor(T,XL))%*%inv%*%g(T,dt)%*%(t(g(T,dt))%*%inv%*%g(T,dt))^-1

      dd3=function(T,X1,X2,XL,inv,dt) {
           ((t(g(T,dt))%*%inv%*%(infor(T,X2)-infor(T,XL))%*%inv%*%(infor(T,X1)-infor(T,XL))%*%
           inv%*%g(T,dt)+t(g(T,dt))%*%inv%*%(infor(T,X1)-infor(T,XL))%*%inv%*%(infor(T,X2)-infor(T,XL))%*%
           inv%*%g(T,dt)%*%t(g(T,dt))%*%inv%*%g(T,dt))-(t(g(T,dt))%*%inv%*%(infor(T,X1)-infor(T,XL))%*%
           inv%*%g(T,dt)%*%t(g(T,dt))%*%inv%*%(infor(T,X2)-infor(T,XL))%*%inv%*%g(T,dt)))%*%
	     (t(g(T,dt))%*%inv%*%g(T,dt))^-2
      }

      c_weight<-function(W,T,X,d) {
           p=length(W)
           k=length(X)
           inv=ginv(upinfor(W,T,X))
           V=g(T,dt)%*%t(g(T,dt))
           M=upinfor(W,T,X)
           f1=rep(0,p)
           f2=matrix(c(rep(f1,p)),nrow=p,ncol=p,byrow=F)
           for (i in 1:p) 
               f1[i]=d3(T,X[i],X[k],inv,dt)
           for(i in 1:p) 
               for(j in 1:p) 
                   f2[i,j]=dd3(T,X[i],X[j],X[k],inv,dt)
           newweight=W-d*(f1%*%ginv(f2))
           newweight
      }
     
      search_weight<-function(X,T) {
           diff=10
           k=length(X)
           W=rep(1/k,k-1)
           while(diff>e1) {
                d=.2
                NW=c_weight(W,T,X,d)
                minW=min(min(NW),1-sum(NW))
                while(minW<0 & d>.0001) {
                      d=d/2
                      NW=c_weight(W,T,X,d)
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
           D
      }

      X=c(LB,LB+(UB-LB)/3,LB+2*(UB-LB)/3,UB)
      W=rep(1/4,3)
      n=1
      p=1
   
      f<-function(T,X) 
           matrix(c(1/(1+exp(T[2]*X+T[3])),(-T[1]*X*exp(T[2]*X+T[3]))/(1+exp(T[2]*X+T[3]))^2,(-T[1]*exp(T[2]*X+T[3]))/(1+exp(T[2]*X+T[3]))^2,1),nrow=4,ncol=1,byrow=F)

      M=upinfor(W,T,X)
      while(n<nit){
           x=seq(LB,UB,gr)
           n1=length(x)
           ds=rep(0,n1)
           inv=ginv(M)
           for (i in 1:n1) 
                ds[i]=t(f(T,x[i]))%*%inv%*%g(T,dt)%*%t(g(T,dt))%*%inv%*%f(T,x[i])
           for (i in 1:n1) 
                if(max(ds)==ds[i])x[i]=x[i] else x[i]=NA
           newX=na.omit(x)
           newX=round(newX[1],2)
           newds=max(ds)
           an=1/(n+1)
           BB=t(g(T,dt))%*%inv%*%g(T,dt)
           p=abs(newds-BB)
           newM=(1-an)*M+an*f(T,newX)%*%t(f(T,newX))
           M=newM
           X=c(X,newX)
           n=n+1
      }
      r=length(X)
      X=unique(X[(r-4):r])
      X=sort(X,decreasing=F)
      p=1
      n=1
      it=1

      while(p>e2){
           x=seq(LB,UB,gr)
           n1=length(x)
           ds=rep(0,n1)
           D=search_weight(X,T)
           X=D[1,]
           k=length(X)
           W=D[2,1:k-1]
           inv=ginv(upinfor(W,T,X))
           for (i in 1:n1) 
                ds[i]=t(f(T,x[i]))%*%inv%*%g(T,dt)%*%t(g(T,dt))%*%inv%*%f(T,x[i])
           for (i in 1:n1) 
                if(max(ds)==ds[i])x[i]=x[i] else x[i]=NA
           newX=na.omit(x)
           newds=max(ds)
           BB=t(g(T,dt))%*%inv%*%g(T,dt)
           X=c(X,newX[1])
           X=sort(X,decreasing=F)
           X=unique(X)
           newp=abs(newds-BB)
           if(abs(newp-p)<.0000000001) newp=10^-20
           if(it>20) newp=10^-20
           p=newp
           it=it+1
      }

      X=D[1,]
      n=length(X)
      W=D[2,1:n-1]
      M=upinfor(W,T,X)
      inv=ginv(M)
      x1=seq(LB,UB,gr)
      n1=length(x1)
      ds=rep(0,n1)
      BB=t(g(T,dt))%*%inv%*%g(T,dt)
      for (i in 1:n1) 
           ds[i]=t(f(T,x1[i]))%*%inv%*%g(T,dt)%*%t(g(T,dt))%*%inv%*%f(T,x1[i])-BB

      plot(x1,ds,cex=.3,main="Verify the c-optimal design for MED",ylab="Sensitivity function",xlab="Dose levels")

      kk=length(dose)
      weight=weight[1:kk-1]
      inv=ginv(upinfor(weight,T,dose))
      invs=ginv(upinfor(W,T,X))
      eff=(t(g(T,dt))%*%invs%*%g(T,dt))/(t(g(T,dt))%*%inv%*%g(T,dt))
      eff=as.vector(eff)
      cat(format("c-optimal design for MED", width=80),"\n")
      print(D)
      cat(format("c-efficiency for MED", width=80),"\n")
      print(eff)     
}
