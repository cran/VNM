ceff1 <-
function(weight,T,dose,nit,LB,UB)
{    
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
     x = seq(LB, UB, 0.01)
     k = length(T)
     nit = nit
     
     infor <- function(T, X) {
        f = matrix(c(1/(1 + exp(T[2] * X + T[3])), (-T[1] * X * 
            exp(T[2] * X + T[3]))/(1 + exp(T[2] * X + T[3]))^2, 
            (-T[1] * exp(T[2] * X + T[3]))/(1 + exp(T[2] * X + 
                T[3]))^2, 1), nrow = 4, ncol = 1, byrow = F)
        f %*% t(f)
     }

     upinfor <- function(W, T, X) {
        k = length(X)
        last_infor = infor(T, X[k])
        infor = (1 - sum(W)) * last_infor
        for (i in 1:(k - 1)) {
            infor = infor + W[i] * infor(T, X[i])
        }
        infor
     }

     g<-function(T) {
        matrix(c(0,-T[3]/(T[2])^2,-1/T[2],0),nrow=4,ncol=1,byrow=F)
     } 

     c_weight<-function(W,T,X,d) {
        p=length(W)
        k=length(X)
        inv=ginv(upinfor(W,T,X))
        V=g(T)%*%t(g(T))
        M=upinfor(W,T,X)
        f1=rep(0,p)
        f2=matrix(c(rep(f1,p)),nrow=p,ncol=p,byrow=F)
        for (i in 1:p) {
             f1[i]=sum(diag(-inv%*%(infor(T,X[i])-infor(T,X[k]))%*%inv%*%V))
        }
        for(i in 1:p) {
             for(j in 1:p) {
                  f2[i,j]=(sum(diag((inv%*%(infor(T,X[j])-infor(T,X[k]))%*%inv%*%(infor(T,X[i])-infor(T,X[k]))%*%inv+inv%*%(infor(T,X[i])-infor(T,X[k]))%*%inv%*%(infor(T,X[j])-infor(T,X[k]))%*%inv)%*%V)))
             }
        }
        newweight=W-d*(f1%*%ginv(f2))
        newweight
     }
     
     search_weight<-function(X,T) {
        diff=10
        k=length(X)
        W=rep(1/k,k-1)
        while(diff>.000000001) {
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
            if (minW<0) {
                 for(i in 1:n) {
                       if (NW[i]==minW)NW[i]=0
                 }
            }
            diff=max(abs(W-NW[1:n-1]))
            D=rbind(X,NW)
            for (i in 1:n) {
                 if (D[2,i]==0) D[,i]=NA
            }
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
   
     f<-function(T,X) {
         matrix(c(1/(1+exp(T[2]*X+T[3])),(-T[1]*X*exp(T[2]*X+T[3]))/(1+exp(T[2]*X+T[3]))^2,(-T[1]*exp(T[2]*X+T[3]))/(1+exp(T[2]*X+T[3]))^2,1),nrow=4,ncol=1,byrow=F)
     }
     M=upinfor(W,T,X)
     while(n<nit){
         x=seq(LB,UB,.01)
         n1=length(x)
         ds=rep(0,n1)
         inv=ginv(M)
         for (i in 1:n1) {
               ds[i]=t(f(T,x[i]))%*%inv%*%g(T)%*%t(g(T))%*%inv%*%f(T,x[i])
         }
         for (i in 1:n1) {
               if(max(ds)==ds[i])x[i]=x[i] else x[i]=NA
         }
         newX=na.omit(x)
         newX=round(newX,2)
         newds=max(ds)
         an=1/(n+1)
         BB=t(g(T))%*%inv%*%g(T)
         p=abs(newds-BB)
         newM=(1-an)*M+an*f(T,newX)%*%t(f(T,newX))
         M=newM
         X=c(X,newX)
         n=n+1
     }
     r=length(X)
     X=unique(X[(r-4):r])
     R=search_weight(X,T)
     X=R[1,]
     k=length(X)
     W=R[2,1:k-1]
     p=1
     n=1
     it=1

     while(p>.0001){
          x=seq(LB,UB,.01)
          n1=length(x)
          ds=rep(0,n1)
          inv=ginv(upinfor(W,T,X))
          for (i in 1:n1) {
               ds[i]=t(f(T,x[i]))%*%inv%*%g(T)%*%t(g(T))%*%inv%*%f(T,x[i])
          }
          for (i in 1:n1) {
               if(max(ds)==ds[i])x[i]=x[i] else x[i]=NA
          }
          newX=na.omit(x)
          newds=max(ds)
          X=c(X,newX)
          X=sort(X,decreasing=F)
          X=unique(X)
          D=search_weight(X,T)
          X=D[1,]
          k=length(X)
          W=D[2,1:k-1]
          BB=t(g(T))%*%inv%*%g(T)
          newp=abs(newds-BB)
          if(abs(newp-p)<.0000000001) newp=.000001
          if(it>20) newp=.000001
          p=newp
          it=it+1
     }

     X=D[1,]
     n=length(X)
     W=D[2,1:n-1]
     M=upinfor(W,T,X)
     inv=ginv(M)
     x1=seq(LB,UB,.01)
     n1=length(x1)
     ds=rep(0,n1)
     BB=t(g(T))%*%inv%*%g(T)
     for (i in 1:n1) {
           ds[i]=t(f(T,x1[i]))%*%inv%*%g(T)%*%t(g(T))%*%inv%*%f(T,x1[i])-BB
     }
     plot(x1,ds,cex=.1,main="Verify the c-optimal design for ED50",ylab="Sensitive function",xlab="Dose levels")

     kk=length(dose)
     weight=weight[1:kk-1]
     inv=ginv(upinfor(weight,T,dose))
     invs=ginv(upinfor(W,T,X))
     eff=(t(g(T))%*%invs%*%g(T))/(t(g(T))%*%inv%*%g(T))
     eff=as.vector(eff)
     cat(format("c-optimal design for ED50", width=80),"\n")
     print(D)
     cat(format("c-efficiency for ED50", width=80),"\n")
     print(eff)     
}
