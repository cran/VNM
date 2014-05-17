\name{MOPT}
\alias{MOPT}
\title{Identify multiple-objective optimal designs for estimating 3 objectives (model parameters, the ED50, the MED).}
\description{
   Identify multiple-objective optimal design (i.e., optimal dose levels and corresponding optimal weights) that maximizes the efficiencies for estimating model parameters, the ED50, and the MED under the 4-parameter logistic model. The 4-parameter logistic model is used to describe the dose-response curve: Y=t1/(1+exp(t2*X+t3))+t4+e, where e~N(0,sigma^2); X=log(dose); t1=the minimum effect-the maximum effect; t2=the slope; t3=log(ED50)^(-t2); and t4=the maximum effect. 
}  
\usage{
MOPT(LB, UB, nit, T, q1, q2, dt)
}
\arguments{
  \item{LB}{ Numeric. Predetermined lower bound of the dose range for the log dose.}
  \item{UB}{ Numeric.  Predetermined upper bound of the dose range for the log dose.}
  \item{nit}{ Numeric. The number of iterations to set an initial design to search multiple-objective optimal designs for estimating the 3 objectives. Default is 1 and needed to be increased (for example, nit=10 or 30 or 50) if the searched multiple-objective optimal design is not a true optimal. It can be verified by the verification plot. }
  \item{T}{ A numeric vector. Predetermined parameter values for the 4-parameter logistic model. T=c(t1,t2,t3,t4).}
  \item{q1}{ Numeric. A weight for estimating model parameters. The weight represents relative importance of estimating the model parameters compared to the other two objectives.}
  \item{q2}{ Numeric. A weight for estimating the ED50. The weight represents relative importance of estimating the ED50 compared to the other two objectives.}
  \item{dt}{ Numeric. Predetermined clinically significant effect to define the MED. The MED is the dose producing the mean response of dt units better than the minimum dose.}
}
\value{ 1. A matrix showing the multiple-objective optimal design for estimating model parameters, the ED50, and the MED. The first row of the matrix represents optimal dose levels to be used and the second row of the matrix represents the optimal weights for the corresponding dose levels. The weight represents the proportional allocation of subjects to the corresponding dose level; 2. A verification plot of the multiple-objective optimal design by the General Equivalence Theorem.}
\author{Seung Won Hyun, Weng Kee Wong, and Yarong Yang }
\references{Seung Won Hyun, Weng Kee Wong, and Yarong Yang (2013), Multiple Objective Optimal Designs to Study the
Interesting Features in a Dose-Response Relationship, forthcoming. 
}
\examples{
   LB1=-3
   UB1=0
   nit1<-10
   T1<-c(-1.563,-0.825,-0.653,1.700)
   q11<-1/2
   q22<-1/3
   dt1<--1
   MOPT(LB=LB1, UB=UB1, nit=nit1, T=T1, q1=q11, q2=q22, dt=dt1)
}

