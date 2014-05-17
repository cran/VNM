\name{ceff2}
\alias{ceff2}
\title{Checking c-efficiency for estimating the MED.}
\description{
  Obtaining c-efficiency for estimating MED for a given design under the 4-parameter logistic model. Here MED is the minimum effective dose level. The 4-parameter logistic model is used to describe the dose-response curve: Y=t1/(1+exp(t2*X+t3))+t4+e, where e~N(0,sigma^2); X=log(dose); t1=the minimum effect-the maximum effect; t2=the slope; t3=log(ED50)^(-t2); and t4=the maximum effect. 
}
\usage{
ceff2(weight,T,dose,nit,LB,UB,dt)
}
\arguments{
  \item{weight}{A numeric vector. Weights for a given design. The weights represent the proportional allocations of subjects to the dose levels in a given design.}
  \item{T}{A numeric vector. Predetermined parameter values for the 4-parameter logistic model. T=c(t1,t2,t3,t4).}
  \item{dose}{ A vector. K dose levels for a given design.}
  \item{nit}{Numeric. The number of iterations to set an initial design to search c-optimal design for estimating the MED. Default is 1 and needed to be increased (for example, nit=10 or 30 or 50) if the searched c-optimal design is not a true optimal. It can be verified by the verification plot. }
  \item{LB}{ Numeric. Predetermined lower bound of the dose range for the log dose. }
  \item{UB}{ Numeric. Predetermined upper bound of the dose range for the log dose.}
  \item{dt}{ Numeric. Predetermined clinically significant effect to define the MED. The MED is the dose producing the mean response of dt units better than the minimum dose.}
}
\value{ 1. c-optimal design for estimating the MED; 2. A verification plot of the c-optimal design for estimating the MED by the General Equivalence Theorem; 3. c-efficiency for estimating the MED for a given design.}
\author{Seung Won Hyun, Weng Kee Wong, and Yarong Yang }
\references{Seung Won Hyun, Weng Kee Wong, and Yarong Yang (2013), Multiple Objective Optimal Designs to Study the
Interesting Features in a Dose-Response Relationship, forthcoming. 
}
\examples{
  T1<-c(-1.563,-0.825,-0.653,1.700)
  dose1<-c(-3,-2.23,-1.14,0)
  weight1<-c(0.309,0.254,0.290)
  weight1<-c(weight1,1-sum(weight1))
  ceff2(weight=weight1,T=T1,dose=dose1,nit=10,LB=-3,UB=0,dt=-1) 
}

