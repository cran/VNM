\name{Deff}
\alias{Deff}
\title{D-efficiency for estimating model parameters. }
\description{
  Obtaining D-efficiency for estimating model parameters for a given design under the 4-parameter logistic model. The 4-parameter logistic model is used to describe the dose-response curve: Y=t1/(1+exp(t2*X+t3))+t4+e, where e~N(0,sigma^2); X=log(dose); t1=the minimum effect-the maximum effect; t2=the slope; t3=log(ED50)^(-t2); and t4=the maximum effect. 
}
\usage{
Deff(weight,T,dose,nit,LB,UB)
}
\arguments{
  \item{weight}{ A numeric vector. Weights for a given design. The weights represent the proportional allocations of subjects to the dose levels in a given design.}
  \item{T}{ A numeric vector. Predetermined parameter values for the 4-parameter logistic model. T=c(t1,t2,t3,t4).}
  \item{dose}{ A vector. K dose levels for a given design.}
  \item{nit}{Numeric. The number of iterations to set an initial design to search D-optimal design for estimating model parameters. Default is 1 and needed to be increased (for example, nit=10 or 30 or 50) if the searched D-optimal design is not a true optimal. It can be verified by the verification plot. }
  \item{LB}{ Numeric. Predetermined lower bound of the dose range for the log dose.}
  \item{UB}{ Numeric. Predetermined upper bound of the dose range for the log dose.}
}
\value{ 1. D-optimal design for estimating model parameters; 2. A verification plot of the D-optimal design by the General Equivalence Theorem; 3. D-efficiency for estimating model parameters for a given design.}
\author{Seung Won Hyun, Weng Kee Wong, and Yarong Yang }
\references{Seung Won Hyun, Weng Kee Wong, and Yarong Yang (2013), Multiple Objective Optimal Designs to Study the
Interesting Features in a Dose-Response Relationship, forthcoming.
}
\examples{
  T1<-c(-1.563,-1.790,-8.442,1.700)
  dose1<-c(-4,-2.23,-1.14,1.00)
  weight1<-c(0.309,0.254,0.290)
  weight1<-c(weight1,1-sum(weight1))
  Deff(weight=weight1,T=T1,dose=dose1,nit=10,LB=-4,UB=1) 
}
