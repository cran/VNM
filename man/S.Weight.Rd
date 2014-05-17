\name{S.Weight}
\alias{S.Weight}
\title{Identify optimal weights for estimating 3 objectives (model parameters, the ED50, the MED) for given dose levels.}
\description{
  Identify optimal weights that maximize the efficiency for estimating 3 objectives (model parameters, the ED50, the MED) for given dose levels under the 4-parameter logistic model. The weights represent the proportional allocations of subjects to given dose levels. The 4-parameter logistic model is used to describe the dose-response curve: Y=t1/(1+exp(t2*X+t3))+t4+e, where e~N(0,sigma^2); X=log(dose); t1=the minimum effect-the maximum effect; t2=the slope; t3=log(ED50)^(-t2); and t4=the maximum effect. 
}
\usage{
S.Weight(X,T,q1,q2,dt)
}
\arguments{
  \item{X}{ A numeric vector. Given dose levels to search the optimal weights.}
  \item{T}{ A numeric vector. Predetermined parameter values for the 4-parameter logistic model. T=c(t1,t2,t3,t4).}
  \item{q1}{ Numeric. A weight for estimating model parameters. The weight represents relative importance of estimating the model parameters compared to the other two objectives.}
  \item{q2}{ Numeric. A weight for estimating the ED50. The weight represents relative importance of estimating the ED50 compared to the other two objectives.}
  \item{dt}{ Numeric. Predetermined clinically significant effect to define the MED. The MED is the dose producing the mean response of dt units better than the minimum dose.}
}
\value{ 1. A matrix showing the optimal weights for the specified X. The first row represents the given dose levels and the second row represents the obtained optimal weights. The matrix only shows the dose levels with non-zero optimal weights (i.e., dose levels with obtained optimal weights are equal to 0 are eliminated from the matrix); 2. The first derivative of the optimality criterion w.r.t weights; 3. The second derivative of the optimality criterion w.r.t weights.}
\author{Seung Won Hyun, Weng Kee Wong, and Yarong Yang }
\references{Seung Won Hyun, Weng Kee Wong, and Yarong Yang (2013), Multiple Objective Optimal Designs to Study the
Interesting Features in a Dose-Response Relationship, forthcoming. 
}
\examples{
  X1<-c(-6.91,-5.26,-4.10,6.91)
  T1<-c(-1.562,-1.03,4.851,1.700)
  q11<-1/3
  q22<-1/3
  dt1<--1
  S.Weight(X=X1,T=T1,q1=q22,q2=q22,dt=dt1)
}

