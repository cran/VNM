\name{Deff}
\alias{Deff}
\title{D-efficiency for 4-parameter logistic model. }
\description{
  Obtaining D-efficiency for estimating model parameters for a given design under the 4-parameter logistic model. 
}
\usage{
Deff(weight,P,dose,LB,UB,r,grid,epsilon,epsilon_w)
}
\arguments{
  \item{weight}{ A numeric vector. K weights for a given design. The weights represent the proportional allocations of subjects to the dose levels in a given design.}
  \item{P}{ A numeric vector. Predetermined parameter values for the 4-parameter logistic model. T=c(t1,t2,t3,t4).}
  \item{dose}{ A vector. K dose levels for a given design.}
  \item{LB}{ Numeric. Predetermined lower bound of the dose range for the log dose.}
  \item{UB}{ Numeric. Predetermined upper bound of the dose range for the log dose.}
  \item{r}{Numeric. The number of iterations to set an initial design to search D-optimal design for estimating the model paramters. Default is 10 and needed to be increased (for example, r=30 or 50) if the searched D-optimal design is not a true optimal. It can be verified by the verification plot.}
  \item{grid}{Numeric. The grid density to discretize the predetermined dose interval. Default is 0.01.}
  \item{epsilon}{Numeric. Stopping criterion for the algorithm to search D-optimal design. Default is 0.001.}
  \item{epsilon_w}{Numeric. Stopping criterion for the Newton Raphson method inside of the algorithm. Default is 10^-6.}
}
\value{ 1. D-optimal design for estimating model parameters; 2. A verification plot of the D-optimal design by the General Equivalence Theorem; 3. D-efficiency for estimating model parameters for a given design.}
\author{Seung Won Hyun, Weng Kee Wong, and Yarong Yang }
\references{Seung Won Hyun, Weng Kee Wong, and Yarong Yang (2014), VNM: An R Package for Finding Multiple-Objective Optimal Designs for the 4-Parameter Logistic Model, submitted to Journal of Statistical Software.

Seung Won Hyun and Weng Kee Wong (2014), Multiple Objective Optimal Designs to Study the Interesting Features in a Dose-Response Relationship, submitted to the International Journal of Biostatistics.
}
\examples{
  Deff(weight=c(.330,.302,.270,.098),P=c(60,340,107.14,1),
dose=c(-6.91, 2.22, 3.75, 4.60),LB=log(.001),UB=log(500),r=30)
}

