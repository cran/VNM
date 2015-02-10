\name{ceff1}
\alias{ceff1}
\title{Checking c-efficiency for estimating the ED50.}
\description{
  Obtaining c-efficiency for estimating the ED50 for a given design under the 4-parameter logistic model. Here the ED50 is the dose producing a 50 percent response between the estimates of the minimum and maximum of the response. The 4-parameter logistic model is used to describe the dose-response curve: Y=t1/(1+exp(t2*X+t3))+t4+e, where e~N(0,sigma^2); X=log(dose); t1=the minimum effect-the maximum effect; t2=the slope; t3=log(ED50)^(-t2); and t4=the maximum effect.
}
\usage{
ceff1(weight,P,dose,LB,UB,r,grid,epsilon,epsilon_w)
}
\arguments{
  \item{weight}{A numeric vector. K weights for a given design. The weights represent the proportional allocations of subjects to the dose levels in a given design. }
  \item{P}{ A numeric vector. Solicited information on nominal values for the vector. P=(p1, p2, p3, p4), where p1 is the lower limit of the response, p2 is Emax, p3 is the ED50 and p4 is the slope at the ED50. For the 4-parameter logistic model, a user needs to specify all 4 nominal values in P: P=(p1, p2, p3, p4). For the 3-parameter logistic model, a user needs to specify only the 3 nominal values, Emax, the ED50, and the slope: P=(p2 ,p3, p4). For the 2-parameter logistic model, a user needs to specify only the 2 nominal values, the ED50 and the slope: P=(p3, p4).  }
  \item{dose}{ A vector. K dose levels for a given design.}
  \item{LB}{ Numeric. Predetermined lower bound of the dose range for the log dose.}
  \item{UB}{ Numeric. Predetermined upper bound of the dose range for the log dose.}
  \item{r}{ Numeric. The number fo iterations to set an initial design to search c-optimal design for estimating the ED50. Default is 10 and needed to be increased (for example, r=30 or 50) if the searched c-optimal design is not a true optimal. It can be verified by the verification plot.} 
  \item{grid}{ Numeric. The grid density to discretize the predetermined dose interval. Default is 0.01.}
  \item{epsilon}{ Numeric. Stopping criterion for the algorithm to search c-optimal design for the ED50. Default is 0.001.}
  \item{epsilon_w}{ Numeric. Stopping criterion for the Newton Raphson method inside of the algorithm. Default is 10^-6.}
}
\value{An object of class OPT.}
\author{Seung Won Hyun, Weng Kee Wong, and Yarong Yang }
\references{Seung Won Hyun, Weng Kee Wong, and Yarong Yang (2014), VNM: An R Package for Finding Multiple-Objective Optimal Designs for the 4-Parameter Logistic Model, submitted to Journal of Statistical Software.

Seung Won Hyun and Weng Kee Wong (2014), Multiple Objective Optimal Designs to Study the Interesting Features in a Dose-Response Relationship, submitted to the International Journal of Biostatistics.
}
\examples{
   ceff1(weight = c(0.330, 0.302, 0.270, 0.098), P = c(22, 16.8, 70, 1),
dose = c(-6.91, 2.22, 3.75, 4.60), LB = log(.001), UB = log(100),r = 30)

}
