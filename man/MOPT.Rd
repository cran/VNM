\name{MOPT}
\alias{MOPT}
\title{Identify multiple-objective optimal designs for estimating model parameters, ED50, MED.}
\description{
   Identify multiple-objective optimal design (i.e., optimal dose levels and corresponding optimal weights) that maximizes the efficiencies for estimating model parameters, the ED50, and the MED under the 4-parameter logistic model. 
}  
\usage{
MOPT(LB, UB, P, lambda, delta, r, grid, epsilon, epsilon_w, verbose)
}
\arguments{
  \item{LB}{ Numeric. Predetermined lower bound of the dose range for the log dose.}
  \item{UB}{ Numeric.  Predetermined upper bound of the dose range for the log dose.}
  \item{P}{ A numeric vector. Solicited information on nominal values for the vector. P=(p1, p2, p3, p4), where p1 is the lower limit of the response, p2 is Emax, p3 is the ED50 and p4 is the slope at the ED50. For the 4-parameter logistic model, a user needs to specify all 4 nominal values in P: P=(p1, p2, p3, p4). For the 3-parameter logistic model, a user needs to specify only the 3 nominal values, Emax, the ED50, and the slope: P=(p2 ,p3, p4). For the 2-parameter logistic model, a user needs to specify only the 2 nominal values, the ED50 and the slope: P=(p3, p4).   }
  \item{lambda}{ A numeric vector. User select weights. lambda=c(q1, q2), where q1, q2 represent weights for estimating model parameter and estimating the ED50 respectively. They are non-negative and q1+q2<=1.}
  \item{delta}{ Numeric. Predetermined clinically significant effect to define the MED. The MED is the dose producing the mean response of delta units better than the minimum dose.}
  \item{r}{ Numeric. The number fo iteritions to set an initial design to search the multiple-objective optimal design. Default is 10 and needed to be increased (for example, r=30 or 50) if the searched multiple-objective optimal design is not a true optimal. It can be verified by the verification plot.}
  \item{grid}{ Numeric. The grid density to discretize the predetermined dose interval. Default is 0.01.}
  \item{epsilon}{ Numeric. Stopping criterion for the algorithm to search the multiple-objective optimal design. Default is 0.001.}
  \item{epsilon_w}{ Numeric. Stopping criterion for the Newton Raphson method inside of the algorithm. Default is 10^-6.}
  \item{verbose}{Set to TRUE if you want more evidence of progress while data is being processed. Set to FALSE if you want your CPU cycles to be used on analysis and not printing messages. Default is FALSE.}
}
\value{An object of class OPT.}
\author{Seung Won Hyun, Weng Kee Wong, and Yarong Yang }
\references{Seung Won Hyun, Weng Kee Wong, and Yarong Yang (2014), VNM: An R Package for Finding Multiple-Objective Optimal Designs for the 4-Parameter Logistic Model, submitted to Journal of Statistical Software.

Seung Won Hyun and Weng Kee Wong (2015), Multiple Objective Optimal Designs to Study the Interesting Features in a Dose-Response Relationship, accepted by the International Journal of Biostatistics.
}
\examples{
 MOPT(LB = -1, UB = 1, P = c(exp(-1/3),3),lambda = c(1/3,1/3), 
delta = 0.2, r = 10,grid=.01, verbose=TRUE)
}

