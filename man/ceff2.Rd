\name{ceff2}
\alias{ceff2}
\title{Checking c-efficiency for estimating the MED.}
\description{
  Obtaining c-efficiency for estimating MED for a given design under the 4-parameter logistic model. Here MED is the minimum effective dose level.  
}
\usage{
ceff2(weight,P,dose,LB,UB,delta,r,grid,epsilon,epsilon_w)
}
\arguments{
  \item{weight}{A numeric vector. K weights for a given design. The weights represent the proportional allocations of subjects to the dose levels in a given design.}
  \item{P}{ A numeric vector. Solicited information on nominal values for the vector. P=(p1, p2, p3, p4), where p1 is the lower limit of the response, p2 is Emax, p3 is the ED50 and p4 is the slope at the ED50. For the 4-parameter logistic model, a user needs to specify all 4 nominal values in P: P=(p1, p2, p3, p4). For the 3-parameter logistic model, a user needs to specify only the 3 nominal values, Emax, the ED50, and the slope: P=(p2 ,p3, p4). For the 2-parameter logistic model, a user needs to specify only the 2 nominal values, the ED50 and the slope: P=(p3, p4).   }
  \item{dose}{ A vector. K dose levels for a given design.}
  \item{LB}{ Numeric. Predetermined lower bound of the dose range for the log dose. }
  \item{UB}{ Numeric. Predetermined upper bound of the dose range for the log dose.}
  \item{delta}{ Numeric. Predetermined clinically significant effect to define the MED. The MED is the dose producing the mean response of delta units better than the minimum dose.}
  \item{r}{Numeric. The number of iteritions to set an initial design to search c-optimal design for estimating the MED. Default is 10 and needed to be increased (for example, r=30 or 50) if the searched c-optimal design is not a true optimal. It can be verified by the verification plot.}
  \item{grid}{Numeric. The grid density to discretize the predetermined dose interval. Default is 0.01.}
  \item{epsilon}{Numeric. Stopping criterion for the algorithm to search c-optimal design for the MED. Default is 0.001.}
  \item{epsilon_w}{Numeric. Stopping criterion for the Newton Raphson method inside of the algorithm. Default is 10^-6.}
}
\value{An object of class OPT.}
\author{Seung Won Hyun, Weng Kee Wong, and Yarong Yang }
\references{Seung Won Hyun, Weng Kee Wong, and Yarong Yang (2014), VNM: An R Package for Finding Multiple-Objective Optimal Designs for the 4-Parameter Logistic Model, submitted to Journal of Statistical Software.

Seung Won Hyun and Weng Kee Wong (2014), Multiple Objective Optimal Designs to Study the Interesting Features in a Dose-Response Relationship, submitted to the International Journal of Biostatistics.
}
\examples{
  ceff2(weight=c(0.25,0.25,0.25,0.25),P=c(exp(-1/3),3),
dose=c(-1,-.5,.5,1),LB=-1,UB=1,delta=.2,r=10)
}

