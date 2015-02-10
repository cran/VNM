\name{S.Weight}
\alias{S.Weight}
\title{Identify optimal weights for given dose levels.}
\description{
  Identify optimal weights that maximize the efficiency for estimating 3 objectives (model parameters, the ED50, the MED) for given dose levels under the 4-parameter logistic model. The weights represent the proportional allocations of subjects to given dose levels.  
}
\usage{
S.Weight(X,P,lambda,delta,epsilon_w)
}
\arguments{
  \item{X}{ A numeric vector. Given dose levels to search the optimal weights.}
  \item{P}{ A numeric vector. Solicited information on nominal values for the vector. P=(p1, p2, p3, p4), where p1 is the lower limit of the response, p2 is Emax, p3 is the ED50 and p4 is the slope at the ED50.}
  \item{lambda}{ A numeric vector. User select weights. lambda=c(q1, q2), where q1,q2 represent weights for estimating model parameter and estimating the ED50 respectively. They are non-negative and q1+q2<=1.}
  \item{delta}{ Numeric. Predetermined clinically significant effect to define the MED. The MED is the dose producing the mean response of dt units better than the minimum dose.}
  \item{epsilon_w}{ Numeric. Stopping criterion for the Newton Raphson method to search the optimal weights for the given dose levels. Default is 10^-6.}
} 
\value{An object of class SW.}
\author{Seung Won Hyun, Weng Kee Wong, and Yarong Yang }
\references{Seung Won Hyun, Weng Kee Wong, and Yarong Yang (2014), VNM: An R Package for Finding Multiple-Objective Optimal Designs for the 4-Parameter Logistic Model, submitted to Journal of Statistical Software.

Seung Won Hyun and Weng Kee Wong (2014), Multiple Objective Optimal Designs to Study the Interesting Features in a Dose-Response Relationship, submitted to the International Journal of Biostatistics.
}
\examples{
   S.Weight(X=c(-6.91, 2.22, 3.75, 4.60),P=c(22,16.8,70,1),lambda=c(1/3,1/3),delta=5)
}

