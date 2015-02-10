\name{VNM-package}
\alias{VNM-package}
\docType{package}
\title{ Using V-algorithm and Newton-Raphson Method to Obtain Multiple-objective Optimal Design.}
\description{Using V-algorithm and Newton-Raphson method to obtain multiple-objective optimal design for estimating the shape of dose-response, the ED50 (the dose producing an effect midway between the expected responses at the extreme doses) and the MED (the minimum effective dose level) 
for the 2,3,4-parameter logistic models and for evaluating its efficiencies for the three objectives.
}
\details{
\tabular{ll}{
Package: \tab VNM\cr
Type: \tab Package\cr
Version: \tab 4.0\cr
Date: \tab 2015-02-09\cr
License: \tab GPL-2\cr
}
 }
\author{ Seung Won Hyun, Weng Kee Wong, and Yarong Yang  }
\references{Seung Won Hyun, Weng Kee Wong, and Yarong Yang (2014), VNM: An R Package for Finding Multiple-Objective Optimal Designs for the 4-Parameter Logistic Model, submitted to Journal of Statistical Software.

Seung Won Hyun and Weng Kee Wong (2014), Multiple Objective Optimal Designs to Study the Interesting Features in a Dose-Response Relationship, submitted to the International Journal of Biostatistics.
}

\keyword{package }

\examples{
MOPT(LB = -1, UB = 1, P = c(exp(-1/3),3),lambda = c(1/3,1/3), 
delta = 0.2, r = 10,grid=.01, verbose=TRUE)
}
