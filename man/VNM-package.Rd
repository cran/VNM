\name{VNM-package}
\alias{VNM-package}
\docType{package}
\title{ Using V-algorithm and Newton-Raphson Method to Obtain Multiple-objective Optimal Design.}
\description{Using V-algorithm and Newton-Raphson method to obtain multiple-objective optimal design for estimating the shape of dose-response,
 the ED50 (the dose producing an effect mid-way between the expected responses at the extreme doses) and the MED (the minimum effective dose level) 
for the 4-parameter logistic model and for evaluating its efficiencies for the three objectives. 
}
\details{
\tabular{ll}{
Package: \tab VNM\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2014-10-31\cr
License: \tab GPL-2\cr
}
 }
\author{ Seung Won Hyun, Weng Kee Wong, and Yarong Yang  }
\references{Seung Won Hyun, Weng Kee Wong, and Yarong Yang (2014), VNM: An R Package for Finding Multiple-Objective Optimal Designs for the 4-Parameter Logistic Model, submitted to Journal of Statistical Software.

Seung Won Hyun and Weng Kee Wong (2014), Multiple Objective Optimal Designs to Study the Interesting Features in a Dose-Response Relationship, submitted to the International Journal of Biostatistics.
}

\keyword{package }

\examples{
    MOPT(LB=-3, UB=0, P=c(0.137,1.563,0.453,-0.825), lambda=c(1/2,1/3), delta=-1)
}
