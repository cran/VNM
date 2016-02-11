\name{PAR}
\Rdversion{1.1}
\docType{class}
\alias{PAR-class}

\title{Class to contain the variables for the verification plots from function MOPT, ceff1, ceff2, and Deff.}
\description{
For function MOPT, ceff1, ceff2, and Deff, the plot of the sensitivity function over the given
dose interval can be used to verify the optimality of the obtained design. Class PAR is used to contain the variables for the plot.
}
\section{Objects from the Class}{
new("PAR",fid="character",LB="numeric",UB="numeric",grid="numeric",ds="vector")
}
\section{Slots}{
  \describe{
    \item{\code{fid}:}{A character showing the name of the function from where the four arguments are.}
    \item{\code{LB}:}{A numeric number showing the predetermined lower bound of the dose range for the log dose.}
    \item{\code{UB}:}{A numeric number showing the predetermined upper bound of the dose range for the log dose.}
    \item{\code{grid}:}{A numeric number showing the grid density to discretize the predetermined dose interval.}
    \item{\code{ds}:}{A numeric vector represents the values of criterion (an absolute value of the difference between the sensitivity function and 1) computed at each iteration. }
  }
}
\references{Seung Won Hyun, Weng Kee Wong, and Yarong Yang (2014), VNM: An R Package for Finding Multiple-Objective Optimal Designs for the 4-Parameter Logistic Model, submitted to Journal of Statistical Software.

Seung Won Hyun and Weng Kee Wong (2015), Multiple Objective Optimal Designs to Study the Interesting Features in a Dose-Response Relationship, accepted by the International Journal of Biostatistics.
}
\author{
Seung Won Hyun, Weng Kee Wong, and Yarong Yang
}
\examples{
showClass("PAR")
}
\keyword{classes}


