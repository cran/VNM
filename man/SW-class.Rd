\name{SW}
\Rdversion{1.1}
\docType{class}
\alias{SW-class}
\alias{summary,SW-method}

\title{Class to contain the results from function S. Weight.}
\description{
 The function S.Weight returns an object of class SW that contains the optimal weights and the first and second derivatives of the criterion.
}
\section{Objects from the Class}{
new("SW",Opt.W=new("matrix"),First.C=new("vector"),Second.C=new("vector"))
}
\section{Slots}{
  \describe{
    \item{\code{Opt.W}:}{A numeric matrix showing the optimal weights.}
    \item{\code{First.C}:}{A numeric vector showing the first derivative of the criterion.}
    \item{\code{Second.C}:}{A numeric vector showing the second derivative of the criterion.}
  }
}
\section{Methods}{
  \describe{
    \item{summary}{\code{signature(object = "SW")}: ... }
	 }
}
\references{Seung Won Hyun, Weng Kee Wong, and Yarong Yang (2014), VNM: An R Package for Finding Multiple-Objective Optimal Designs for the 4-Parameter Logistic Model, submitted to Journal of Statistical Software.

Seung Won Hyun and Weng Kee Wong (2015), Multiple Objective Optimal Designs to Study the Interesting Features in a Dose-Response Relationship, accepted by the International Journal of Biostatistics.
}
\author{
Seung Won Hyun, Weng Kee Wong, and Yarong Yang
}
\examples{
showClass("SW")
}
\keyword{classes}
