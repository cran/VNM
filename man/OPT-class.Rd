\name{OPT}
\Rdversion{1.1}
\docType{class}
\alias{OPT-class}
\alias{plot,OPT-method}
\alias{summary,OPT-method}

\title{Class to contain the results from function MOPT, ceff1, ceff2, and Deff.}
\description{
The function MOPT, ceff1, ceff2, Deff return object of class OPT that contains the variables for the verification plot, the optimal matrix, and the efficiency.
}
\section{Objects from the Class}{
new("OPT",Par="PAR",OPT=new("matrix"),Eff=new("numeric"))
}
\section{Slots}{
  \describe{
    \item{\code{Par}:}{An object of class PAR containing the variables for the verification plot.}
    \item{\code{Opt}:}{A numeric matrix showing the multiple-objective optimal design for estimating model parameters, the ED50, and the MED (for function MOPT), or c-optimal design for ED50 (for function ceff1),
   or c-optimal design for MED (for function ceff2), or D-optimal design (for function Deff).}
    \item{\code{Eff}:}{A numeric number showing the c-efficiency for ED50 (for function ceff1), or c-efficiency for MED (for function ceff2), or D-efficiency (for function Deff).}
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "OPT")}: ... }
    \item{summary}{\code{signature(object = "OPT")}: ... }
	 }
}
\references{Seung Won Hyun, Weng Kee Wong, and Yarong Yang (2014), VNM: An R Package for Finding Multiple-Objective Optimal Designs for the 4-Parameter Logistic Model, submitted to Journal of Statistical Software.

Seung Won Hyun and Weng Kee Wong (2014), Multiple Objective Optimal Designs to Study the Interesting Features in a Dose-Response Relationship, submitted to the International Journal of Biostatistics.
}
\author{
Seung Won Hyun, Weng Kee Wong, and Yarong Yang
}
\examples{
showClass("OPT")
}
\keyword{classes}
