\name{VNM-package}
\alias{VNM-package}
\docType{package}
\title{ Obtaining multiple objective optimal design.}
\description{ This package is for obtaining multiple objective optimal design for 
estimating three objectives: estimating model paramters; estimating ED50; estimating 
the minimum effective dose (MED).
}
\details{
\tabular{ll}{
Package: \tab VNM\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-05-15\cr
License: \tab GPL-2\cr
}
 }
\author{ Seung Won Hyun, Weng Kee Wong, and Yarong Yang  }
\references{ Multiple Objective Optimal Designs to Study the
Interesting Features in a Dose-Response Relationship, forthcoming.  }

 \keyword{package }

\examples{
   LB1=-3
   UB1=0
   nit1<-10
   T1<-c(-1.563,-0.825,-0.653,1.700)
   q11<-1/2
   q22<-1/3
   dt1<--1
   MOPT(LB=LB1, UB=UB1, nit=nit1, T=T1, q1=q11, q2=q22, dt=dt1)
}
