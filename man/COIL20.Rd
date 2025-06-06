% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{COIL20}
\alias{COIL20}
\title{COIL-20 Object Image Dataset}
\format{
A data frame with 1,440 rows and 1,025 columns:
\describe{
\item{X1, X2, ..., X1024}{Grayscale pixel values (numeric).}
\item{class}{Factor indicating the object label (1–20).}
}
}
\usage{
COIL20
}
\description{
The COIL-20 dataset consists of grayscale images of 20 distinct objects. Each object
was placed on a turntable and imaged at pose intervals of 5 degrees, resulting in 72 images per object (totaling 1,440 images).
}
\details{
Each row represents an image sample, while each column (except the last) corresponds to a pixel value.
The final column \code{Object} indicates the object's identity (from 1 to 20).
}
\references{
Nene, S.A., Nayar, S.K., & Murase, H. (1996). Columbia Object Image Library (COIL-20). Technical Report CUCS-005-96, Department of Computer Science, Columbia University.
}
\keyword{datasets}
