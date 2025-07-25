% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{qGrab}
\alias{qGrab}
\title{Read and delete a .qs file}
\usage{
qGrab(qsFile)
}
\arguments{
\item{qsFile}{Name of .qs file with path.}
}
\value{
The content of the .qs file.
}
\description{
This functions reads a .qs file, deletes it, and returns its content.
}
\examples{
library(qs)
qsave(c(1, 2, 3), 'temp.qs')
qGrab('temp.qs')

}
