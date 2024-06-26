% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{bin_profile}
\alias{bin_profile}
\title{Creating a Class by Attribute Matrix}
\usage{
bin_profile(natt)
}
\arguments{
\item{natt}{An integer containing the number of assessed attributes.}
}
\value{
`profiles` A tibble containing a class by attribute matrix listing
which attributes are mastered by each latent class.
}
\description{
Automating the creation of Class by Attribute Matrix
}
\examples{
bin_profile(natt = 3)
}
