% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cpUtils.R
\name{getHomHetRatio}
\alias{getHomHetRatio}
\title{Calculate Log2(Hom:Het) Ratios on Scaffold Windows}
\usage{
getHomHetRatio(cD, gamety)
}
\arguments{
\item{cD}{A dataframe of normalized read counts obtained from a GR object}

\item{gamety}{The vector from the original SE object indicating homo- or
heterogamety.}
}
\value{
A vector of log2(Hom:Het) coverage ratios
}
\description{
This function is used to get the log2(Hom:Het) ratios for a specific
scaffold, where normalized read count values are extracted as GR
from the master SE object. This used internally in other functions.
cD is a DF of normalized read counts obtained from a GR object
gamety is the vector from the original SE object.
}
\details{
This is primarily a utility function used by \code{runChangpointOnGR} and should
not normally be called directly in routine usage of ZWYX.
}
