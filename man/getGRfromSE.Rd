% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getNormDataByScaff.R
\name{getGRfromSE}
\alias{getGRfromSE}
\title{Make \code{GenomicRanges} from \code{SummarizedExperiment}}
\usage{
getGRfromSE(se, select_assay = 1)
}
\arguments{
\item{se}{A \code{SummarizedExperiment} object, typically prepared for
analysis in ZWYX}

\item{select_assay}{Indicates which assay from the input \code{SummarizedExperiment}
object will be extracted into a \code{GenomicRanges} object. Input value can
be either a single integer or character value corresponding to named assay; this
value is used for subsetting a single assay from the input \code{SummarizedExperiment}.}
}
\value{
A \code{GenomicRanges} object corresponding to the selected assay.
}
\description{
The \code{getGRfromSE} function is a simple utility to extract a single
assay matrix from a \code{SummarizedExperiment} object and return it as
a corresponding stand-alone \code{GenomicRanges} object.
}
\details{
The \code{getGRfromSE} function is typically used inside the \code{getNormDataByScaff}
function and should not normally be called direction in routine usage of ZWYX.
}
