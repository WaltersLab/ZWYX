#' Make \code{GenomicRanges} from \code{SummarizedExperiment}
#'
#' The \code{getGRfromSE} function is a simple utility to extract a single
#' assay matrix from a \code{SummarizedExperiment} object and return it as
#' a corresponding stand-alone \code{GenomicRanges} object.
#'
#' The \code{getGRfromSE} function is typically used inside the \code{getNormDataByScaff}
#' function and should not normally be called direction in routine usage of ZWYX.
#'
#'
#' @param se A \code{SummarizedExperiment} object, typically prepared for
#' analysis in ZWYX
#' @param select_assay Indicates which assay from the input \code{SummarizedExperiment}
#' object will be extracted into a \code{GenomicRanges} object. Input value can
#' be either a single integer or character value corresponding to named assay; this
#' value is used for subsetting a single assay from the input \code{SummarizedExperiment}.
#'
#' @return A \code{GenomicRanges} object corresponding to the selected assay.
#' @export
#'
#' @examples
getGRfromSE <- function(se, select_assay = 1) {
   assay_data <- SummarizedExperiment::assays(se)[[select_assay]]
   myRanges <- SummarizedExperiment::rowRanges(se)
   SummarizedExperiment::mcols(myRanges) <- assay_data
   return(myRanges)
}



# Function to make GRanges list of windows from each scaffold.
# `split` sorts output, so calling unique and rearranging list
# will keep the order consistent with input.
#' Split \code{GenomicRanges} by Sequence
#'
#' The \code{splitGRangesByScaffold} function takes a \code{GenomicRanges}
#' object and splits it into a list of \code{GenomicRanges}, one per unique
#' sequence identified in the input, which typically corresponds to the contig
#' or scaffold of the genome assembly being analyze with ZWYX.
#'
#' @param myGR A \code{GenomicRanges} object.
#'
#' @return A list of \code{GenomicRanges} objects, one per unique sequence
#' @export
#'
#' @examples
splitGRangesByScaffold <- function(myGR) {
   scaffolds <- unique(GenomicRanges::seqnames(myGR))
   myGRlist <- split(x=myGR, f=GenomicRanges::seqnames(myGR))
   myGRlist <- myGRlist[scaffolds]
   return(myGRlist)
}




# function to get normalized coverage values in sample per scaffold
# Must run `normalizeCounts()` first to get `normdata` assay
# Extracts normalized windows as metadata in GRanges object
# Then splits GRanges object into list by scaffold.
# Then generates a single a single normalized depth value per scaffold.
# If normalized by total reads, then the sum of normalized values across windows is used
# If normalized by window median, then the median across normalized values is used.
#' Obtain a normalized value for each scaffold
#'
#' The \code{getNormDataByScaff} function generates a single representative
#' value for normalized coverage depth for each scaffold (or contig) for each
#' sample.
#'
#' @param se A \code{SummarizedExperiment} object that contains normalized
#' values generated via \code{\link{normalizeCounts}}
#'
#' @return A dataframe with one row per scaffold (or contig) containing a single
#' representative value for coverage depth in each sample.
#' @export
#'
#' @examples
getNormDataByScaff  <- function(se) {
   if (! any (names(SummarizedExperiment::assays(se)) == "normdata")) {
      stop("Missing 'normdata' assay from SummarizedExperiment input.
           You must run `normalizeCounts()` first to generate normalized values.")
   }
   if (! (all(SummarizedExperiment::colData(se)$NormMethod == "totalreads") ||
          all(SummarizedExperiment::colData(se)$NormMethod == "windowmedian")) ) {
      stop("Cannot determine method of normalization performed on read counts per window.
           `colData(se)$NormMethod` must all be either 'totalreads' or 'windowmedian'
           Run `normalizeCounts()` first to generate normalized values and indicate method.")
   }

   normmeth <- SummarizedExperiment::colData(se)$NormMethod[1]

   gr_counts <- getGRfromSE(se, select_assay =  "normdata")
   grListByScaff <- splitGRangesByScaffold(gr_counts)

   if (normmeth == "totalreads") {
      scaffoldData <- lapply(grListByScaff, function(x) {
         apply( X = GenomicRanges::mcols(x), MARGIN = 2, FUN = sum )
      }
      )
   }
   if (normmeth == "windowmedian") {
      scaffoldData <- lapply(grListByScaff, function(x) {
         apply( X = GenomicRanges::mcols(x), MARGIN = 2, FUN = stats::median )
      }
      )
   }
   scaffTR_df <- do.call(what = rbind, args = scaffoldData)
   row.names(scaffTR_df) <- names(grListByScaff)
   return(scaffTR_df)
}


