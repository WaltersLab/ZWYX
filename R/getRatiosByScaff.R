#' Get lengths of each sequence
#'
#' A simple utility function that extracts the largest end coordinate for
#' each sequence in a GenomicRangesList object into a numeric vector.
#'
#' The working assumption here inside ZWYX is that the GenomicRanges objects
#' have been prepared such that the largest end coordinate for each unique
#' sequence is the final base and there indicates the length of that sequence.
#'
#' @param myGRlist A GenomicRanges List object
#'
#' @return A vector
#' @export
#'
#' @examples
getSeqLengths <- function(myGRlist) {
   max_ends <- vapply( myGRlist,
                       FUN = function(gr) {
                          max(GenomicRanges::end(GenomicRanges::ranges(gr)))
                       },
                       FUN.VALUE = numeric(1)
   )
   return(max_ends)
}






# function to generate hom:het ratios for each scaffold
# and build dataframe of values and results
# for later plotting.
#' Generate Hom:Het ratios per scaffold
#'
#' The \code{getRatiosByScaff} function generates a point estimate for each
#' scaffold (or contig) of the normalized \emph{Homogametic:Heterogametic}
#' coverage ratio. It also classifies each scaffold as being autosomal, X|Z, or
#' Y|W based on specified thresholds for the \emph{Hom:Het} ratio.
#'
#' @param se A \code{SummarizedExperiment} object that contains normalized
#' values generated via \code{\link{normalizeCounts}}
#' @param XZthr A numeric value corresponding to the log2(Hom:Het) ratio that
#' determines whether a scaffold is classified as being X|Z linked. Scaffolds
#' with log2(Hom:Het) > \code{XZthr} are considered to be located on the X or Z
#' chromosome. Default value is 0.5.
#' @param YWthr A numeric value corresponding to the log2(Hom:Het) ratio that
#' determines whether a scaffold is classified as being Y|W linked. Scaffolds
#' with log2(Hom:Het) < \code{XZthr} are considered to be located on the Y or W
#' chromosome. Default value is -0.5.
#'
#' @return A dataframe with four columns corresponding to the sequence (scaffold)
#' identifier, estimated log2(Hom:Het) ratio, assigned linkage, and sequence
#' (scaffold) length.
#' @export
#'
#' @examples
getRatiosByScaff <- function(se, XZthr = 0.5, YWthr= -0.5) {
   scaffTR_df <- getNormDataByScaff(se = se)
   gamety <- SummarizedExperiment::colData(se)$gamety
   # using drop = F allows single column matrix, not conversion to vector
   # so rowMeans will still work on it, even though it is meaningless.
   dat_het <- rowMeans(scaffTR_df[, gamety == "Het" , drop = FALSE])
   dat_hom <- rowMeans(scaffTR_df[, gamety == "Hom" , drop = FALSE])
   homhet_ratio <- log2(dat_hom/dat_het)
   linkage <- as.character(cut(x = homhet_ratio, breaks = c(-Inf,YWthr, XZthr, Inf), labels = c("YW", "Auto", "XZ")))
   se_GRList <- splitGRangesByScaffold(myGR = SummarizedExperiment::rowRanges(se))
   scaffL <- getSeqLengths(se_GRList)

   outframe <- data.frame("scaffold" = row.names(scaffTR_df),
                          "log2homhet" = homhet_ratio,
                          "Linkage" = linkage,
                          "scaffLen" = scaffL)
   return(outframe)
}
