
# method, q, penalty, and minseglength are all changepoint.mean args; ... is for other args.
#' Run Changepoint on GenomicRanges
#'
#' This function runs changepoint on counts from a single scaffold,
#' provided as GenomicRanges object. It returns a list that contains the input
#' as well as the result Changepoint results, and any filtered windows.
#'
#' This is primarily a utility function used by \code{runChangepointOnSE} and should
#' not normally be called directly in routine usage of ZWYX.
#'
#' @param gr A GenomicRanges object for a single scaffold, with associated
#' metadata containg readcounts per window for each sample.
#' @param gamety The "gamety" vector from a ZWYX-prepared SummarizedExperiment object.
#' @param method Corresponds to argument for \code{cpt.mean} function from Changepoint package
#' @param penalty Corresponds to argument for \code{cpt.mean} function from Changepoint package
#' @param minseglen Corresponds to argument for \code{cpt.mean} function from Changepoint package
#' @param XZthr A numeric value corresponding to the log2(Hom:Het) ratio that
#' determines whether a scaffold is classified as being X|Z linked. Scaffolds
#' with log2(Hom:Het) > \code{XZthr} are considered to be located on the X or Z
#' chromosome. Default value is 0.5.
#' @param YWthr A numeric value corresponding to the log2(Hom:Het) ratio that
#' determines whether a scaffold is classified as being Y|W linked. Scaffolds
#' with log2(Hom:Het) < \code{XZthr} are considered to be located on the Y or W
#' chromosome. Default value is -0.5.
#' @param ... Other arguments to be provided to \code{cpt.mean} function.
#'
#' @return A list of length 4: 1) The original GR (with ratios in column),
#' 2) GR of segments detected by Changepoint analysis, 3) the Changepoint
#' object resulting from analysis, 4) a vector of indices indicating which
#' windows had Inf coverage ratios and thus were removed prior to Changepoint
#' analysis.
#' @export
#'
#' @examples
runChangpointOnGR <- function(gr, gamety,  method = "PELT",  penalty= "Hannan-Quinn",
                              minseglen=5, XZthr = 0.5, YWthr= -0.5, ...) {
   if (XZthr <= YWthr) {
      stop("Threshold for Z/X chrom must be greater than threshold for Y/W chromosome")
   }
   countsData <- GenomicRanges::mcols(gr)  # extract metadata columns from GR object
   homhet_ratio <- getHomHetRatio(cD = countsData, gamety = gamety) # get vector of ratios
   GenomicRanges::mcols(gr)$log2homhet <- homhet_ratio
   finite.ratio <- is.finite(homhet_ratio)
   gr.finite <- gr[finite.ratio]  # subset input GR to only those windows with finite log2ratio
   out.cp <- changepoint::cpt.mean(homhet_ratio[finite.ratio], method = method,  penalty= penalty, minseglen = minseglen, ...)
   grfromcp <- makeGRfromCP(gr = gr.finite, cpt = out.cp, XZthr = XZthr, YWthr= YWthr)
   outlist <- list("InputGR" = gr, "SegmentsGR" = grfromcp, "Changepoint" = out.cp, "FilteredWindows" = which(!finite.ratio))
   return(outlist)
}


# This function takes an SE object containing read counts
# that have been normalized (or if not, performs normalization),
# decomposes SE object into GRlist by scaffold, runs changepoint on it.
# Still need to figure out passthrough arguments with default values for cpt.mean,
# but this website indicates a possible solution.
# https://stackoverflow.com/questions/35587265/pass-arguments-in-nested-function-to-update-default-arguments/35587633
#' Run Changepoint on All Scaffolds
#'
#' This function takes a ZWYX-prepared SummarizedExperiment object containing
#' normalized read counts and runs changepoint on each scaffold (contig).
#'
#' @param se A ZWYX-prepared SummarizedExperiment object containing
#' normalized read counts.
#' @param XZthr A numeric value corresponding to the log2(Hom:Het) ratio that
#' determines whether a scaffold is classified as being X|Z linked. Scaffolds
#' with log2(Hom:Het) > \code{XZthr} are considered to be located on the X or Z
#' chromosome. Default value is 0.5.
#' @param YWthr A numeric value corresponding to the log2(Hom:Het) ratio that
#' determines whether a scaffold is classified as being Y|W linked. Scaffolds
#' with log2(Hom:Het) < \code{XZthr} are considered to be located on the Y or W
#' chromosome. Default value is -0.5.
#'
#' @return A list with one element for each scaffold. Each element is itself a
#' list of length four -- the output of \code{runChangpointOnGR}.
#' For each scaffold the list contains:
#' 1) The original GR (with ratios in column),
#' 2) GR of segments detected by Changepoint analysis,
#' 3) the Changepoint
#' object resulting from analysis,
#' 4) a vector of indices indicating which
#' windows had Inf coverage ratios and thus were removed prior to Changepoint
#' analysis.
#' @export
#'
#' @examples
runChangepointOnSE <- function(se, XZthr = 0.5, YWthr= -0.5 ) {
   # Make sure SE contains normalized values.
   if (! any (names(SummarizedExperiment::assays(se)) == "normdata")) {
      stop("Missing 'normdata' assay from SummarizedExperiment input.
           You must run `normalizeCounts()` first to generate normalized values.")
   }

   # Parse SE into GRList of scaffolds
   grNorm <- getGRfromSE(se, select_assay = "normdata")
   grList <- splitGRangesByScaffold(grNorm)
   gamety.se <- SummarizedExperiment::colData(se)$gamety
   # Run Changepoint on each scaffold
   changeList <- lapply(grList, FUN = function(x) { runChangpointOnGR( gr=x, gamety = gamety.se, XZthr = XZthr, YWthr= YWthr )})
   return(changeList)
}

