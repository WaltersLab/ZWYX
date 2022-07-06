
#' Calculate Log2(Hom:Het) Ratios on Scaffold Windows
#'
#' This function is used to get the log2(Hom:Het) ratios for a specific
#' scaffold, where normalized read count values are extracted as GR
#' from the master SE object. This used internally in other functions.
#' cD is a DF of normalized read counts obtained from a GR object
#' gamety is the vector from the original SE object.
#'
#' This is primarily a utility function used by \code{runChangpointOnGR} and should
#' not normally be called directly in routine usage of ZWYX.
#'
#' @param cD A dataframe of normalized read counts obtained from a GR object
#' @param gamety The vector from the original SE object indicating homo- or
#' heterogamety.
#'
#' @return A vector of log2(Hom:Het) coverage ratios
#' @export
#'
#' @examples
getHomHetRatio <- function(cD, gamety) {
   # using drop = F allows single column matrix, not conversion to vector
   # so rowMeans will still work on it, even though it is meaningless.
   cD <- as.matrix(cD) # convert from S4 DataFrame to matrix to allow rowMeans
   dat_het <- rowMeans(cD[, gamety == "Het" , drop = FALSE])
   dat_hom <- rowMeans(cD[, gamety == "Hom" , drop = FALSE])
   homhet_ratio <- log2(dat_hom/dat_het)
   return(homhet_ratio)
}



#' Make GenomicRanges from Changepoint
#'
#' This function takes as input a genomic ranges object that was used to run
#' changepoint and returns a genomic ranges object giving the ranges of CP
#' segments and their associated means.
#'
#' This is primarily a utility function used by \code{runChangpointOnGR} and should
#' not normally be called directly in routine usage of ZWYX.
#'
#' @param gr A \code{GenomicRanges} object that was used as input for running
#' changepoint analysis.
#' @param cpt A \code{changepoint} object.
#' @param XZthr A numeric value corresponding to the log2(Hom:Het) ratio that
#' determines whether a scaffold is classified as being X|Z linked. Scaffolds
#' with log2(Hom:Het) > \code{XZthr} are considered to be located on the X or Z
#' chromosome. Default value is 0.5.
#' @param YWthr A numeric value corresponding to the log2(Hom:Het) ratio that
#' determines whether a scaffold is classified as being Y|W linked. Scaffolds
#' with log2(Hom:Het) < \code{XZthr} are considered to be located on the Y or W
#' chromosome. Default value is -0.5.
#'
#' @return A \code{GenomicRanges} object corresponding to segments detected
#' by Changepoint analysis.
#' @export
#'
#' @examples
makeGRfromCP <- function(gr, cpt, XZthr = 0.5, YWthr= -0.5){
   if (XZthr <= YWthr) {
      stop("Threshold for Z/X chrom must be greater than threshold for Y/W chromosome")
   }
   rng <- GenomicRanges::ranges(gr) # get ranges
   rng.start <- GenomicRanges::start(rng) # get start coords

   cpoints <- changepoint::cpts(cpt) # get changepoint positions
   cpcoords <- rng.start[cpoints] # get seq coords of window positions
   cpstarts <- c(1,cpcoords) # define start of segments
   cpends <- c(cpcoords-1, utils::tail(GenomicRanges::end(rng), n=1) ) # define end of segments
   cpmeans <- cpt@param.est$mean # extract means
   linkage <- as.character(cut(x = cpmeans, breaks = c(-Inf,YWthr, XZthr, Inf), labels = c("YW", "Auto", "XZ")))
   cptdf <- data.frame(seqname=GenomicRanges::seqnames(gr)[1], "start"=cpstarts, "end"=cpends, "MeanLog2HetHom" = cpmeans, "Linkage" = linkage)
   cptgr <- GenomicRanges::makeGRangesFromDataFrame(df = cptdf, keep.extra.columns = T )
   return(cptgr)
}
