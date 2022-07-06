
#' Plot Chromosomal Segments
#'
#' This function provides a visualiztion of the distinct segments detected by
#' Changepoint analysis of male:female sequencing coverage ratios.
#'
#' @param InputGR The GenomicRanges object input to Changepoint analysis,
#' containing log2(Hom:Het) coverage ratios in metadata
#' @param SegmentGR The GenomicRanges object representing the results of Changepoint
#' analysis, with coordinates of segment boundaries and means as metadata
#' @param bpunit A character value, one of "bp", "kbp", or "mbp". This defines
#' the scale of the plot's X-axis coordinates to be individual base pairs,
#' kilobases, or megabases, accordingly.
#' @param addLines Logical. TRUE will include line segments drawn at the mean
#' log2(Hom:Het) value for each sequence segment detected by Changepoint
#' @param XZLineCol The color for lines indicating Z|X linkage
#' @param YWLineCol The color for lines indicating W|Y linkage
#' @param AutoLineCol The color for lines indicating autosomal linkage
#' @param seglwd Numeric, the line width parameter for lines drawn when
#'  \code{addLines} is TRUE.
#' @param ... Other standard graphical parameter arguments provided to the \code{plot}
#' function.
#'
#' @return A base-graphics plot of log2(Hom:Het) ratios in each window, along with
#' line segments indicating the sequence segments detected by Changepoint.
#' @export
#'
#' @examples
zwyxPlot <- function(InputGR, SegmentGR, bpunit = "kbp", addLines = TRUE, XZLineCol = "dodger blue", YWLineCol = "deep pink", AutoLineCol="yellow2", seglwd=1, ... ) {
   if (! any(bpunit == c("bp", "kbp", "mbp")) ) {
      stop("Value for argument `bpunit` must be one of: 'bp', 'kbp', 'mbp'")
   }
   inputdf <- as.data.frame(InputGR)
   segmentdf <- as.data.frame(SegmentGR)
   # Rescale coordinates
   if (! bpunit == "bp") {
      if (bpunit == "kbp") { bpscale <- 1000 }
      if (bpunit == "mbp") { bpscale <- 1e6 }
      inputdf[,c("start","end")] <- inputdf[,c("start","end")]/bpscale
      segmentdf[,c("start","end")] <- segmentdf[,c("start","end")]/bpscale
   }
   plot(x=inputdf$start, y = inputdf$log2homhet, ...)
   if (addLines) {
      segcols <- rep(x = AutoLineCol, times = nrow(segmentdf))
      segcols[segmentdf$Linkage == "XZ"] <- XZLineCol
      segcols[segmentdf$Linkage == "YW"] <- YWLineCol
      graphics::segments(x0=segmentdf$start, x1=segmentdf$end, y0=segmentdf$MeanLog2HetHom, col=segcols, lwd = seglwd)
   }
}
