## ---- include = FALSE----------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------------------------------
library(ZWYX)
data("dmel_counts")
head(dmel_counts)

## ------------------------------------------------------------------------------------------------
suppressPackageStartupMessages( library(SummarizedExperiment) )
counts_ranges.dmel <- makeGRangesFromDataFrame(df = dmel_counts, seqnames.field = "scaffold")
counts_ranges.dmel

## ------------------------------------------------------------------------------------------------
sample.dmel <- data.frame("sampleID" = names(dmel_counts)[-(1:3)],
                          "gamety" = c("Het","Het","Hom","Hom","Het","Hom")
               )

## ------------------------------------------------------------------------------------------------
count.se <- SummarizedExperiment(assays = list("counts" = dmel_counts[,-(1:3)]), 
                                 rowRanges = counts_ranges.dmel,
                                 colData = sample.dmel)
count.se

## ------------------------------------------------------------------------------------------------
norm.tot <- normalizeCounts(se = count.se, method = "totalreads")
norm.tot

## ------------------------------------------------------------------------------------------------
head(assays(norm.tot)[["normdata"]])

## ------------------------------------------------------------------------------------------------
scaffoldRatiosTot <- getRatiosByScaff(norm.tot)
head(scaffoldRatiosTot)

## ---- fig.width=6--------------------------------------------------------------------------------
library(ggplot2)
ggplot(data = scaffoldRatiosTot, aes(x=log10(scaffLen), y = log2homhet, col = Linkage)) + geom_point() +
          labs(y = "Log2(Female:Male)", x = "Log10(Fragment Length)") + 
   scale_color_manual(values = c("black","red", "blue")) + theme_light()

## ---- fig.width = 6------------------------------------------------------------------------------

scaffoldRatiosTot$chimera <- grepl(x = scaffoldRatiosTot$scaffold, pattern="chimera")
scaffoldRatiosTot$chrom <- gsub(x=scaffoldRatiosTot$scaffold, 
                                pattern = "(\\w)\\w?\\..*", rep="\\1", perl = T)

ggtot <- ggplot(data = scaffoldRatiosTot, aes(x=log10(scaffLen), y = log2homhet, 
               col = chimera, pch = chrom)) + geom_point() + 
               scale_color_manual(values = c("black","red"))  + 
               labs(y = "Log2(Female:Male)", x = "Log10(Fragment Length)") + 
               theme_light()
print(ggtot) 

## ------------------------------------------------------------------------------------------------
cp.out <- runChangepointOnSE(se = norm.tot)

## ------------------------------------------------------------------------------------------------
cp.out[["2L.2_5878429_9317088_chimeraAX"]]

## ------------------------------------------------------------------------------------------------
outsegs.df <- mergeSegments(changeList = cp.out[1:3], asDF = TRUE)
head(outsegs.df)

## ---- fig.width=7--------------------------------------------------------------------------------
# Extract a known chimeric scaffold for demonstration
cp <- cp.out[[2]]
zwyxPlot(InputGR = cp[[1]], SegmentGR = cp[[2]], bpunit = "mbp", 
         pch = 19, cex = 0.2, seglwd = 3, main = unique(seqnames(cp[[1]])) ,
         xlab = "Position (Mbp)", ylab = "Log2(Het:Hom")

## ---- fig.height=10, fig.width=7-----------------------------------------------------------------
chimera.ind <- which(scaffoldRatiosTot$chimera)


par(mfrow = c(4,2))
for (nscaff in chimera.ind) {
   zwyxPlot(InputGR = cp.out[[nscaff]][[1]], SegmentGR = cp.out[[nscaff]][[2]], 
            bpunit="mbp", pch = 19, cex = .2, seglwd = 3, 
            main = names(cp.out [nscaff]),  xlab = "Position (Mbp)", ylab = "Log2(Het:Hom" )
}

