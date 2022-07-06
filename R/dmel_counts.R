#' Read count data for 6 Drosophila samples
#'
#' A dataset containing read counts in 2 Kbp windows from 6 different Illumina
#' DNA sequencing experiments of \emph{Drosophila melanogaster}. Paired end
#' reads were aligned to a modified version of the Dmel flybase version r.6.43
#' chromosomal assembly. The assembly has been altered in two ways. First, each
#' chromosome (or arm) was divided into several smaller chunks. Additionally,
#' specific translocations between autosomes and sex chromosomes were introduced
#' to simulate chimeric autosome-sexchrom misassemblies.
#'
#' @format A dataframe with 9 rows and 68818 columns.
#' \describe{
#' \item{scaffold}{The identifier for each modified chromosome segment created
#' from the original \emph{D. melanogaster} r.6.43 assembly}
#' \item{start end}{The start and end coordinates of the windows analyzed for
#' read counts}
#' \item{SRR#####}{Six columns corresponding to the six sequencing experiments
#' included in this data set.}
#' }
#'
#' @section Sample Sexes:
#' \describe{
#' \item{Male Samples}{SRR5217228, SRR5217251, SRR5217184}
#' \item{Female Samples}{SRR5217254, SRR5217250, SRR5217239}
#' }
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/bioproject/PRJNA369048}
"dmel_counts"
