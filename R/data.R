#' Preprocessed data for vignettes and examples 
#'
#' The following data objects are prepared for use
#' in the alpine vignette and examples pages,
#' as the preparation of these objects requires
#' either long running time or a large amount of disk
#' space.
#'
#' \itemize{
#'   \item \strong{ebt.fit} - the GRangesList prepared in the vignette
#' for fitting the bias models
#'   \item \strong{fitpar} - the fitted parameters, similar to those
#' made in the vignette, but using \code{minsize=80} and \code{maxsize=350}
#'   \item \strong{fitpar.small} - the fitted parameters from the
#' vignette, returned by fitBiasModels
#'   \item \strong{res} - the results object from the vignette,
#' returned by estimateAbundance
#'   \item \strong{ebt.theta} - the GRangesList prepared in the vignette
#' for running estimateAbundance
#'   \item \strong{genes.theta} - the names of genes used in the vignette
#' for running estimateAbundance
#'   \item \strong{txdf.theta} - the DataFrame of gene and transcript
#' information used in the vignette for running estimateAbundance
#' }
#'
#' @docType data
#' @name preprocessedData
#' @aliases ebt.fit ebt.theta fitpar fitpar.small genes.theta res txdf.theta
#' 
#' 
#' @format \code{ebt.fit} and \code{ebt.theta} are GRangesList.
#' \code{fitpar}, \code{fitpar.small}, \code{res} are lists created
#' by alpine functions. \code{genes.theta} is a character vector.
#' \code{txdf.theta} is a DataFrame.
#' @source See vignette for details of object construction.
#' The alignments come from alpineData (4 samples from GEUVADIS project),
#' the Ensembl gene annotations come from \code{Homo_sapiens.GRCh38.84.gtf},
#' and the genome is \code{BSgenome.Hsapiens.NCBI.GRCh38}.
NULL
