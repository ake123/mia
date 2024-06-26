#' Estimate divergence
#'
#' Estimate divergence against a given reference sample.
#' 
#' @param x a \code{\link{SummarizedExperiment}} object.
#'
#' @param assay.type the name of the assay used for calculation of the
#'   sample-wise estimates.
#'   
#' @param assay_name a single \code{character} value for specifying which
#'   assay to use for calculation.
#'   (Please use \code{assay.type} instead. At some point \code{assay_name}
#'   will be disabled.)
#'
#' @param name a name for the column of the colData the results should be
#'   stored in. By default, \code{name} is \code{"divergence"}.
#'   
#' @param reference a numeric vector that has length equal to number of
#'   features, or a non-empty character value; either 'median' or 'mean'.
#'   \code{reference} specifies the reference that is used to calculate
#'   \code{divergence}. by default, \code{reference} is  \code{"median"}.
#'   
#' @param FUN a \code{function} for distance calculation. The function must
#'   expect the input matrix as its first argument. With rows as samples 
#'   and columns as features. By default, \code{FUN} is
#'   \code{vegan::vegdist}.
#'   
#' @param method a method that is used to calculate the distance. Method is
#'   passed to the function that is specified by \code{FUN}. By default,
#'   \code{method} is \code{"bray"}.
#'
#' @param ... optional arguments
#' 
#' @return \code{x} with additional \code{\link{colData}} named \code{*name*}
#' 
#' @details
#'
#' Microbiota divergence (heterogeneity / spread) within a given sample
#' set can be quantified by the average sample dissimilarity or beta
#' diversity with respect to a given reference sample.
#'
#' This measure is sensitive to sample size.
#' Subsampling or bootstrapping can be applied to equalize sample sizes
#' between comparisons.
#' 
#' @seealso
#' \code{\link[scater:plotColData]{plotColData}}
#' \itemize{
#'   \item{\code{\link[mia:estimateRichness]{estimateRichness}}}
#'   \item{\code{\link[mia:estimateEvenness]{estimateEvenness}}}
#'   \item{\code{\link[mia:estimateDominance]{estimateDominance}}}
#' }
#' 
#' @name addDivergence
#' @export
#'
#' @author Leo Lahti and Tuomas Borman. Contact: \url{microbiome.github.io}
#' 
#' @examples
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#' 
#' # By default, reference is median of all samples. The name of column where results
#' # is "divergence" by default, but it can be specified. 
#' tse <- addDivergence(tse)
#' 
#' # The method that are used to calculate distance in divergence and 
#' # reference can be specified. Here, euclidean distance and dist function from 
#' # stats package are used. Reference is the first sample.
#' tse <- addDivergence(tse, name = "divergence_first_sample", 
#'                           reference = assays(tse)$counts[,1], 
#'                           FUN = stats::dist, method = "euclidean")
#' 
#' # Reference can also be median or mean of all samples. 
#' # By default, divergence is calculated by using median. Here, mean is used.
#' tse <- addDivergence(tse, name = "divergence_average", reference = "mean")
#' 
#' # All three divergence results are stored in colData.
#' colData(tse)
#' 
NULL

#' @rdname addDivergence
#' @export
setGeneric("addDivergence",signature = c("x"),
           function(x, assay.type = assay_name, assay_name = "counts", 
                    name = "divergence", reference = "median", 
                    FUN = vegan::vegdist, method = "bray", ...)
             standardGeneric("addDivergence"))

#' @rdname addDivergence
#' @export
setMethod("addDivergence", signature = c(x="SummarizedExperiment"),
    function(x, assay.type = assay_name, assay_name = "counts", 
             name = "divergence", reference = "median", 
             FUN = vegan::vegdist, method = "bray", ...){
        
        ################### Input check ###############
        # Check assay.type
        .check_assay_present(assay.type, x)
        # Check name
        if(!.is_non_empty_character(name) || length(name) != 1L){
            stop("'name' must be a non-empty character value.",
                 call. = FALSE)
        }
        # Check reference
        # If "reference" is not right: 
        # it is not numeric or character
        # its length does not equal to number of samples when it's numeric,
        # reference is not "median" or "mean"
        reference_stop_msg <- 
            paste0("'reference' must be a numeric vector that has lenght equal",
                   " to number of features, or 'reference' must be either",
                   " 'median' or 'mean'.")
        if( !(is.numeric(reference) || is.character(reference)) ){
            stop(reference_stop_msg, call. = FALSE)
        } else {
            if( is.numeric(reference) && length(reference) != nrow(x) ){
                stop(reference_stop_msg, call. = FALSE)
            }
            if( is.character(reference) && length(reference) != 1L && 
               !any(c("median","mean") %in% reference) ){
                stop(reference_stop_msg, call. = FALSE)
            }
        }

        ################# Input check end #############
        divergence <- .calc_reference_dist(mat = assay(x, assay.type),
                                       reference = reference, 
                                       FUN = FUN,
                                       method = method, ...)

        .add_values_to_colData(x, list(divergence), name)

    }
)

############################## HELP FUNCTIONS ##############################


.calc_reference_dist <- function(mat, reference, FUN = stats::dist, method, ...){

    # Calculates median or mean if that is specified
    if (is.character(reference)) {
        if( "median" %in% reference || "mean" %in% reference ){
            reference <- apply(mat, 1, reference)
        } else if( !reference %in% colnames(mat) ) {
            stop("Reference ", reference, " not recognized.", call. = FALSE)
        }
    }

    # Distance between all samples against one reference sample
    # FIXME: could be be optimzed with sweep / parallelization
    v <- seq_len(ncol(mat))
    vapply(v, function (i) {FUN(rbind(mat[,i], reference), method=method, ...)},FUN.VALUE = numeric(1))
}

