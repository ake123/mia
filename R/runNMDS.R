#' Perform non-metric MDS on sample-level data
#'
#' Perform non-metric multi-dimensional scaling (nMDS) on samples, based on the
#' data in a \code{SingleCellExperiment} object.
#'
#' @param x For \code{getNMDS}, a numeric matrix of expression values
#'   where rows are features and columns are cells.
#'   Alternatively, a \code{TreeSummarizedExperiment} containing such a matrix.
#'
#'   For \code{addNMDS} a \linkS4class{SingleCellExperiment}
#'
#' @param ncomponents Numeric scalar indicating the number of NMDS dimensions
#'   to obtain.
#'
#' @param ntop Numeric scalar specifying the number of features with the highest
#'   variances to use for dimensionality reduction.
#'
#' @param subset.row Vector specifying the subset of features to use for
#'   dimensionality reduction. This can be a character vector of row names, an
#'   integer vector of row indices or a logical vector.
#' 
#' @param subset_row Deprecated. Use \code{subset.row} instead.
#'
#' @param scale Logical scalar, should the expression values be standardized?
#'
#' @param keep.dist Logical scalar indicating whether the \code{dist} object
#'   calculated by \code{FUN} should be stored as \sQuote{dist} attribute of
#'   the matrix returned/stored by \code{getNMDS}/ \code{addNMDS}.
#' 
#' @param keep_dist Deprecated. Use \code{keep.dist} instead.
#'
#' @param transposed Logical scalar, is x transposed with cells in rows?
#'
#' @param assay.type a single \code{character} value for specifying which
#'   assay to use for calculation.
#'   
#' @param exprs_values a single \code{character} value for specifying which
#'   assay to use for calculation.
#'   (Please use \code{assay.type} instead.)
#' 
#' @param assay_name a single \code{character} value for specifying which
#'   assay to use for calculation.
#'   (Please use \code{assay.type} instead. At some point \code{assay_name}
#'   will be disabled.)
#'
#' @param FUN a \code{function} or \code{character} value with a function
#'   name returning a \code{\link[stats:dist]{dist}} object
#'
#' @param nmds.fun a \code{character} value to choose the scaling
#'   implementation, either \dQuote{isoMDS} for
#'   \code{\link[MASS:isoMDS]{MASS::isoMDS}} or \dQuote{monoMDS} for
#'   \code{\link[vegan:monoMDS]{vegan::monoMDS}}
#'   
#' @param nmdsFUN Deprecated. Use \code{nmds.fun} instead.
#'
#' @param ... additional arguments to pass to \code{FUN} and
#'   \code{nmds.fun}.
#'
#' @param dimred String or integer scalar specifying the existing dimensionality
#'   reduction results to use.
#'
#' @param ndimred Integer scalar or vector specifying the dimensions to use if
#'   dimred is specified.
#' 
#' @param n_dimred Deprecated. Use \code{ndimred} instead.
#'
#' @param altexp String or integer scalar specifying an alternative experiment
#'   containing the input data.
#'
#' @param name String specifying the name to be used to store the result in the
#'   reducedDims of the output.
#'
#' @return For \code{getNMDS}, a matrix is returned containing the MDS
#'   coordinates for each sample (row) and dimension (column).
#'
#' @details
#' Either \code{\link[MASS:isoMDS]{MASS::isoMDS}} or
#' \code{\link[vegan:monoMDS]{vegan::monoMDS}} are used internally to compute
#' the NMDS components. If you supply a custom \code{FUN}, make sure that
#' the arguments of \code{FUN} and \code{nmds.fun} do not collide.
#'
#' @name runNMDS
#'
#' @seealso
#' \code{\link[MASS:isoMDS]{MASS::isoMDS}},
#' \code{\link[vegan:monoMDS]{vegan::monoMDS}}
#' for NMDS component calculation.
#'
#' \code{\link[scater:plotReducedDim]{plotMDS}}, to quickly visualize the
#' results.
#'
#' @author
#' Felix Ernst
#'
#' @examples
#' # generate some example data
#' mat <- matrix(1:60, nrow = 6)
#' df <- DataFrame(n = c(1:6))
#' tse <- TreeSummarizedExperiment(assays = list(counts = mat),
#'                                 rowData = df)
#' #
#' getNMDS(tse)
#'
#' #
#' data(esophagus)
#' esophagus <- addNMDS(esophagus, FUN = vegan::vegdist, name = "BC")
#' esophagus <- addNMDS(esophagus, FUN = vegan::vegdist, name = "euclidean",
#'                      method = "euclidean")
#' reducedDims(esophagus)
NULL


#' @rdname runNMDS
#' @export
setGeneric("getNMDS", function(x, ...) standardGeneric("getNMDS"))


.format_nmds_isoMDS <- function(nmds){
    ans <- nmds$points
    attr(ans, "Stress") <- nmds$stress
    ans
}

#' @importFrom vegan scores
.format_nmds_monoMDS <- function(nmds){
    ans <- scores(nmds)
    colnames(ans) <- NULL
    attr(ans, "Stress") <- nmds$stress
    ans
}

.format_nmds <- function(nmds, nmds.fun, sample_names){
    ans <- switch(nmds.fun,
                  "isoMDS" = .format_nmds_isoMDS(nmds),
                  "monoMDS" = .format_nmds_monoMDS(nmds))
    rownames(ans) <- sample_names
    ans
}

.get_nmds_args_isoMDS <- function(args){
    args[c("maxit","trace","tol","p")]
}

.get_nmds_args_monoMDS <- function(args){
    args[c("model","threshold","maxit","weakties","stress","scaling","pc",
           "smin","sfgrmin","sratmax")]
}

.get_nmds_args <- function(nmds.fun, ...){
    args <- list(...)
    args <- switch(nmds.fun,
                   "isoMDS" = .get_nmds_args_isoMDS(args),
                   "monoMDS" = .get_nmds_args_monoMDS(args))
    args <- args[!vapply(args,is.null,logical(1))]
    args
}

#' @importFrom MASS isoMDS
#' @importFrom stats cmdscale
#' @importFrom vegan vegdist monoMDS
.calculate_nmds <- function(x, FUN = vegdist, 
                            nmds.fun = nmdsFUN,
                            nmdsFUN = c("isoMDS","monoMDS"),
                            ncomponents = 2, ntop = 500, subset.row = subset_row, 
                            subset_row = NULL, scale = FALSE, transposed = FALSE,
                            keep.dist = keep_dist,
                            keep_dist = FALSE, ...){
    nmds.fun <- match.arg(nmds.fun)
    nmdsArgs <- .get_nmds_args(nmds.fun, ...)
    if(!transposed) {
        x <- .get_mat_for_reddim(x, subset_row = subset.row, ntop = ntop,
                                 scale = scale)
    }
    x <- as.matrix(x)
    sample_names <- rownames(x)
    sample_dist <- do.call(FUN,
                           c(list(x),
                             list(...)))
    attributes(sample_dist) <- attributes(sample_dist)[c("class","Size")]
    y <- cmdscale(sample_dist, k = ncomponents)
    ans <- do.call(nmds.fun,
                   c(list(sample_dist, y = y, k = ncomponents),
                     nmdsArgs))
    ans <- .format_nmds(ans, nmds.fun, sample_names)
    if (keep.dist) {
        attr(ans,"dist") <- sample_dist
    }
    ans
}

#' @rdname runNMDS
#' @export
setMethod("getNMDS", "ANY", .calculate_nmds)

#' @rdname runNMDS
#' @importFrom SummarizedExperiment assay
#' @export
setMethod("getNMDS", "SummarizedExperiment",
    function(x, ..., assay.type = assay_name, assay_name = exprs_values, 
             exprs_values = "counts", FUN = vegdist) {
        .calculate_nmds(assay(x, assay.type), FUN = FUN, ...)
    }
)

#' @rdname runNMDS
#' @export
setMethod("getNMDS", "SingleCellExperiment",
    function(x, ..., assay.type = assay_name, assay_name = exprs_values, 
            exprs_values = "counts", dimred = NULL, ndimred = n_dimred, 
            n_dimred = NULL, FUN = vegdist){
        mat <- .get_mat_from_sce(x, exprs_values = assay.type,
                                 dimred = dimred, n_dimred = ndimred)
        getNMDS(mat, transposed = !is.null(dimred), FUN = FUN,...)
    }
)

#' @rdname runNMDS
#' @export
#' @aliases getNMDS
calculateNMDS <- function(x,...){
    getNMDS(x,...)
}

#' @rdname runNMDS
#' @export
#' @importFrom SingleCellExperiment reducedDim<- altExp
addNMDS <- function(x, ..., altexp = NULL, name = "NMDS") {
    if (!is.null(altexp)) {
        y <- altExp(x, altexp)
    } else {
        y <- x
    }
    reducedDim(x, name) <- getNMDS(y, ...)
    x
}

#' @rdname runNMDS
#' @export
#' @aliases addNMDS
runNMDS <- function(x,...){
    addNMDS(x,...)
}
