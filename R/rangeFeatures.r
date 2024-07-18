#' @import data.table
#' @import Biostrings
#' @importFrom dplyr %>%

#' @export
rangeFeatures <- function(ov, queryRanges, subjectRanges, featureType, allowMiss = FALSE) {

    frangedt <- switch(
        featureType,
        categorical = rangeCatFeature(ov, queryRanges, subjectRanges, allowMiss),
        stop("Invalid 'featureType'")
    )

    return(frangedt)

}

getNoOverlapQuery <- function(ov, queryRanges) {

    inov <- setdiff(1:length(queryRanges), S4Vectors::queryHits(ov))
    novdt <- data.table::data.table(
        seqnames = as.character(GenomicRanges::seqnames(queryRanges)[inov]),
        start = GenomicRanges::start(queryRanges)[inov],
        end = GenomicRanges::end(queryRanges)[inov]
    )

    return(novdt)

}

rangeCatFeature <- function(ov, queryRanges, subjectRanges, allowMiss) {

    # make a data table with the overlap results
    ovdt <- data.table::data.table(
        seqnames = as.character(GenomicRanges::seqnames(queryRanges)[S4Vectors::queryHits(ov)]),
        start = GenomicRanges::start(queryRanges)[S4Vectors::queryHits(ov)],
        end = GenomicRanges::end(queryRanges)[S4Vectors::queryHits(ov)],
        feature = subjectRanges$feature[S4Vectors::subjectHits(ov)]
    )

    # identify query ranges with multiple hits
    qovdt <- ovdt[, .SD, .SDcols = -c("feature")]
    udupdt <- unique(qovdt[duplicated(qovdt)])

    # join feature values for queries with multiple hits
    udupdt[, "duplicate" := TRUE]
    ovdt[udupdt, "duplicate" := i.duplicate, on = c("seqnames", "start", "end")]
    dupdt <- ovdt[order(feature)][ # order is so that joined labels always appear in the same order (permutations are not important, just combinations)
        duplicate == TRUE,
        list("feature" = stringi::stri_join(feature, collapse = "")),
        by = c("seqnames", "start", "end")
    ]

    # handle query ranges with no overlap with the subject
    novdt <- getNoOverlapQuery(ov, queryRanges)
    novdt[, "feature" := NA_character_]
    if (!allowMiss) novdt[, "feature" := "notAssignable"]

    frangedt <- rbind(ovdt[is.na(duplicate), .SD, .SDcols = -c("duplicate")], dupdt, novdt)
    return(frangedt)

}

#' @export
txAnnotate <- function(.seqnames, .start, .end, codingStrand, genome) {

    rangeSites <- mapply(':', .start, .end, SIMPLIFY = FALSE)
    siteIdxs <- unlist(rangeSites, recursive = FALSE, use.names = FALSE)

    n <- sapply(rangeSites, length)
    siteRanges <- GenomicRanges::GRanges(rep(.seqnames, n), IRanges::IRanges(siteIdxs, siteIdxs))
    txdt <- data.table::data.table(
        seqnames = rep(.seqnames, n),
        start = siteIdxs,
        end = siteIdxs,
        ref = as.character(genome[siteRanges], use.names = FALSE),
        codingStrand = rep(codingStrand, n)
    )

    txdt[, ':=' ("pyriStrand" = "+", "txPyri" = "T")]
    txdt[ref %in% c("G", "A"), "pyriStrand" := "-"]
    txdt[pyriStrand == codingStrand, "txPyri" := "U"]

    return(txdt)

}

