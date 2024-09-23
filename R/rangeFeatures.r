#' @import data.table
#' @import Biostrings
#' @importFrom dplyr %>%

#' @export
rangeFeatures <- function(ov, queryRanges, subjectRanges, featureType, allowMiss = FALSE) {

    frangedt <- switch(
        featureType,
        categorical = rangeCatFeature(ov, queryRanges, subjectRanges, allowMiss),
        numerical = rangeNumFeature(ov, queryRanges, subjectRanges, allowMiss),
        stop("Invalid 'featureType'")
    )

    return(frangedt)

}

getNoOverlapQuery <- function(ov, queryRanges) {

    inov <- setdiff(1:length(queryRanges), S4Vectors::queryHits(ov))
    novdt <- data.table::data.table(
        seqnames = as.character(GenomicRanges::seqnames(queryRanges)[inov]),
        start = GenomicRanges::start(queryRanges)[inov],
        end = GenomicRanges::end(queryRanges)[inov],
        iquery = inov
    )

    return(novdt)

}

rangeNumFeature <- function(ov, queryRanges, subjectRanges, allowMiss) {

    # get intersections of overlaps to know the widths
    intRanges <- GenomicRanges::pintersect(
        queryRanges[S4Vectors::queryHits(ov)],
        subjectRanges[S4Vectors::subjectHits(ov)]
    )

    # gather overlap data
    ovdt <- data.table::data.table(
        iquery = S4Vectors::queryHits(ov),
        ovWidth = GenomicRanges::width(intRanges),
        feature = subjectRanges$feature[S4Vectors::subjectHits(ov)]
    )
    
    # aggregate signals by query range weighted by overlap size
    qovdt <- ovdt[, list("wfeature" = sum(feature * ovWidth), "nsitesQuery" = sum(ovWidth)), by = "iquery"]

    # get weighted average of signal per query range
    if(!allowMiss) qovdt[, "nsitesQuery" := GenomicRanges::width(queryRanges)[iquery]]
    qovdt[, "feature" := wfeature / nsitesQuery]

    # format data correctly
    qovdt[, "wfeature" := NULL]
    qovdt[
        ,
        ':=' (
            "seqnames" = as.character(GenomicRanges::seqnames(queryRanges)[iquery]),
            "start" = GenomicRanges::start(queryRanges)[iquery],
            "end" = GenomicRanges::end(queryRanges)[iquery]
        )
    ]

    # handle query ranges with no overlap
    novdt <- getNoOverlapQuery(ov, queryRanges)
    novdt[, ':=' ("feature" = NA_real_, "nsitesQuery" = 0L)]
    if (!allowMiss) novdt[, ':=' ("feature" = 0L, "nsitesQuery" = GenomicRanges::width(queryRanges)[iquery])]

    frangedt <- rbind(qovdt[, .SD, .SDcols = names(novdt)], novdt)
    frangedt[, "iquery" := NULL]
    return(frangedt)

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

