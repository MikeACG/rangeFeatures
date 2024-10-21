#' @import data.table
#' @import Biostrings
#' @importFrom dplyr %>%

# https://www.biostars.org/p/489350/
withFeatureSdiff <- function(gr1, gr2) {

    # all ranges in .D will overlap exactly once with one and only one range in gr1 if gr1 has non-overlapping ranges
    adgr <- GenomicRanges::disjoin(c(gr1, gr2))
    .D <- adgr[!IRanges::overlapsAny(adgr, gr2)]

    dov <- GenomicRanges::findOverlaps(.D, gr1)
    .D <- .D[S4Vectors::queryHits(dov)]
    .D$feature <- gr1$feature[S4Vectors::subjectHits(dov)]
    .D$n <- gr1$n[S4Vectors::subjectHits(dov)]

    return(.D)

}

# does this produce the same output as the thing above?
# withFeatureSdiff <- function(gr1, gr2) {

#     # compute difference between ranges
#     .D <- setdiff(gr1, gr2)

#     # intersections of gr1 with .D should be the parts of the ranges in gr1 that don't intersect with gr2
#     # hence we can directly get their features from gr1 
#     .ov <- GenomicRanges::findOverlaps(gr1, .D)
#     dr <- GenomicRanges::pintersect(gr1[S4Vectors::queryHits(.ov)], .D[S4Vectors::subjectHits(.ov)])
#     dr$feature <- gr1$feature[S4Vectors::queryHits(.ov)]
#     dr$n <- gr1$n[S4Vectors::subjectHits(.ov)]

#     return(dr)

# }

#' @export
sumFeature <- function(gr) {

    dgr <- GenomicRanges::disjoin(gr, with.revmap = TRUE)

    dgr$n <- sapply(dgr$revmap, length)
    sdt <- data.table::data.table(
        disIdx = rep(1:length(dgr), dgr$n),
        fidx = unlist(dgr$revmap, recursive = FALSE, use.names = FALSE)
    )
    sdt[, "feature" := gr$feature[fidx]]

    dgr$featureSum <- sdt[, list("featureSum" = sum(feature)), by = "disIdx"]$featureSum

    return(dgr)

}

# # assumes disjoint ranges within each of accgr and newgr 
# #' @export
# sumFeature <- function(ov, accgr, newgr) {

#     # intersections between overlapping ranges and their sum of feature signal
#     .I <- GenomicRanges::pintersect(accgr[S4Vectors::queryHits(ov)], newgr[S4Vectors::subjectHits(ov)])
#     .I$feature <- accgr$feature[S4Vectors::queryHits(ov)] + newgr$feature[S4Vectors::subjectHits(ov)]
#     .I$n <- accgr$n[S4Vectors::queryHits(ov)] + 1L # 

#     # get parts that don't intersect from both overlapping & non-overlapping ranges 
#     d1 <- withFeatureSdiff(accgr, .I)
#     d2 <- withFeatureSdiff(newgr, .I)

#     .I$hit <- NULL
#     return(c(.I, d1, d2))

# }

#' @export
frangesLoad <- function(.seq, .dir) {

    .files <- list.files(.dir)
    chrFile <- .files[grep(paste0("\\b", .seq, "\\b"), .files)]
    rangesdt <- fread(paste0(.dir, chrFile), nThread = 1)[, .SD, .SDcols = 2:4]

    ranges <- GRanges(
        rep(.seq, nrow(rangesdt)),
        IRanges(rangesdt[[1]], rangesdt[[2]]),
        feature = rangesdt[[3]]
    )

    return(ranges)

}

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

    frangedt <- rbind(
        ovdt[is.na(duplicate), .SD, .SDcols = -c("duplicate")],
        dupdt,
        novdt[, .SD, .SDcols = -c("iquery")]
    )
    return(frangedt)

}

#' @export
pyriAnnotate <- function(.seqnames, .start, .end, featureStrand, genome) {

    rangeSites <- mapply(':', .start, .end, SIMPLIFY = FALSE)
    siteIdxs <- unlist(rangeSites, recursive = FALSE, use.names = FALSE)

    n <- sapply(rangeSites, length)
    siteRanges <- GenomicRanges::GRanges(rep(.seqnames, n), IRanges::IRanges(siteIdxs, siteIdxs))
    sitedt <- data.table::data.table(
        seqnames = rep(.seqnames, n),
        start = siteIdxs,
        end = siteIdxs,
        ref = as.character(genome[siteRanges], use.names = FALSE),
        featureStrand = rep(featureStrand, n)
    )

    sitedt[, "pyriStrand" := "+"]
    sitedt[ref %in% c("G", "A"), "pyriStrand" := "-"]

    return(sitedt)

}

txAnnotate <- function(txdt) {

    txdt[, ':=' ("pyriStrand" = "+", "txPyri" = "T")]
    txdt[ref %in% c("G", "A"), "pyriStrand" := "-"]
    txdt[pyriStrand == codingStrand, "txPyri" := "U"]

    return()

}



