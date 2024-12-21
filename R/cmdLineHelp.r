#' @export
parseColnames <- function(.string, pcolnames) {

    .colnames <- unlist(strsplit(.string, ","))
    r <- .colnames
    if (.colnames[1] == "none") r <- NULL
    if (.colnames[1] == "all") r <- setdiff(pcolnames, c("seqname", "start", "end")) 

    return(r)

}
