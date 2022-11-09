##################################
######## SummarizedMotifs ######
##################################

methods::setClass("SummarizedMotifs", contains = c("SummarizedExperiment"))

SummarizedMotifs <- function(granges, seqs, taxon = "4932", db = JASPAR2020::JASPAR2020, ncores = 1) {

    # prepare colData
    motifs <- TFBSTools::getMatrixSet(db, opts = list(species = taxon))[c(1:3, 76)]
    colData <- motifs2colData(motifs)

    # check motifs
    s <- scanGRangesForJASPARMotifs(granges, seqs, taxon, db, ncores)
    mat <- as.matrix(mcols(s)[, colData$motif])
    colnames(mat) <- NULL

    # Create an `AggregatedCoverage` object
    SM <- methods::new(
        "SummarizedMotifs",
        SummarizedExperiment::SummarizedExperiment(
            rowRanges = granges,
            colData = colData,
            assays = list(
                TFBS = mat
            )
        )
    )
    return(SM)
}

##################################
######## Odds ######
##################################

methods::setClass("Odds", contains = c("SummarizedExperiment"))
