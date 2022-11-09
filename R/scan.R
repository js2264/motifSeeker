#' scanGRangesForJASPARMotifs
#'
#' Scan a set of GRanges for JASPAR motifs
#'
#' @param granges
#' @param seqs
#' @param taxon
#' @param db
#' @param ncores
#' @return GRanges input with extra columns. Each column indicates whether a TF
#' binding site is found in the GRanges object
#'
#' @importFrom JASPAR2020 JASPAR2020
#' @import GenomicRanges
#' @import TFBSTools
#' @import BiocParallel
#' @export
#'
#' @examples
#' ATAC_peaks <- readRDS(here::here('../../../Share/Day1/ATAC_peaks.rds'))
#' granges <- ATAC_peaks[ATAC_peaks$annot != "non-DA"]
#' seqlevelsStyle(granges) <- "UCSC"
#' seqs <- getSeq(BSgenome.Scerevisiae.UCSC.sacCer3::BSgenome.Scerevisiae.UCSC.sacCer3)
#' x <- scanGRangesForJASPARMotifs(granges, seqs)

scanGRangesForJASPARMotifs <- function(
        granges,
        seqs,
        taxon = "4932",
        db = JASPAR2020::JASPAR2020,
        ncores = 1
) {

    motifs <- TFBSTools::getMatrixSet(db, opts = list(species = taxon))[c(1:3, 76)]
    motif_names <- TFBSTools::name(motifs)

    motif_overalps <- BiocParallel::bplapply(
        BPPARAM = BiocParallel::MulticoreParam(workers = ncores),
        motif_names,
        function(motif_name) {
            message(motif_name)
            motif <- TFBSTools::toPWM(motifs[[which(motif_names == motif_name)]])
            hits <- GenomicRanges::GRanges(
                TFBSTools::searchSeq(
                    motif,
                    seqs,
                    strand = '*',
                    min.score = '60%'
                )
            )
            hits <- hits[hits$relScore >= 0.90]
            if (length(hits) > 1) {
                res <- data.frame(motif_name = IRanges::overlapsAny(granges, hits))
            }
            else {
                res <- data.frame(motif_name = rep(FALSE, length(granges)))
            }
            colnames(res) <- motif_name
            res
        }
    )
    motif_overalps <- do.call(cbind, motif_overalps)
    GenomicRanges::mcols(granges) <- cbind(GenomicRanges::mcols(granges), motif_overalps)

    return(granges)
}
