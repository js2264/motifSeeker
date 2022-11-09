motifs2colData <- function(motifs) {
    colData <- data.frame(
        motif = TFBSTools::name(motifs),
        ID = TFBSTools::ID(motifs),
        org = TFBSTools::tags(motifs)[[1]]$species,
        taxon = names(TFBSTools::tags(motifs)[[1]]$species),
        type = unlist(purrr::map(TFBSTools::tags(motifs), 'type'))
    )
    rownames(colData) <- colData$motif
    colData
}
