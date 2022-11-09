## SummarizedMotifs-specific Odds method

methods::setMethod("Odds", c("SummarizedMotifs", "factor"), function(x, group_by) {

    mat <- SummarizedExperiment::assay(x, 'TFBS')
    motifs <- colData(x)$motif
    n_gr <- nrow(x)
    groups <- levels(group_by)
    colData(x)$n <- 0
    rowData <- data.frame(group = groups)
    rownames(rowData) <- groups
    rowData$n <- table(group_by)

    res <- lapply(motifs, function(motif) {
        message(motif)
        m <- which(colData(x)$motif == motif)
        tests <- lapply(groups, function(group) {
            tab <- data.frame(
                in_group = group_by == group,
                has_motif = mat[, m]
            ) |>
                table()
            if (!"TRUE" %in% rownames(tab) | !"FALSE" %in% rownames(tab)) {
                tab <- rbind(tab, c(0, 0))
            }
            fisher.test(tab)
        })
        odds <- purrr::map(tests, 'estimate') |> unlist()
        pvals <- purrr::map(tests, 'p.value') |> unlist()
        list(
            n_total = sum(mat[, m]),
            n = purrr::map(groups, ~sum(mat[group_by == .x, m])),
            odds = odds,
            pvals = pvals
        )
    })
    colData(x)$n <- purrr::map(res, 'n_total') |> unlist()
    n <- purrr::map(res, 'n') %>% do.call(cbind, .)
    rownames(n) <- groups
    colnames(n) <- motifs
    odds <- purrr::map(res, 'odds') %>% do.call(cbind, .)
    rownames(odds) <- groups
    colnames(odds) <- motifs
    pvals <- purrr::map(res, 'pvals') %>% do.call(cbind, .)
    rownames(pvals) <- groups
    colnames(pvals) <- motifs

    # Create an `Odds` object
    Odds <- methods::new(
        "Odds",
        SummarizedExperiment::SummarizedExperiment(
            rowData = rowData,
            colData = colData(x),
            assays = list(
                hits = n,
                odds = odds,
                p.value = pvals
            )
        )
    )
    return(Odds)

})
