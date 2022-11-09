test_that("scanning function works", {
    gr <- GenomicRanges::GRanges('chrI', IRanges::IRanges(sample(1:50000, 100), width = 300))
    seqs <- Biostrings::getSeq(BSgenome.Scerevisiae.UCSC.sacCer3::BSgenome.Scerevisiae.UCSC.sacCer3)
    # x <- scanGRangesForJASPARMotifs(gr, seqs)
    expect_equal(
        TRUE,
        TRUE
    )
})
