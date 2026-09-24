test_that("color schemes return named colors", {
  plasmid <- ggplasmid_colors("plasmid")
  phage <- ggplasmid_colors("phage")

  expect_true(length(plasmid) > 0)
  expect_true(length(phage) > 0)
  expect_named(plasmid)
  expect_named(phage)
  expect_true(all(vapply(plasmid, is.character, logical(1))))
  expect_true(all(vapply(phage, is.character, logical(1))))
})

test_that("GC skew calculation returns the documented columns", {
  result <- compute_gc_skew(
    sequence = "GCGCGCGC",
    genome_length = 8,
    window = 4,
    step = 2,
    circular = FALSE
  )

  expect_s3_class(result, "data.frame")
  expect_named(result, c("position", "gc_skew", "gc_content", "window", "step"))
  expect_equal(result$window, rep(4L, nrow(result)))
  expect_equal(result$step, rep(2L, nrow(result)))
})

test_that("packaged example files can be read", {
  gbk <- system.file("extdata", "pSGNDM_5.gbk", package = "ggplasmidZY")
  fasta <- system.file("extdata", "pSGNDM_5.fasta", package = "ggplasmidZY")

  expect_true(file.exists(gbk))
  expect_true(file.exists(fasta))
  expect_gt(nrow(read_gbk(gbk)), 0)
  expect_true(nzchar(read_plasmid_fasta(fasta)$sequence[[1]]))
})
