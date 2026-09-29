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

test_that("corner legends are outside and align to their selected edges", {
  expect_equal(legend_position_for_theme("left_top"), "left")
  expect_equal(legend_position_for_theme("right_top"), "right")
  expect_equal(legend_position_for_theme("left_bottom"), "left")
  expect_equal(legend_position_for_theme("right_bottom"), "right")
  expect_null(legend_position_inside_for_theme("right_top"))
  expect_equal(legend_justification_for_theme("left_top"), c(1, 1))
  expect_equal(legend_justification_for_theme("right_top"), c(0, 1))
  expect_equal(legend_justification_for_theme("left_bottom"), c(1, 0))
  expect_equal(legend_justification_for_theme("right_bottom"), c(0, 0))
  expect_equal(legend_box_justification_for_theme("right_top"), "left")
  expect_equal(legend_box_justification_for_theme("left"), "left")
  expect_equal(legend_box_justification_for_theme("right"), "left")
  expect_equal(legend_justification_for_theme("left"), c(1, 0.5))
  expect_equal(legend_justification_for_theme("right"), c(0, 0.5))
  expect_equal(legend_justification_for_theme("top"), c(0.5, 0))
  expect_equal(legend_justification_for_theme("bottom"), c(0.5, 1))
  expect_equal(validate_legend_plot_spacing(position = "right_top"), 3)
})

test_that("gene palette and manual fill resolve by category", {
  resolved <- ggplasmid_resolve_colors(
    scheme = "phage",
    gene_palette = "npg",
    gene_manual_fill = c(Tail = "#123456")
  )

  expect_equal(unname(resolved[["Tail"]]), "#123456")
  expect_true("Lysis" %in% names(resolved))
  expect_equal(unname(resolved[["GC skew+"]]), "#008000")
})

test_that("automatic BLAST rings stack outside GC and before genes", {
  geometry <- blast_ring_geometry(
    n_rings = 3,
    gene_height = 0.085,
    gc_skew_radius = 0.78,
    gc_skew_height = 0.10,
    gc_content_radius = 0.67,
    gc_content_height = 0.045,
    show_gc_skew = TRUE
  )

  expect_gte(geometry$inner_edge, 0.83)
  expect_gt(geometry$gene_radius - 0.085 / 2, geometry$outer_edge)
  expect_gt(geometry$radius, geometry$inner_edge)
  expect_equal(geometry$spacing, 0.01)

  no_blast <- blast_ring_geometry(
    n_rings = 0,
    gene_height = 0.085,
    gc_skew_radius = 0.78,
    gc_skew_height = 0.10,
    gc_content_radius = 0.67,
    gc_content_height = 0.045,
    show_gc_skew = TRUE
  )
  expect_null(no_blast$radius)
  expect_equal(no_blast$gene_radius - 0.085 / 2, 0.84)
})

test_that("BLAST rings use independent gradient fill scales", {
  features <- data.frame(
    start = c(10, 150),
    end = c(90, 240),
    strand = c("+", "-"),
    product = c("replication protein", "tail protein"),
    category = c("Other functions", "Other functions")
  )
  hits <- data.frame(
    qseqid = "query",
    sseqid = "subject",
    pident = c(96, 88),
    length = c(50, 60),
    qstart = c(20, 220),
    qend = c(69, 279),
    sstart = c(1, 70),
    send = c(50, 129),
    evalue = c(1e-20, 1e-12),
    bitscore = c(100, 80)
  )
  plot <- ggplasmid(
    annotation = features,
    genome_length = 500,
    name = "synthetic genome",
    blast_rings = list(
      first = list(label = "first comparison", hits = hits),
      second = list(label = "second comparison", hits = hits, min_identity = 85)
    ),
    show_gc_skew = FALSE,
    show_labels = FALSE
  )

  expect_no_error(ggplot2::ggplotGrob(plot))
  expect_no_error(ggplot2::ggplotGrob(
    plot + gene_highlight(`Other functions` = "#FF0000")
  ))
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

test_that("pairwise BLAST input accepts in-memory sequences", {
  query <- c(query_gene = "ACGTN", query_other = "TTTT")
  subject <- data.frame(
    name = "reference",
    sequence = "ACGTNACGT",
    stringsAsFactors = FALSE
  )

  work_dir <- tempfile("blast-input-test-")
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

  query_file <- blast_pairwise_sequence_input(query, "query", work_dir)
  subject_file <- blast_pairwise_sequence_input(subject, "subject", work_dir)

  expect_true(file.exists(query_file$path))
  expect_true(file.exists(subject_file$path))
  expect_equal(
    readLines(query_file$path),
    c(">query_gene", "ACGTN", ">query_other", "TTTT")
  )
  expect_equal(readLines(subject_file$path), c(">reference", "ACGTNACGT"))
})

test_that("pairwise BLAST input accepts Biostrings DNAStringSet", {
  skip_if_not_installed("Biostrings")
  sequences <- Biostrings::DNAStringSet(
    c(reference_a = "ACGTN", reference_b = "TTTT")
  )
  work_dir <- tempfile("blast-xstring-test-")
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

  input_file <- blast_pairwise_sequence_input(sequences, "query", work_dir)

  expect_equal(
    readLines(input_file$path),
    c(">reference_a", "ACGTN", ">reference_b", "TTTT")
  )
})

test_that("pairwise BLAST output is normalized for circular plotting", {
  work_dir <- tempfile("blast-output-test-")
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)
  output <- file.path(work_dir, "hits.tsv")
  writeLines(
    c(
      "query_1\tsubject_1\t99.5\t100\t0\t0\t20\t119\t900\t801\t1e-40\t150",
      "query_1\tsubject_1\t97\t50\t1\t0\t200\t249\t100\t149\t2e-12\t80"
    ),
    output
  )

  hits <- blast_pairwise_read_output(output)
  expect_equal(nrow(hits), 2)
  expect_equal(hits$query_start, c(20, 200))
  expect_equal(hits$query_end, c(119, 249))
  expect_equal(hits$subject_start, c(801, 100))
  expect_equal(hits$subject_end, c(900, 149))
  expect_equal(hits$subject_strand, c("-", "+"))
})

test_that("pairwise BLAST output format is passed as one argument", {
  quoted <- blast_pairwise_quote_arg("6 qseqid sseqid pident")
  expect_true(grepl("qseqid", quoted, fixed = TRUE))
  expect_true(grepl("pident", quoted, fixed = TRUE))
  expect_gt(nchar(quoted), nchar("6 qseqid sseqid pident"))
})

test_that("BLAST plotting modes prepare reference segments", {
  hits <- data.frame(
    qseqid = "query_1",
    sseqid = "subject_1",
    pident = c(95, 85, 99),
    length = c(100, 100, 50),
    qstart = c(100, 150, 400),
    qend = c(199, 249, 449),
    sstart = c(1, 1, 1),
    send = c(100, 100, 50),
    bitscore = c(100, 80, 120),
    stringsAsFactors = FALSE
  )

  hsp <- blast_prepare_plot_data(
    hits,
    mode = "hsp",
    reference = "query",
    genome_length = 500,
    min_identity = 90
  )
  expect_equal(nrow(hsp$data), 2)
  expect_true(all(hsp$data$fill_key %in% hsp$breaks))

  region <- blast_prepare_plot_data(
    hits,
    mode = "region",
    reference = "query",
    genome_length = 500,
    min_identity = 80
  )
  expect_equal(nrow(region$data), 4)
  expect_true(all(region$data$start <= region$data$end))
})

test_that("BLAST gene mode colors sufficiently covered features", {
  hits <- data.frame(
    pident = c(92, 75),
    length = c(90, 100),
    query_start = c(101, 300),
    query_end = c(190, 399),
    bitscore = c(100, 90),
    stringsAsFactors = FALSE
  )
  features <- data.frame(
    feature_id = c("gene_a", "gene_b"),
    start = c(100, 300),
    end = c(199, 399),
    wraps_origin = FALSE,
    stringsAsFactors = FALSE
  )

  gene <- blast_prepare_plot_data(
    hits,
    mode = "gene",
    reference = "query",
    genome_length = 500,
    features = features,
    min_identity = 80,
    min_gene_coverage = 0.8
  )
  expect_equal(nrow(gene$data), 1)
  expect_equal(gene$data$feature_id, "gene_a")
  expect_equal(gene$data$start, 100)
  expect_equal(gene$data$end, 199)
})

test_that("BLAST gene mode ranks multiple candidate hits without length mismatch", {
  hits <- data.frame(
    qseqid = "query",
    sseqid = "subject",
    pident = c(92, 92),
    length = c(51, 100),
    qstart = c(100, 100),
    qend = c(150, 199),
    sstart = c(1, 52),
    send = c(51, 151),
    evalue = c(1e-10, 1e-20),
    bitscore = c(80, 80)
  )
  features <- data.frame(
    feature_id = "gene_a",
    start = 100,
    end = 199,
    wraps_origin = FALSE
  )

  gene <- blast_prepare_plot_data(
    hits,
    mode = "gene",
    reference = "query",
    genome_length = 500,
    features = features,
    min_identity = 80,
    min_gene_coverage = 0.80
  )

  expect_equal(nrow(gene$data), 1)
  expect_equal(gene$data$length, 100)
})
