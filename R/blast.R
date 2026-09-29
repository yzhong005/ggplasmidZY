# Standard BLAST tabular fields needed by the pairwise wrapper and the map
# overlay planned for a later stage.
blast_pairwise_fields <- c(
  "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
  "qstart", "qend", "sstart", "send", "evalue", "bitscore"
)

blast_pairwise_empty_hits <- function() {
  out <- data.frame(
    qseqid = character(),
    sseqid = character(),
    pident = numeric(),
    length = numeric(),
    mismatch = numeric(),
    gapopen = numeric(),
    qstart = numeric(),
    qend = numeric(),
    sstart = numeric(),
    send = numeric(),
    evalue = numeric(),
    bitscore = numeric(),
    stringsAsFactors = FALSE
  )
  out$query_start <- numeric()
  out$query_end <- numeric()
  out$subject_start <- numeric()
  out$subject_end <- numeric()
  out$subject_strand <- character()
  out
}

blast_pairwise_resolve_executable <- function(blastn = NULL) {
  if (is.null(blastn)) {
    blastn <- Sys.which("blastn")
  }
  if (!is.character(blastn) || length(blastn) != 1L || !nzchar(blastn)) {
    stop(
      "Could not find `blastn`. Install NCBI BLAST+ or supply its path with `blastn`.",
      call. = FALSE
    )
  }
  if (file.exists(blastn)) {
    return(normalizePath(blastn, mustWork = TRUE))
  }
  resolved <- Sys.which(blastn)
  if (!nzchar(resolved)) {
    stop("The supplied `blastn` executable was not found: ", blastn, call. = FALSE)
  }
  normalizePath(resolved, mustWork = TRUE)
}

blast_pairwise_validate_sequence <- function(sequence, label, index) {
  sequence <- toupper(gsub("\\s+", "", as.character(sequence)))
  sequence <- chartr("U", "T", sequence)
  if (!nzchar(sequence)) {
    stop(label, " sequence ", index, " is empty.", call. = FALSE)
  }
  if (!grepl("^[ACGTRYSWKMBDHVN.-]+$", sequence)) {
    stop(
      label, " sequence ", index,
      " contains characters that are not valid nucleotide symbols.",
      call. = FALSE
    )
  }
  sequence
}

blast_pairwise_sequence_input <- function(input, label, work_dir) {
  if (is.character(input) && length(input) == 1L && file.exists(input)) {
    return(list(
      path = normalizePath(input, mustWork = TRUE),
      temporary = FALSE
    ))
  }

  if (inherits(input, "XStringSet")) {
    sequences <- as.character(input)
    ids <- names(input)
  } else if (is.data.frame(input)) {
    fasta <- read_plasmid_fasta(input)
    ids <- clean_text(fasta$name)
    sequences <- fasta$sequence
  } else if (is.character(input)) {
    sequences <- vapply(
      seq_along(input),
      function(i) blast_pairwise_validate_sequence(input[[i]], label, i),
      character(1L)
    )
    ids <- names(input)
    if (is.null(ids)) {
      ids <- paste0(tolower(label), "_", seq_along(sequences))
    }
    ids <- clean_text(ids)
    ids[!nzchar(ids)] <- paste0(tolower(label), "_", which(!nzchar(ids)))
  } else {
    stop(
      "`", label, "` must be a FASTA path, a Biostrings XStringSet, a data " ,
      "frame with a `sequence` column, or a character vector of nucleotide " ,
      "sequences.",
      call. = FALSE
    )
  }

  if (!is.null(ids)) {
    ids <- clean_text(ids)
  }
  ids <- gsub("[^A-Za-z0-9_.|:-]+", "_", ids)
  ids[!nzchar(ids)] <- paste0(tolower(label), "_", which(!nzchar(ids)))
  ids <- make.unique(ids, sep = "_")

  sequences <- vapply(
    seq_along(sequences),
    function(i) blast_pairwise_validate_sequence(sequences[[i]], label, i),
    character(1L)
  )
  if (length(ids) != length(sequences)) {
    ids <- paste0(tolower(label), "_", seq_along(sequences))
  }

  fasta_path <- file.path(work_dir, paste0(tolower(label), ".fasta"))
  fasta_lines <- unlist(
    lapply(seq_along(sequences), function(i) c(
      paste0(">", ids[[i]]),
      sequences[[i]]
    )),
    use.names = FALSE
  )
  writeLines(fasta_lines, fasta_path, useBytes = TRUE)
  list(path = fasta_path, temporary = TRUE)
}

blast_pairwise_numeric_arg <- function(args, flag, value, minimum = 0,
                                       maximum = Inf) {
  if (is.null(value)) {
    return(args)
  }
  value <- as.numeric(value[[1L]])
  if (!is.finite(value) || value < minimum || value > maximum) {
    stop("`", sub("^-", "", flag), "` must be a valid number.", call. = FALSE)
  }
  c(args, flag, format(value, scientific = FALSE, trim = TRUE))
}

blast_pairwise_integer_arg <- function(args, flag, value, minimum = 1L) {
  if (is.null(value)) {
    return(args)
  }
  value <- as.numeric(value[[1L]])
  if (!is.finite(value) || value != floor(value) || value < minimum) {
    stop("`", sub("^-", "", flag), "` must be a positive integer.", call. = FALSE)
  }
  c(args, flag, as.character(as.integer(value)))
}

blast_pairwise_quote_arg <- function(value) {
  quote_type <- if (.Platform$OS.type == "windows") "cmd" else "sh"
  shQuote(value, type = quote_type)
}

blast_pairwise_read_output <- function(path) {
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size == 0) {
    return(blast_pairwise_empty_hits())
  }

  out <- tryCatch(
    utils::read.delim(
      path,
      header = FALSE,
      sep = "\t",
      quote = "",
      comment.char = "",
      check.names = FALSE,
      stringsAsFactors = FALSE,
      col.names = blast_pairwise_fields
    ),
    error = function(error) {
      stop("Could not parse BLAST output: ", conditionMessage(error), call. = FALSE)
    }
  )
  if (!nrow(out)) {
    return(blast_pairwise_empty_hits())
  }

  numeric_columns <- setdiff(blast_pairwise_fields, c("qseqid", "sseqid"))
  for (column in numeric_columns) {
    out[[column]] <- suppressWarnings(as.numeric(out[[column]]))
  }
  out$query_start <- pmin(out$qstart, out$qend)
  out$query_end <- pmax(out$qstart, out$qend)
  out$subject_start <- pmin(out$sstart, out$send)
  out$subject_end <- pmax(out$sstart, out$send)
  out$subject_strand <- ifelse(out$sstart <= out$send, "+", "-")
  out
}

blast_plot_empty <- function() {
  data.frame(
    start = numeric(),
    end = numeric(),
    pident = numeric(),
    length = numeric(),
    bitscore = numeric(),
    fill_key = character(),
    feature_id = character(),
    ring_key = character(),
    ring_label = character(),
    ring_order = numeric(),
    radius = numeric(),
    stringsAsFactors = FALSE
  )
}

blast_identity_breaks <- function() {
  c(
    "BLAST <50%", "BLAST 50-70%", "BLAST 70-80%",
    "BLAST 80-90%", "BLAST 90-95%", "BLAST 95-100%"
  )
}

blast_identity_fill_key <- function(pident) {
  cut(
    as.numeric(pident),
    breaks = c(-Inf, 50, 70, 80, 90, 95, Inf),
    labels = blast_identity_breaks(),
    right = FALSE,
    include.lowest = TRUE
  )
}

blast_identity_colors <- function(colour = "#2C7FB8") {
  colour <- as.character(colour[[1L]])
  if (!nzchar(colour)) {
    stop("`blast_colour` must be a non-empty colour value.", call. = FALSE)
  }
  colors <- grDevices::colorRampPalette(c("#EEEEEE", colour))(6)
  stats::setNames(colors, blast_identity_breaks())
}

blast_ring_geometry <- function(n_rings, radius = NULL, height = NULL,
                                spacing = NULL, gene_radius = NULL,
                                gene_height = 0.10,
                                gc_skew_radius = 0.78,
                                gc_skew_height = 0.10,
                                gc_content_radius = 0.67,
                                gc_content_height = 0.045,
                                show_gc_skew = TRUE,
                                ruler_radius = NULL,
                                ruler_major_tick = 0.035,
                                clearance = 0.01) {
  n_rings <- as.integer(n_rings[[1L]])
  if (is.na(n_rings) || n_rings < 0L) {
    stop("`n_rings` must be a non-negative integer.", call. = FALSE)
  }

  if (is.null(spacing)) spacing <- 0.01
  spacing <- as.numeric(spacing[[1L]])
  if (!is.finite(spacing) || spacing < 0) {
    stop("`blast_spacing` must be a non-negative number.", call. = FALSE)
  }

  ruler_radius <- ruler_radius %||% (gc_content_radius - 0.08)
  lower_tracks <- c(ruler_radius + ruler_major_tick)
  if (isTRUE(show_gc_skew)) {
    lower_tracks <- c(
      lower_tracks,
      gc_skew_radius + gc_skew_height / 2,
      gc_content_radius + gc_content_height
    )
  }
  lower_edge <- max(lower_tracks) + clearance

  if (n_rings == 0L) {
    if (!is.null(gene_radius)) {
      gene_radius <- as.numeric(gene_radius[[1L]])
      if (!is.finite(gene_radius) || gene_radius <= 0) {
        stop("`gene_radius` must be a positive number.", call. = FALSE)
      }
    } else {
      gene_radius <- lower_edge + gene_height / 2
    }
    return(list(
      radius = NULL, height = NULL, spacing = spacing,
      inner_edge = lower_edge, outer_edge = lower_edge,
      gene_radius = gene_radius
    ))
  }

  if (!is.null(gene_radius)) {
    gene_radius <- as.numeric(gene_radius[[1L]])
    if (!is.finite(gene_radius) || gene_radius <= 0) {
      stop("`gene_radius` must be a positive number.", call. = FALSE)
    }
    upper_edge <- gene_radius - gene_height / 2 - clearance
    available <- upper_edge - lower_edge - (n_rings - 1L) * spacing
  } else {
    available <- Inf
  }

  if (is.null(height)) {
    height <- if (is.finite(available)) min(0.08, available / n_rings) else 0.08
  }
  height <- as.numeric(height[[1L]])
  if (!is.finite(height) || height <= 0) {
    stop("`blast_height` must be a positive number.", call. = FALSE)
  }
  if (is.finite(available) &&
      n_rings * height > available + sqrt(.Machine$double.eps)) {
    max_height <- available / n_rings
    stop(
      "BLAST rings do not fit between the GC tracks and gene ring with the requested `blast_height`/`blast_spacing`. ",
      "For ", n_rings, " rings, `blast_height` must be at most ",
      format(max_height, digits = 3),
      " with this spacing. Set `blast_height = NULL` for automatic sizing, or reduce `blast_spacing`.",
      call. = FALSE
    )
  }

  if (is.null(radius)) {
    radius <- lower_edge + height / 2
  }
  radius <- as.numeric(radius[[1L]])
  if (!is.finite(radius) || radius <= 0) {
    stop("`blast_radius` must be a positive number.", call. = FALSE)
  }
  final_outer_edge <- radius + (n_rings - 1L) * (height + spacing) + height / 2
  if (radius - height / 2 < lower_edge - 1e-8 ||
      (is.finite(available) && final_outer_edge > upper_edge + 1e-8)) {
    stop(
      "The requested BLAST ring placement overlaps a GC track or the gene ring. ",
      "Set `blast_radius = NULL` for automatic placement.",
      call. = FALSE
    )
  }

  if (is.null(gene_radius)) {
    gene_radius <- final_outer_edge + clearance + gene_height / 2
  }

  list(radius = radius, height = height, spacing = spacing,
       inner_edge = radius - height / 2,
       outer_edge = final_outer_edge,
       gene_radius = gene_radius)
}

blast_normalize_plot_hits <- function(hits, reference = c("query", "subject"),
                                      min_identity = 0,
                                      min_alignment_length = 1,
                                      coordinate_offset = 0) {
  reference <- match.arg(reference)
  if (!is.data.frame(hits)) {
    stop("`blast_hits` must be a data frame returned by `blast_pairwise()`.",
         call. = FALSE)
  }
  required <- c("pident", "length")
  missing <- setdiff(required, names(hits))
  if (length(missing)) {
    stop(
      "`blast_hits` is missing required column(s): ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }
  coordinate_names <- if (reference == "query") {
    c("query_start", "query_end")
  } else {
    c("subject_start", "subject_end")
  }
  fallback_names <- if (reference == "query") c("qstart", "qend") else c("sstart", "send")
  if (!all(coordinate_names %in% names(hits))) {
    if (!all(fallback_names %in% names(hits))) {
      stop(
        "`blast_hits` must contain normalized coordinates `",
        paste(coordinate_names, collapse = "`, `"), "` (or BLAST columns `",
        paste(fallback_names, collapse = "`, `"), "`).", call. = FALSE
      )
    }
    hits[[coordinate_names[[1L]]]] <- pmin(hits[[fallback_names[[1L]]]], hits[[fallback_names[[2L]]]])
    hits[[coordinate_names[[2L]]]] <- pmax(hits[[fallback_names[[1L]]]], hits[[fallback_names[[2L]]]])
  }

  out <- data.frame(
    start = suppressWarnings(as.numeric(hits[[coordinate_names[[1L]]]])) - coordinate_offset,
    end = suppressWarnings(as.numeric(hits[[coordinate_names[[2L]]]])) - coordinate_offset,
    pident = suppressWarnings(as.numeric(hits$pident)),
    length = suppressWarnings(as.numeric(hits$length)),
    bitscore = if ("bitscore" %in% names(hits)) {
      suppressWarnings(as.numeric(hits$bitscore))
    } else {
      0
    },
    stringsAsFactors = FALSE
  )
  keep <- is.finite(out$start) & is.finite(out$end) &
    is.finite(out$pident) & is.finite(out$length) &
    out$end >= out$start & out$pident >= min_identity &
    out$length >= min_alignment_length
  out <- out[keep, , drop = FALSE]
  if (!nrow(out)) {
    return(out)
  }
  starts <- out$start
  ends <- out$end
  out$start <- pmin(starts, ends)
  out$end <- pmax(starts, ends)
  out
}

blast_region_segments <- function(hits) {
  if (!nrow(hits)) {
    return(blast_plot_empty())
  }
  breakpoints <- sort(unique(c(hits$start, hits$end + 1)))
  pieces <- vector("list", max(length(breakpoints) - 1L, 0L))
  piece_i <- 0L
  for (i in seq_len(max(length(breakpoints) - 1L, 0L))) {
    start <- breakpoints[[i]]
    end <- breakpoints[[i + 1L]] - 1
    active <- which(hits$start <= start & hits$end >= start)
    if (!length(active) || end < start) {
      next
    }
    rank <- order(
      -hits$pident[active],
      -hits$bitscore[active],
      -hits$length[active]
    )
    winner <- active[[rank[[1L]]]]
    piece_i <- piece_i + 1L
    pieces[[piece_i]] <- data.frame(
      start = start,
      end = end,
      pident = hits$pident[[winner]],
      length = hits$length[[winner]],
      bitscore = hits$bitscore[[winner]],
      fill_key = NA_character_,
      feature_id = NA_character_,
      stringsAsFactors = FALSE
    )
  }
  if (!piece_i) {
    return(blast_plot_empty())
  }
  do.call(rbind, pieces[seq_len(piece_i)])
}

blast_interval_union_width <- function(starts, ends) {
  if (!length(starts)) {
    return(0)
  }
  order_i <- order(starts, ends)
  starts <- starts[order_i]
  ends <- ends[order_i]
  total <- 0
  current_start <- starts[[1L]]
  current_end <- ends[[1L]]
  if (length(starts) > 1L) {
    for (i in 2:length(starts)) {
      if (starts[[i]] <= current_end + 1) {
        current_end <- max(current_end, ends[[i]])
      } else {
        total <- total + current_end - current_start + 1
        current_start <- starts[[i]]
        current_end <- ends[[i]]
      }
    }
  }
  total + current_end - current_start + 1
}

blast_gene_segments <- function(hits, features, genome_length,
                                min_coverage = 0.80) {
  if (!nrow(hits) || is.null(features) || !nrow(features)) {
    return(blast_plot_empty())
  }
  if (!all(c("start", "end") %in% names(features))) {
    stop("Gene-level BLAST plotting requires feature `start` and `end` columns.",
         call. = FALSE)
  }
  pieces <- list()
  piece_i <- 0L
  for (i in seq_len(nrow(features))) {
    feature <- features[i, , drop = FALSE]
    intervals <- feature_segments(feature, genome_length, circular = TRUE)
    if (!nrow(intervals)) {
      next
    }
    gene_start <- intervals$segment_start
    gene_end <- intervals$segment_end
    overlap_start <- outer(
      hits$start, gene_start,
      function(hit_start, feature_start) pmax(hit_start, feature_start)
    )
    overlap_end <- outer(
      hits$end, gene_end,
      function(hit_end, feature_end) pmin(hit_end, feature_end)
    )
    # A hit is a candidate if it overlaps at least one part of the feature.
    candidate <- which(rowSums(overlap_end >= overlap_start) > 0)
    if (!length(candidate)) {
      next
    }
    covered_start <- numeric()
    covered_end <- numeric()
    candidate_coverage <- numeric(length(candidate))
    feature_width <- sum(gene_end - gene_start + 1)
    for (candidate_i in seq_along(candidate)) {
      hit_i <- candidate[[candidate_i]]
      hit_start <- pmax(hits$start[[hit_i]], gene_start)
      hit_end <- pmin(hits$end[[hit_i]], gene_end)
      overlaps <- hit_end >= hit_start
      candidate_coverage[[candidate_i]] <- blast_interval_union_width(
        hit_start[overlaps],
        hit_end[overlaps]
      ) / feature_width
      for (interval_i in seq_len(nrow(intervals))) {
        if (overlaps[[interval_i]]) {
          covered_start <- c(covered_start, hit_start[[interval_i]])
          covered_end <- c(covered_end, hit_end[[interval_i]])
        }
      }
    }
    coverage <- blast_interval_union_width(covered_start, covered_end) /
      feature_width
    if (!is.finite(coverage) || coverage < min_coverage) {
      next
    }
    rank <- order(
      -hits$pident[candidate],
      -candidate_coverage,
      -hits$bitscore[candidate],
      -hits$length[candidate],
      na.last = TRUE
    )
    winner <- candidate[[rank[[1L]]]]
    feature_id <- if ("feature_id" %in% names(feature)) {
      as.character(feature$feature_id[[1L]])
    } else {
      as.character(i)
    }
    for (interval_i in seq_len(nrow(intervals))) {
      piece_i <- piece_i + 1L
      pieces[[piece_i]] <- data.frame(
        start = gene_start[[interval_i]],
        end = gene_end[[interval_i]],
        pident = hits$pident[[winner]],
        length = hits$length[[winner]],
        bitscore = hits$bitscore[[winner]],
        fill_key = NA_character_,
        feature_id = feature_id,
        stringsAsFactors = FALSE
      )
    }
  }
  if (!piece_i) {
    return(blast_plot_empty())
  }
  do.call(rbind, pieces[seq_len(piece_i)])
}

blast_prepare_plot_data <- function(hits, mode = c("hsp", "region", "gene"),
                                    reference = c("query", "subject"),
                                    genome_length, features = NULL,
                                    min_identity = NULL,
                                    min_alignment_length = 1,
                                    min_gene_coverage = 0.80,
                                    colour = "#2C7FB8",
                                    coordinate_offset = 0) {
  mode <- match.arg(mode)
  reference <- match.arg(reference)
  mode_default_identity <- c(hsp = 0, region = 90, gene = 80)[[mode]]
  min_identity <- if (is.null(min_identity)) mode_default_identity else as.numeric(min_identity[[1L]])
  min_alignment_length <- as.numeric(min_alignment_length[[1L]])
  min_gene_coverage <- as.numeric(min_gene_coverage[[1L]])
  genome_length <- as.numeric(genome_length[[1L]])
  if (!is.finite(min_identity) || min_identity < 0 || min_identity > 100) {
    stop("`blast_min_identity` must be between 0 and 100.", call. = FALSE)
  }
  if (!is.finite(min_alignment_length) || min_alignment_length < 1) {
    stop("`blast_min_alignment_length` must be at least 1.", call. = FALSE)
  }
  if (!is.finite(min_gene_coverage) || min_gene_coverage < 0 || min_gene_coverage > 1) {
    stop("`blast_min_gene_coverage` must be between 0 and 1.", call. = FALSE)
  }
  if (!is.finite(genome_length) || genome_length < 1) {
    stop("A valid `genome_length` is required for BLAST plotting.", call. = FALSE)
  }
  normalized <- blast_normalize_plot_hits(
    hits,
    reference = reference,
    min_identity = min_identity,
    min_alignment_length = min_alignment_length,
    coordinate_offset = coordinate_offset
  )
  if (nrow(normalized)) {
    normalized$start <- pmax(normalized$start, 1)
    normalized$end <- pmin(normalized$end, genome_length)
    normalized <- normalized[normalized$end >= normalized$start, , drop = FALSE]
  }
  data <- switch(
    mode,
    hsp = normalized,
    region = blast_region_segments(normalized),
    gene = {
      if (is.null(features) || !nrow(features)) {
        stop(
          "`blast_plotting_mode = \"gene\"` requires annotation features.",
          call. = FALSE
        )
      }
      blast_gene_segments(
        normalized,
        features = features,
        genome_length = genome_length,
        min_coverage = min_gene_coverage
      )
    }
  )
  if (nrow(data)) {
    data$fill_key <- as.character(blast_identity_fill_key(data$pident))
  }
  list(
    data = data,
    colors = blast_identity_colors(colour),
    breaks = blast_identity_breaks(),
    mode = mode,
    reference = reference,
    min_identity = min_identity
  )
}

blast_normalize_ring_specs <- function(blast_rings) {
  if (is.null(blast_rings)) {
    return(list())
  }
  if (is.data.frame(blast_rings)) {
    if (!"label" %in% names(blast_rings) ||
        !any(c("hits", "sequence") %in% names(blast_rings))) {
      stop(
        "A blast_rings data frame must have label and either a hits or sequence list-column.",
        call. = FALSE
      )
    }
    for (column_name in intersect(c("hits", "sequence"), names(blast_rings))) {
      if (!is.list(blast_rings[[column_name]])) {
        stop(
          "blast_rings$", column_name,
          " must be a list-column.",
          call. = FALSE
        )
      }
    }
    specs <- lapply(seq_len(nrow(blast_rings)), function(i) {
      stats::setNames(
        lapply(blast_rings, function(column) {
          if (is.list(column)) column[[i]] else column[[i]]
        }),
        names(blast_rings)
      )
    })
  } else if (is.list(blast_rings)) {
    specs <- blast_rings
    spec_names <- names(specs)
    if (is.null(spec_names) || any(!nzchar(spec_names))) {
      stop(
        "List-form blast_rings must name entries ring1, ring2, and so on.",
        call. = FALSE
      )
    }
  } else {
    stop("blast_rings must be a list or a data frame.", call. = FALSE)
  }
  if (!is.null(names(specs)) && anyDuplicated(names(specs))) {
    stop("BLAST ring list entry names must be unique.", call. = FALSE)
  }
  if (!length(specs)) {
    return(specs)
  }
  for (i in seq_along(specs)) {
    spec <- specs[[i]]
    if (!is.list(spec) || is.data.frame(spec)) {
      stop("Each blast_rings entry must be a list of ring settings.",
           call. = FALSE)
    }
    if (is.null(spec$label) || length(spec$label) != 1L ||
        is.na(spec$label) || !nzchar(trimws(as.character(spec$label)))) {
      stop("Each BLAST ring needs a non-empty label.", call. = FALSE)
    }
    has_hits <- !is.null(spec$hits) && is.data.frame(spec$hits)
    has_sequence <- !is.null(spec$sequence) &&
      !(length(spec$sequence) == 1L && is.atomic(spec$sequence) &&
        is.na(spec$sequence))
    if (identical(has_hits, has_sequence)) {
      stop(
        "Each BLAST ring must supply exactly one of hits or sequence.",
        call. = FALSE
      )
    }
    spec$label <- trimws(as.character(spec$label))
    specs[[i]] <- spec
  }
  labels <- vapply(specs, function(spec) spec$label, character(1L))
  if (anyDuplicated(labels)) {
    stop("Each BLAST ring label must be unique.", call. = FALSE)
  }
  if (is.null(names(specs)) || any(!nzchar(names(specs)))) {
    names(specs) <- paste0("ring", seq_along(specs))
  }
  specs
}

blast_prepare_rings_plot_data <- function(
    blast_rings, genome_length, features, reference = "query",
    mode = "region", min_identity = NULL, min_alignment_length = 1,
    min_gene_coverage = 0.80, palette = "npg", radius = NULL,
    height = 0.08, spacing = 0.025, coordinate_offset = 0,
    backbone_sequence = NULL, blastn = NULL, evalue = 1e-10,
    blast_args = character(), color_manual = NULL) {
  specs <- blast_normalize_ring_specs(blast_rings)
  if (!length(specs)) {
    return(NULL)
  }
  available_palettes <- ggplasmid_palette_names()
  if (length(palette) != 1L || is.na(palette) ||
      !palette %in% available_palettes) {
    stop(
      "blast_palette must be one of: ",
      paste(available_palettes, collapse = ", "),
      call. = FALSE
    )
  }
  ring_order <- vapply(seq_along(specs), function(i) {
    value <- specs[[i]]$order %||% i
    value <- suppressWarnings(as.numeric(value[[1L]]))
    if (!is.finite(value) || value < 1 || value != floor(value)) {
      stop("Each ring order must be a positive integer.", call. = FALSE)
    }
    value
  }, numeric(1L))
  if (anyDuplicated(ring_order)) {
    stop("BLAST ring order values must be unique.", call. = FALSE)
  }
  draw_order <- order(ring_order)
  palette_fun <- getExportedValue("ggsci", paste0("pal_", palette))
  palette_colors <- suppressWarnings(palette_fun()(max(length(specs), 1L)))
  palette_colors <- palette_colors[!is.na(palette_colors) & nzchar(palette_colors)]
  if (!length(palette_colors)) {
    stop("The selected ggsci palette did not provide any colors.", call. = FALSE)
  }
  if (length(palette_colors) < length(specs)) {
    palette_colors <- grDevices::colorRampPalette(palette_colors)(length(specs))
  }
  base_colors <- stats::setNames(rep(NA_character_, length(specs)), names(specs))
  base_colors[draw_order] <- palette_colors[seq_along(specs)]
  if (!is.null(color_manual)) {
    if (!is.character(color_manual) || length(color_manual) != length(specs)) {
      stop(
        "blast_color_manual must provide exactly one color for each BLAST ring.",
        call. = FALSE
      )
    }
    if (!is.null(names(color_manual)) && all(nzchar(names(color_manual)))) {
      ring_labels <- vapply(specs, function(spec) spec$label, character(1L))
      if (!setequal(names(color_manual), ring_labels)) {
        stop(
          "Named blast_color_manual colors must be named for every ring label.",
          call. = FALSE
        )
      }
      color_manual <- color_manual[match(ring_labels, names(color_manual))]
    }
    for (i in seq_along(color_manual)) {
      tryCatch(
        grDevices::col2rgb(color_manual[[i]], alpha = TRUE),
        error = function(error) {
          stop("Invalid color in blast_color_manual.", call. = FALSE)
        }
      )
    }
    base_colors[] <- unname(color_manual)
  }
  for (i in seq_along(specs)) {
    if (!is.null(specs[[i]]$colour) &&
        length(specs[[i]]$colour) == 1L &&
        !is.na(specs[[i]]$colour)) {
      base_colors[[i]] <- as.character(specs[[i]]$colour[[1L]])
      tryCatch(
        grDevices::col2rgb(base_colors[[i]], alpha = TRUE),
        error = function(error) {
          stop("Invalid BLAST ring colour for ", specs[[i]]$label, ".",
               call. = FALSE)
        }
      )
    }
  }
  if (is.null(radius)) {
    stop("A base radius is required to draw BLAST rings.", call. = FALSE)
  }
  radius <- as.numeric(radius[[1L]])
  height <- as.numeric(height[[1L]])
  spacing <- as.numeric(spacing[[1L]])
  if (!is.finite(radius) || radius <= 0 ||
      !is.finite(height) || height <= 0 ||
      !is.finite(spacing) || spacing < 0) {
    stop("BLAST ring radius/height must be positive and spacing non-negative.",
         call. = FALSE)
  }
  ring_positions <- stats::setNames(rep(NA_real_, length(specs)), names(specs))
  ring_positions[draw_order] <- radius +
    (seq_along(specs) - 1L) * (height + spacing)

  all_data <- vector("list", length(specs))
  all_colors <- character()
  all_breaks <- character()
  identity_thresholds <- stats::setNames(rep(NA_real_, length(specs)), names(specs))
  for (i in draw_order) {
    spec <- specs[[i]]
    ring_mode <- match.arg(spec$mode %||% mode, c("hsp", "region", "gene"))
    ring_identity <- spec$min_identity %||% min_identity
    ring_alignment <- spec$min_alignment_length %||% min_alignment_length
    ring_coverage <- spec$min_gene_coverage %||% min_gene_coverage
    ring_reference <- match.arg(spec$reference %||% reference, c("query", "subject"))
    hits <- spec$hits
    if (is.null(hits)) {
      if (is.null(backbone_sequence) || length(backbone_sequence) != 1L ||
          !nzchar(as.character(backbone_sequence))) {
        stop(
          "Automatic BLAST requires a backbone sequence. Supply it through a GenBank ORIGIN sequence or FASTA.",
          call. = FALSE
        )
      }
      search_identity <- ring_identity
      if (is.null(search_identity)) {
        search_identity <- c(hsp = 0, region = 90, gene = 80)[[ring_mode]]
      }
      hits <- blast_pairwise(
        query = as.character(backbone_sequence),
        subject = spec$sequence,
        blastn = blastn,
        perc_identity = search_identity,
        evalue = evalue,
        blast_args = blast_args
      )
      specs[[i]]$hits <- hits
    }
    prepared <- blast_prepare_plot_data(
      hits = hits,
      mode = ring_mode,
      reference = if (is.null(spec$hits)) "query" else ring_reference,
      genome_length = genome_length,
      features = features,
      min_identity = ring_identity,
      min_alignment_length = ring_alignment,
      min_gene_coverage = ring_coverage,
      colour = base_colors[[i]],
      coordinate_offset = coordinate_offset
    )
    data <- prepared$data
    identity_thresholds[[i]] <- prepared$min_identity
    if (!"feature_id" %in% names(data)) {
      data$feature_id <- rep(NA_character_, nrow(data))
    }
    if (nrow(data)) {
      data$ring_key <- names(specs)[[i]]
      data$ring_label <- spec$label
      data$ring_order <- ring_order[[i]]
      data$radius <- ring_positions[[i]]
      shades <- grDevices::colorRampPalette(
        c("#E2E2E2", base_colors[[i]])
      )(length(blast_identity_breaks()))
      names(shades) <- blast_identity_breaks()
      identity_bin <- data$fill_key
      data$fill_key <- paste(data$ring_key, identity_bin, sep = "::")
      observed_bins <- rev(blast_identity_breaks())[
        rev(blast_identity_breaks()) %in% unique(identity_bin)
      ]
      all_colors <- c(
        all_colors,
        stats::setNames(
          shades[observed_bins],
          paste(names(specs)[[i]], observed_bins, sep = "::")
        )
      )
      all_breaks <- c(
        all_breaks,
        paste(names(specs)[[i]], observed_bins, sep = "::")
      )
    }
    all_data[[i]] <- data
  }
  all_data <- Filter(function(x) !is.null(x) && nrow(x), all_data)
  combined_data <- if (length(all_data)) do.call(rbind, all_data) else blast_plot_empty()
  row.names(combined_data) <- NULL
  labels <- vapply(all_breaks, function(key) {
    parts <- strsplit(key, "::", fixed = TRUE)[[1L]]
    ring <- specs[[parts[[1L]]]]$label
    paste0(ring, " — ", parts[[2L]])
  }, character(1L))
  list(
    data = combined_data,
    colors = all_colors,
    breaks = all_breaks,
    labels = labels,
    specs = specs,
    order = ring_order,
    base_colors = base_colors,
    identity_thresholds = identity_thresholds,
    radius = radius,
    height = height,
    spacing = spacing
  )
}

#' Run a pairwise nucleotide BLAST search
#'
#' This is a small wrapper around the local NCBI BLAST+ `blastn` executable.
#' It uses BLAST's direct `-query`/`-subject` mode, so a persistent BLAST
#' database is not required. Query and subject may be FASTA paths, Biostrings
#' `XStringSet` objects, data frames with a `sequence` column, or character
#' vectors of nucleotide sequences.
#'
#' @param query Query FASTA path, Biostrings `XStringSet`, sequence data frame,
#'   or character vector.
#' @param subject Subject FASTA path, Biostrings `XStringSet`, sequence data
#'   frame, or character vector.
#' @param blastn Optional path to the `blastn` executable. By default it is
#'   located using `Sys.which("blastn")`.
#' @param perc_identity Optional minimum percent identity passed to BLAST.
#' @param evalue Optional maximum E-value passed to BLAST.
#' @param word_size Optional BLAST word size.
#' @param num_threads Optional number of BLAST threads.
#' @param max_target_seqs Optional maximum number of target sequences per query.
#' @param blast_args Additional BLAST command-line arguments as a character
#'   vector, for example `c("-task", "blastn-short")`.
#' @param output Optional path for the raw tabular BLAST output. If `NULL`, a
#'   temporary output file is used.
#' @param keep_files Whether temporary FASTA and output files should be kept.
#' @param verbose Whether to print the executed command and BLAST messages.
#' @return A data frame with standard BLAST tabular fields plus normalized
#'   `query_start`, `query_end`, `subject_start`, `subject_end`, and
#'   `subject_strand` columns.
#' @export
blast_pairwise <- function(query, subject, blastn = NULL,
                            perc_identity = NULL, evalue = NULL,
                            word_size = NULL, num_threads = NULL,
                            max_target_seqs = NULL, blast_args = character(),
                            output = NULL, keep_files = FALSE,
                            verbose = FALSE) {
  if (!is.character(blast_args)) {
    stop("`blast_args` must be a character vector.", call. = FALSE)
  }
  reserved_args <- c("-query", "-subject", "-out", "-outfmt")
  if (any(blast_args %in% reserved_args)) {
    stop(
      "`blast_args` cannot override the query, subject, output, or output-format arguments.",
      call. = FALSE
    )
  }
  if (length(keep_files) != 1L || is.na(keep_files)) {
    stop("`keep_files` must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }

  blastn <- blast_pairwise_resolve_executable(blastn)
  work_dir <- tempfile("ggplasmid-blast-")
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit({
    if (!isTRUE(keep_files)) {
      unlink(work_dir, recursive = TRUE, force = TRUE)
    }
  }, add = TRUE)

  query_file <- blast_pairwise_sequence_input(query, "query", work_dir)
  subject_file <- blast_pairwise_sequence_input(subject, "subject", work_dir)
  output_file <- if (is.null(output)) {
    file.path(work_dir, "blast.tsv")
  } else {
    if (!is.character(output) || length(output) != 1L || !nzchar(output)) {
      stop("`output` must be a single non-empty file path.", call. = FALSE)
    }
    output <- normalizePath(output, mustWork = FALSE)
    parent <- dirname(output)
    if (!dir.exists(parent)) {
      stop("The directory for `output` does not exist: ", parent, call. = FALSE)
    }
    output
  }

  args <- c(
    "-query", query_file$path,
    "-subject", subject_file$path,
    "-out", output_file,
    "-outfmt",
    blast_pairwise_quote_arg(paste(c("6", blast_pairwise_fields), collapse = " "))
  )
  args <- blast_pairwise_numeric_arg(
    args, "-perc_identity", perc_identity, 0, maximum = 100
  )
  args <- blast_pairwise_numeric_arg(args, "-evalue", evalue, 0)
  args <- blast_pairwise_integer_arg(args, "-word_size", word_size, 1L)
  args <- blast_pairwise_integer_arg(args, "-num_threads", num_threads, 1L)
  args <- blast_pairwise_integer_arg(args, "-max_target_seqs", max_target_seqs, 1L)
  args <- c(args, blast_args)

  if (isTRUE(verbose)) {
    message("Running: ", paste(c(shQuote(blastn), args), collapse = " "))
  }
  status <- system2(
    command = blastn,
    args = args,
    stdout = if (isTRUE(verbose)) "" else TRUE,
    stderr = if (isTRUE(verbose)) "" else TRUE
  )
  exit_status <- attr(status, "status", exact = TRUE)
  if (is.null(exit_status)) {
    exit_status <- if (is.numeric(status) && length(status) == 1L) {
      status[[1L]]
    } else {
      0L
    }
  }
  if (!identical(as.integer(exit_status[[1L]]), 0L)) {
    message_text <- if (is.character(status)) paste(status, collapse = "\n") else ""
    stop(
      "blastn failed with exit status ", exit_status[[1L]],
      if (nzchar(message_text)) paste0(":\n", message_text) else "",
      call. = FALSE
    )
  }
  if (!file.exists(output_file)) {
    stop("blastn completed but did not create an output file.", call. = FALSE)
  }

  hits <- blast_pairwise_read_output(output_file)
  attr(hits, "blast_command") <- c(blastn, args)
  attr(hits, "blast_query") <- query_file$path
  attr(hits, "blast_subject") <- subject_file$path
  attr(hits, "blast_output") <- output_file
  hits
}
