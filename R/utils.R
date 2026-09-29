`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

utils::globalVariables(c(
  "anchor_x", "anchor_y", "angle", "base_y", "category", "hjust",
  "content_y", "label", "leader_x", "leader_y", "neg_y", "npcx", "npcy", "panel_x", "panel_y",
  "polygon_id", "pos_y", "position", "row_index", "row_y", "type",
  "track", "text_x", "text_y", "vjust", "x", "x0", "x1", "x2", "xend",
  "xmid_linear", "y", "y0", "y1", "y2", "yend"
))

first_or <- function(x, y) {
  if (is.null(x) || !length(x) || is.na(x[[1]])) y else x[[1]]
}

clean_text <- function(x) {
  x <- ifelse(is.na(x), "", as.character(x))
  trimws(gsub("\\s+", " ", x))
}

normalise_key <- function(x) {
  tolower(gsub("[^a-z0-9]+", "", clean_text(x)))
}

read_delim_auto <- function(path) {
  if (is.data.frame(path)) {
    return(as.data.frame(path, stringsAsFactors = FALSE))
  }
  first <- readLines(path, n = 1, warn = FALSE)
  if (!length(first)) {
    stop("Input table is empty: ", path, call. = FALSE)
  }
  tab_count <- lengths(regmatches(first, gregexpr("\t", first, fixed = TRUE)))
  comma_count <- lengths(regmatches(first, gregexpr(",", first, fixed = TRUE)))
  sep <- if (tab_count >= comma_count) "\t" else ","

  utils::read.delim(
    path,
    sep = sep,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "\"",
    comment.char = "",
    row.names = NULL,
    fill = TRUE
  )
}

find_column <- function(data, candidates, required = TRUE, label = "column") {
  lookup <- stats::setNames(names(data), normalise_key(names(data)))
  for (candidate in candidates) {
    key <- normalise_key(candidate)
    if (key %in% names(lookup)) {
      return(unname(lookup[[key]]))
    }
  }
  if (required) {
    stop(
      "Missing ", label, "; expected one of: ",
      paste(candidates, collapse = ", "),
      call. = FALSE
    )
  }
  NA_character_
}

as_number <- function(x) {
  suppressWarnings(as.numeric(gsub(",", "", as.character(x))))
}

strand_direction <- function(x) {
  text <- tolower(clean_text(x))
  ifelse(text %in% c("-", "-1", "minus", "reverse", "rev", "complement"), -1L, 1L)
}

legend_position_is_vertical <- function(position) {
  position %in% c(
    "right", "left",
    "left_top", "right_top", "left_bottom", "right_bottom"
  )
}

legend_position_for_theme <- function(position) {
  switch(
    position,
    left_top = "left",
    left_bottom = "left",
    right_top = "right",
    right_bottom = "right",
    position
  )
}

legend_position_inside_for_theme <- function(position, label_bounds = NULL) {
  # Corner positions use an outside side position with corner-specific
  # justification; no legend is anchored inside the map panel.
  NULL
}

legend_justification_for_theme <- function(position) {
  switch(
    position,
    left_top = c(1, 1),
    right_top = c(0, 1),
    left_bottom = c(1, 0),
    right_bottom = c(0, 0),
    left = c(1, 0.5),
    right = c(0, 0.5),
    top = c(0.5, 0),
    bottom = c(0.5, 1),
    "center"
  )
}

legend_box_justification_for_theme <- function(position) {
  "left"
}

validate_legend_columns <- function(columns = NULL, position = "bottom") {
  if (is.null(columns)) {
    return(if (legend_position_is_vertical(position)) 1L else 3L)
  }
  columns <- as.integer(columns)
  if (length(columns) != 1L || !is.finite(columns) || columns < 1L) {
    stop("`legend_columns` must be a positive integer.", call. = FALSE)
  }
  columns
}

validate_gc_legend_columns <- function(columns = NULL, position = "bottom") {
  if (is.null(columns)) {
    return(if (legend_position_is_vertical(position)) 1L else 3L)
  }
  columns <- as.integer(columns)
  if (length(columns) != 1L || !is.finite(columns) || columns < 1L) {
    stop("`gc_legend_columns` must be a positive integer.", call. = FALSE)
  }
  columns
}

validate_legend_plot_spacing <- function(spacing = NULL, position = "bottom") {
  if (is.null(spacing)) {
    return(if (legend_position_is_vertical(position)) 3.0 else 0.35)
  }
  spacing <- as.numeric(spacing)
  if (length(spacing) != 1L || !is.finite(spacing) || spacing < 0) {
    stop("`legend_plot_spacing` must be a non-negative number.", call. = FALSE)
  }
  spacing
}

gene_legend_breaks <- function(colors) {
  setdiff(names(colors), c("GC skew+", "GC skew-"))
}

is_unknown_feature_label <- function(...) {
  parts <- lapply(list(...), clean_text)
  text <- tolower(do.call(paste, c(parts, sep = " ")))
  grepl("hypothetical|unknown|uncharacterized|uncharacterised", text)
}

gc_track_colors <- function(fill_colors) {
  c(
    "GC content" = "grey20",
    "GC skew+" = unname(fill_colors[["GC skew+"]]),
    "GC skew-" = unname(fill_colors[["GC skew-"]])
  )
}

gc_skew_legend_data <- function() {
  data.frame(
    x = 0,
    xend = 0,
    y = 0,
    yend = 0,
    track = c("GC skew+", "GC skew-"),
    stringsAsFactors = FALSE
  )
}

wrap_label <- function(text, width = 22, max_lines = 3) {
  text <- clean_text(text)
  if (!nzchar(text)) {
    return("")
  }
  lines <- strwrap(text, width = width)
  lines <- trimws(lines)
  if (length(lines) > max_lines) {
    lines <- lines[seq_len(max_lines)]
    lines[[max_lines]] <- paste0(sub("[[:space:].]+$", "", lines[[max_lines]]), "...")
  }
  paste(lines, collapse = "\n")
}

shorten_feature_label <- function(text) {
  text <- clean_text(text)
  text <- gsub("\\s*\\(EC [^)]+\\)", "", text)

  if (grepl("@", text, fixed = TRUE)) {
    text <- sub("^.*@\\s*", "", text)
  }

  if (grepl("=>", text, fixed = TRUE)) {
    text <- sub(".*=>\\s*", "", text)
  }

  text <- gsub("^Chromosome \\(plasmid\\) ", "", text)
  text <- gsub("^IncF plasmid conjugative transfer ", "", text)
  text <- gsub("^Type I restriction-modification system, ", "Type I R-M ", text)
  text <- gsub("^Small multidrug resistance \\(SMR\\) efflux transporter", "SMR efflux transporter", text)
  text <- gsub("^Dihydropteroate synthase type-2", "Sulfonamide resistance protein", text)
  text <- gsub("^Subclass B1 beta-lactamase", "NDM beta-lactamase", text)

  if (nchar(text, type = "chars") > 42 && grepl(",", text, fixed = TRUE)) {
    text <- sub(",.*$", "", text)
  }
  clean_text(text)
}

format_bp <- function(x) {
  format(round(as.numeric(x)), big.mark = ",", scientific = FALSE, trim = TRUE)
}

read_fasta_sequence <- function(fasta) {
  fasta_data <- read_plasmid_fasta(fasta)
  sequence <- fasta_data$sequence[[1]]
  if (!nzchar(sequence)) {
    stop("No DNA sequence found in the first FASTA record.", call. = FALSE)
  }
  sequence
}

feature_midpoint <- function(start, end, genome_length) {
  start <- as.numeric(start)
  end <- as.numeric(end)
  ifelse(
    is.finite(genome_length) & start > end,
    ((start + end + genome_length) / 2 - 1) %% genome_length + 1,
    (pmin(start, end) + pmax(start, end)) / 2
  )
}

radial_coord <- function(theta = "x", start = 0, inner.radius = 0.35,
                         clip = "off", expand = FALSE) {
  coord_radial <- get0("coord_radial", envir = asNamespace("ggplot2"), inherits = FALSE)
  if (is.function(coord_radial)) {
    coord_radial(
      theta = theta,
      start = start,
      inner.radius = inner.radius,
      clip = clip,
      expand = expand
    )
  } else {
    ggplot2::coord_polar(theta = theta, start = start, clip = clip)
  }
}
