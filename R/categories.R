#' Return ggplasmid category colors
#'
#' @param scheme Either `"plasmid"` or `"phage"`.
#' @return A named character vector of hex colors.
#' @export
ggplasmid_colors <- function(scheme = c("plasmid", "phage")) {
  scheme <- match.arg(scheme)
  if (scheme == "phage") {
    c(
      "Connectors" = "#D6A800",
      "DNA, RNA and nucleotide metabolism" = "#FF0000",
      "Head and packaging" = "#000000",
      "Lysis" = "#FF00FF",
      "Moron, auxiliary metabolic gene and host takeover" = "#FF9900",
      "Other functions" = "#009999",
      "Tail" = "#002060",
      "Transcription regulation" = "#008000",
      "tRNA" = "#7030A0",
      "Unknown function, hypothetical protein" = "#AEAAAA",
      "GC skew+" = "#008000",
      "GC skew-" = "#7030A0"
    )
  } else {
    c(
      "Replication" = "#0072B2",
      "Conjugation and transfer" = "#D55E00",
      "Partition and maintenance" = "#009E73",
      "Antimicrobial resistance" = "#CC79A7",
      "Virulence" = "#E69F00",
      "Mobile element" = "#8B4513",
      "Metabolism" = "#56B4E9",
      "Regulation" = "#F0E442",
      "Other functions" = "#999999",
      "tRNA" = "#7030A0",
      "Unknown function, hypothetical protein" = "#C7C7C7",
      "GC skew+" = "#008000",
      "GC skew-" = "#7030A0"
    )
  }
}

#' Define category color highlights
#'
#' Use this helper with the `gene_manual_fill` argument in [ggplasmid()], or
#' add it to an existing plot with `+`, to override one or more category fills.
#' `gene_highlight` remains a compatibility alias in plotting functions.
#'
#' @param ... Named color values, for example
#'   `"Antimicrobial resistance" = "#B2182B"`, or a data frame.
#' @param scheme Either `"plasmid"` or `"phage"`; used to validate category
#'   names.
#' @param highlight Optional named character vector or data frame with
#'   category/color columns.
#' @return A named character vector of category colors. When added to a
#'   ggplasmid plot, the plot is returned with an updated fill scale.
#' @export
gene_highlight <- function(..., scheme = c("plasmid", "phage"), highlight = NULL) {
  scheme <- match.arg(scheme)
  dots <- list(...)
  if (length(dots) == 1L && is.data.frame(dots[[1L]])) {
    values <- dots[[1L]]
  } else {
    values <- c(...)
    dot_names <- names(dots)
    if (length(dot_names) && all(nzchar(dot_names)) &&
        (is.null(names(values)) || any(!nzchar(names(values))))) {
      values <- stats::setNames(unlist(dots, use.names = FALSE), dot_names)
    }
  }
  if (!is.null(highlight)) {
    if (is.data.frame(highlight) && !length(values)) {
      values <- highlight
    } else {
      values <- c(values, highlight)
    }
  }

  if (is.data.frame(values)) {
    category_col <- find_column(values, c("category", "name", "feature_category"))
    color_col <- find_column(values, c("color", "colour", "hex", "fill"))
    values <- stats::setNames(as.character(values[[color_col]]), values[[category_col]])
  }
  if (is.list(values) && !is.data.frame(values)) {
    values <- unlist(values, use.names = TRUE)
  }
  value_names <- names(values)
  values <- as.character(values)
  names(values) <- value_names
  if (!length(values)) {
    return(stats::setNames(character(), character()))
  }
  if (is.null(names(values)) || any(!nzchar(names(values)))) {
    stop("`gene_highlight()` values must be named by category.", call. = FALSE)
  }

  available <- names(ggplasmid_colors(scheme))
  unknown <- setdiff(names(values), available)
  if (length(unknown)) {
    stop(
      "Unknown category in `gene_highlight`: ",
      paste(unknown, collapse = ", "),
      ". Available categories are: ",
      paste(available, collapse = ", "),
      call. = FALSE
    )
  }
  class(values) <- c("ggplasmid_gene_highlight", class(values))
  values
}

#' @export
ggplot_add.ggplasmid_gene_highlight <- function(object, plot, object_name) {
  style <- attr(plot, "ggplasmid_style", exact = TRUE)
  if (is.null(style)) {
    stop(
      "`gene_highlight()` can be added only to a plot returned by ",
      "`ggplasmid()`, `plot_plasmid_map()`, or `plot_phage_map()`.",
      call. = FALSE
    )
  }

  plot_data <- attr(plot, "ggplasmid_data", exact = TRUE)
  layout <- attr(plot, "ggplasmid_layout", exact = TRUE)
  scheme <- style$category_scheme %||% "plasmid"
  base <- style$fill_colors %||% ggplasmid_resolve_colors(
    scheme = scheme,
    gene_palette = style$gene_palette %||% "npg",
    gene_manual_fill = style$gene_manual_fill
  )
  highlight <- gene_highlight(highlight = unclass(object), scheme = scheme)
  base[names(highlight)] <- unname(highlight)

  fill_scale_indices <- which(vapply(
    plot$scales$scales,
    function(scale) any(grepl("^fill($|_)", scale$aesthetics)),
    logical(1L)
  ))
  gene_fill_aesthetic <- if (length(fill_scale_indices)) {
    plot$scales$scales[[fill_scale_indices[[1L]]]]$aesthetics
  } else {
    "fill"
  }
  if (length(fill_scale_indices)) {
    plot$scales$scales <- plot$scales$scales[-fill_scale_indices[[1L]]]
  }
  plot <- plot +
    ggplot2::scale_fill_manual(
      aesthetics = gene_fill_aesthetic,
      values = base,
      breaks = gene_legend_breaks(base),
      drop = TRUE,
      na.value = "#999999",
      guide = ggplot2::guide_legend(
        order = 2,
        ncol = style$legend_columns %||% 1L,
        byrow = TRUE,
        keyheight = grid::unit(0.34, "cm"),
        keywidth = grid::unit(0.45, "cm"),
        override.aes = list(colour = NA)
      )
    )

  style$fill_colors <- base
  style$fill_breaks <- gene_legend_breaks(base)
  attr(plot, "ggplasmid_style") <- style
  attr(plot, "ggplasmid_data") <- plot_data
  attr(plot, "ggplasmid_layout") <- layout
  plot
}

ggplasmid_palette_names <- function() {
  exports <- getNamespaceExports("ggsci")
  sort(sub("^pal_", "", exports[grepl("^pal_[[:alnum:]]+$", exports)]))
}

ggplasmid_resolve_colors <- function(scheme = c("plasmid", "phage"),
                                     gene_palette = "npg",
                                     gene_manual_fill = NULL) {
  scheme <- match.arg(scheme)
  gene_palette <- match.arg(gene_palette, choices = ggplasmid_palette_names())
  base <- ggplasmid_colors(scheme)

  pal_fun <- if (identical(gene_palette, "npg")) {
    ggsci::pal_npg
  } else {
    getExportedValue("ggsci", paste0("pal_", gene_palette))
  }
  generated <- suppressWarnings(pal_fun()(length(base)))
  available <- generated[!is.na(generated) & nzchar(generated)]
  if (length(available) < length(base)) {
    generated <- if (length(available) == 1L) {
      rep(available, length(base))
    } else {
      grDevices::colorRampPalette(available)(length(base))
    }
  }
  base <- stats::setNames(generated, names(base))
  base[c("GC skew+", "GC skew-")] <- c("#008000", "#7030A0")

  if (!is.null(gene_manual_fill)) {
    highlight <- do.call(
      "gene_highlight",
      list(highlight = gene_manual_fill, scheme = scheme)
    )
    base[names(highlight)] <- unname(highlight)
  }
  base
}

normalise_supplied_category <- function(value, scheme) {
  text <- tolower(clean_text(value))
  if (!nzchar(text)) {
    return(NA_character_)
  }

  if (grepl("\\btrna\\b|\\brrna\\b|tmrna", text)) {
    return("tRNA")
  }

  if (scheme == "phage") {
    if (grepl("hypothetical|unknown|uncharacterized|uncharacterised", text)) {
      return("Unknown function, hypothetical protein")
    }
    if (grepl("connector|neck|adaptor|adapter", text)) return("Connectors")
    if (grepl("dna|rna|nucleotide", text)) return("DNA, RNA and nucleotide metabolism")
    if (grepl("head|packaging|capsid|terminase|portal", text)) return("Head and packaging")
    if (grepl("lysis|lysin|holin|spanin", text)) return("Lysis")
    if (grepl("moron|auxiliary|host takeover|host", text)) {
      return("Moron, auxiliary metabolic gene and host takeover")
    }
    if (grepl("tail|baseplate|fiber|fibre|sheath|spike|tube|tape measure", text)) return("Tail")
    if (grepl("transcription|regulation|regulator|repressor|sigma", text)) {
      return("Transcription regulation")
    }
    if (grepl("other", text)) return("Other functions")
  } else {
    if (grepl("hypothetical|unknown|uncharacterized|uncharacterised", text)) {
      return("Unknown function, hypothetical protein")
    }
    if (grepl("replication|replicase|rep\\b|origin|\\bori\\b|primase", text)) return("Replication")
    if (grepl("transpos|integrase|insertion|\\bis\\b|mobile|recombinase", text)) {
      return("Mobile element")
    }
    if (grepl("conjug|transfer|mobil|relaxase|\\btra[a-z]?\\b|\\btrb[a-z]?\\b|type iv", text)) {
      return("Conjugation and transfer")
    }
    if (grepl("partition|maintenance|para|parb|toxin|antitoxin|addiction", text)) {
      return("Partition and maintenance")
    }
    if (grepl("resistance|antimicrobial|antibiotic|amr|beta-lactamase|betalactamase", text)) {
      return("Antimicrobial resistance")
    }
    if (grepl("virulence|vfdb|pathogen|toxin", text)) return("Virulence")
    if (grepl("metabolism|metabolic|dehydrogenase|hydrolase|oxidase|reductase", text)) {
      return("Metabolism")
    }
    if (grepl("regulation|regulator|repressor|transcription|sigma|response regulator", text)) {
      return("Regulation")
    }
    if (grepl("other", text)) return("Other functions")
  }

  NA_character_
}

classify_one_feature <- function(category_source, label, type, scheme) {
  supplied <- normalise_supplied_category(category_source, scheme)
  if (!is.na(supplied)) {
    return(supplied)
  }

  text <- tolower(paste(clean_text(type), clean_text(label), sep = " "))
  normalise_supplied_category(text, scheme) %||% "Other functions"
}

classify_features <- function(data, category_col = NA_character_,
                              label_col = "label", type_col = "type",
                              scheme = "plasmid") {
  category_source <- rep("", nrow(data))
  if (!is.na(category_col) && category_col %in% names(data)) {
    category_source <- data[[category_col]]
  }
  vapply(
    seq_len(nrow(data)),
    function(i) {
      classify_one_feature(
        category_source[[i]],
        data[[label_col]][[i]],
        data[[type_col]][[i]],
        scheme
      )
    },
    character(1)
  )
}
