utils::globalVariables(c(
  "colour", "label_x", "label_y", "line_colour", "line_x1", "line_y0", "line_y1",
  "linear_label_gene_y"
))

feature_segments <- function(features, genome_length, circular = TRUE) {
  pieces <- vector("list", nrow(features))
  for (i in seq_len(nrow(features))) {
    row <- features[i, , drop = FALSE]
    wraps_origin <- if ("wraps_origin" %in% names(row)) {
      isTRUE(row$wraps_origin[[1]])
    } else {
      row$start > row$end
    }
    if (circular && wraps_origin) {
      pieces[[i]] <- rbind(
        transform(row, segment_start = row$start, segment_end = genome_length),
        transform(row, segment_start = 1, segment_end = row$end)
      )
    } else {
      pieces[[i]] <- transform(
        row,
        segment_start = pmin(row$start, row$end),
        segment_end = pmax(row$start, row$end)
      )
    }
  }
  out <- do.call(rbind, pieces)
  row.names(out) <- NULL
  out$segment_id <- seq_len(nrow(out))
  out
}

gene_polygon_one <- function(row, y, gene_height, arrow_head_bp,
                             arrow_head_fraction = 0.55) {
  x1 <- as.numeric(row$segment_start)
  x2 <- as.numeric(row$segment_end)
  width <- max(abs(x2 - x1), 1)
  head <- min(width * arrow_head_fraction, arrow_head_bp)
  y_low <- y - gene_height / 2
  y_high <- y + gene_height / 2

  if (row$direction >= 0L) {
    x <- c(x1, x2 - head, x2, x2 - head, x1)
  } else {
    x <- c(x2, x1 + head, x1, x1 + head, x2)
  }
  data.frame(
    feature_id = row$feature_id,
    polygon_id = row$segment_id %||% row$feature_id,
    x = x,
    y = c(y_low, y_low, y, y_high, y_high),
    category = row$category,
    stringsAsFactors = FALSE
  )
}

circular_gene_polygons <- function(features, genome_length,
                                    radius = 1, gene_height = 0.10,
                                    arrow_head_bp = NULL,
                                    arrow_head_fraction = 0.55,
                                    gene_gap_bp = 0) {
  if (is.null(arrow_head_bp)) {
    arrow_head_bp <- max(80, genome_length / 240)
  }
  segments <- feature_segments(features, genome_length, circular = TRUE)
  gene_gap_bp <- as.numeric(gene_gap_bp %||% 0)
  if (is.finite(gene_gap_bp) && gene_gap_bp > 0) {
    width <- pmax(segments$segment_end - segments$segment_start, 0)
    shrink <- pmin(gene_gap_bp / 2, pmax((width - 1) / 2, 0))
    segments$segment_start <- segments$segment_start + shrink
    segments$segment_end <- segments$segment_end - shrink
  }
  polys <- lapply(seq_len(nrow(segments)), function(i) {
    gene_polygon_one(
      segments[i, , drop = FALSE],
      radius,
      gene_height,
      arrow_head_bp,
      arrow_head_fraction = arrow_head_fraction
    )
  })
  do.call(rbind, polys)
}

linear_gene_polygons <- function(features, genome_length, rows = 4,
                                 gene_height = 0.18, arrow_head_bp = NULL,
                                 arrow_head_fraction = 0.35,
                                 row_spacing = 3.40) {
  rows <- as.integer(rows[[1L]])
  if (!is.finite(rows) || rows < 1L) {
    stop("`rows` must be a positive integer.", call. = FALSE)
  }
  row_spacing <- as.numeric(row_spacing[[1L]])
  if (!is.finite(row_spacing) || row_spacing <= 0) {
    stop("`linear_row_spacing` must be a positive number.", call. = FALSE)
  }
  segment <- ceiling(genome_length / rows)
  if (!is.finite(segment) || segment <= 0) {
    segment <- genome_length
  }
  if (is.null(arrow_head_bp)) {
    arrow_head_bp <- max(35, segment / 180)
  }
  row_start <- seq(0, by = segment, length.out = rows)
  row_end <- pmin(row_start + segment, genome_length)
  row_y <- 1 + (rows - seq_len(rows)) * row_spacing
  xmid <- feature_midpoint(features$start, features$end, genome_length)
  row_index <- pmin(floor((xmid - 1) / segment) + 1L, rows)

  segments <- feature_segments(features, genome_length, circular = TRUE)
  segment_mid <- (segments$segment_start + segments$segment_end) / 2
  segments$row_index <- pmin(floor((segment_mid - 1) / segment) + 1L, rows)
  segments$row_y <- row_y[segments$row_index]
  segments$segment_start <- pmax(segments$segment_start - row_start[segments$row_index], 0)
  segments$segment_end <- pmin(segments$segment_end - row_start[segments$row_index], segment)

  polys <- lapply(seq_len(nrow(segments)), function(i) {
    gene_polygon_one(
      segments[i, , drop = FALSE],
      y = segments$row_y[[i]],
      gene_height = gene_height,
      arrow_head_bp = arrow_head_bp,
      arrow_head_fraction = arrow_head_fraction
    )
  })
  gene_polys <- do.call(rbind, polys)
  gene_polys$row_index <- segments$row_index[match(gene_polys$feature_id, segments$feature_id)]

  list(
    polygons = gene_polys,
    segment = segment,
    rows = rows,
    row_spacing = row_spacing,
    row_start = row_start,
    row_end = row_end,
    row_y = row_y,
    row_index = row_index
  )
}

linear_label_effective_line_angle <- function(line_angle, text_angle) {
  line_angle <- as.numeric(line_angle[[1L]])
  text_angle <- as.numeric(text_angle[[1L]])
  if (abs(text_angle %% 180) > 1e-10) {
    text_angle
  } else {
    line_angle
  }
}

linear_label_leader_endpoints <- function(label_data, segment, y_span,
                                          line_angle, text_angle) {
  if (!nrow(label_data)) {
    return(label_data)
  }
  rotated_text <- abs(as.numeric(text_angle[[1L]]) %% 180) > 1e-10
  if (rotated_text && all(c("xmid_linear", "line_y0") %in% names(label_data))) {
    gene_y <- label_data$linear_label_gene_y %||% label_data$line_y0
    direction_x <- (label_data$label_x - label_data$xmid_linear) / segment
    direction_y <- (label_data$label_y - gene_y) / y_span
    line_angle <- atan2(abs(direction_y), pmax(abs(direction_x), 1e-12)) * 180 / pi
  } else {
    line_angle <- linear_label_effective_line_angle(line_angle, text_angle)
  }
  line_angle <- abs(line_angle %% 180)
  line_angle <- ifelse(line_angle > 90, 180 - line_angle, line_angle)
  line_angle <- pmin(pmax(line_angle, 5), 90)
  line_theta <- line_angle * pi / 180
  text_theta <- abs(as.numeric(text_angle[[1L]]) %% 180) * pi / 180

  # Find the first intersection of the connector ray with the text rectangle.
  # A projection of the rectangle is too long for diagonal leaders and leaves
  # a visible gap below horizontal labels. The intersection also works when
  # the text itself is rotated because the ray is evaluated in text coordinates.
  text_half_width_norm <- label_data$label_text_width / segment / 2
  text_half_height_norm <- label_data$label_text_height / y_span / 2
  ray_in_text_x <- abs(cos(line_theta - text_theta))
  ray_in_text_y <- abs(sin(line_theta - text_theta))
  edge_distance <- pmin(
    text_half_width_norm / pmax(ray_in_text_x, 1e-12),
    text_half_height_norm / pmax(ray_in_text_y, 1e-12)
  )
  edge_x <- edge_distance * cos(line_theta) * segment
  edge_y <- edge_distance * sin(line_theta) * y_span
  has_leader <- if ("linear_label_has_leader" %in% names(label_data)) {
    as.logical(label_data$linear_label_has_leader)
  } else {
    label_data$label_lane > 0L
  }
  has_leader[is.na(has_leader)] <- FALSE
  label_side <- label_data$linear_label_side %||% rep(1, nrow(label_data))
  vertical_side <- label_data$linear_label_vertical_sign %||%
    rep(1, nrow(label_data))
  vertical_side <- ifelse(vertical_side < 0, -1, 1)
  label_data$line_x1 <- ifelse(
    has_leader,
    label_data$label_x - label_side * edge_x,
    NA_real_
  )
  label_data$line_y1 <- ifelse(
    has_leader,
    label_data$label_y - vertical_side * edge_y,
    NA_real_
  )
  label_data
}

linear_label_text_anchor <- function(label_data, segment, y_span, text_angle) {
  if (!nrow(label_data)) {
    return(label_data)
  }
  text_angle_value <- as.numeric(text_angle[[1L]])
  text_theta <- (text_angle_value %% 360) * pi / 180
  # Keep the old centered anchor for horizontal labels. For rotated labels,
  # place the text's left edge on the calculated gene-to-label ray.
  if (abs(text_angle_value %% 180) < 1e-10) {
    label_data$text_x <- label_data$label_x
    label_data$text_y <- label_data$label_y
  } else {
    label_data$text_x <- label_data$label_x -
      0.5 * label_data$label_text_width * cos(text_theta)
    label_data$text_y <- label_data$label_y -
      0.5 * label_data$label_text_width * sin(text_theta) * y_span / segment
  }
  label_data
}

linear_label_boxes_overlap <- function(candidate_x, candidate_y,
                                       candidate_width, candidate_height,
                                       candidate_text_width,
                                       candidate_text_height, placed,
                                       segment, y_span, text_angle,
                                       padding_x, padding_y) {
  if (!nrow(placed)) {
    return(FALSE)
  }
  text_angle_value <- as.numeric(text_angle[[1L]])
  rotated <- abs(text_angle_value %% 180) > 1e-10
  if (!rotated) {
    frame <- data.frame(
      xmin = candidate_x - candidate_width / 2,
      xmax = candidate_x + candidate_width / 2,
      ymin = candidate_y - candidate_height / 2 - padding_y,
      ymax = candidate_y + candidate_height / 2 + padding_y
    )
    placed_xmin <- placed$center_x - placed$box_width / 2
    placed_xmax <- placed$center_x + placed$box_width / 2
    placed_ymin <- placed$center_y - placed$box_height / 2 - padding_y
    placed_ymax <- placed$center_y + placed$box_height / 2 + padding_y
    return(any(
      frame$xmin < placed_xmax + padding_x &
        frame$xmax > placed_xmin - padding_x &
        frame$ymin < placed_ymax &
        frame$ymax > placed_ymin
    ))
  }

  # All linear labels share one text angle, so rotated rectangles can be
  # checked on the text axis and its perpendicular instead of their larger
  # axis-aligned bounding boxes.
  theta <- abs(text_angle_value %% 180) * pi / 180
  cos_theta <- cos(theta)
  sin_theta <- sin(theta)
  candidate_x <- candidate_x / segment
  candidate_y <- candidate_y / y_span
  placed_x <- placed$center_x / segment
  placed_y <- placed$center_y / y_span
  dx <- placed_x - candidate_x
  dy <- placed_y - candidate_y
  along <- abs(dx * cos_theta + dy * sin_theta)
  across <- abs(-dx * sin_theta + dy * cos_theta)
  candidate_half_width <- candidate_text_width / segment / 2
  candidate_half_height <- candidate_text_height / y_span / 2
  placed_half_width <- placed$text_width / segment / 2
  placed_half_height <- placed$text_height / y_span / 2
  padding <- max(padding_x / segment, padding_y / y_span)
  any(
    along < candidate_half_width + placed_half_width + padding &
      across < candidate_half_height + placed_half_height + padding
  )
}

linear_label_layout <- function(label_data, segment, rows = 4,
                                row_spacing = 3.4, gene_height = 0.24,
                                label_text_size = 3.3,
                                label_wrap_width = 18L,
                                label_max_lines = 2L,
                                label_offset = 0.14,
                                label_lane_gap = NULL,
                                label_lane_step = NULL,
                                line_angle = 90,
                                text_angle = 0,
                                min_gap_fraction = 0.012,
                                allow_gene_line_crossing = FALSE,
                                row_upper_limits = NULL,
                                row_lower_limits = NULL,
                                label_side = c("top", "bottom", "both")) {
  if (!nrow(label_data)) {
    return(label_data)
  }
  rows <- as.integer(rows[[1L]])
  row_spacing <- as.numeric(row_spacing[[1L]])
  if (!is.finite(rows) || rows < 1L || !is.finite(row_spacing) || row_spacing <= 0) {
    stop("Linear label layout requires positive `rows` and `row_spacing`.", call. = FALSE)
  }
  label_offset <- as.numeric(label_offset[[1L]])
  if (!is.finite(label_offset) || label_offset < 0) {
    stop("`linear_label_offset` must be a non-negative number.", call. = FALSE)
  }
  allow_gene_line_crossing <- as.logical(allow_gene_line_crossing[[1L]])
  if (is.na(allow_gene_line_crossing)) {
    stop("`linear_label_allow_gene_line_crossing` must be TRUE or FALSE.", call. = FALSE)
  }
  if (is.null(row_upper_limits)) {
    row_upper_limits <- rep(Inf, rows)
  } else {
    row_upper_limits <- as.numeric(row_upper_limits)
    if (length(row_upper_limits) != rows) {
      stop("`row_upper_limits` must have one value per linear row.", call. = FALSE)
    }
    row_upper_limits[!is.finite(row_upper_limits)] <- Inf
  }
  if (is.null(row_lower_limits)) {
    row_lower_limits <- rep(Inf, rows)
  } else {
    row_lower_limits <- as.numeric(row_lower_limits)
    if (length(row_lower_limits) != rows) {
      stop("`row_lower_limits` must have one value per linear row.", call. = FALSE)
    }
    row_lower_limits[!is.finite(row_lower_limits)] <- Inf
  }
  label_side_mode <- match.arg(label_side)

  label_text <- wrap_label_text_max_lines(
    label_data$label_display %||% label_data$label,
    width = label_wrap_width,
    max_lines = label_max_lines
  )
  dims <- label_box_dimensions(label_text, label_text_size)

  text_theta <- abs(as.numeric(text_angle[[1L]]) %% 180) * pi / 180
  rotated_width <- dims$width * abs(cos(text_theta)) + dims$height * abs(sin(text_theta))
  rotated_height <- dims$width * abs(sin(text_theta)) + dims$height * abs(cos(text_theta))
  y_span <- max(row_spacing * max(rows - 1L, 1L), 1)
  width_bp <- pmax(rotated_width * segment * 1.05, segment * 0.014)
  height_y <- pmax(rotated_height * y_span * 0.65, 0.30)
  # The padding is part of the collision box and keeps wrapped text clear of
  # the arrow edge, including when a label has two lines.
  label_padding_y <- max(0.03, gene_height * 0.08)

  if (!is.null(label_lane_step)) {
    label_lane_gap <- as.numeric(label_lane_step[[1L]])
  } else if (is.null(label_lane_gap)) {
    label_lane_gap <- max(0.25, stats::median(height_y, na.rm = TRUE) + 0.02)
  }
  label_lane_gap <- as.numeric(label_lane_gap[[1L]])
  if (!is.finite(label_lane_gap) || label_lane_gap <= 0) {
    stop("`linear_label_lane_step` must be a positive number.", call. = FALSE)
  }

  min_gap_bp <- segment * min_gap_fraction
  label_padding_x <- min_gap_bp
  line_angle <- linear_label_effective_line_angle(line_angle, text_angle)
  line_angle <- abs(line_angle %% 180)
  line_angle <- if (line_angle > 90) 180 - line_angle else line_angle
  line_angle <- if (line_angle >= 90) 90 else max(line_angle, 5)
  line_theta <- line_angle * pi / 180
  # `label_lane_gap` is the distance travelled along the connector. Its
  # vertical component changes with the requested angle, so text layers do not
  # become artificially far apart when the connector is diagonal.
  lane_y_step <- if (line_angle == 90) {
    label_lane_gap
  } else {
    label_lane_gap * sin(line_theta)
  }
  x_per_y <- segment / y_span
  # Linear leaders use one consistent horizontal direction: to the right.
  # The vertical sign is selected independently when both label sides are used.
  horizontal_side <- rep(1, nrow(label_data))

  # Convert a vertical rise in plot coordinates into the horizontal distance
  # needed to keep each connector at the requested angle.
  label_x_for_y <- function(i, y, side = 1) {
    if (line_angle == 90) {
      return(label_data$xmid_linear[[i]])
    }
    gene_y <- if (side < 0) {
      label_data$line_y_bottom[[i]]
    } else {
      label_data$line_y0[[i]]
    }
    rise <- abs(y - gene_y)
    label_data$xmid_linear[[i]] + horizontal_side[[i]] *
      rise * x_per_y / tan(line_theta)
  }

  label_data$label <- label_text
  label_data$label_lane <- 0L
  label_data$line_y0 <- label_data$row_y + gene_height / 2
  label_data$line_y_bottom <- label_data$row_y - gene_height / 2
  label_data$label_text_width <- dims$width * segment * 1.05
  label_data$label_text_height <- dims$height * y_span * 0.65
  label_data$label_box_width <- width_bp
  label_data$label_box_height <- height_y
  initial_vertical_sign <- if (identical(label_side_mode, "bottom")) -1 else 1
  label_data$linear_label_vertical_sign <- rep(
    initial_vertical_sign,
    nrow(label_data)
  )
  label_data$label_y <- if (initial_vertical_sign > 0) {
    label_data$line_y0 + label_offset +
      height_y / 2 + label_padding_y
  } else {
    label_data$line_y_bottom - label_offset - height_y / 2 - label_padding_y
  }
  label_data$label_x <- vapply(
    seq_len(nrow(label_data)),
    function(i) label_x_for_y(
      i,
      label_data$label_y[[i]],
      initial_vertical_sign
    ),
    numeric(1L)
  )
  label_data$linear_label_initial_x <- label_data$label_x
  label_data$linear_label_initial_y <- label_data$label_y
  label_data$line_y1 <- label_data$label_y
  removed_indices <- integer()
  rotated_text <- abs(as.numeric(text_angle[[1L]]) %% 180) > 1e-10
  if (rotated_text && !identical(label_side_mode, "top")) {
    stop(
      "`linear_label_side = 'bottom'` or `'both'` currently requires `linear_label_text_angle = 0`.",
      call. = FALSE
    )
  }
  horizontal_limit <- segment * 1.025 + max(width_bp, na.rm = TRUE)

  placed_global <- data.frame(
    center_x = numeric(), center_y = numeric(),
    box_width = numeric(), box_height = numeric(),
    text_width = numeric(), text_height = numeric()
  )
  row_order <- sort(unique(label_data$row_index))
  for (row in row_order) {
    placed <- placed_global
    idx <- which(label_data$row_index == row)
    current_row_y <- label_data$row_y[idx[[1L]]]
    # Place labels in the same left-to-right order as the genes in this row.
    # The original row order is only a deterministic tie-break for identical centers.
    idx <- idx[order(label_data$xmid_linear[idx], idx)]
    upper_limit <- if (row > 1L) {
      above_row_y <- label_data$row_y[label_data$row_index == row - 1L]
      if (length(above_row_y)) {
        above_row_y[[1L]] - gene_height / 2 - label_padding_y
      } else {
        Inf
      }
    } else {
      Inf
    }
    if (!allow_gene_line_crossing && is.finite(row_upper_limits[[row]])) {
      upper_limit <- min(upper_limit, row_upper_limits[[row]])
    }
    if (identical(label_side_mode, "both") && row > 1L) {
      above_row_y <- label_data$row_y[label_data$row_index == row - 1L]
      if (length(above_row_y)) {
        # Keep the top labels of this row below the midpoint between rows.
        # The matching bottom corridor is applied below, so connectors from
        # neighboring rows cannot cross in the same vertical space.
        upper_limit <- min(
          upper_limit,
          (above_row_y[[1L]] + current_row_y) / 2
        )
      }
    }

    if (rotated_text) {
      # Rotated labels are packed horizontally first. Only after the current
      # height is full do they move to a higher layer.
      vertical_step <- label_lane_gap
      for (i in idx) {
        base_y <- label_data$line_y0[[i]] + label_offset
        initial_y <- base_y + height_y[[i]] / 2 + label_padding_y
        text_theta <- (as.numeric(text_angle[[1L]]) %% 360) * pi / 180
        initial_x <- label_data$xmid_linear[[i]] +
          0.5 * label_data$label_text_width[[i]] * cos(text_theta)
        label_data$linear_label_initial_x[[i]] <- initial_x
        label_data$linear_label_initial_y[[i]] <- initial_y
        found <- FALSE
        lane <- 0L
        candidate_x <- initial_x
        candidate_y <- initial_y
        frame <- NULL

        try_horizontal <- function(y, x) {
          repeat {
            candidate_frame <- data.frame(
              xmin = x - width_bp[[i]] / 2,
              xmax = x + width_bp[[i]] / 2,
              ymin = y - height_y[[i]] / 2 - label_padding_y,
              ymax = y + height_y[[i]] / 2 + label_padding_y
            )
            overlap <- linear_label_boxes_overlap(
              candidate_x = x,
              candidate_y = y,
              candidate_width = width_bp[[i]],
              candidate_height = height_y[[i]],
              candidate_text_width = label_data$label_text_width[[i]],
              candidate_text_height = label_data$label_text_height[[i]],
              placed = placed,
              segment = segment,
              y_span = y_span,
              text_angle = text_angle,
              padding_x = label_padding_x,
              padding_y = label_padding_y
            )
            if (!overlap) {
              if (candidate_frame$ymax <= upper_limit + 1e-8 &&
                  candidate_frame$xmax <= horizontal_limit + 1e-8) {
                return(list(found = TRUE, x = x, y = y, frame = candidate_frame))
              }
              return(list(found = FALSE))
            }
            placed_xmin <- placed$center_x - placed$box_width / 2
            placed_xmax <- placed$center_x + placed$box_width / 2
            placed_ymin <- placed$center_y - placed$box_height / 2 - label_padding_y
            placed_ymax <- placed$center_y + placed$box_height / 2 + label_padding_y
            overlapping <- placed_xmin < x + width_bp[[i]] / 2 + label_padding_x &
              placed_xmax > x - width_bp[[i]] / 2 - label_padding_x &
              placed_ymin < y + height_y[[i]] / 2 + label_padding_y &
              placed_ymax > y - height_y[[i]] / 2 - label_padding_y
            overlapping_right <- if (any(overlapping)) {
              max(placed$center_x[overlapping] + placed$box_width[overlapping] / 2)
            } else {
              max(placed$center_x + placed$box_width / 2)
            }
            x_next <- max(
              x + label_padding_x,
              overlapping_right + label_padding_x + width_bp[[i]] / 2
            )
            if (!is.finite(x_next) || x_next <= x + 1e-10 ||
                x_next + width_bp[[i]] / 2 > horizontal_limit + 1e-8) {
              return(list(found = FALSE))
            }
            x <- x_next
          }
        }

        horizontal_result <- try_horizontal(candidate_y, candidate_x)
        if (isTRUE(horizontal_result$found)) {
          found <- TRUE
          candidate_x <- horizontal_result$x
          candidate_y <- horizontal_result$y
          frame <- horizontal_result$frame
        }

        if (!found) {
          lane <- 1L
          repeat {
            if (lane > max(100L, nrow(label_data) * 2L + 10L)) {
              removed_indices <- c(removed_indices, i)
              break
            }
            candidate_y <- initial_y + lane * vertical_step
            horizontal_result <- try_horizontal(candidate_y, initial_x)
            if (isTRUE(horizontal_result$found)) {
              found <- TRUE
              candidate_x <- horizontal_result$x
              frame <- horizontal_result$frame
              break
            }
            if (is.finite(upper_limit) &&
                candidate_y + height_y[[i]] / 2 + label_padding_y >
                  upper_limit + 1e-8 && !allow_gene_line_crossing) {
              removed_indices <- c(removed_indices, i)
              break
            }
            if (is.finite(upper_limit) &&
                candidate_y + height_y[[i]] / 2 + label_padding_y >
                  upper_limit + 1e-8 && allow_gene_line_crossing) {
              # Continue above the previous gene row only when explicitly allowed.
              lane <- lane + 1L
              next
            }
            lane <- lane + 1L
          }
        }

        if (!found) {
          next
        }
        label_data$label_lane[[i]] <- lane
        label_data$label_x[[i]] <- candidate_x
        label_data$label_y[[i]] <- candidate_y
        placed <- rbind(
          placed,
          data.frame(
            center_x = candidate_x,
            center_y = candidate_y,
            box_width = width_bp[[i]],
            box_height = height_y[[i]],
            text_width = label_data$label_text_width[[i]],
            text_height = label_data$label_text_height[[i]]
          )
        )
      }
      placed_global <- placed
      next
    }

    side_counts <- c(top = 0L, bottom = 0L)
    side_options <- switch(
      label_side_mode,
      top = 1L,
      bottom = -1L,
      both = c(1L, -1L)
    )
    candidate_for_side <- function(i, side) {
      top_base_y <- label_data$line_y0[[i]] + label_offset
      bottom_base_y <- label_data$line_y_bottom[[i]] - label_offset
      bottom_ceiling <- if (is.finite(row_lower_limits[[row]])) {
        row_lower_limits[[row]] - label_padding_y
      } else {
        label_data$row_y[[i]] - gene_height / 2 - label_padding_y
      }
      bottom_floor <- if (row < rows) {
        next_row_y <- label_data$row_y[label_data$row_index == row + 1L]
        if (length(next_row_y)) {
          next_gene_floor <- next_row_y[[1L]] + gene_height / 2 + label_padding_y
          if (identical(label_side_mode, "both")) {
            max(
              next_gene_floor,
              (current_row_y + next_row_y[[1L]]) / 2
            )
          } else {
            next_gene_floor
          }
        } else {
          -Inf
        }
      } else {
        -Inf
      }
      if (side > 0) {
        max_lane <- if (is.finite(upper_limit)) {
          floor((upper_limit - top_base_y - height_y[[i]] / 2 - label_padding_y) /
            lane_y_step)
        } else {
          max(100L, nrow(label_data) * 2L + 10L)
        }
        candidate_lanes <- if (max_lane >= 0) {
          seq.int(0L, max(0L, max_lane))
        } else {
          integer()
        }
      } else {
        min_lane <- max(
          0L,
          ceiling((bottom_base_y + height_y[[i]] / 2 + label_padding_y -
            bottom_ceiling) / lane_y_step)
        )
        max_lane <- if (is.finite(bottom_floor)) {
          floor((bottom_base_y - height_y[[i]] / 2 - label_padding_y -
            bottom_floor) / lane_y_step)
        } else {
          max(100L, nrow(label_data) * 2L + 10L)
        }
        candidate_lanes <- if (max_lane >= min_lane) {
          seq.int(min_lane, max_lane)
        } else {
          integer()
        }
      }
      for (lane in candidate_lanes) {
        candidate_y <- if (side > 0) {
          top_base_y + height_y[[i]] / 2 + label_padding_y +
            lane * lane_y_step
        } else {
          bottom_base_y - height_y[[i]] / 2 - label_padding_y -
            lane * lane_y_step
        }
        candidate_x <- label_x_for_y(i, candidate_y, side)
        frame <- data.frame(
          xmin = candidate_x - width_bp[[i]] / 2,
          xmax = candidate_x + width_bp[[i]] / 2,
          ymin = candidate_y - height_y[[i]] / 2 - label_padding_y,
          ymax = candidate_y + height_y[[i]] / 2 + label_padding_y
        )
        overlap <- linear_label_boxes_overlap(
          candidate_x = candidate_x,
          candidate_y = candidate_y,
          candidate_width = width_bp[[i]],
          candidate_height = height_y[[i]],
          candidate_text_width = label_data$label_text_width[[i]],
          candidate_text_height = label_data$label_text_height[[i]],
          placed = placed,
          segment = segment,
          y_span = y_span,
          text_angle = text_angle,
          padding_x = label_padding_x,
          padding_y = label_padding_y
        )
        within_bounds <- if (side > 0) {
          frame$ymax <= upper_limit + 1e-8
        } else {
          frame$ymax <= bottom_ceiling + 1e-8 &&
            frame$ymin >= bottom_floor - 1e-8
        }
        if (!overlap && within_bounds) {
          return(list(
            found = TRUE,
            side = side,
            lane = lane,
            x = candidate_x,
            y = candidate_y,
            frame = frame,
            base_x = label_x_for_y(
              i,
              if (side > 0) {
                top_base_y + height_y[[i]] / 2 + label_padding_y
              } else {
                bottom_base_y - height_y[[i]] / 2 - label_padding_y
              },
              side
            ),
            base_y = if (side > 0) {
              top_base_y + height_y[[i]] / 2 + label_padding_y
            } else {
              bottom_base_y - height_y[[i]] / 2 - label_padding_y
            }
          ))
        }
      }
      list(found = FALSE)
    }

    for (i in idx) {
      candidates <- lapply(side_options, function(side) {
        candidate_for_side(i, side)
      })
      found_candidates <- vapply(candidates, function(candidate) {
        isTRUE(candidate$found)
      }, logical(1L))
      if (!any(found_candidates)) {
        if (!allow_gene_line_crossing && identical(label_side_mode, "top")) {
          removed_indices <- c(removed_indices, i)
        }
        next
      }
      candidate_indices <- which(found_candidates)
      if (length(candidate_indices) > 1L) {
        candidate_score <- vapply(candidates[candidate_indices], function(candidate) {
          side_name <- if (candidate$side > 0) "top" else "bottom"
          candidate$lane + side_counts[[side_name]] * 0.001
        }, numeric(1L))
        chosen_index <- candidate_indices[which.min(candidate_score)]
      } else {
        chosen_index <- candidate_indices[[1L]]
      }
      chosen <- candidates[[chosen_index]]
      side_name <- if (chosen$side > 0) "top" else "bottom"
      side_counts[[side_name]] <- side_counts[[side_name]] + 1L
      label_data$label_lane[[i]] <- chosen$lane
      label_data$label_x[[i]] <- chosen$x
      label_data$label_y[[i]] <- chosen$y
      label_data$linear_label_initial_x[[i]] <- chosen$base_x
      label_data$linear_label_initial_y[[i]] <- chosen$base_y
      label_data$linear_label_vertical_sign[[i]] <- chosen$side
      placed <- rbind(
        placed,
        data.frame(
          center_x = chosen$x,
          center_y = chosen$y,
          box_width = width_bp[[i]],
          box_height = height_y[[i]],
          text_width = label_data$label_text_width[[i]],
          text_height = label_data$label_text_height[[i]]
        )
      )
    }
    placed_global <- placed
  }

  removed_labels <- if (length(removed_indices)) {
    as.character((label_data$label_display %||% label_data$label)[removed_indices])
  } else {
    character()
  }
  kept_indices <- setdiff(seq_len(nrow(label_data)), removed_indices)
  initial_x <- label_data$linear_label_initial_x
  initial_y <- label_data$linear_label_initial_y
  label_data <- label_data[kept_indices, , drop = FALSE]
  label_data$linear_label_base_x <- initial_x[kept_indices]
  label_data$linear_label_base_y <- initial_y[kept_indices]
  label_data$linear_label_side <- horizontal_side[kept_indices]
  label_data$linear_label_gene_y <- ifelse(
    label_data$linear_label_vertical_sign < 0,
    label_data$line_y_bottom,
    label_data$line_y0
  )
  label_data$linear_label_has_leader <- label_data$label_lane > 0L
  label_data <- linear_label_leader_endpoints(
    label_data,
    segment = segment,
    y_span = y_span,
    line_angle = line_angle,
    text_angle = text_angle
  )
  label_data <- linear_label_text_anchor(
    label_data,
    segment = segment,
    y_span = y_span,
    text_angle = text_angle
  )
  attr(label_data, "linear_label_summary") <- list(
    n_requested = nrow(label_data) + length(removed_indices),
    n_placed = nrow(label_data),
    n_removed = length(removed_indices),
    n_removed_cross_gene_line = length(removed_indices),
    removed_labels = removed_labels,
    allow_gene_line_crossing = allow_gene_line_crossing,
    side_requested = label_side_mode,
    n_top = sum(label_data$linear_label_vertical_sign > 0),
    n_bottom = sum(label_data$linear_label_vertical_sign < 0)
  )
  label_data
}

prepare_gc_skew_for_plot <- function(skew_table = NULL, fasta = NULL, sequence = NULL,
                                     genome_length, circular = TRUE,
                                     window = 100L, step = 1L) {
  if (!is.null(skew_table)) {
    skew <- read_gc_skew_table(skew_table)
  } else if (!is.null(fasta) || (is.character(sequence) && nzchar(sequence))) {
    skew <- compute_gc_skew(
      fasta = fasta,
      sequence = sequence,
      genome_length = genome_length,
      window = window,
      step = step,
      circular = circular
    )
  } else {
    return(NULL)
  }
  skew <- skew[is.finite(skew$position) & is.finite(skew$gc_skew), , drop = FALSE]
  skew <- skew[skew$position >= 1 & skew$position <= genome_length, , drop = FALSE]
  if (!nrow(skew)) {
    return(NULL)
  }
  max_abs <- max(abs(skew$gc_skew), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs <= 0) {
    max_abs <- 1
  }
  skew$scaled_score <- skew$gc_skew / max_abs
  if ("gc_content" %in% names(skew)) {
    gc_content <- as.numeric(skew$gc_content)
    mean_gc <- mean(gc_content, na.rm = TRUE)
    max_gc_delta <- max(abs(gc_content - mean_gc), na.rm = TRUE)
    if (!is.finite(max_gc_delta) || max_gc_delta <= 0) {
      max_gc_delta <- 1
    }
    skew$scaled_gc_content <- (gc_content - mean_gc) / max_gc_delta
  } else {
    skew$scaled_gc_content <- NA_real_
  }
  skew
}

nice_tick_step <- function(genome_length, target_ticks = 8) {
  raw <- genome_length / target_ticks
  base <- 10 ^ floor(log10(raw))
  candidates <- c(1, 2, 5, 10) * base
  candidates[[which.min(abs(candidates - raw))]]
}

circular_ruler_data <- function(genome_length, major_bp = NULL, minor_bp = NULL,
                                ruler_radius = 0.59,
                                ruler_minor_tick = 0.015,
                                ruler_major_tick = 0.035,
                                ruler_label_radius = 0.47) {
  if (is.null(major_bp) || !is.finite(major_bp)) {
    major_bp <- nice_tick_step(genome_length, target_ticks = 8)
  }
  if (is.null(minor_bp) || !is.finite(minor_bp)) {
    minor_bp <- major_bp / 5
  }

  major <- seq(0, genome_length - 1, by = major_bp)
  minor <- seq(0, genome_length - 1, by = minor_bp)
  minor <- setdiff(round(minor), round(major))

  major_labels <- if (major_bp >= 1000) {
    paste0(format(round(major / 1000, 1), trim = TRUE), " kb")
  } else {
    paste0(format(round(major), trim = TRUE), " bp")
  }
  major_labels[[1]] <- "0"

  list(
    ring = data.frame(
      x = seq(0, genome_length, length.out = 721),
      y = ruler_radius
    ),
    ticks = rbind(
      data.frame(
        x = minor,
        y = ruler_radius - ruler_minor_tick,
        yend = ruler_radius + ruler_minor_tick,
        type = "minor"
      ),
      data.frame(
        x = major,
        y = ruler_radius - ruler_major_tick,
        yend = ruler_radius + ruler_major_tick,
        type = "major"
      )
    ),
    labels = data.frame(
      x = major,
      y = ruler_label_radius,
      label = major_labels,
      theta = (major / genome_length) * 360
    )
  )
}

default_label_pattern <- function(category_scheme) {
  terms <- if (category_scheme == "phage") {
    c(
      "terminase", "portal", "capsid", "head", "tail", "fiber", "fibre",
      "baseplate", "holin", "lysin", "spanin", "integrase", "recombinase",
      "polymerase", "helicase", "nuclease", "tRNA"
    )
  } else {
    c(
      "NDM", "beta-lactamase", "CTX", "OXA", "TEM", "Tet\\(", "Qac",
      "AAC", "ANT", "AadA", "sulfonamide", "dihydrofolate", "Rep",
      "replication", "Tra", "conjugative transfer", "relaxase", "FinO",
      "ParA", "ParB", "partition", "integrase", "transposase",
      "\\bIS[0-9A-Za-z-]*\\b", "toxin", "antitoxin"
    )
  }
  paste(terms, collapse = "|")
}

select_label_data <- function(features, label_unknown = FALSE,
                              label_exclude_categories = NULL,
                              label_pattern = NULL, max_labels = Inf,
                              min_feature_bp = 1,
                              category_scheme = "plasmid",
                              genome_length = NULL,
                              label_min_gap_deg = 0) {
  label_data <- features[nzchar(features$label), , drop = FALSE]
  if (!label_unknown) {
    unknown_label <- label_data$category == "Unknown function, hypothetical protein" |
      is_unknown_feature_label(
        label_data$label_raw,
        label_data$label_display,
        label_data$gene_id
      )
    label_data <- label_data[
      !unknown_label,
      ,
      drop = FALSE
    ]
  }
  if (!is.null(label_exclude_categories) && length(label_exclude_categories)) {
    label_exclude_categories <- as.character(label_exclude_categories)
    label_exclude_categories <- label_exclude_categories[nzchar(label_exclude_categories)]
    if (length(label_exclude_categories)) {
      label_data <- label_data[
        !(label_data$category %in% label_exclude_categories),
        ,
        drop = FALSE
      ]
    }
  }
  if (!nrow(label_data)) {
    return(label_data)
  }

  label_data$bp <- abs(label_data$end - label_data$start) + 1
  label_data <- label_data[label_data$bp >= min_feature_bp, , drop = FALSE]
  if (!nrow(label_data)) {
    return(label_data)
  }

  pattern <- label_pattern %||% default_label_pattern(category_scheme)
  label_display <- label_data$label_display %||% label_data$label_raw
  label_data$pattern_hit <- if (is.null(pattern) || !nzchar(pattern)) {
    FALSE
  } else {
    grepl(pattern, label_data$label_raw, ignore.case = TRUE, perl = TRUE) |
      grepl(pattern, label_display, ignore.case = TRUE, perl = TRUE) |
      grepl(pattern, label_data$gene_id, ignore.case = TRUE, perl = TRUE)
  }

  priority <- c(
    "Antimicrobial resistance" = 1,
    "Replication" = 2,
    "Conjugation and transfer" = 3,
    "Partition and maintenance" = 4,
    "Mobile element" = 5,
    "Virulence" = 6,
    "Regulation" = 7,
    "Metabolism" = 8,
    "Other functions" = 9,
    "Unknown function, hypothetical protein" = 10,
    "DNA, RNA and nucleotide metabolism" = 3,
    "Head and packaging" = 3,
    "Tail" = 3,
    "Lysis" = 3,
    "Connectors" = 4,
    "Transcription regulation" = 7,
    "Moron, auxiliary metabolic gene and host takeover" = 8,
    "tRNA" = 8
  )
  label_data$priority <- unname(priority[label_data$category])
  label_data$priority[is.na(label_data$priority)] <- 9
  label_data$priority <- label_data$priority - ifelse(label_data$pattern_hit, 4, 0)
  label_data <- label_data[order(label_data$priority, -label_data$bp, label_data$xmid), , drop = FALSE]

  if (is.finite(label_min_gap_deg) && label_min_gap_deg > 0 &&
      is.finite(genome_length) && nrow(label_data) > 1L) {
    label_data$theta <- (label_data$xmid %% genome_length) / genome_length * 360
    keep <- rep(FALSE, nrow(label_data))
    kept_theta <- numeric()
    for (i in seq_len(nrow(label_data))) {
      distance <- if (!length(kept_theta)) {
        Inf
      } else {
        abs(((label_data$theta[[i]] - kept_theta + 180) %% 360) - 180)
      }
      if (all(distance >= label_min_gap_deg)) {
        keep[[i]] <- TRUE
        kept_theta <- c(kept_theta, label_data$theta[[i]])
      }
    }
    label_data <- label_data[keep, , drop = FALSE]
    label_data$theta <- NULL
  }

  if (is.finite(max_labels) && nrow(label_data) > max_labels) {
    label_data <- label_data[seq_len(max_labels), , drop = FALSE]
  }

  label_data[order(label_data$xmid), , drop = FALSE]
}

radial_npc <- function(x, y, genome_length, y_max = 1.95,
                       inner_radius = 0.38) {
  theta <- 2 * pi * as.numeric(x) / genome_length
  radius <- 0.4 * (
    inner_radius + (as.numeric(y) / y_max) * (1 - inner_radius)
  )
  data.frame(
    npcx = 0.5 + radius * sin(theta),
    npcy = 0.5 + radius * cos(theta)
  )
}

spread_label_positions <- function(y, lower = 0.08, upper = 0.92,
                                   min_gap = 0.055) {
  if (!length(y)) {
    return(numeric())
  }
  if (length(y) == 1L) {
    return(pmin(pmax(y, lower), upper))
  }

  ordering <- order(y)
  sorted <- pmin(pmax(y[ordering], lower), upper)
  gap <- min(min_gap, (upper - lower) / (length(sorted) - 1L))
  for (i in 2:length(sorted)) {
    sorted[[i]] <- max(sorted[[i]], sorted[[i - 1L]] + gap)
  }
  if (sorted[[length(sorted)]] > upper) {
    sorted <- sorted - (sorted[[length(sorted)]] - upper)
  }
  for (i in (length(sorted) - 1L):1L) {
    sorted[[i]] <- min(sorted[[i]], sorted[[i + 1L]] - gap)
  }
  if (sorted[[1]] < lower) {
    sorted <- sorted + (lower - sorted[[1]])
  }

  out <- numeric(length(y))
  out[ordering] <- sorted
  out
}

spread_label_boxes <- function(y, height, lower = 0.04, upper = 0.96,
                               padding = 0.008) {
  if (!length(y)) {
    return(numeric())
  }
  height <- pmax(as.numeric(height), 0.01)
  if (length(y) == 1L) {
    return(pmin(pmax(y, lower + height / 2), upper - height / 2))
  }

  ordering <- order(y)
  sorted <- y[ordering]
  sorted_height <- height[ordering]
  total_height <- sum(sorted_height) + padding * (length(sorted) - 1L)
  if (total_height >= upper - lower) {
    centers <- seq(lower, upper, length.out = length(sorted))
    out <- numeric(length(y))
    out[ordering] <- centers
    return(out)
  }

  sorted <- pmin(pmax(sorted, lower + sorted_height / 2), upper - sorted_height / 2)
  for (i in 2:length(sorted)) {
    sorted[[i]] <- max(
      sorted[[i]],
      sorted[[i - 1L]] + sorted_height[[i - 1L]] / 2 + sorted_height[[i]] / 2 + padding
    )
  }
  overflow <- sorted[[length(sorted)]] + sorted_height[[length(sorted)]] / 2 - upper
  if (overflow > 0) {
    sorted <- sorted - overflow
  }
  for (i in (length(sorted) - 1L):1L) {
    sorted[[i]] <- min(
      sorted[[i]],
      sorted[[i + 1L]] - sorted_height[[i + 1L]] / 2 - sorted_height[[i]] / 2 - padding
    )
  }
  underflow <- lower - (sorted[[1L]] - sorted_height[[1L]] / 2)
  if (underflow > 0) {
    sorted <- sorted + underflow
  }

  out <- numeric(length(y))
  out[ordering] <- sorted
  out
}

wrap_label_text <- function(text, width = 22L) {
  text <- clean_text(text)
  text <- gsub("[\r\n]+", " ", text)
  if (is.null(width) || !is.finite(width) || width <= 0) {
    return(text)
  }
  vapply(
    text,
    function(one) {
      paste(strwrap(one, width = width), collapse = "\n")
    },
    character(1)
  )
}

balanced_two_line_label <- function(text, width = 18L) {
  text <- clean_text(text)
  text <- gsub("[\r\n]+", " ", text)
  text <- gsub("[[:space:]]+", " ", text)
  text <- trimws(text)
  if (is.null(width) || !length(width) || !is.finite(width[[1L]]) || width[[1L]] <= 0) {
    return(text)
  }
  width <- as.integer(floor(as.numeric(width[[1L]])))
  if (nchar(text) <= width) {
    return(text)
  }

  characters <- strsplit(text, "", fixed = TRUE)[[1L]]
  split_positions <- seq_len(length(characters) - 1L)
  left <- vapply(
    split_positions,
    function(position) trimws(paste(characters[seq_len(position)], collapse = "")),
    character(1L)
  )
  right <- vapply(
    split_positions,
    function(position) paste(characters[(position + 1L):length(characters)], collapse = ""),
    character(1L)
  )
  right <- trimws(right)
  usable <- nzchar(left) & nzchar(right)
  if (!any(usable)) {
    return(text)
  }

  left_lengths <- nchar(left)
  right_lengths <- nchar(right)
  # Prefer a whitespace boundary when it can satisfy the requested width.
  whitespace_boundary <- vapply(
    split_positions,
    function(position) {
      characters[[position]] == " " || characters[[position + 1L]] == " "
    },
    logical(1L)
  )
  candidates <- if (any(usable & whitespace_boundary)) {
    # Keep ordinary phrases intact even when the requested width is too
    # narrow for both lines; only unspaced words use a character split.
    which(usable & whitespace_boundary)
  } else {
    # A single word can be longer than the requested width. Split it by
    # characters near the midpoint rather than leaving an unbreakable line.
    which(usable)
  }
  score <- abs(left_lengths - right_lengths)
  split_at <- candidates[[which.min(score[candidates])]]
  paste(left[[split_at]], right[[split_at]], sep = "\n")
}

wrap_label_text_max_lines <- function(text, width = 18L, max_lines = 2L) {
  text <- clean_text(text)
  max_lines <- as.integer(max_lines[[1L]])
  if (!is.finite(max_lines) || max_lines < 1L) {
    stop("`linear_label_max_lines` must be a positive integer.", call. = FALSE)
  }
  if (max_lines == 2L) {
    return(vapply(text, balanced_two_line_label, character(1), width = width))
  }
  text <- wrap_label_text(text, width = width)
  vapply(
    strsplit(text, "\n", fixed = TRUE),
    function(lines) {
      if (length(lines) <= max_lines) {
        return(paste(lines, collapse = "\n"))
      }
      if (max_lines == 1L) {
        return(paste(lines, collapse = " "))
      }
      paste(
        c(lines[seq_len(max_lines - 1L)], paste(lines[max_lines:length(lines)], collapse = " ")),
        collapse = "\n"
      )
    },
    character(1)
  )
}

max_label_line_chars <- function(text) {
  vapply(
    strsplit(as.character(text), "\n", fixed = TRUE),
    function(lines) max(nchar(lines, type = "chars"), 1),
    numeric(1)
  )
}

label_box_dimensions <- function(text, text_size, lineheight = 0.9) {
  line_count <- lengths(strsplit(as.character(text), "\n", fixed = TRUE))
  data.frame(
    width = 0.012 + max_label_line_chars(text) * text_size * 0.0024,
    height = 0.006 + pmax(line_count, 1L) * text_size * 0.0065 * lineheight
  )
}

label_box_overlaps <- function(candidate, placed, padding = 0.006) {
  if (!nrow(placed)) {
    return(FALSE)
  }
  any(
    candidate$xmin < placed$xmax + padding &
      candidate$xmax > placed$xmin - padding &
      candidate$ymin < placed$ymax + padding &
      candidate$ymax > placed$ymin - padding
  )
}

label_box_frame <- function(x, y, width, height, right = TRUE, hjust = NULL) {
  if (is.null(hjust)) {
    hjust <- if (isTRUE(right)) 0 else 1
  }
  data.frame(
    xmin = x - width * hjust,
    xmax = x + width * (1 - hjust),
    ymin = y - height / 2,
    ymax = y + height / 2
  )
}

label_adjust_column <- function(data, candidates) {
  matched <- candidates[candidates %in% names(data)]
  if (length(matched)) matched[[1L]] else NA_character_
}

normalise_label_adjust <- function(adjustment = NULL) {
  empty <- data.frame(
    label = character(),
    hjust = numeric(),
    vjust = numeric(),
    stringsAsFactors = FALSE
  )
  if (is.null(adjustment)) {
    return(empty)
  }

  if (inherits(adjustment, "ggplasmid_label_adjust")) {
    adjustment <- unclass(adjustment)
  }

  if (is.data.frame(adjustment)) {
    label_col <- label_adjust_column(
      adjustment,
      c("label", "label_display", "label_raw", "product", "gene_id", "gene", "name")
    )
    if (is.na(label_col)) {
      stop(
        "`label_adjust` data frames need a label/name column.",
        call. = FALSE
      )
    }
    h_col <- label_adjust_column(adjustment, c("hjust", "x", "dx", "x_shift"))
    v_col <- label_adjust_column(adjustment, c("vjust", "y", "dy", "y_shift"))
    out <- data.frame(
      label = as.character(adjustment[[label_col]]),
      hjust = if (is.na(h_col)) 0 else as.numeric(adjustment[[h_col]]),
      vjust = if (is.na(v_col)) 0 else as.numeric(adjustment[[v_col]]),
      stringsAsFactors = FALSE
    )
  } else if (is.list(adjustment)) {
    single_label_name <- intersect(
      c("label", "label_display", "label_raw", "product", "gene_id", "gene", "name"),
      names(adjustment)
    )
    if (length(single_label_name)) {
      out <- data.frame(
        label = as.character(adjustment[[single_label_name[[1L]]]]),
        hjust = as.numeric(adjustment[["hjust"]] %||% 0),
        vjust = as.numeric(adjustment[["vjust"]] %||% 0),
        stringsAsFactors = FALSE
      )
    } else if (!is.null(names(adjustment)) && all(nzchar(names(adjustment)))) {
      pieces <- lapply(names(adjustment), function(label) {
        value <- adjustment[[label]]
        if (is.list(value)) {
          data.frame(
            label = label,
            hjust = as.numeric(value[["hjust"]] %||% value[["x"]] %||% 0),
            vjust = as.numeric(value[["vjust"]] %||% value[["y"]] %||% 0),
            stringsAsFactors = FALSE
          )
        } else {
          value <- unlist(value)
          data.frame(
            label = label,
            hjust = as.numeric(value[["hjust"]] %||% value[["x"]] %||% 0),
            vjust = as.numeric(value[["vjust"]] %||% value[["y"]] %||% 0),
            stringsAsFactors = FALSE
          )
        }
      })
      out <- do.call(rbind, pieces)
    } else {
      pieces <- lapply(adjustment, normalise_label_adjust)
      pieces <- pieces[vapply(pieces, nrow, integer(1L)) > 0L]
      out <- if (length(pieces)) do.call(rbind, pieces) else empty
    }
  } else {
    stop(
      "`label_adjust` must be a data frame, list, or label_adjust() result.",
      call. = FALSE
    )
  }

  out$label <- trimws(out$label)
  out$hjust[!is.finite(out$hjust)] <- 0
  out$vjust[!is.finite(out$vjust)] <- 0
  out[nzchar(out$label), , drop = FALSE]
}

#' Manually shift selected circular labels
#'
#' @param label Label text to match. The match is checked against the plotted
#'   label, raw label, display label, and wrapped label text.
#' @param hjust Horizontal shift in plot-panel units. Positive values move the
#'   label to the right.
#' @param vjust Vertical shift in plot-panel units. Positive values move the
#'   label upward.
#' @return A data frame usable with the `label_adjust` argument in
#'   [ggplasmid()] or [plot_phage_map()].
#' @export
label_adjust <- function(label, hjust = 0, vjust = 0) {
  if (is.data.frame(label) || is.list(label)) {
    out <- normalise_label_adjust(label)
  } else {
    out <- normalise_label_adjust(
      data.frame(
        label = label,
        hjust = hjust,
        vjust = vjust,
        stringsAsFactors = FALSE
      )
    )
  }
  class(out) <- c("ggplasmid_label_adjust", class(out))
  out
}

normalise_label_key <- function(x) {
  x <- gsub("\\s+", " ", as.character(x))
  trimws(x)
}

apply_label_adjust <- function(label_position, label_data, text,
                               adjustment = NULL) {
  adjustments <- normalise_label_adjust(adjustment)
  if (!nrow(adjustments) || !nrow(label_position)) {
    return(label_position)
  }

  candidates <- data.frame(
    label_display = normalise_label_key(label_data$label_display %||% ""),
    label_raw = normalise_label_key(label_data$label_raw %||% ""),
    label = normalise_label_key(label_data$label %||% ""),
    gene_id = normalise_label_key(label_data$gene_id %||% ""),
    gene = normalise_label_key(label_data$gene %||% ""),
    text = normalise_label_key(text),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(adjustments))) {
    target <- normalise_label_key(adjustments$label[[i]])
    match_index <- which(
      candidates$label_display == target |
      candidates$label_raw == target |
      candidates$label == target |
      candidates$gene_id == target |
      candidates$gene == target |
      candidates$text == target
    )
    if (!length(match_index)) {
      next
    }
    label_position$npcx[match_index] <-
      label_position$npcx[match_index] + adjustments$hjust[[i]]
    label_position$npcy[match_index] <-
      label_position$npcy[match_index] + adjustments$vjust[[i]]
  }

  label_position
}

label_overlay_npc_bounds <- function(label_overlay, padding = 0.006) {
  if (is.null(label_overlay) || !nrow(label_overlay)) {
    return(NULL)
  }
  keep <- is.finite(label_overlay$npcx) & is.finite(label_overlay$npcy)
  if (!any(keep)) {
    return(NULL)
  }
  label_overlay <- label_overlay[keep, , drop = FALSE]
  dims <- label_box_dimensions(label_overlay$label, label_overlay$size)
  boxes <- label_box_frame(
    label_overlay$npcx,
    label_overlay$npcy,
    dims$width,
    dims$height,
    hjust = label_overlay$hjust
  )
  c(
    xmin = min(boxes$xmin, na.rm = TRUE) - padding,
    xmax = max(boxes$xmax, na.rm = TRUE) + padding,
    ymin = min(boxes$ymin, na.rm = TRUE) - padding,
    ymax = max(boxes$ymax, na.rm = TRUE) + padding
  )
}

auto_legend_plot_spacing <- function(spacing, position, label_bounds,
                                     cm_per_npc = 18, padding_cm = 0.8) {
  if (is.null(label_bounds) || !length(label_bounds)) {
    return(spacing)
  }
  overflow <- switch(
    position,
    right = max(0, unname(label_bounds[["xmax"]]) - 1),
    left = max(0, -unname(label_bounds[["xmin"]])),
    top = max(0, unname(label_bounds[["ymax"]]) - 1),
    left_top = max(0, -unname(label_bounds[["xmin"]])),
    right_top = max(0, unname(label_bounds[["xmax"]]) - 1),
    bottom = max(0, -unname(label_bounds[["ymin"]])),
    left_bottom = max(0, -unname(label_bounds[["xmin"]])),
    right_bottom = max(0, unname(label_bounds[["xmax"]]) - 1),
    0
  )
  if (!is.finite(overflow) || overflow <= 0) {
    return(spacing)
  }
  max(spacing, overflow * cm_per_npc + padding_cm)
}

circular_plot_margin <- function(label_bounds = NULL,
                                 base = c(top = 16, right = 28, bottom = 16, left = 28),
                                 overflow_scale = 420) {
  margin_values <- base
  if (!is.null(label_bounds) && length(label_bounds)) {
    overflow <- c(
      top = max(0, unname(label_bounds[["ymax"]]) - 1),
      right = max(0, unname(label_bounds[["xmax"]]) - 1),
      bottom = max(0, -unname(label_bounds[["ymin"]])),
      left = max(0, -unname(label_bounds[["xmin"]]))
    )
    margin_values <- margin_values + overflow * overflow_scale
  }

  ggplot2::margin(
    t = unname(margin_values[["top"]]),
    r = unname(margin_values[["right"]]),
    b = unname(margin_values[["bottom"]]),
    l = unname(margin_values[["left"]])
  )
}

pack_side_label_boxes <- function(y, height, lower = 0.015, upper = 0.985,
                                  padding = 0.006) {
  n <- length(y)
  if (!n) {
    return(list(y = numeric(), lane = integer()))
  }

  y <- as.numeric(y)
  y[!is.finite(y)] <- 0.5
  height <- pmax(as.numeric(height), 0.01)
  available <- upper - lower
  if (!is.finite(available) || available <= 0) {
    return(list(y = y, lane = integer(n)))
  }

  total_height <- function(idx) {
    if (!length(idx)) {
      return(0)
    }
    sum(height[idx]) + padding * max(length(idx) - 1L, 0L)
  }

  ordering <- order(y)
  column_count <- max(
    1L,
    ceiling(total_height(seq_len(n)) / available)
  )
  column_count <- min(column_count, n)

  repeat {
    lane_sorted <- (seq_along(ordering) - 1L) %% column_count
    overfull <- FALSE
    for (lane in 0:(column_count - 1L)) {
      idx <- ordering[lane_sorted == lane]
      if (total_height(idx) > available && column_count < n) {
        overfull <- TRUE
        break
      }
    }
    if (!overfull || column_count >= n) {
      break
    }
    column_count <- column_count + 1L
  }

  y_out <- numeric(n)
  lanes <- integer(n)
  lane_sorted <- (seq_along(ordering) - 1L) %% column_count
  for (lane in 0:(column_count - 1L)) {
    idx <- ordering[lane_sorted == lane]
    if (!length(idx)) {
      next
    }
    y_out[idx] <- spread_label_boxes(
      y[idx],
      height[idx],
      lower = lower,
      upper = upper,
      padding = padding
    )
    lanes[idx] <- lane
  }

  list(y = y_out, lane = lanes)
}

place_radial_label_boxes <- function(x, genome_length, text, right,
                                     inner_radius = 0.38,
                                     label_anchor_radius = 1.07,
                                     label_line_length = 0.02,
                                     label_lane_spacing = 0.24,
                                     label_lane_step = NULL,
                                     label_text_size = 3.4,
                                     priority = NULL,
                                     lower = 0.04, upper = 0.96,
                                     y_max = 1.95,
                                     gene_outer_radius = NULL) {
  n <- length(x)
  if (!n) {
    return(data.frame(npcx = numeric(), npcy = numeric(), lane = integer()))
  }

  dims <- label_box_dimensions(text, label_text_size)
  data_per_npc_radius <- y_max / (0.4 * (1 - inner_radius))
  lane_height <- stats::quantile(
    dims$height[is.finite(dims$height)],
    0.75,
    names = FALSE,
    na.rm = TRUE
  )
  if (!is.finite(lane_height)) {
    lane_height <- 0.02
  }
  if (is.null(label_lane_step)) {
    lane_step <- max(
      0.026,
      lane_height * data_per_npc_radius * 0.55,
      label_lane_spacing * 0.12
    )
  } else {
    lane_step <- as.numeric(label_lane_step)
    if (length(lane_step) != 1L || !is.finite(lane_step) || lane_step <= 0) {
      stop("`label_lane_step` must be a positive number.", call. = FALSE)
    }
  }
  padding <- 0.006

  right <- as.logical(right)
  right[is.na(right)] <- TRUE
  theta <- 2 * pi * as.numeric(x) / genome_length
  vertical_label <- abs(sin(theta)) < 0.20
  label_hjust <- ifelse(vertical_label, 0.5, ifelse(right, 0, 1))
  vertical_edge_gap <- 0.008
  npcx <- rep(NA_real_, n)
  npcy <- rep(NA_real_, n)
  lanes <- rep(NA_integer_, n)
  placed <- data.frame(
    xmin = numeric(), xmax = numeric(),
    ymin = numeric(), ymax = numeric()
  )

  label_order <- order(
    dims$width,
    dims$height,
    max_label_line_chars(text),
    x
  )
  remaining <- label_order
  lane <- 0L
  max_lanes <- max(32L, n * 4L)

  move_box_outside_gene_ring <- function(text_x, text_y, i) {
    if (is.null(gene_outer_radius) ||
        length(gene_outer_radius) != 1L ||
        !is.finite(gene_outer_radius)) {
      return(c(x = text_x, y = text_y))
    }
    ring <- radial_npc(
      x[[i]],
      y = gene_outer_radius,
      genome_length = genome_length,
      y_max = y_max,
      inner_radius = inner_radius
    )
    ring_radius <- sqrt(
      (ring$npcx[[1L]] - 0.5)^2 + (ring$npcy[[1L]] - 0.5)^2
    )
    theta_i <- theta[[i]]
    outward <- c(sin(theta_i), cos(theta_i))
    hjust_i <- label_hjust[[i]]
    box_distance <- function(delta) {
      box <- label_box_frame(
        text_x + outward[[1L]] * delta,
        text_y + outward[[2L]] * delta,
        dims$width[[i]],
        dims$height[[i]],
        hjust = hjust_i
      )
      corners <- rbind(
        c(box$xmin, box$ymin),
        c(box$xmin, box$ymax),
        c(box$xmax, box$ymin),
        c(box$xmax, box$ymax)
      )
      min(sqrt((corners[, 1] - 0.5)^2 + (corners[, 2] - 0.5)^2))
    }
    target <- ring_radius + 0.004
    if (box_distance(0) >= target) {
      return(c(x = text_x, y = text_y))
    }
    upper_delta <- 0.01
    while (box_distance(upper_delta) < target && upper_delta < 0.50) {
      upper_delta <- upper_delta * 2
    }
    if (box_distance(upper_delta) < target) {
      delta <- upper_delta
    } else {
      delta <- stats::uniroot(
        function(value) box_distance(value) - target,
        interval = c(0, upper_delta)
      )$root
    }
    c(
      x = text_x + outward[[1L]] * delta,
      y = text_y + outward[[2L]] * delta
    )
  }

  candidate_box <- function(i, lane) {
    position <- radial_npc(
      x[[i]],
      y = label_anchor_radius + label_line_length + lane * lane_step,
      genome_length = genome_length,
      y_max = y_max,
      inner_radius = inner_radius
    )
    vertical_shift <- if (vertical_label[[i]]) {
      dims$height[[i]] * 0.65 + vertical_edge_gap
    } else {
      0
    }
    text_x <- position$npcx[[1L]]
    text_y <- position$npcy[[1L]]
    if (vertical_shift > 0) {
      text_x <- text_x + sin(theta[[i]]) * vertical_shift
      text_y <- text_y + cos(theta[[i]]) * vertical_shift
    }
    shifted <- move_box_outside_gene_ring(text_x, text_y, i)
    text_x <- shifted[["x"]]
    text_y <- shifted[["y"]]
    list(
      x = text_x,
      y = text_y,
      box = label_box_frame(
        text_x,
        text_y,
        dims$width[[i]],
        dims$height[[i]],
        hjust = label_hjust[[i]]
      )
    )
  }

  while (length(remaining) && lane <= max_lanes) {
    next_remaining <- integer()
    for (i in remaining[order(dims$width[remaining], dims$height[remaining], x[remaining])]) {
      candidate <- candidate_box(i, lane)
      if (!label_box_overlaps(candidate$box, placed, padding = padding)) {
        npcx[[i]] <- candidate$x
        npcy[[i]] <- candidate$y
        lanes[[i]] <- lane
        placed <- rbind(placed, candidate$box)
      } else {
        next_remaining <- c(next_remaining, i)
      }
    }
    remaining <- next_remaining
    lane <- lane + 1L
  }

  for (i in remaining) {
    repeat {
      candidate <- candidate_box(i, lane)
      if (!label_box_overlaps(candidate$box, placed, padding = padding) ||
          lane > max_lanes + n * 8L) {
        npcx[[i]] <- candidate$x
        npcy[[i]] <- candidate$y
        lanes[[i]] <- lane
        placed <- rbind(placed, candidate$box)
        lane <- lane + 1L
        break
      }
      lane <- lane + 1L
    }
  }

  data.frame(npcx = npcx, npcy = npcy, lane = lanes, hjust = label_hjust)
}

radial_label_lanes <- function(theta, text, min_gap_deg = 5,
                               char_width_deg = 0.35, max_lanes = Inf) {
  if (!length(theta)) {
    return(integer())
  }
  if (length(theta) == 1L) {
    return(0L)
  }

  theta <- theta %% 360
  sorted <- sort(theta)
  gaps <- diff(c(sorted, sorted[[1L]] + 360))
  start_angle <- sorted[(which.max(gaps) %% length(sorted)) + 1L]
  unwrapped <- ifelse(theta < start_angle, theta + 360, theta)

  line_count <- lengths(strsplit(as.character(text), "\n", fixed = TRUE))
  span <- pmax(
    5,
    pmin(70, (max_label_line_chars(text) + pmax(line_count - 1L, 0L) * 6) * char_width_deg)
  )
  label_start <- unwrapped - span / 2
  label_end <- unwrapped + span / 2

  if (!is.finite(max_lanes)) {
    max_lanes <- length(theta)
  }
  max_lanes <- max(1L, as.integer(max_lanes))
  lanes <- integer(length(theta))
  lane_end <- rep(-Inf, max_lanes)
  for (i in order(unwrapped)) {
    available <- which(label_start[[i]] > lane_end + min_gap_deg)
    lane <- if (length(available)) {
      available[[1L]]
    } else {
      which.min(lane_end)
    }
    lanes[[i]] <- lane - 1L
    lane_end[[lane]] <- label_end[[i]]
  }
  lanes
}

screen_label_lanes <- function(y, height, right, padding = 0.006,
                               max_lanes = 3L) {
  if (!length(y)) {
    return(integer())
  }
  lanes <- integer(length(y))
  max_lanes <- max(1L, as.integer(max_lanes))

  for (side in c(FALSE, TRUE)) {
    index <- which(right == side)
    if (!length(index)) {
      next
    }
    ordering <- index[order(y[index])]
    lane_end <- rep(-Inf, max_lanes)
    for (i in ordering) {
      label_start <- y[[i]] - height[[i]] / 2
      label_end <- y[[i]] + height[[i]] / 2
      available <- which(label_start > lane_end + padding)
      lane <- if (length(available)) {
        available[[1L]]
      } else {
        which.min(lane_end)
      }
      lanes[[i]] <- lane - 1L
      lane_end[[lane]] <- max(lane_end[[lane]], label_end)
    }
  }

  lanes
}

GeomPlasmidTextOverlay <- ggplot2::ggproto(
  "GeomPlasmidTextOverlay",
  ggplot2::Geom,
  required_aes = c("label"),
  default_aes = ggplot2::aes(
    npcx = NA_real_,
    npcy = NA_real_,
    panel_x = NA_real_,
    panel_y = NA_real_,
    angle = 0,
    hjust = 0.5,
    vjust = 0.5,
    anchor_x = NA_real_,
    anchor_y = NA_real_,
    leader_x = NA_real_,
    leader_y = NA_real_,
    x0 = NA_real_,
    y0 = NA_real_,
    x1 = NA_real_,
    y1 = NA_real_,
    x2 = NA_real_,
    y2 = NA_real_,
    line_colour = NA_character_,
    label_colour = NA_character_
  ),
  draw_key = function(...) grid::nullGrob(),
  draw_panel = function(data, panel_params, coord, family = "sans",
                        text_colour = "black", text_size = 3,
                        fontface = "plain", lineheight = 0.9,
                        segment_colour = "grey45",
                        segment_linewidth = 0.25,
                        segment_linetype = "solid", na.rm = FALSE) {
    grobs <- list()

    text_x <- data$npcx
    text_y <- data$npcy
    has_panel_text <- is.finite(data$panel_x) & is.finite(data$panel_y)
    if (any(has_panel_text)) {
      panel_text <- coord$transform(
        data.frame(
          x = data$panel_x[has_panel_text],
          y = data$panel_y[has_panel_text]
        ),
        panel_params
      )
      text_x[has_panel_text] <- panel_text$x
      text_y[has_panel_text] <- panel_text$y
    }

    line_x0 <- data$x0
    line_y0 <- data$y0
    has_anchor <- is.finite(data$anchor_x) & is.finite(data$anchor_y)
    if (any(has_anchor)) {
      anchor <- coord$transform(
        data.frame(
          x = data$anchor_x[has_anchor],
          y = data$anchor_y[has_anchor]
        ),
        panel_params
      )
      line_x0[has_anchor] <- anchor$x
      line_y0[has_anchor] <- anchor$y
    }

    line_x1 <- data$x1
    line_y1 <- data$y1
    line_x2 <- data$x2
    line_y2 <- data$y2
    has_leader <- is.finite(data$leader_x) & is.finite(data$leader_y)
    if (any(has_leader)) {
      leader <- coord$transform(
        data.frame(
          x = data$leader_x[has_leader],
          y = data$leader_y[has_leader]
        ),
        panel_params
      )
      line_x1[has_leader] <- leader$x
      line_y1[has_leader] <- leader$y
      line_x2[has_leader] <- leader$x
      line_y2[has_leader] <- leader$y
    }
    use_text_endpoint <- has_anchor & is.finite(text_x) & is.finite(text_y) &
      !has_leader &
      (!is.finite(line_x1) | !is.finite(line_y1) |
        !is.finite(line_x2) | !is.finite(line_y2))
    if (any(use_text_endpoint)) {
      gap <- 0.008
      side_gap <- ifelse(
        data$hjust <= 0.25,
        gap,
        ifelse(data$hjust >= 0.75, -gap, 0)
      )
      line_x1[use_text_endpoint] <- text_x[use_text_endpoint] - side_gap[use_text_endpoint]
      line_y1[use_text_endpoint] <- text_y[use_text_endpoint]
      line_x2[use_text_endpoint] <- line_x1[use_text_endpoint]
      line_y2[use_text_endpoint] <- text_y[use_text_endpoint]
    }

    has_line <- is.finite(line_x0) & is.finite(line_y0) &
      is.finite(line_x1) & is.finite(line_y1) &
      is.finite(line_x2) & is.finite(line_y2)
    if (any(has_line)) {
      line_data <- data[has_line, , drop = FALSE]
      line_data$x0 <- line_x0[has_line]
      line_data$y0 <- line_y0[has_line]
      line_data$x1 <- line_x1[has_line]
      line_data$y1 <- line_y1[has_line]
      line_data$x2 <- line_x2[has_line]
      line_data$y2 <- line_y2[has_line]
      for (j in seq_len(nrow(line_data))) {
        line_colour <- if ("line_colour" %in% names(line_data) &&
          is.character(line_data$line_colour) &&
          !is.na(line_data$line_colour[[j]]) &&
          nzchar(line_data$line_colour[[j]])) {
          line_data$line_colour[[j]]
        } else {
          segment_colour
        }
        grobs[[length(grobs) + 1L]] <- grid::polylineGrob(
          x = grid::unit(
            c(line_data$x0[[j]], line_data$x1[[j]], line_data$x2[[j]]),
            "npc"
          ),
          y = grid::unit(
            c(line_data$y0[[j]], line_data$y1[[j]], line_data$y2[[j]]),
            "npc"
          ),
          gp = grid::gpar(
            col = line_colour,
            lwd = segment_linewidth * 2.845276,
            lty = segment_linetype,
            lineend = "round",
            linejoin = "round"
          )
        )
      }
    }

    for (i in seq_len(nrow(data))) {
      label_colour <- if ("label_colour" %in% names(data) &&
        is.character(data$label_colour) &&
        !is.na(data$label_colour[[i]]) &&
        nzchar(data$label_colour[[i]])) {
        data$label_colour[[i]]
      } else {
        text_colour
      }
      grobs[[length(grobs) + 1L]] <- grid::textGrob(
        label = data$label[[i]],
        x = grid::unit(text_x[[i]], "npc"),
        y = grid::unit(text_y[[i]], "npc"),
        hjust = data$hjust[[i]],
        vjust = data$vjust[[i]],
        rot = data$angle[[i]],
        gp = grid::gpar(
          col = label_colour,
          fontsize = text_size * 2.845276,
          fontfamily = family,
          fontface = fontface,
          lineheight = lineheight
        )
      )
    }
    if (!length(grobs)) {
      return(grid::nullGrob())
    }
    do.call(grid::grobTree, grobs)
  }
)

geom_plasmid_text_overlay <- function(data, family = "sans",
                                      text_colour = "black",
                                      text_size = 3,
                                      fontface = "plain",
                                      lineheight = 0.9,
                                      segment_colour = "grey45",
                                      segment_linewidth = 0.25,
                                      segment_linetype = "solid") {
  if (!"line_colour" %in% names(data)) {
    data$line_colour <- NA_character_
  }
  ggplot2::layer(
    data = data,
    mapping = ggplot2::aes(
      npcx = npcx,
      npcy = npcy,
      panel_x = panel_x,
      panel_y = panel_y,
      label = label,
      angle = angle,
      hjust = hjust,
      vjust = vjust,
      anchor_x = anchor_x,
      anchor_y = anchor_y,
      leader_x = leader_x,
      leader_y = leader_y,
      x0 = x0,
      y0 = y0,
      x1 = x1,
      y1 = y1,
      x2 = x2,
      y2 = y2,
      line_colour = line_colour,
      label_colour = colour
    ),
    stat = "identity",
    geom = GeomPlasmidTextOverlay,
    position = "identity",
    inherit.aes = FALSE,
    show.legend = FALSE,
    params = list(
      family = family,
      text_colour = text_colour,
      text_size = text_size,
      fontface = fontface,
      lineheight = lineheight,
      segment_colour = segment_colour,
      segment_linewidth = segment_linewidth,
      segment_linetype = segment_linetype,
      na.rm = TRUE
    )
  )
}

resolve_label_text_colours <- function(label_data, feature_colors,
                                       label_text_colour = "black") {
  n <- nrow(label_data)
  if (!n) {
    return(character())
  }
  if (is.null(label_text_colour)) {
    return(rep("black", n))
  }

  column_name <- NULL
  if (length(label_text_colour) == 1L && is.character(label_text_colour)) {
    requested_name <- label_text_colour[[1L]]
    exact_match <- requested_name %in% names(label_data)
    normalised_matches <- which(
      normalise_key(names(label_data)) == normalise_key(requested_name)
    )
    if (exact_match) {
      column_name <- requested_name
    } else if (length(normalised_matches) == 1L) {
      column_name <- names(label_data)[[normalised_matches[[1L]]]]
    }
  }

  if (!is.null(column_name)) {
    values <- as.character(label_data[[column_name]])
    if (normalise_key(column_name) == normalise_key("category")) {
      colours <- unname(feature_colors[values])
      colours[is.na(colours) | !nzchar(colours)] <- "black"
      source <- "the selected category palette"
    } else {
      colours <- values
      source <- paste0("annotation column `", column_name, "`")
    }
  } else if (length(label_text_colour) == n) {
    colours <- as.character(label_text_colour)
    source <- "colour vector"
  } else if (length(label_text_colour) == 1L) {
    colours <- rep(as.character(label_text_colour[[1L]]), n)
    source <- "fixed colour"
  } else {
    stop(
      "`label_text_colour` must be a fixed colour, a column name in the " ,
      "annotation table, or a vector with one colour per label.",
      call. = FALSE
    )
  }

  missing <- is.na(colours) | !nzchar(trimws(colours))
  colours[missing] <- "black"
  valid <- vapply(
    colours,
    function(value) {
      tryCatch({
        grDevices::col2rgb(value)
        TRUE
      }, error = function(error) FALSE)
    },
    logical(1L)
  )
  if (any(!valid)) {
    stop(
      "Values from ", source, " are not valid R colours: ",
      paste(unique(colours[!valid]), collapse = ", "),
      call. = FALSE
    )
  }
  colours
}

circular_feature_label_overlay <- function(label_data, genome_length,
                                           inner_radius = 0.38,
                                           label_line_length = 0.02,
                                           label_wrap_width = 28L,
                                           label_anchor_radius = 1.07,
                                           label_line_anchor_radius = label_anchor_radius,
                                           label_lane_spacing = 0.24,
                                           label_lane_step = NULL,
                                           label_text_size = 3.4,
                                           label_text_colour = "black",
                                           label_line_colour = "grey70",
                                           label_adjust = NULL,
                                           label_y_max = 1.95,
                                           gene_outer_radius = NULL,
                                           feature_colors = NULL) {
  theta <- (label_data$xmid %% genome_length) / genome_length * 360
  right <- theta < 180

  text <- wrap_label_text(
    label_data$label_display %||% label_data$label,
    width = label_wrap_width
  )
  label_colours <- resolve_label_text_colours(
    label_data,
    feature_colors %||% stats::setNames(character(), character()),
    label_text_colour = label_text_colour
  )
  line_colours <- resolve_label_text_colours(
    label_data,
    feature_colors %||% stats::setNames(character(), character()),
    label_text_colour = label_line_colour
  )

  label_position <- place_radial_label_boxes(
    x = label_data$xmid,
    genome_length = genome_length,
    text = text,
    right = right,
    inner_radius = inner_radius,
    label_anchor_radius = label_anchor_radius,
    label_line_length = label_line_length,
    label_lane_spacing = label_lane_spacing,
    label_lane_step = label_lane_step,
    label_text_size = label_text_size,
    priority = label_data$priority,
    y_max = label_y_max,
    gene_outer_radius = gene_outer_radius
  )
  keep <- is.finite(label_position$npcx) & is.finite(label_position$npcy)
  label_position <- label_position[keep, , drop = FALSE]
  label_data <- label_data[keep, , drop = FALSE]
  text <- text[keep]
  label_colours <- label_colours[keep]
  line_colours <- line_colours[keep]
  theta <- theta[keep]
  label_position <- apply_label_adjust(
    label_position,
    label_data,
    text,
    adjustment = label_adjust
  )
  if (!nrow(label_position)) {
    return(data.frame(
      npcx = numeric(),
      npcy = numeric(),
      panel_x = numeric(),
      panel_y = numeric(),
      label = character(),
      angle = numeric(),
      hjust = numeric(),
      vjust = numeric(),
      colour = character(),
      size = numeric(),
      fontface = character(),
      lineheight = numeric(),
      anchor_x = numeric(),
      anchor_y = numeric(),
      leader_x = numeric(),
      leader_y = numeric(),
      x0 = numeric(),
      y0 = numeric(),
      x1 = numeric(),
      y1 = numeric(),
      x2 = numeric(),
      y2 = numeric(),
      line_colour = character(),
      stringsAsFactors = FALSE
    ))
  }

  label_dims <- label_box_dimensions(text, label_text_size)
  line_x <- label_position$npcx
  line_y <- label_position$npcy
  theta_rad <- theta * pi / 180
  top_bottom <- abs(sin(theta_rad)) < 0.20
  edge_gap <- 0.008
  edge_shift <- label_dims$height * 0.50 + edge_gap
  line_x[top_bottom] <- label_position$npcx[top_bottom] -
    sin(theta_rad[top_bottom]) * edge_shift[top_bottom]
  line_y[top_bottom] <- label_position$npcy[top_bottom] -
    cos(theta_rad[top_bottom]) * edge_shift[top_bottom]

  data.frame(
    npcx = label_position$npcx,
    npcy = label_position$npcy,
    panel_x = NA_real_,
    panel_y = NA_real_,
    label = text,
    angle = 0,
    hjust = label_position$hjust,
    vjust = 0.5,
    colour = label_colours,
    size = label_text_size,
    fontface = "bold",
    lineheight = 0.9,
    anchor_x = label_data$xmid,
    anchor_y = label_line_anchor_radius,
    leader_x = NA_real_,
    leader_y = NA_real_,
    x0 = NA_real_,
    y0 = NA_real_,
    x1 = line_x,
    y1 = line_y,
    x2 = line_x,
    y2 = line_y,
    line_colour = line_colours,
    stringsAsFactors = FALSE
  )
}

circular_ruler_label_overlay <- function(ruler, genome_length,
                                         inner_radius = 0.38) {
  theta <- ruler$labels$theta
  angle <- -theta
  flip <- theta > 90 & theta < 270
  angle[flip] <- angle[flip] + 180
  data.frame(
    npcx = NA_real_,
    npcy = NA_real_,
    panel_x = ruler$labels$x,
    panel_y = ruler$labels$y,
    label = ruler$labels$label,
    angle = angle,
    hjust = 0.5,
    vjust = 0.5,
    colour = "grey25",
    size = 2.25,
    fontface = "plain",
    lineheight = 0.9,
    anchor_x = NA_real_,
    anchor_y = NA_real_,
    leader_x = NA_real_,
    leader_y = NA_real_,
    x0 = NA_real_,
    y0 = NA_real_,
    x1 = NA_real_,
    y1 = NA_real_,
    x2 = NA_real_,
    y2 = NA_real_,
    stringsAsFactors = FALSE
  )
}

center_text_overlay <- function(name, genome_length) {
  data.frame(
    npcx = 0.5,
    npcy = 0.5,
    panel_x = NA_real_,
    panel_y = NA_real_,
    label = paste(c(name, paste0(format_bp(genome_length), " bp")), collapse = "\n"),
    angle = 0,
    hjust = 0.5,
    vjust = 0.5,
    colour = "black",
    size = 4.3,
    fontface = "bold",
    lineheight = 0.9,
    anchor_x = NA_real_,
    anchor_y = NA_real_,
    leader_x = NA_real_,
    leader_y = NA_real_,
    x0 = NA_real_,
    y0 = NA_real_,
    x1 = NA_real_,
    y1 = NA_real_,
    x2 = NA_real_,
    y2 = NA_real_,
    stringsAsFactors = FALSE
  )
}

plot_circular_plasmid <- function(features, genome_length, name = NULL,
                                  phage_topology = c("circular", "linear"),
                                  category_scheme = "plasmid",
                                  gc_skew = NULL, show_labels = TRUE,
                                  label_unknown = FALSE,
                                  label_exclude_categories = NULL,
                                  label_pattern = NULL, max_labels = Inf,
                                  label_min_gap_deg = 0,
                                  label_line_length = 0.02,
                                  mini_label_line_length = NULL,
                                  label_wrap_width = 28L,
                                  label_anchor_radius = 1.07,
                                  label_lane_spacing = 0.24,
                                  label_lane_step = NULL,
                                  label_line_colour = "grey70",
                                  label_linewidth = 0.30,
                                  label_line_linetype = "dashed",
                                  label_adjust = NULL,
                                  min_feature_bp = 1,
                                  show_ruler = TRUE,
                                  ruler_major_bp = NULL,
                                  ruler_minor_bp = NULL,
                                  gene_palette = "npg",
                                  gene_manual_fill = NULL,
                                  blast_plot_data = NULL,
                                  blast_radius = NULL,
                                  blast_height = 0.08,
                                  gene_radius = 1,
                                  gene_height = 0.10,
                                  gene_color = "black",
                                  gene_linewidth = 0.22,
                                  gene_border_linewidth = NULL,
                                  gene_arrow_head_bp = NULL,
                                  gene_arrow_head_fraction = 0.55,
                                  gene_gap_bp = 0,
                                  gc_skew_radius = 0.78,
                                  gc_skew_height = 0.10,
                                  gc_content_radius = 0.67,
                                  gc_content_height = 0.045,
                                  gc_content_linewidth = 0.28,
                                  gc_legend_linewidth = 1.8,
                                  ruler_radius = 0.59,
                                  ruler_minor_tick = 0.015,
                                  ruler_major_tick = 0.035,
                                  ruler_label_radius = 0.47,
                                  ruler_linewidth = 0.28,
                                  ruler_minor_linewidth = 0.18,
                                  ruler_major_linewidth = 0.38,
                                  label_text_size = 3.4,
                                  label_text_colour = "black",
                                  label_fontface = "bold",
                                  ruler_text_size = 2.25,
                                  center_text_size = 4.3,
                                  center_text_colour = "black",
                                  legend_position = "bottom",
                                  legend_text_size = 8.5,
                                  legend_columns = NULL,
                                  gc_legend_columns = NULL,
                                  legend_plot_spacing = NULL,
                                  font_family = "sans",
                                  label_font_family = NULL,
                                  center_text_font_family = NULL,
                                  legend_font_family = NULL,
                                  ruler_fontface = "plain",
                                  inner_radius = 0.38) {
  phage_topology <- match.arg(phage_topology)
  label_fontface <- match.arg(
    label_fontface,
    c("plain", "bold", "italic", "bold.italic")
  )
  label_font_family <- label_font_family %||% font_family
  center_text_font_family <- center_text_font_family %||% font_family
  legend_font_family <- legend_font_family %||% font_family
  ruler_fontface <- match.arg(
    ruler_fontface,
    c("plain", "bold", "italic", "bold.italic")
  )
  if (!is.null(mini_label_line_length)) {
    label_line_length <- mini_label_line_length
  }
  gene_linewidth <- gene_border_linewidth %||% gene_linewidth
  legend_columns <- validate_legend_columns(legend_columns, legend_position)
  gc_legend_columns <- validate_gc_legend_columns(gc_legend_columns, legend_position)
  legend_plot_spacing_auto <- is.null(legend_plot_spacing)
  legend_plot_spacing <- validate_legend_plot_spacing(legend_plot_spacing, legend_position)
  corner_legend <- legend_position %in% c(
    "left_top", "right_top", "left_bottom", "right_bottom"
  )
  effective_legend_text_size <- if (corner_legend) {
    min(legend_text_size, 10)
  } else {
    legend_text_size
  }
  effective_legend_key_height <- if (corner_legend) 0.22 else 0.34
  effective_legend_key_width <- if (corner_legend) 0.36 else 0.45
  effective_gc_key_width <- if (corner_legend) 0.55 else 0.70
  effective_legend_spacing_y <- if (corner_legend) 0.16 else 0.35
  plot_colors <- ggplasmid_resolve_colors(
    category_scheme,
    gene_palette = gene_palette,
    gene_manual_fill = gene_manual_fill
  )
  feature_colors <- plot_colors[names(plot_colors) %in% unique(c(features$category, "GC skew+", "GC skew-"))]
  feature_breaks <- gene_legend_breaks(feature_colors)
  blast_data <- if (is.null(blast_plot_data)) blast_plot_empty() else blast_plot_data$data
  blast_colors <- if (is.null(blast_plot_data)) {
    stats::setNames(character(), character())
  } else {
    blast_plot_data$colors
  }
  blast_breaks <- if (is.null(blast_plot_data)) character() else blast_plot_data$breaks
  fill_colors <- feature_colors
  gc_colors <- gc_track_colors(plot_colors)
  gene_polys <- circular_gene_polygons(
    features,
    genome_length,
    radius = gene_radius,
    gene_height = gene_height,
    arrow_head_bp = gene_arrow_head_bp,
    arrow_head_fraction = gene_arrow_head_fraction,
    gene_gap_bp = gene_gap_bp
  )
  gene_polys$category <- factor(gene_polys$category, levels = names(plot_colors))
  ruler <- NULL

  p <- ggplot2::ggplot()

  if (!is.null(gc_skew)) {
    skew <- gc_skew
    skew$base_y <- gc_skew_radius
    skew$pos_y <- skew$base_y + pmax(skew$scaled_score, 0) * gc_skew_height
    skew$neg_y <- skew$base_y + pmin(skew$scaled_score, 0) * gc_skew_height
    skew$content_y <- gc_content_radius + skew$scaled_gc_content * gc_content_height
    p <- p +
      ggplot2::geom_ribbon(
        data = skew,
        ggplot2::aes(x = position, ymin = base_y, ymax = pos_y, fill = "GC skew+"),
        inherit.aes = FALSE,
        colour = NA,
        alpha = 0.9
      ) +
      ggplot2::geom_ribbon(
        data = skew,
        ggplot2::aes(x = position, ymin = base_y, ymax = neg_y, fill = "GC skew-"),
        inherit.aes = FALSE,
        colour = NA,
        alpha = 0.9
      ) +
      ggplot2::geom_segment(
        data = gc_skew_legend_data(),
        ggplot2::aes(x = x, xend = xend, y = y, yend = yend, colour = track),
        inherit.aes = FALSE,
        alpha = 0,
        linewidth = gc_content_linewidth,
        show.legend = TRUE
      )
    if (any(is.finite(skew$scaled_gc_content))) {
      p <- p +
        ggplot2::geom_path(
          data = skew,
          ggplot2::aes(x = position, y = content_y, colour = "GC content"),
          inherit.aes = FALSE,
          linewidth = gc_content_linewidth
        )
    }
  }

  p <- p +
    ggplot2::geom_polygon(
      data = gene_polys,
      ggplot2::aes(x = x, y = y, group = polygon_id, fill = category),
      colour = gene_color,
      linewidth = gene_linewidth
    )

  p <- p + ggplot2::scale_fill_manual(
    values = fill_colors,
    breaks = feature_breaks,
    name = NULL,
    drop = TRUE,
    na.value = "#999999",
    guide = ggplot2::guide_legend(
      order = 2,
      ncol = legend_columns,
      byrow = TRUE,
      keyheight = grid::unit(effective_legend_key_height, "cm"),
      keywidth = grid::unit(effective_legend_key_width, "cm"),
      override.aes = list(colour = NA)
    )
  )

  if (!is.null(blast_plot_data) && length(blast_plot_data$specs)) {
    blast_ring_keys <- names(blast_plot_data$specs)[
      order(blast_plot_data$order)
    ]
    p <- p + ggnewscale::new_scale_fill()
    for (ring_index in seq_along(blast_ring_keys)) {
      if (ring_index > 1L) {
        p <- p + ggnewscale::new_scale_fill()
      }
      ring_key <- blast_ring_keys[[ring_index]]
      ring_spec <- blast_plot_data$specs[[ring_key]]
      ring_data <- blast_data[
        !is.na(blast_data$ring_key) & blast_data$ring_key == ring_key,
        , drop = FALSE
      ]
      identity_min <- blast_plot_data$identity_thresholds[[ring_key]]
      gradient_min <- min(identity_min, 99.999)
      identity_breaks <- if (identity_min >= 100) {
        100
      } else {
        c(identity_min, 100)
      }
      identity_labels <- paste0(format(identity_breaks, trim = TRUE), "%")
      p <- p +
        ggplot2::geom_rect(
          data = ring_data,
          ggplot2::aes(
            xmin = start,
            xmax = end,
            ymin = radius - blast_height / 2,
            ymax = radius + blast_height / 2,
            fill = pident
          ),
          inherit.aes = FALSE,
          colour = NA,
          alpha = 0.95,
          show.legend = FALSE
        ) +
        ggplot2::scale_fill_gradient(
          low = "#E2E2E2",
          high = blast_plot_data$base_colors[[ring_key]],
          limits = c(gradient_min, 100),
          breaks = identity_breaks,
          labels = identity_labels,
          name = ring_spec$label,
          guide = ggplot2::guide_colorbar(
            order = 2L + ring_index,
            direction = "horizontal",
            reverse = TRUE,
            title.position = "top",
            label.position = "bottom",
            barwidth = grid::unit(1.8, "cm"),
            barheight = grid::unit(0.45, "cm"),
            frame.colour = "grey30",
            frame.linewidth = 0.4
          )
        )
    }
  }

  if (isTRUE(show_ruler)) {
    ruler <- circular_ruler_data(
      genome_length,
      major_bp = ruler_major_bp,
      minor_bp = ruler_minor_bp,
      ruler_radius = ruler_radius,
      ruler_minor_tick = ruler_minor_tick,
      ruler_major_tick = ruler_major_tick,
      ruler_label_radius = ruler_label_radius
    )

    p <- p +
      ggplot2::geom_path(
        data = ruler$ring,
        ggplot2::aes(x = x, y = y),
        inherit.aes = FALSE,
        colour = "grey35",
        linewidth = ruler_linewidth
      ) +
      ggplot2::geom_segment(
        data = ruler$ticks,
        ggplot2::aes(x = x, xend = x, y = y, yend = yend, linewidth = type),
        inherit.aes = FALSE,
        colour = "grey35"
      ) +
      ggplot2::scale_linewidth_manual(
        values = c(minor = ruler_minor_linewidth, major = ruler_major_linewidth),
        guide = "none"
      )
  }

  if (phage_topology == "linear") {
    p <- p +
      ggplot2::geom_segment(
        data = data.frame(x = 1, xend = 1, y = 0.70, yend = 1.10),
        ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
        inherit.aes = FALSE,
        linewidth = 0.75,
        colour = "black"
      )
  }

  label_data <- select_label_data(
    features,
    label_unknown = label_unknown,
    label_exclude_categories = label_exclude_categories,
    label_pattern = label_pattern,
    max_labels = max_labels,
    label_min_gap_deg = label_min_gap_deg,
    min_feature_bp = min_feature_bp,
    category_scheme = category_scheme,
    genome_length = genome_length
  )
  y_limit <- max(
    1.50,
    gene_radius + gene_height / 2 + 0.45,
    if (nrow(blast_data)) max(blast_data$radius) + blast_height / 2 + 0.15 else 0,
    gc_skew_radius + gc_skew_height / 2 + 0.24,
    gc_content_radius + gc_content_height / 2 + 0.22,
    ruler_radius + ruler_major_tick + 0.18
  )
  label_bounds <- NULL
  if (show_labels && nrow(label_data)) {
    label_count <- nrow(label_data)
    effective_label_text_size <- if (label_count > 40L) {
      label_size_floor <- 2.2
      max(label_size_floor, min(label_text_size, label_text_size * sqrt(40 / label_count)))
    } else {
      label_text_size
    }
    effective_label_wrap_width <- if (label_count > 18L) {
      if (is.null(label_wrap_width) || !is.finite(label_wrap_width)) {
        18L
      } else {
        min(as.integer(label_wrap_width), 18L)
      }
    } else {
      label_wrap_width
    }
    effective_label_lane_spacing <- label_lane_spacing
    effective_label_anchor_radius <- max(
      label_anchor_radius,
      gene_radius + gene_height / 2,
      na.rm = TRUE
    )
    label_overlay <- circular_feature_label_overlay(
      label_data,
      genome_length = genome_length,
      inner_radius = inner_radius,
      label_line_length = label_line_length,
      label_wrap_width = effective_label_wrap_width,
      label_anchor_radius = effective_label_anchor_radius,
      label_line_anchor_radius = gene_radius,
      label_lane_spacing = effective_label_lane_spacing,
      label_lane_step = label_lane_step,
      label_text_size = effective_label_text_size,
      label_text_colour = label_text_colour,
      label_line_colour = label_line_colour,
      label_adjust = label_adjust,
      label_y_max = y_limit,
      gene_outer_radius = gene_radius + gene_height / 2,
      feature_colors = plot_colors
    )
    label_bounds <- label_overlay_npc_bounds(label_overlay)
    if (isTRUE(legend_plot_spacing_auto)) {
      legend_plot_spacing <- auto_legend_plot_spacing(
        legend_plot_spacing,
        legend_position,
        label_bounds
      )
    }
    if (any(is.finite(label_overlay$panel_y))) {
      y_limit <- max(y_limit, max(label_overlay$panel_y, na.rm = TRUE) + 0.12)
    }
    p <- p + geom_plasmid_text_overlay(
      label_overlay,
      family = label_font_family,
      text_colour = "black",
      text_size = effective_label_text_size,
      fontface = label_fontface,
      lineheight = 0.9,
      segment_colour = label_line_colour,
      segment_linewidth = label_linewidth,
      segment_linetype = label_line_linetype
    )
  }

  if (!is.null(ruler)) {
    p <- p + geom_plasmid_text_overlay(
      circular_ruler_label_overlay(
        ruler,
        genome_length = genome_length,
        inner_radius = inner_radius
      ),
      family = font_family,
      text_colour = "grey25",
      text_size = ruler_text_size,
      fontface = ruler_fontface
    )
  }

  gc_colour_scale <- if (!is.null(gc_skew)) {
    ggplot2::scale_colour_manual(
      values = gc_colors,
      breaks = names(gc_colors),
      drop = TRUE
    )
  } else {
    NULL
  }

  p <- p +
    geom_plasmid_text_overlay(
      center_text_overlay(name, genome_length),
      family = center_text_font_family,
      text_colour = center_text_colour,
      text_size = center_text_size,
      fontface = "bold",
      lineheight = 0.9
    ) +
    ggplot2::scale_x_continuous(limits = c(0, genome_length), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, y_limit), expand = c(0, 0)) +
    gc_colour_scale +
    ggplot2::labs(colour = NULL) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        order = 1,
        ncol = gc_legend_columns,
        byrow = TRUE,
        keyheight = grid::unit(effective_legend_key_height, "cm"),
        keywidth = grid::unit(effective_gc_key_width, "cm"),
        override.aes = list(alpha = 1, linewidth = gc_legend_linewidth)
      )
    ) +
    radial_coord(theta = "x", start = 0, inner.radius = inner_radius, clip = "off", expand = FALSE) +
    ggplot2::theme_void(base_family = font_family) +
    ggplot2::theme(
      legend.position = legend_position_for_theme(legend_position),
      legend.position.inside = legend_position_inside_for_theme(
        legend_position,
        label_bounds
      ),
      legend.justification = legend_justification_for_theme(legend_position),
      legend.box = "vertical",
      legend.box.just = legend_box_justification_for_theme(legend_position),
      legend.box.spacing = grid::unit(legend_plot_spacing, "cm"),
      legend.spacing.y = grid::unit(effective_legend_spacing_y, "cm"),
      legend.text = ggplot2::element_text(
        family = legend_font_family,
        face = "bold", size = effective_legend_text_size, colour = "black"
      ),
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.key.height = grid::unit(effective_legend_key_height, "cm"),
      legend.key.width = grid::unit(effective_legend_key_width, "cm"),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = circular_plot_margin(label_bounds)
    )

  attr(p, "ggplasmid_layout") <- "circular"
  attr(p, "ggplasmid_style") <- list(
    category_scheme = category_scheme,
    gene_palette = gene_palette,
    fill_colors = plot_colors,
    fill_breaks = feature_breaks,
    blast_colors = blast_colors,
    blast_breaks = blast_breaks,
    blast_labels = if (is.null(blast_plot_data)) character() else blast_plot_data$labels,
    legend_columns = legend_columns,
    gc_legend_columns = gc_legend_columns,
    legend_plot_spacing = legend_plot_spacing
  )
  p
}

plot_linear_plasmid <- function(features, genome_length, name = NULL, rows = 4,
                                category_scheme = "plasmid", gc_skew = NULL,
                                show_labels = TRUE, label_unknown = FALSE,
                                label_exclude_categories = NULL,
                                label_pattern = NULL, max_labels = Inf,
                                label_min_gap_deg = 0,
                                min_feature_bp = 1,
                                gene_palette = "npg",
                                gene_manual_fill = NULL,
                                gene_height = 0.38,
                                gene_color = "black",
                                gene_linewidth = 0.24,
                                gene_border_linewidth = NULL,
                                gene_arrow_head_bp = NULL,
                                gene_arrow_head_fraction = 0.35,
                                linear_row_spacing = 4.00,
                                gc_skew_height = 0.16,
                                gc_content_height = 0.065,
                                gc_content_linewidth = 0.45,
                                gc_legend_linewidth = 1.8,
                                sequence_linewidth = 0.55,
                                label_text_size = 4.0,
                                label_text_colour = "black",
                                label_fontface = "bold",
                                label_wrap_width = 28L,
                                linear_label_wrap_width = 18L,
                                linear_label_max_lines = 2L,
                                label_line_colour = "grey70",
                                label_linewidth = 0.30,
                                label_line_linetype = "dashed",
                                label_adjust = NULL,
                                linear_label_offset = 0.14,
                                linear_label_lane_step = NULL,
                                linear_label_allow_gene_line_crossing = FALSE,
                                linear_label_line_angle = 90,
                                linear_label_text_angle = 0,
                                linear_label_side = "top",
                                row_label_text_size = 3.4,
                                legend_position = "right",
                                legend_text_size = 8.5,
                                legend_columns = NULL,
                                gc_legend_columns = NULL,
                                legend_plot_spacing = NULL,
                                font_family = "sans",
                                label_font_family = NULL,
                                center_text_font_family = NULL,
                                legend_font_family = NULL,
                                ruler_fontface = "plain") {
  label_font_family <- label_font_family %||% font_family
  legend_font_family <- legend_font_family %||% font_family
  label_fontface <- match.arg(
    label_fontface,
    c("plain", "bold", "italic", "bold.italic")
  )
  gene_linewidth <- gene_border_linewidth %||% gene_linewidth
  legend_columns <- validate_legend_columns(legend_columns, legend_position)
  gc_legend_columns <- validate_gc_legend_columns(gc_legend_columns, legend_position)
  legend_plot_spacing <- validate_legend_plot_spacing(legend_plot_spacing, legend_position)
  linear_row_spacing <- as.numeric(linear_row_spacing[[1L]])
  if (!is.finite(linear_row_spacing) || linear_row_spacing <= 0) {
    stop("`linear_row_spacing` must be a positive number.", call. = FALSE)
  }
  linear_label_text_angle <- as.numeric(linear_label_text_angle[[1L]])
  if (!is.finite(linear_label_text_angle)) {
    stop("`linear_label_text_angle` must be a finite number.", call. = FALSE)
  }
  linear_label_line_angle <- as.numeric(linear_label_line_angle[[1L]])
  if (!is.finite(linear_label_line_angle) ||
      linear_label_line_angle < 0 || linear_label_line_angle > 180) {
    stop("`linear_label_line_angle` must be between 0 and 180 degrees.", call. = FALSE)
  }
  linear_label_side <- match.arg(
    linear_label_side,
    c("top", "both")
  )
  plot_colors <- ggplasmid_resolve_colors(
    category_scheme,
    gene_palette = gene_palette,
    gene_manual_fill = gene_manual_fill
  )
  feature_breaks <- gene_legend_breaks(plot_colors)
  gc_colors <- gc_track_colors(plot_colors)
  linear <- linear_gene_polygons(
    features,
    genome_length,
    rows = rows,
    gene_height = gene_height,
    arrow_head_bp = gene_arrow_head_bp,
    arrow_head_fraction = gene_arrow_head_fraction,
    row_spacing = linear_row_spacing
  )
  rows <- linear$rows
  gene_polys <- linear$polygons
  gene_polys$category <- factor(gene_polys$category, levels = names(plot_colors))
  row_id <- sprintf(
    "%s: %d-%d kbp",
    name %||% "sequence",
    floor(linear$row_start / 1000),
    floor(linear$row_end / 1000)
  )
  row_label_x <- -linear$segment * 0.006
  plot_x_min <- -linear$segment * 0.07
  plot_x_max <- linear$segment * 1.025
  seq_lines <- data.frame(
    x = 0,
    xend = pmax(linear$row_end - linear$row_start, 1),
    y = linear$row_y,
    yend = linear$row_y,
    label = row_id,
    row_label_x = row_label_x
  )

  p <- ggplot2::ggplot()
  label_y_limit <- -Inf
  label_y_min <- Inf
  gc_track_y_min <- -Inf
  gc_row_upper_limits <- rep(Inf, rows)
  gc_row_lower_limits <- rep(Inf, rows)
  linear_label_summary <- list(
    n_requested = 0L,
    n_placed = 0L,
    n_removed = 0L,
    n_removed_cross_gene_line = 0L,
    removed_labels = character(),
    allow_gene_line_crossing = isTRUE(linear_label_allow_gene_line_crossing)
  )

  if (!is.null(gc_skew)) {
    skew <- gc_skew
    skew$row_index <- pmin(floor((skew$position - 1L) / linear$segment) + 1L, rows)
    skew$row_y <- linear$row_y[skew$row_index]
    skew$x <- skew$position - linear$row_start[skew$row_index]
    gc_track_gap <- 0.08
    skew$base_y <- skew$row_y - gene_height / 2 - gc_track_gap - gc_skew_height
    skew$pos_y <- skew$base_y + pmax(skew$scaled_score, 0) * gc_skew_height
    skew$neg_y <- skew$base_y + pmin(skew$scaled_score, 0) * gc_skew_height
    skew$content_y <- skew$base_y - 0.10 - gc_content_height +
      skew$scaled_gc_content * gc_content_height
    gc_track_y_min <- min(
      c(skew$neg_y, skew$content_y),
      na.rm = TRUE
    )
    gc_row_lower <- vapply(
      seq_len(rows),
      function(row) {
        idx <- which(skew$row_index == row)
        if (!length(idx)) {
          return(Inf)
        }
        min(
          c(skew$base_y[idx], skew$pos_y[idx], skew$neg_y[idx], skew$content_y[idx]),
          na.rm = TRUE
        )
      },
      numeric(1L)
    )
    gc_row_lower_limits <- gc_row_lower
    if (rows > 1L) {
      gc_row_upper_limits[seq.int(2L, rows)] <-
        gc_row_lower[seq.int(1L, rows - 1L)] - 0.02
    }
    p <- p +
      ggplot2::geom_ribbon(
        data = skew,
        ggplot2::aes(x = x, ymin = base_y, ymax = pos_y, group = row_index, fill = "GC skew+"),
        inherit.aes = FALSE,
        colour = NA,
        alpha = 0.9
      ) +
      ggplot2::geom_ribbon(
        data = skew,
        ggplot2::aes(x = x, ymin = base_y, ymax = neg_y, group = row_index, fill = "GC skew-"),
        inherit.aes = FALSE,
        colour = NA,
        alpha = 0.9
      ) +
      ggplot2::geom_segment(
        data = gc_skew_legend_data(),
        ggplot2::aes(x = x, xend = xend, y = y, yend = yend, colour = track),
        inherit.aes = FALSE,
        alpha = 0,
        linewidth = gc_content_linewidth,
        show.legend = TRUE
      )
    if (any(is.finite(skew$scaled_gc_content))) {
      p <- p +
        ggplot2::geom_path(
          data = skew,
          ggplot2::aes(x = x, y = content_y, group = row_index, colour = "GC content"),
          inherit.aes = FALSE,
          linewidth = gc_content_linewidth
        )
    }
  }

  p <- p +
    ggplot2::geom_segment(
      data = seq_lines,
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE,
      linewidth = sequence_linewidth,
      colour = "grey35"
    ) +
    ggplot2::geom_polygon(
      data = gene_polys,
      ggplot2::aes(x = x, y = y, group = polygon_id, fill = category),
      colour = gene_color,
      linewidth = gene_linewidth
    )

  label_data <- select_label_data(
    features,
    label_unknown = label_unknown,
    label_exclude_categories = label_exclude_categories,
    label_pattern = label_pattern,
    max_labels = max_labels,
    label_min_gap_deg = label_min_gap_deg,
    min_feature_bp = min_feature_bp,
    category_scheme = category_scheme,
    genome_length = genome_length
  )
  if (show_labels && nrow(label_data)) {
    label_data$row_index <- linear$row_index[match(label_data$feature_id, features$feature_id)]
    label_data$row_y <- linear$row_y[label_data$row_index]
    label_data$xmid_linear <- label_data$xmid - linear$row_start[label_data$row_index]
    label_data$xmid_linear <- pmin(pmax(label_data$xmid_linear, 0), linear$segment)
    label_data$label_colour <- resolve_label_text_colours(
      label_data,
      plot_colors,
      label_text_colour = label_text_colour
    )
    label_data$line_colour <- resolve_label_text_colours(
      label_data,
      plot_colors,
      label_text_colour = label_line_colour
    )
    layout_input <- label_data
    layout_for_angle <- function(angle) {
      linear_label_layout(
        layout_input,
        segment = linear$segment,
        gene_height = gene_height,
        label_text_size = label_text_size,
        label_wrap_width = linear_label_wrap_width,
        label_max_lines = linear_label_max_lines,
        label_offset = linear_label_offset,
        label_lane_step = linear_label_lane_step,
        allow_gene_line_crossing = linear_label_allow_gene_line_crossing,
        rows = rows,
        row_spacing = linear$row_spacing,
        line_angle = linear_label_line_angle,
        text_angle = angle,
        row_upper_limits = gc_row_upper_limits,
        row_lower_limits = gc_row_lower_limits,
        label_side = linear_label_side
      )
    }
    used_text_angle <- linear_label_text_angle
    label_data <- layout_for_angle(used_text_angle)
    first_summary <- attr(label_data, "linear_label_summary", exact = TRUE) %||%
      linear_label_summary
    if (abs(used_text_angle %% 180) > 1e-10 &&
        first_summary$n_removed > 0L) {
      fallback_angles <- c(0, 90)
      fallback_layouts <- lapply(fallback_angles, layout_for_angle)
      fallback_placed <- vapply(
        fallback_layouts,
        function(candidate) {
          candidate_summary <- attr(candidate, "linear_label_summary", exact = TRUE)
          candidate_summary$n_placed %||% nrow(candidate)
        },
        numeric(1L)
      )
      best_fallback <- which.max(fallback_placed)
      if (length(best_fallback) &&
          fallback_placed[[best_fallback]] > first_summary$n_placed) {
        used_text_angle <- fallback_angles[[best_fallback]]
        label_data <- fallback_layouts[[best_fallback]]
      }
    }
    linear_label_summary <- attr(label_data, "linear_label_summary", exact = TRUE) %||%
      linear_label_summary
    linear_label_summary$text_angle_requested <- linear_label_text_angle
    linear_label_summary$text_angle_used <- used_text_angle

    if (nrow(label_data)) {
      adjusted_position <- apply_label_adjust(
        data.frame(
          npcx = label_data$label_x,
          npcy = label_data$label_y,
          stringsAsFactors = FALSE
        ),
        label_data,
        label_data$label,
        adjustment = label_adjust
      )
      label_data$label_x <- adjusted_position$npcx
      label_data$label_y <- adjusted_position$npcy
      # A connector is needed only after a label leaves its first lane. Manual
      # movement also creates a connector so the adjustment remains visible.
      label_data$linear_label_has_leader <- label_data$label_lane > 0L |
        abs(label_data$label_x - label_data$linear_label_base_x) > 1e-10 |
        abs(label_data$label_y - label_data$linear_label_base_y) > 1e-10
      label_data <- linear_label_leader_endpoints(
        label_data,
        segment = linear$segment,
        y_span = max(linear$row_spacing * max(rows - 1L, 1L), 1),
        line_angle = linear_label_line_angle,
        text_angle = used_text_angle
      )
      label_data <- linear_label_text_anchor(
        label_data,
        segment = linear$segment,
        y_span = max(linear$row_spacing * max(rows - 1L, 1L), 1),
        text_angle = used_text_angle
      )

      label_y_limit <- max(
        label_y_limit,
        max(label_data$label_y + label_data$label_box_height / 2, na.rm = TRUE),
        na.rm = TRUE
      )
      label_y_min <- min(
        label_y_min,
        min(label_data$label_y - label_data$label_box_height / 2, na.rm = TRUE),
        na.rm = TRUE
      )
      plot_x_min <- min(
        plot_x_min,
        min(label_data$label_x - label_data$label_box_width / 2, na.rm = TRUE)
      )
      plot_x_max <- max(
        plot_x_max,
        max(label_data$label_x + label_data$label_box_width / 2, na.rm = TRUE)
      )
      leader_data <- label_data[label_data$linear_label_has_leader, , drop = FALSE]
      if (nrow(leader_data)) {
        for (leader_colour in unique(leader_data$line_colour)) {
          colour_data <- leader_data[
            leader_data$line_colour == leader_colour,
            ,
            drop = FALSE
          ]
          p <- p +
            ggplot2::geom_segment(
              data = colour_data,
              ggplot2::aes(
                x = xmid_linear,
                xend = line_x1,
                y = linear_label_gene_y,
                yend = line_y1
              ),
              inherit.aes = FALSE,
              colour = I(leader_colour),
              linewidth = label_linewidth,
              linetype = label_line_linetype
            )
        }
      }
      for (label_colour in unique(label_data$label_colour)) {
        colour_data <- label_data[
          label_data$label_colour == label_colour,
          ,
          drop = FALSE
        ]
        p <- p +
          ggplot2::geom_text(
            data = colour_data,
            ggplot2::aes(x = text_x, y = text_y, label = label),
            inherit.aes = FALSE,
            angle = used_text_angle,
            hjust = if (abs(used_text_angle %% 180) < 1e-10) 0.5 else 0,
            vjust = 0.5,
            colour = I(label_colour),
            size = label_text_size,
            fontface = label_fontface,
            family = label_font_family,
            lineheight = 0.90
          )
      }
    }
  }

  label_colours_for_scale <- if (nrow(label_data) &&
    all(c("label_colour", "line_colour") %in% names(label_data))) {
    unique(c(label_data$label_colour, label_data$line_colour))
  } else {
    character()
  }
  label_colours_for_scale <- label_colours_for_scale[
    nzchar(label_colours_for_scale) &
      !label_colours_for_scale %in% names(gc_colors)
  ]
  colour_scale_values <- c(
    gc_colors,
    stats::setNames(label_colours_for_scale, label_colours_for_scale)
  )
  gc_colour_scale <- if (!is.null(gc_skew)) {
    ggplot2::scale_colour_manual(
      values = colour_scale_values,
      breaks = names(gc_colors),
      drop = TRUE
    )
  } else {
    NULL
  }

  p <- p +
    ggplot2::geom_text(
      data = seq_lines,
      ggplot2::aes(x = row_label_x, y = y - 0.42, label = label),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 1,
      fontface = "bold",
      size = row_label_text_size
    ) +
    ggplot2::coord_cartesian(
      xlim = c(plot_x_min, plot_x_max),
      ylim = c(
          min(
            min(linear$row_y, na.rm = TRUE) - 0.55,
          if (is.finite(gc_track_y_min)) gc_track_y_min - 0.10 else Inf,
          if (is.finite(label_y_min)) label_y_min - 0.25 else Inf
          ),
        max(
          max(linear$row_y, na.rm = TRUE) + 0.75,
          if (is.finite(label_y_limit)) {
            label_y_limit + 0.25
          } else {
            -Inf
          }
        )
      ),
      clip = "off"
    ) +
    ggplot2::scale_fill_manual(
      values = plot_colors,
      breaks = feature_breaks,
      drop = TRUE,
      na.value = "#999999"
    ) +
    gc_colour_scale +
    ggplot2::labs(fill = NULL, colour = NULL) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        order = 2,
        ncol = legend_columns,
        keyheight = grid::unit(0.32, "cm"),
        keywidth = grid::unit(0.45, "cm")
      ),
      colour = ggplot2::guide_legend(
        order = 1,
        ncol = gc_legend_columns,
        byrow = TRUE,
        keyheight = grid::unit(0.34, "cm"),
        keywidth = grid::unit(0.70, "cm"),
        override.aes = list(alpha = 1, linewidth = gc_legend_linewidth)
      )
    ) +
    ggplot2::theme_void(base_family = font_family) +
    ggplot2::theme(
      legend.position = legend_position_for_theme(legend_position),
      legend.position.inside = legend_position_inside_for_theme(legend_position),
      legend.justification = legend_justification_for_theme(legend_position),
      legend.box = "vertical",
      legend.box.spacing = grid::unit(legend_plot_spacing, "cm"),
      legend.spacing.y = grid::unit(0.35, "cm"),
      legend.text = ggplot2::element_text(
        family = legend_font_family,
        face = "bold", size = legend_text_size, colour = "black"
      ),
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.key.height = grid::unit(0.30, "cm"),
      legend.key.width = grid::unit(0.45, "cm"),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(12, 10, 12, 48)
    )

  attr(p, "ggplasmid_layout") <- "linear"
  attr(p, "ggplasmid_linear_label_summary") <- linear_label_summary
  attr(p, "ggplasmid_style") <- list(
    category_scheme = category_scheme,
    gene_palette = gene_palette,
    fill_colors = plot_colors,
    fill_breaks = feature_breaks,
    legend_columns = legend_columns,
    gc_legend_columns = gc_legend_columns,
    legend_plot_spacing = legend_plot_spacing
  )
  p
}

#' Draw a plasmid map
#'
#' @param annotation A data frame or path to a CSV/TSV annotation table.
#' @param gbk Optional GenBank file. Its `ORIGIN` sequence is used directly for
#'   genome length and GC-track calculation when present.
#' @param fasta Optional FASTA file for annotation-table input, or as an
#'   alternative sequence source when the GenBank record has no `ORIGIN`.
#' @param skew_table Optional precomputed GC-skew table with `position` and `gc_skew`.
#' @param output Optional output path. When supplied, the plot is saved and
#'   invisibly returned.
#' @param layout Either `"circular"` or `"linear"`.
#' @param phage_topology Either `"circular"` or `"linear"`. Only used when
#'   `layout = "circular"`; `"linear"` adds a terminal boundary line.
#' @param name Display name.
#' @param genome_length Optional plasmid/genome length.
#' @param region_start,region_end Optional 1-based inclusive coordinates for
#'   the sequence window to display. Features are filtered and rebased to the
#'   window; sequence-derived GC tracks and FASTA records are sliced to match.
#' @param feature_types GenBank feature types to keep.
#' @param label_mode `"auto"` prefers compact gene labels, `"product"` uses
#'   product names, and `"gene"` prefers gene IDs.
#' @param category_scheme Either `"plasmid"` or `"phage"`.
#' @param gene_palette Any palette name exposed by an exported `ggsci::pal_*`
#'   function. Examples include `"npg"`, `"aaas"`, `"lancet"`, `"jco"`,
#'   `"locuszoom"`, and `"futurama"`.
#' @param gene_manual_fill Optional named character vector, or data frame with
#'   category/color columns, used to override selected gene-category fill
#'   colors after applying `gene_palette`.
#' @param gene_highlight Deprecated compatibility alias for
#'   `gene_manual_fill`.
#' @param blast_rings Optional list or data frame describing one or more BLAST
#'   comparison rings. List entries use names such as ring1 and each contain
#'   label and either a hits data frame or a sequence to BLAST against the
#'   reference sequence from the map's GenBank/FASTA input. The data-frame form
#'   has one row per ring and a list-column named hits or sequence. Per-ring
#'   settings can override the shared BLAST mode and thresholds. Ring order
#'   follows input order unless an order value is supplied. Each entry may
#'   also set colour.
#' @param blast_plotting_mode Shared BLAST display mode: hsp draws each
#'   qualifying high-scoring pair, region consolidates overlapping hits into
#'   identity-colored reference regions, and gene colors annotated features
#'   whose BLAST coverage passes the threshold.
#' @param blast_reference Which BLAST coordinate system is the map backbone:
#'   "query" (the default, suitable when P0413 is supplied as query) or
#'   "subject".
#' @param blast_min_identity Shared minimum percent identity. For automatic
#'   sequence comparisons it is passed directly to BLAST as `perc_identity`,
#'   and it is also used to filter the hits before plotting. For supplied hit
#'   tables it filters the existing hits only. If NULL, each ring uses a
#'   mode-specific default: 0 for hsp, 90 for region, and 80 for gene. A ring
#'   can override it with min_identity.
#' @param blast_min_alignment_length Shared minimum aligned length in bp. A
#'   ring can override it with min_alignment_length.
#' @param blast_min_gene_coverage Shared minimum fraction of an annotated
#'   feature covered in gene mode. A ring can override it with
#'   min_gene_coverage.
#' @param blast_palette One ggsci palette name used to assign distinct base
#'   colors to rings in order. Each ring can override its color with colour.
#'   Identity is shown as a continuous grey-to-ring-color gradient, with a
#'   separate horizontal legend bar titled by the ring label. The bar spans
#'   that ring's effective minimum identity through 100%.
#' @param blast_color_manual Optional character vector with one manual color
#'   per ring, in input order, or named by ring label. A ring-level colour
#'   setting overrides this vector.
#' @param blastn Optional path to the local NCBI BLAST+ blastn executable used
#'   for automatic ring searches.
#' @param blast_evalue,blast_args E-value and extra command-line arguments
#'   used for each automatic BLAST search. The identity threshold is passed
#'   from the shared or per-ring min_identity setting.
#' @param blast_radius,blast_height,blast_spacing Base radius, thickness, and
#'   spacing for the circular BLAST rings. By default, ring placement and
#'   thickness are calculated to keep all BLAST rings between the GC tracks and
#'   gene ring. Explicit dimensions that overlap those tracks produce an error.
#' @param rows Number of rows for linear layout.
#' @param genome_line_num,plot_line_num Aliases for `rows`. When supplied,
#'   `genome_line_num` overrides `rows`, and `plot_line_num` overrides both.
#' @param show_gc_skew Whether to draw GC-skew and GC-content tracks when
#'   sequence/skew data are available.
#' @param show_labels Whether to draw feature labels.
#' @param label_unknown Whether to label unknown/hypothetical features.
#' @param label_exclude_categories Optional category names whose features should
#'   still be drawn but not labeled.
#' @param label_pattern Optional regular expression. Matching features are
#'   prioritized for labels.
#' @param max_labels Maximum number of labels to consider. The default is `Inf`,
#'   so all eligible labels are considered. Circular labels are placed from
#'   shorter to longer text on farther outer lanes when needed. Use a finite
#'   value to label only the top eligible features.
#' @param label_min_gap_deg Minimum angular gap between circular labels. Larger
#'   values omit more crowded labels. Set to `0` to disable angular thinning.
#' @param label_line_length Distance from the feature anchor to radial labels,
#'   in circular plot radius units.
#' @param mini_label_line_length Clearer alias for `label_line_length`. If
#'   supplied, it overrides `label_line_length`.
#' @param label_wrap_width Number of characters per label line before wrapping.
#'   Set to `Inf` or `NULL` to disable wrapping.
#' @param label_lane_spacing Extra spacing for circular labels when nearby
#'   labels need a small outward offset.
#' @param label_lane_step Optional direct radial step between circular label
#'   lanes, in circular plot radius units. Smaller values pack lanes tighter.
#'   When `NULL`, the step is calculated from label height.
#' @param linear_label_offset Vertical gap between linear gene arrows and the
#'   first label lane.
#' @param linear_label_lane_step Optional tuning value for linear label spacing.
#'   This is the distance moved along each connector direction; larger values
#'   leave more space between layers. Linear labels are checked against other
#'   labels on the same genome row and remain below the gene line of the row
#'   above.
#' @param linear_label_allow_gene_line_crossing Whether a linear label may be
#'   pushed above the gene line of the row above when no free lane remains.
#'   The default `FALSE` removes that label and records it in
#'   `plasmid_data(plot)$linear_label_summary`.
#' @param linear_label_wrap_width Number of characters per line for linear
#'   labels before the maximum-line rule is applied. The default is `18` to
#'   keep linear labels compact.
#' @param linear_label_max_lines Maximum number of text rows for a linear gene
#'   label. The default is `2`; excess wrapped text is folded into the final
#'   row.
#' @param linear_label_line_angle Connector-line angle from the gene center to
#'   the label center, in degrees. `90` is vertical and `45` is diagonal. It is
#'   used when `linear_label_text_angle = 0`. For rotated text, labels are
#'   packed rightward first and then upward, so the connector follows the final
#'   label position instead of forcing this angle.
#' @param linear_label_text_angle Rotation angle for linear labels, in degrees.
#'   A nonzero value rotates the text and activates right-first, upward-second
#'   label packing. If that angle cannot place all labels, horizontal and
#'   vertical fallback angles are tested and the best-fitting angle is used.
#'   Use `0` to keep text horizontal and use `linear_label_line_angle` for the
#'   connector.
#' @param linear_label_side Which side of each linear gene may receive labels:
#'   `"top"` (the default) or `"both"`. In `"both"` mode labels are packed
#'   independently above and below each genome row.
#' @param label_line_angle,label_text_angle Short aliases for the two linear
#'   label-angle parameters. When supplied, they override the corresponding
#'   `linear_` values.
#' @param linear_row_spacing Vertical distance between adjacent linear genome
#'   rows.
#' @param label_line_colour,label_linewidth,label_line_linetype Label leader
#'   line colour, width, and linetype. `label_line_colour` may be a fixed R
#'   colour, a vector with one colour per label, or the name of any annotation
#'   column containing valid R colours.
#' @param label_adjust Manual label shifts for either layout. Use [label_adjust()], a
#'   data frame with `label`, `hjust`, and `vjust` columns, or a list with those
#'   entries. A data frame may use `gene` or `gene_id` instead of `label`.
#'   Positive `hjust` moves a label right and positive `vjust` moves it up;
#'   shifts are applied after automatic placement.
#' @param min_feature_bp Minimum feature length to label.
#' @param show_ruler Whether to draw the circular clock/ruler track. Only used
#'   when `layout = "circular"`.
#' @param ruler_major_bp,ruler_minor_bp Optional major/minor ruler tick spacing
#'   in bp. When omitted, sensible values are chosen from genome length.
#' @param gene_radius,gene_height Circular gene-ring radius and gene arrow
#'   height. When `gene_radius = NULL`, the circular gene ring is placed just
#'   outside the GC tracks or the outermost BLAST ring. For linear maps,
#'   `gene_height` controls arrow height.
#' @param gene_color Fixed outline color for gene arrows, independent of their
#'   category fill colors and the BLAST ring colors.
#' @param gene_linewidth,gene_border_linewidth Outline width for gene arrows.
#'   `gene_border_linewidth` is a clearer alias; if supplied, it overrides
#'   `gene_linewidth`.
#' @param gene_arrow_head_bp Optional gene arrow-head length in bp. The default
#'   scales with genome length.
#' @param gene_arrow_head_fraction Maximum fraction of each feature that can be
#'   used as the arrow head. When `NULL`, circular maps use `0.55` and linear
#'   maps use `0.35`.
#' @param gene_gap_bp Gap removed from adjacent circular feature arrows, in bp.
#'   Use `0` for no gap.
#' @param gc_skew_radius,gc_skew_height Circular GC-skew track radius and
#'   height. For linear maps, `gc_skew_height` controls the skew-track
#'   amplitude.
#' @param gc_content_radius,gc_content_height Circular GC-content track radius
#'   and height. For linear maps, `gc_content_height` controls line amplitude.
#' @param gc_content_linewidth Line width for the GC-content track. When `NULL`,
#'   the default is `0.28` for circular maps and `0.45` for linear maps.
#' @param gc_legend_linewidth Line width for GC content/skew entries in the
#'   legend.
#' @param ruler_minor_tick,ruler_major_tick,ruler_label_radius Circular ruler
#'   tick lengths and label radius.
#' @param ruler_linewidth,ruler_minor_linewidth,ruler_major_linewidth Circular
#'   ruler line widths.
#' @param sequence_linewidth Linear layout baseline width.
#' @param label_text_size,ruler_text_size,center_text_size,legend_text_size Text
#'   sizes.
#' @param label_fontface Font face for feature labels: `"plain"`, `"bold"`,
#'   `"italic"`, or `"bold.italic"`. The default is `"bold"`.
#' @param center_text_colour Colour for the circular map sample/genome label.
#' @param label_text_colour Label text colour. Use `"category"` to use the
#'   selected palette, a fixed R colour such as `"black"`, a vector with one
#'   colour per label, or the name of an annotation column containing valid R
#'   colours.
#' @param label_anchor_radius Starting radius for outside circular labels.
#'   Increase it to push labels farther from the gene ring.
#' @param row_label_text_size Text size for linear row labels.
#' @param legend_position Legend position passed to ggplot2, such as
#'   `"bottom"`, `"right"`, `"left"`, `"top"`, or `"none"`. Corner positions
#'   `"left_top"`, `"right_top"`, `"left_bottom"`, and `"right_bottom"`
#'   use outside side positions, aligned to the selected top or bottom corner.
#'   Side and corner legends align their guide contents to the left.
#'   Corner legends reserve outside plot space, with additional spacing for
#'   labels extending beyond the panel. Circular BLAST maps default to `"right"`.
#' @param legend_columns Number of columns in the feature-category legend. When
#'   omitted, side and corner legends use one column, while top and bottom
#'   legends use three columns.
#' @param gc_legend_columns Number of columns in the GC content/skew legend.
#'   When omitted, side and corner legends use one column, while top and bottom
#'   legends use three columns so the GC entries fit on one row.
#' @param legend_plot_spacing Space in cm between the plot panel and legend.
#'   When omitted, side and corner legends are moved outward based on the
#'   measured outer label range.
#' @param inner_radius Inner radius passed to `coord_radial()` for circular
#'   maps. Larger values make a larger center hole and push tracks outward.
#'   Use `NULL` to choose a layout-aware default.
#' @param font_family Font family passed to ggplot2 text themes. The default
#'   `"sans"` avoids Windows font database warnings.
#' @param label_font_family Optional font family for feature labels. Defaults
#'   to `font_family`.
#' @param center_text_font_family Optional font family for the circular center
#'   text. Defaults to `font_family`.
#' @param legend_font_family Optional font family for legend text. Defaults to
#'   `font_family`.
#' @param ruler_fontface Font face for circular ruler labels: `"plain"`,
#'   `"bold"`, `"italic"`, or `"bold.italic"`.
#' @param ruler_radius Radius of the circular ruler ring. The default `NULL`
#'   places it just inside the GC-content track.
#' @param gc_window GC-skew window.
#' @param gc_step GC-skew step. Use `NULL` to choose an automatic step that is
#'   faster for long genomes.
#' @param width,height,dpi Output size for `output`.
#' @return A ggplot object. Use [plasmid_data()] to retrieve the parsed
#'   annotation, FASTA, normalized features, GC-skew data, and the linear label
#'   placement summary used by the plot.
#' @examples
#' features <- data.frame(
#'   start = c(100, 1800, 3200),
#'   end = c(900, 2700, 4200),
#'   strand = c("+", "-", "+"),
#'   product = c("replication protein", "mobilization protein", "beta-lactamase"),
#'   category = c("Replication", "Mobile element", "Antimicrobial resistance")
#' )
#' ggplasmid(
#'   annotation = features,
#'   genome_length = 5000,
#'   name = "pExample",
#'   layout = "linear",
#'   plot_line_num = 2,
#'   show_gc_skew = FALSE,
#'   legend_position = "none"
#' )
#' @export
ggplasmid <- function(annotation = NULL, gbk = NULL, fasta = NULL, skew_table = NULL,
                      output = NULL, layout = c("circular", "linear"),
                      phage_topology = c("circular", "linear"),
                      name = NULL, genome_length = NULL,
                      region_start = NULL, region_end = NULL,
                      feature_types = c(
                        "CDS", "tRNA", "rRNA", "tmRNA", "ncRNA",
                        "rep_origin", "mobile_element"
                      ),
                      label_mode = c("auto", "product", "gene"),
                      category_scheme = c("plasmid", "phage"),
                      gene_palette = "npg", gene_manual_fill = NULL,
                      gene_highlight = NULL,
                      blast_rings = NULL,
                      blast_plotting_mode = c("hsp", "region", "gene"),
                      blast_reference = c("query", "subject"),
                      blast_min_identity = NULL,
                      blast_min_alignment_length = 1,
                      blast_min_gene_coverage = 0.80,
                      blast_palette = "npg",
                      blast_color_manual = NULL,
                      blastn = NULL,
                      blast_evalue = 1e-10,
                      blast_args = character(),
                      blast_radius = NULL,
                      blast_height = NULL,
                      blast_spacing = 0.01,
                      rows = 4, genome_line_num = NULL, plot_line_num = NULL,
                      show_gc_skew = TRUE, show_labels = TRUE,
                      label_unknown = FALSE, label_pattern = NULL,
                      label_exclude_categories = NULL,
                      max_labels = NULL, label_min_gap_deg = NULL,
                      label_line_length = 0.02,
                      mini_label_line_length = NULL,
                      label_wrap_width = 28L,
                      label_anchor_radius = 1.07,
                      label_lane_spacing = 0.24,
                      label_lane_step = NULL,
                      linear_label_offset = 0.14,
                      linear_label_lane_step = NULL,
                      linear_label_allow_gene_line_crossing = FALSE,
                      linear_label_line_angle = 90,
                      linear_label_text_angle = 0,
                      linear_label_side = "top",
                      label_line_angle = NULL,
                      label_text_angle = NULL,
                      label_line_colour = "grey70",
                      label_linewidth = 0.30,
                      label_line_linetype = "dashed",
                      label_adjust = NULL,
                      min_feature_bp = 1,
                      show_ruler = TRUE,
                      ruler_major_bp = NULL, ruler_minor_bp = NULL,
                      gene_radius = NULL,
                      gene_height = NULL,
                      gene_color = "black",
                      gene_linewidth = NULL,
                      gene_border_linewidth = NULL,
                      gene_arrow_head_bp = NULL,
                      gene_arrow_head_fraction = NULL,
                      linear_row_spacing = 4.00,
                      gene_gap_bp = NULL,
                      gc_skew_radius = 0.78,
                      gc_skew_height = NULL,
                      gc_content_radius = 0.67,
                      gc_content_height = NULL,
                      gc_content_linewidth = NULL,
                      gc_legend_linewidth = 1.8,
                      ruler_minor_tick = 0.015,
                      ruler_major_tick = 0.035,
                      ruler_label_radius = 0.47,
                      ruler_linewidth = 0.28,
                      ruler_minor_linewidth = 0.18,
                      ruler_major_linewidth = 0.38,
                      sequence_linewidth = 0.55,
                      label_text_size = NULL,
                      label_text_colour = "black",
                      label_fontface = "bold",
                      linear_label_wrap_width = 18L,
                      linear_label_max_lines = 2L,
                      row_label_text_size = 3.4,
                      ruler_text_size = 2.25,
                      center_text_size = 4.3,
                      center_text_colour = "black",
                      legend_position = NULL,
                      legend_text_size = 8.5,
                      legend_columns = NULL,
                      gc_legend_columns = NULL,
                      legend_plot_spacing = NULL,
                      inner_radius = NULL,
                      font_family = "sans",
                      label_font_family = NULL,
                      center_text_font_family = NULL,
                      legend_font_family = NULL,
                      ruler_fontface = "plain",
                      ruler_radius = NULL,
                      gc_window = 100L, gc_step = NULL,
                      width = NULL, height = NULL, dpi = 300) {
  layout <- match.arg(layout)
  phage_topology <- match.arg(phage_topology)
  label_mode <- match.arg(label_mode)
  category_scheme <- match.arg(category_scheme)
  blast_plotting_mode <- match.arg(blast_plotting_mode)
  blast_reference <- match.arg(blast_reference)
  if (!is.null(gene_highlight)) {
    if (!is.null(gene_manual_fill)) {
      stop("Use only one of `gene_manual_fill` and the legacy `gene_highlight` argument.",
           call. = FALSE)
    }
    warning(
      "`gene_highlight` is retained for compatibility; use `gene_manual_fill` instead.",
      call. = FALSE
    )
    gene_manual_fill <- gene_highlight
  }
  if (!is.null(blast_rings) && layout != "circular") {
    stop("BLAST plotting is currently supported for circular layout only.",
         call. = FALSE)
  }
  if (!is.null(mini_label_line_length)) {
    label_line_length <- mini_label_line_length
  }
  if (!is.null(label_line_angle)) {
    linear_label_line_angle <- label_line_angle
  }
  if (!is.null(label_text_angle)) {
    linear_label_text_angle <- label_text_angle
  }
  if (!is.null(genome_line_num)) {
    rows <- genome_line_num
  }
  if (!is.null(plot_line_num)) {
    rows <- plot_line_num
  }

  annotation_data <- read_plasmid_annotation(
    annotation = annotation,
    gbk = gbk,
    feature_types = feature_types
  )
  annotation_sequence <- attr(annotation_data, "sequence", exact = TRUE) %||% ""
  has_annotation_sequence <- is.character(annotation_sequence) &&
    length(annotation_sequence) == 1L && nzchar(annotation_sequence)
  fasta_data <- if (!is.null(fasta) && !has_annotation_sequence) {
    read_plasmid_fasta(fasta)
  } else {
    NULL
  }

  features <- prepare_plasmid_features(
    annotation = annotation_data,
    fasta = fasta_data,
    genome_length = genome_length,
    feature_types = feature_types,
    label_mode = label_mode,
    category_scheme = category_scheme
  )
  backbone_sequence <- attr(features, "sequence", exact = TRUE) %||% ""
  if ((!is.character(backbone_sequence) || !length(backbone_sequence) ||
       !nzchar(backbone_sequence[[1L]])) &&
      !is.null(fasta_data) && nrow(fasta_data)) {
    backbone_sequence <- as.character(fasta_data$sequence[[1L]])
  }
  full_genome_length <- attr(features, "genome_length", exact = TRUE)
  region <- resolve_plot_region(
    region_start = region_start,
    region_end = region_end,
    genome_length = full_genome_length
  )
  if (isTRUE(region$active)) {
    annotation_data <- subset_annotation_to_region(annotation_data, region)
    fasta_data <- subset_fasta_to_region(fasta_data, region)
    if (!is.null(skew_table)) {
      skew_table <- subset_gc_table_to_region(skew_table, region)
    }
    features <- subset_features_to_region(
      features,
      region = region,
      full_genome_length = full_genome_length
    )
  }
  genome_length <- attr(features, "genome_length", exact = TRUE)
  name <- name %||% attr(features, "name", exact = TRUE) %||% "plasmid"
  sequence <- attr(features, "sequence", exact = TRUE) %||% ""
  if (is.null(max_labels)) {
    max_labels <- Inf
  }
  if (is.null(label_min_gap_deg)) {
    label_min_gap_deg <- 0
  }
  if (is.null(inner_radius)) {
    inner_radius <- if (layout == "circular" && category_scheme == "phage") {
      0.14
    } else {
      0.38
    }
  }
  if (is.null(gene_height)) {
    gene_height <- if (layout == "circular" && category_scheme == "phage") {
      0.085
    } else if (layout == "circular") {
      0.10
    } else {
      0.38
    }
  }
  if (is.null(gene_linewidth)) {
    gene_linewidth <- if (layout == "circular" && category_scheme == "phage") {
      0.03
    } else if (layout == "circular") {
      0.22
    } else {
      0.24
    }
  }
  if (!is.null(gene_border_linewidth)) {
    gene_linewidth <- gene_border_linewidth
  }
  if (is.null(gene_arrow_head_bp) && layout == "circular" && category_scheme == "phage") {
    gene_arrow_head_bp <- max(60, genome_length / 700)
  }
  if (is.null(gene_arrow_head_fraction)) {
    gene_arrow_head_fraction <- if (layout == "linear") 0.35 else 0.55
  }
  if (is.null(gene_gap_bp)) {
    gene_gap_bp <- if (layout == "circular" && category_scheme == "phage") {
      max(25, genome_length / 5000)
    } else {
      0
    }
  }
  if (is.null(gc_skew_height)) {
    gc_skew_height <- if (layout == "circular") 0.10 else 0.16
  }
  if (is.null(gc_content_height)) {
    gc_content_height <- if (layout == "circular") 0.045 else 0.065
  }
  ring_count <- if (is.data.frame(blast_rings)) {
    nrow(blast_rings)
  } else if (is.list(blast_rings)) {
    length(blast_rings)
  } else {
    0L
  }
  if (layout == "circular" && (ring_count > 0L || is.null(gene_radius))) {
    ring_geometry <- blast_ring_geometry(
      n_rings = ring_count,
      radius = blast_radius,
      height = blast_height,
      spacing = blast_spacing,
      gene_radius = gene_radius,
      gene_height = gene_height,
      gc_skew_radius = gc_skew_radius,
      gc_skew_height = gc_skew_height,
      gc_content_radius = gc_content_radius,
      gc_content_height = gc_content_height,
      show_gc_skew = show_gc_skew,
      ruler_radius = ruler_radius,
      ruler_major_tick = ruler_major_tick
    )
    gene_radius <- ring_geometry$gene_radius
    if (ring_count > 0L) {
      blast_radius <- ring_geometry$radius
      blast_height <- ring_geometry$height
      blast_spacing <- ring_geometry$spacing
    }
  }
  if (is.null(gc_content_linewidth)) {
    gc_content_linewidth <- if (layout == "circular") 0.28 else 0.45
  }
  if (is.null(label_text_size)) {
    label_text_size <- if (layout == "circular") 3.4 else 4.0
  }
  if (is.null(legend_position)) {
    legend_position <- if (layout == "circular" && !is.null(blast_rings)) {
      "right"
    } else if (layout == "circular") {
      "bottom"
    } else {
      "right"
    }
  }
  legend_columns <- validate_legend_columns(legend_columns, legend_position)
  gc_legend_columns <- validate_gc_legend_columns(gc_legend_columns, legend_position)

  gc_skew <- NULL
  if (isTRUE(show_gc_skew)) {
    gc_skew <- prepare_gc_skew_for_plot(
      skew_table = skew_table,
      fasta = fasta_data,
      sequence = sequence,
      genome_length = genome_length,
      circular = layout == "circular",
      window = gc_window,
      step = gc_step
    )
  }

  blast_data_for_plot <- NULL
  if (!is.null(blast_rings)) {
    blast_data_for_plot <- blast_prepare_rings_plot_data(
      blast_rings = blast_rings,
      mode = blast_plotting_mode,
      reference = blast_reference,
      genome_length = genome_length,
      features = features,
      min_identity = blast_min_identity,
      min_alignment_length = blast_min_alignment_length,
      min_gene_coverage = blast_min_gene_coverage,
      palette = blast_palette,
      radius = blast_radius,
      height = blast_height,
      spacing = blast_spacing,
      coordinate_offset = if (isTRUE(region$active)) region$start - 1 else 0,
      backbone_sequence = backbone_sequence,
      blastn = blastn,
      evalue = blast_evalue,
      blast_args = blast_args,
      color_manual = blast_color_manual
    )
  }

  p <- if (layout == "circular") {
    plot_circular_plasmid(
      features = features,
      genome_length = genome_length,
      name = name,
      phage_topology = phage_topology,
      category_scheme = category_scheme,
      gc_skew = gc_skew,
      show_labels = show_labels,
      label_unknown = label_unknown,
      label_exclude_categories = label_exclude_categories,
      label_pattern = label_pattern,
      max_labels = max_labels,
      label_min_gap_deg = label_min_gap_deg,
      label_line_length = label_line_length,
      mini_label_line_length = mini_label_line_length,
      label_wrap_width = label_wrap_width,
      label_anchor_radius = label_anchor_radius,
      label_lane_spacing = label_lane_spacing,
      label_lane_step = label_lane_step,
      label_line_colour = label_line_colour,
      label_linewidth = label_linewidth,
      label_line_linetype = label_line_linetype,
      label_adjust = label_adjust,
      min_feature_bp = min_feature_bp,
      show_ruler = show_ruler,
      ruler_major_bp = ruler_major_bp,
      ruler_minor_bp = ruler_minor_bp,
      gene_palette = gene_palette,
      gene_manual_fill = gene_manual_fill,
      blast_plot_data = blast_data_for_plot,
      blast_radius = blast_radius,
      blast_height = blast_height,
      gene_radius = gene_radius,
      gene_height = gene_height,
      gene_color = gene_color,
      gene_linewidth = gene_linewidth,
      gene_border_linewidth = gene_border_linewidth,
      gene_arrow_head_bp = gene_arrow_head_bp,
      gene_arrow_head_fraction = gene_arrow_head_fraction,
      gene_gap_bp = gene_gap_bp,
      gc_skew_radius = gc_skew_radius,
      gc_skew_height = gc_skew_height,
      gc_content_radius = gc_content_radius,
      gc_content_height = gc_content_height,
      gc_content_linewidth = gc_content_linewidth,
      gc_legend_linewidth = gc_legend_linewidth,
      ruler_radius = ruler_radius %||% (gc_content_radius - 0.08),
      ruler_minor_tick = ruler_minor_tick,
      ruler_major_tick = ruler_major_tick,
      ruler_label_radius = ruler_label_radius,
      ruler_linewidth = ruler_linewidth,
      ruler_minor_linewidth = ruler_minor_linewidth,
      ruler_major_linewidth = ruler_major_linewidth,
      label_text_size = label_text_size,
      label_text_colour = label_text_colour,
      label_fontface = label_fontface,
      ruler_text_size = ruler_text_size,
      center_text_size = center_text_size,
      center_text_colour = center_text_colour,
      legend_position = legend_position,
      legend_text_size = legend_text_size,
      legend_columns = legend_columns,
      gc_legend_columns = gc_legend_columns,
      legend_plot_spacing = legend_plot_spacing,
      inner_radius = inner_radius,
      font_family = font_family,
      label_font_family = label_font_family,
      center_text_font_family = center_text_font_family,
      legend_font_family = legend_font_family,
      ruler_fontface = ruler_fontface
    )
  } else {
    plot_linear_plasmid(
      features = features,
      genome_length = genome_length,
      name = name,
      rows = rows,
      category_scheme = category_scheme,
      gc_skew = gc_skew,
      show_labels = show_labels,
      label_unknown = label_unknown,
      label_exclude_categories = label_exclude_categories,
      label_pattern = label_pattern,
      max_labels = max_labels,
      label_min_gap_deg = label_min_gap_deg,
      min_feature_bp = min_feature_bp,
      gene_palette = gene_palette,
      gene_manual_fill = gene_manual_fill,
      gene_color = gene_color,
      gene_height = gene_height,
      gene_linewidth = gene_linewidth,
      gene_border_linewidth = gene_border_linewidth,
      gene_arrow_head_bp = gene_arrow_head_bp,
      gene_arrow_head_fraction = gene_arrow_head_fraction,
      linear_row_spacing = linear_row_spacing,
      gc_skew_height = gc_skew_height,
      gc_content_height = gc_content_height,
      gc_content_linewidth = gc_content_linewidth,
      gc_legend_linewidth = gc_legend_linewidth,
      sequence_linewidth = sequence_linewidth,
      label_text_size = label_text_size,
      label_text_colour = label_text_colour,
      label_fontface = label_fontface,
      label_wrap_width = label_wrap_width,
      linear_label_wrap_width = linear_label_wrap_width,
      linear_label_max_lines = linear_label_max_lines,
      label_line_colour = label_line_colour,
      label_linewidth = label_linewidth,
      label_line_linetype = label_line_linetype,
      label_adjust = label_adjust,
      linear_label_offset = linear_label_offset,
      linear_label_lane_step = linear_label_lane_step,
      linear_label_allow_gene_line_crossing = linear_label_allow_gene_line_crossing,
      linear_label_line_angle = linear_label_line_angle,
      linear_label_text_angle = linear_label_text_angle,
      linear_label_side = linear_label_side,
      row_label_text_size = row_label_text_size,
      legend_position = legend_position,
      legend_text_size = legend_text_size,
      legend_columns = legend_columns,
      gc_legend_columns = gc_legend_columns,
      legend_plot_spacing = legend_plot_spacing,
      font_family = font_family,
      label_font_family = label_font_family,
      center_text_font_family = center_text_font_family,
      legend_font_family = legend_font_family,
      ruler_fontface = ruler_fontface
    )
  }

  if (is.null(fasta_data) && nzchar(sequence)) {
    fasta_data <- data.frame(
      name = name,
      description = name,
      sequence = sequence,
      length = nchar(sequence),
      stringsAsFactors = FALSE
    )
  }
  attr(p, "ggplasmid_data") <- list(
    annotation = annotation_data,
    fasta = fasta_data,
    features = features,
    gc_skew = gc_skew,
    genome_length = genome_length,
    name = name,
    region = region,
    blast = blast_data_for_plot,
    blast_plotting_mode = if (is.null(blast_data_for_plot)) NULL else blast_plotting_mode,
    blast_reference = if (is.null(blast_data_for_plot)) NULL else blast_reference,
    linear_label_summary = attr(p, "ggplasmid_linear_label_summary", exact = TRUE)
  )

  if (!is.null(output)) {
    save_plasmid_map(p, output, width = width, height = height, dpi = dpi)
    return(invisible(p))
  }
  p
}

resolve_plot_region <- function(region_start = NULL, region_end = NULL,
                                genome_length) {
  genome_length <- as.numeric(genome_length[[1L]])
  if (!is.finite(genome_length) || genome_length < 1) {
    stop("A finite genome length is required before selecting a sequence region.",
         call. = FALSE)
  }
  start <- if (is.null(region_start)) 1 else as.numeric(region_start[[1L]])
  end <- if (is.null(region_end)) genome_length else as.numeric(region_end[[1L]])
  if (!is.finite(start) || !is.finite(end) || start != floor(start) ||
      end != floor(end) || start < 1 || end < start || end > genome_length) {
    stop(
      "`region_start` and `region_end` must be integer coordinates within the " ,
      "genome, with `region_start <= region_end`.",
      call. = FALSE
    )
  }
  start <- as.integer(start)
  end <- as.integer(end)
  list(
    start = start,
    end = end,
    length = end - start + 1L,
    active = start != 1L || end != as.integer(genome_length)
  )
}

subset_annotation_to_region <- function(annotation, region) {
  if (is.null(annotation) || !nrow(annotation)) {
    return(annotation)
  }
  start_col <- find_column(
    annotation,
    c("start", "begin", "from", "left"),
    required = FALSE
  )
  end_col <- find_column(
    annotation,
    c("end", "stop", "to", "right"),
    required = FALSE
  )
  if (is.na(start_col) || is.na(end_col)) {
    return(annotation)
  }
  start <- as_number(annotation[[start_col]])
  end <- as_number(annotation[[end_col]])
  keep <- is.finite(start) & is.finite(end) &
    pmax(start, end) >= region$start & pmin(start, end) <= region$end
  annotation[keep, , drop = FALSE]
}

subset_fasta_to_region <- function(fasta_data, region) {
  if (is.null(fasta_data) || !nrow(fasta_data)) {
    return(fasta_data)
  }
  out <- as.data.frame(fasta_data, stringsAsFactors = FALSE)
  out$sequence <- vapply(
    out$sequence,
    function(sequence) substring(as.character(sequence), region$start, region$end),
    character(1L)
  )
  out$length <- nchar(out$sequence)
  if (!"name" %in% names(out)) {
    out$name <- paste0("sequence_", seq_len(nrow(out)))
  }
  if (!"description" %in% names(out)) {
    out$description <- out$name
  }
  out$fasta <- format_plasmid_fasta(out[c("name", "description", "sequence", "length")])
  out
}

subset_gc_table_to_region <- function(skew_table, region) {
  if (is.null(skew_table)) {
    return(skew_table)
  }
  skew <- read_gc_skew_table(skew_table)
  skew <- skew[
    is.finite(skew$position) &
      skew$position >= region$start & skew$position <= region$end,
    ,
    drop = FALSE
  ]
  skew$position <- skew$position - region$start + 1
  skew
}

subset_features_to_region <- function(features, region, full_genome_length) {
  if (!isTRUE(region$active)) {
    return(features)
  }
  pieces <- list()
  for (i in seq_len(nrow(features))) {
    start <- as.numeric(features$start[[i]])
    end <- as.numeric(features$end[[i]])
    if (!is.finite(start) || !is.finite(end)) {
      next
    }
    wraps_origin <- "wraps_origin" %in% names(features) &&
      isTRUE(features$wraps_origin[[i]]) && start > end
    intervals <- if (wraps_origin) {
      list(c(start, full_genome_length), c(1, end))
    } else {
      list(c(min(start, end), max(start, end)))
    }
    for (interval in intervals) {
      if (interval[[2L]] < region$start || interval[[1L]] > region$end) {
        next
      }
      row <- features[i, , drop = FALSE]
      row$start <- max(interval[[1L]], region$start) - region$start + 1
      row$end <- min(interval[[2L]], region$end) - region$start + 1
      row$wraps_origin <- FALSE
      pieces[[length(pieces) + 1L]] <- row
    }
  }
  if (!length(pieces)) {
    stop("No annotation features overlap the selected sequence region.", call. = FALSE)
  }
  out <- do.call(rbind, pieces)
  row.names(out) <- NULL
  out$feature_id <- seq_len(nrow(out))
  out$xmid <- feature_midpoint(out$start, out$end, region$length)
  sequence <- attr(features, "sequence", exact = TRUE) %||% ""
  attr(out, "genome_length") <- region$length
  attr(out, "name") <- attr(features, "name", exact = TRUE)
  attr(out, "sequence") <- if (nzchar(sequence)) {
    substring(sequence, region$start, region$end)
  } else {
    ""
  }
  out
}

#' Draw a plasmid map
#'
#' Alias for [ggplasmid()].
#' @inheritParams ggplasmid
#' @return A ggplot object.
#' @export
plot_plasmid_map <- ggplasmid

#' Retrieve data used by a ggplasmid plot
#'
#' @param plot A plot returned by [ggplasmid()] or [plot_phage_map()].
#' @return A list containing `annotation`, `fasta`, `features`, `gc_skew`,
#'   `genome_length`, `name`, and the selected `region`.
#' @export
plasmid_data <- function(plot) {
  data <- attr(plot, "ggplasmid_data", exact = TRUE)
  if (is.null(data)) {
    stop("`plot` does not contain ggplasmid data.", call. = FALSE)
  }
  data
}

#' Draw a phage map
#'
#' A light wrapper around [ggplasmid()] that uses the phage category scheme.
#' `phage_topology = "linear"` only adds a terminal boundary line when
#' `layout = "circular"`.
#'
#' @inheritParams ggplasmid
#' @return A ggplot object.
#' @export
plot_phage_map <- function(annotation = NULL, gbk = NULL, fasta = NULL,
                           skew_table = NULL, output = NULL,
                           layout = c("circular", "linear"),
                           phage_topology = c("circular", "linear"),
                           name = NULL, genome_length = NULL,
                           region_start = NULL, region_end = NULL,
                           feature_types = c(
                             "CDS", "tRNA", "rRNA", "tmRNA", "ncRNA",
                             "rep_origin", "mobile_element"
                           ),
                           label_mode = c("auto", "product", "gene"),
                           gene_palette = "npg", gene_manual_fill = NULL,
                           gene_highlight = NULL,
                           blast_rings = NULL,
                           blast_plotting_mode = c("hsp", "region", "gene"),
                           blast_reference = c("query", "subject"),
                           blast_min_identity = NULL,
                           blast_min_alignment_length = 1,
                           blast_min_gene_coverage = 0.80,
                           blast_palette = "npg",
                           blast_color_manual = NULL,
                           blastn = NULL,
                           blast_evalue = 1e-10,
                           blast_args = character(),
                           blast_radius = NULL,
                           blast_height = NULL,
                           blast_spacing = 0.01,
                           rows = 4, genome_line_num = NULL, plot_line_num = NULL,
                           show_gc_skew = TRUE,
                           show_labels = TRUE,
                           label_unknown = FALSE,
                           label_exclude_categories = NULL,
                           label_pattern = NULL, max_labels = NULL,
                           label_min_gap_deg = NULL,
                           label_line_length = 0.02,
                           mini_label_line_length = NULL,
                           min_feature_bp = 1,
                           label_wrap_width = 28L,
                           label_anchor_radius = 1.07,
                           label_lane_spacing = 0.24,
                           label_lane_step = NULL,
                           linear_label_offset = 0.14,
                           linear_label_lane_step = NULL,
                           linear_label_allow_gene_line_crossing = FALSE,
                           linear_label_line_angle = 90,
                           linear_label_text_angle = 0,
                           linear_label_side = "top",
                           label_line_angle = NULL,
                           label_text_angle = NULL,
                           label_line_colour = "grey70",
                           label_linewidth = 0.30,
                           label_line_linetype = "dashed",
                           label_adjust = NULL,
                           show_ruler = TRUE,
                           ruler_major_bp = NULL, ruler_minor_bp = NULL,
                           gene_radius = NULL,
                           gene_height = NULL,
                           gene_color = "black",
                           gene_linewidth = NULL,
                           gene_border_linewidth = NULL,
                           gene_arrow_head_bp = NULL,
                           gene_arrow_head_fraction = NULL,
                           linear_row_spacing = 4.00,
                           gene_gap_bp = NULL,
                           gc_skew_radius = 0.78,
                           gc_skew_height = NULL,
                           gc_content_radius = 0.67,
                           gc_content_height = NULL,
                           gc_content_linewidth = NULL,
                           gc_legend_linewidth = 1.8,
                           ruler_minor_tick = 0.015,
                           ruler_major_tick = 0.035,
                           ruler_label_radius = 0.47,
                           ruler_linewidth = 0.28,
                           ruler_minor_linewidth = 0.18,
                           ruler_major_linewidth = 0.38,
                           sequence_linewidth = 0.55,
                           label_text_size = NULL,
                           label_text_colour = "black",
                           label_fontface = "bold",
                           linear_label_wrap_width = 18L,
                           linear_label_max_lines = 2L,
                           row_label_text_size = 3.4,
                           ruler_text_size = 2.25,
                           center_text_size = 4.3,
                           center_text_colour = "black",
                           legend_position = NULL,
                           legend_text_size = 8.5,
                           legend_columns = NULL,
                           gc_legend_columns = NULL,
                           legend_plot_spacing = NULL,
                           inner_radius = NULL,
                           font_family = "sans",
                           label_font_family = NULL,
                           center_text_font_family = NULL,
                           legend_font_family = NULL,
                           ruler_fontface = "plain",
                           ruler_radius = NULL,
                           gc_window = 100L, gc_step = NULL,
                           width = NULL, height = NULL, dpi = 300) {
  ggplasmid(
    annotation = annotation,
    gbk = gbk,
    fasta = fasta,
    skew_table = skew_table,
    output = output,
    layout = match.arg(layout),
    phage_topology = match.arg(phage_topology),
    name = name,
    genome_length = genome_length,
    region_start = region_start,
    region_end = region_end,
    feature_types = feature_types,
    label_mode = match.arg(label_mode),
    category_scheme = "phage",
    gene_palette = gene_palette,
    gene_manual_fill = gene_manual_fill,
    gene_highlight = gene_highlight,
    blast_rings = blast_rings,
    blast_plotting_mode = match.arg(blast_plotting_mode),
    blast_reference = match.arg(blast_reference),
    blast_min_identity = blast_min_identity,
    blast_min_alignment_length = blast_min_alignment_length,
    blast_min_gene_coverage = blast_min_gene_coverage,
    blast_palette = blast_palette,
    blast_color_manual = blast_color_manual,
    blastn = blastn,
    blast_evalue = blast_evalue,
    blast_args = blast_args,
    blast_radius = blast_radius,
    blast_height = blast_height,
    blast_spacing = blast_spacing,
    rows = rows,
    genome_line_num = genome_line_num,
    plot_line_num = plot_line_num,
    show_gc_skew = show_gc_skew,
    show_labels = show_labels,
    label_unknown = label_unknown,
    label_exclude_categories = label_exclude_categories,
    label_pattern = label_pattern,
    max_labels = max_labels,
    label_min_gap_deg = label_min_gap_deg,
    label_line_length = label_line_length,
    mini_label_line_length = mini_label_line_length,
    label_wrap_width = label_wrap_width,
    linear_label_wrap_width = linear_label_wrap_width,
    linear_label_max_lines = linear_label_max_lines,
    label_anchor_radius = label_anchor_radius,
    label_lane_spacing = label_lane_spacing,
    label_lane_step = label_lane_step,
    linear_label_offset = linear_label_offset,
    linear_label_lane_step = linear_label_lane_step,
        linear_label_allow_gene_line_crossing = linear_label_allow_gene_line_crossing,
        linear_label_line_angle = linear_label_line_angle,
        linear_label_text_angle = linear_label_text_angle,
        linear_label_side = linear_label_side,
        label_line_angle = label_line_angle,
    label_text_angle = label_text_angle,
    label_line_colour = label_line_colour,
    label_linewidth = label_linewidth,
    label_line_linetype = label_line_linetype,
    label_adjust = label_adjust,
    min_feature_bp = min_feature_bp,
    show_ruler = show_ruler,
    ruler_major_bp = ruler_major_bp,
    ruler_minor_bp = ruler_minor_bp,
    gene_radius = gene_radius,
    gene_height = gene_height,
    gene_color = gene_color,
    gene_linewidth = gene_linewidth,
    gene_border_linewidth = gene_border_linewidth,
    gene_arrow_head_bp = gene_arrow_head_bp,
    gene_arrow_head_fraction = gene_arrow_head_fraction,
    linear_row_spacing = linear_row_spacing,
    gene_gap_bp = gene_gap_bp,
    gc_skew_radius = gc_skew_radius,
    gc_skew_height = gc_skew_height,
    gc_content_radius = gc_content_radius,
    gc_content_height = gc_content_height,
    gc_content_linewidth = gc_content_linewidth,
    gc_legend_linewidth = gc_legend_linewidth,
    ruler_minor_tick = ruler_minor_tick,
    ruler_major_tick = ruler_major_tick,
    ruler_label_radius = ruler_label_radius,
    ruler_linewidth = ruler_linewidth,
    ruler_minor_linewidth = ruler_minor_linewidth,
    ruler_major_linewidth = ruler_major_linewidth,
    sequence_linewidth = sequence_linewidth,
    label_text_size = label_text_size,
    label_text_colour = label_text_colour,
    label_fontface = label_fontface,
    row_label_text_size = row_label_text_size,
    ruler_text_size = ruler_text_size,
    center_text_size = center_text_size,
    center_text_colour = center_text_colour,
    legend_position = legend_position,
    legend_text_size = legend_text_size,
    legend_columns = legend_columns,
    gc_legend_columns = gc_legend_columns,
    legend_plot_spacing = legend_plot_spacing,
    inner_radius = inner_radius,
    font_family = font_family,
    label_font_family = label_font_family,
    center_text_font_family = center_text_font_family,
    legend_font_family = legend_font_family,
    ruler_fontface = ruler_fontface,
    ruler_radius = ruler_radius,
    gc_window = gc_window,
    gc_step = gc_step,
    width = width,
    height = height,
    dpi = dpi
  )
}

#' Save a plasmid map
#'
#' @param plot A ggplot object returned by [ggplasmid()].
#' @param filename Output filename.
#' @param width,height Plot size in inches. Defaults depend on layout.
#' @param dpi Output resolution.
#' @param bg Background color.
#' @return The normalized output path, invisibly.
#' @export
save_plasmid_map <- function(plot, filename, width = NULL, height = NULL,
                             dpi = 300, bg = "white") {
  layout <- attr(plot, "ggplasmid_layout", exact = TRUE) %||% "circular"
  if (is.null(width)) {
    width <- if (layout == "linear") 22 else 9
  }
  if (is.null(height)) {
    height <- if (layout == "linear") 11 else 9
  }
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = bg,
    limitsize = FALSE
  )
  invisible(normalizePath(filename, mustWork = FALSE))
}
