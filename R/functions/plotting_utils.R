#' @title ERA Plotting Utilities
#' @description
#' Functions for visualising seasonal raster outputs and building stacks.

NULL

# --- palettes ----

get_circ_pal <- function(name = "phenology") {

  if (name == "phenology") {
    cols <- c(
      "#2c7bb6", "#00a6ca", "#00ccbc", "#90eb9d",
      "#ffff8c", "#f9d057", "#f29e2e", "#e76818",
      "#d7191c", "#7f3b08", "#5e4fa2"
    )
  } else {
    stop("Unknown circular palette")
  }

  cols
}

get_seq_pal <- function(name = "magma") {
  viridisLite::viridis(100, option = name)
}

# --- layer detection ----

is_circular_layer <- function(nm) {
  grepl("sos|eos", nm, ignore.case = TRUE)
}

# --- plotting ----

plot_season_stack <- function(r,
                              season_name = NULL,
                              circ_palette = "phenology",
                              seq_palette = "magma") {

  nl <- terra::nlyr(r)
  ncol <- 3
  nrow <- ceiling(nl / ncol)

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  par(mfrow = c(nrow, ncol), mar = c(2,2,2,4))

  circ_cols <- get_circ_pal(circ_palette)
  seq_cols  <- get_seq_pal(seq_palette)

  for (i in seq_len(nl)) {

    lyr <- r[[i]]
    nm  <- names(r)[i]

    if (is_circular_layer(nm)) {

      terra::plot(
        lyr,
        col = circ_cols,
        main = nm,
        axes = FALSE,
        legend = TRUE
      )

    } else {

      terra::plot(
        lyr,
        col = seq_cols,
        main = nm,
        axes = FALSE,
        legend = TRUE
      )
    }
  }

  if (!is.null(season_name)) {
    mtext(season_name, outer = TRUE, cex = 1.2, line = -1)
  }
}

# --- stack builder ----

build_season_stack <- function(raster_list,
                               vars = c("sos", "eos", "slen", "rain"),
                               include_n = TRUE) {

  season_ids   <- as.integer(sub("^season_", "", names(raster_list)))
  season_names <- names(raster_list)[order(season_ids)]

  layer_list <- list()

  for (snm in season_names) {
    r <- raster_list[[snm]]
    season_num <- sub("^season_", "", snm)

    for (v in vars) {
      idx <- grep(paste0("_", v, "_"), names(r), ignore.case = TRUE)

      if (length(idx) == 0L) {
        warning(sprintf("No layers found for var '%s' in %s", v, snm))
        next
      }

      layer_list <- c(layer_list, as.list(r[[idx]]))
    }

    if (include_n) {
      n_name <- paste0("season", season_num, "_n")
      if (n_name %in% names(r)) {
        layer_list <- c(layer_list, list(r[[n_name]]))
      } else {
        warning(sprintf("Layer '%s' not found in %s", n_name, snm))
      }
    }
  }

  if (length(layer_list) == 0L) {
    stop("No layers found when building season stack.")
  }

  x <- rast(layer_list)

  nm <- names(x)
  is_n <- grepl("_n$", nm)
  x <- x[[c(which(!is_n), which(is_n))]]

  x
}

# --- masking ----

mask_raster_list_by_nprop <- function(raster_list, n_prop_min = 0.1) {

  out <- lapply(names(raster_list), function(snm) {

    r <- raster_list[[snm]]
    season_num <- sub("^season_", "", snm)

    nprop_name <- paste0("season", season_num, "_n_prop")

    if (!(nprop_name %in% names(r))) return(r)

    mask_layer <- r[[nprop_name]] >= n_prop_min

    target_idx <- grep("_sos_|_eos_|_rain_|_slen_", names(r), ignore.case = TRUE)

    r<-terra::mask(r,mask_layer, maskvalues = 0)

    r
  })

  names(out) <- names(raster_list)
  out
}

plot_season_vars <- function(
    raster_list,
    vars = c("sos"),
    include_n = TRUE,
    season_name = NULL,
    circ_palette = "phenology",
    seq_palette = "ylgnbu"
) {

  x <- build_season_stack(
    raster_list = raster_list,
    vars = vars,
    include_n = include_n
  )

  plot_season_stack(
    x,
    season_name = season_name,
    circ_palette = circ_palette,
    seq_palette = seq_palette
  )

  invisible(x)
}


save_plot_season_vars <- function(raster_list,
                                  output_dir_plot,
                                  filename,
                                  vars = c("sos"),
                                  include_n = TRUE,
                                  ncol = 3,
                                  circ_palette = "phenology",
                                  seq_palette = "ylgnbu",
                                  panel_width = 4,
                                  panel_height = 4,
                                  res = 300) {

  if (!dir.exists(output_dir_plot)) {
    dir.create(output_dir_plot, recursive = TRUE)
  }

  # Build stack first so we know how many panels there will be
  x <- build_season_stack(
    raster_list = raster_list,
    vars = vars,
    include_n = include_n
  )

  n_panels <- terra::nlyr(x)
  nrow <- ceiling(n_panels / ncol)

  outfile <- file.path(output_dir_plot, filename)

  png(
    filename = outfile,
    width = ncol * panel_width,
    height = nrow * panel_height,
    units = "in",
    res = res,
    bg = "white"
  )

  plot_season_stack(
    x,
    season_name = NULL,
    circ_palette = circ_palette,
    seq_palette = seq_palette
  )

  dev.off()

  invisible(outfile)
}

save_plot_season_stack <- function(r,
                                   output_dir_plot,
                                   filename,
                                   season_name = NULL,
                                   circ_palette = "phenology",
                                   seq_palette = "magma",
                                   ncol = 3,
                                   panel_width = 4,
                                   panel_height = 4,
                                   res = 300) {

  if (!dir.exists(output_dir_plot)) {
    dir.create(output_dir_plot, recursive = TRUE)
  }

  n_panels <- terra::nlyr(r)
  nrow <- ceiling(n_panels / ncol)

  outfile <- file.path(output_dir_plot, filename)

  png(
    filename = outfile,
    width = ncol * panel_width,
    height = nrow * panel_height,
    units = "in",
    res = res,
    bg = "white"
  )

  plot_season_stack(
    r,
    season_name = season_name,
    circ_palette = circ_palette,
    seq_palette = seq_palette
  )

  dev.off()

  invisible(outfile)
}
