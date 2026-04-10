# map_utils.R
#
# Functions for creating geographic map visualisations of the ddsynth datasets.
#
# Requires: rnaturalearth, rnaturalearthdata, sf (all in Suggests)
#           ggplot2, dplyr, tibble, cli (all in Imports)


# --- Pathogen groups and colour palette --------------------------------------

# Six epidemiological groups matching the WHO priority pathogen classification.
# Colours are shaded within each group: darker = larger / better-known pathogen.

.pathogen_group_map <- c(
  # Viral haemorrhagic fevers — red/coral family
  "EVD"      = "Viral haemorrhagic fevers",
  "MVD"      = "Viral haemorrhagic fevers",
  "Lassa"    = "Viral haemorrhagic fevers",
  "CCHF"     = "Viral haemorrhagic fevers",
  # Pandemic respiratory — blue family
  "COVID-19" = "Pandemic respiratory",
  "SARS"     = "Pandemic respiratory",
  "MERS"     = "Pandemic respiratory",
  "Flu"      = "Pandemic respiratory",
  # Arboviral / vector-borne — orange/amber family
  "Dengue"   = "Arboviral / vector-borne",
  "Zika"     = "Arboviral / vector-borne",
  "RVF"      = "Arboviral / vector-borne",
  "YFV"      = "Arboviral / vector-borne",
  # Bat-reservoir zoonoses — purple
  "Nipah"    = "Bat-reservoir zoonoses",
  # Human-to-human viral — green family
  "Mpox"     = "Human-to-human viral",
  "Measles"  = "Human-to-human viral",
  # Environmental / zoonotic bacterial — teal
  "Cholera"  = "Environmental / zoonotic bacterial"
)

# Named colour vector — one hex per pathogen, shaded within group
.pathogen_colours <- c(
  # Viral haemorrhagic fevers: dark → light red/coral
  "EVD"      = "#922B21",
  "MVD"      = "#C0392B",
  "CCHF"     = "#E74C3C",
  "Lassa"    = "#F1948A",
  # Pandemic respiratory: dark → light blue
  "COVID-19" = "#1A3C5E",
  "SARS"     = "#2471A3",
  "MERS"     = "#5DADE2",
  "Flu"      = "#AED6F1",
  # Arboviral / vector-borne: dark → light orange/amber
  "Dengue"   = "#7E5109",
  "RVF"      = "#CA6F1E",
  "Zika"     = "#E59866",
  "YFV"      = "#FAD7A0",
  # Bat-reservoir zoonoses: purple
  "Nipah"    = "#6C3483",
  # Human-to-human viral: dark → light green
  "Mpox"     = "#1D6A39",
  "Measles"  = "#52BE80",
  # Environmental / zoonotic bacterial: teal
  "Cholera"  = "#0E6655"
)

# Ordered list of pathogens for the legend (no fake header entries).
.legend_breaks_full <- c(
  "EVD",      "MVD",    "CCHF",   "Lassa",
  "COVID-19", "SARS",   "MERS",   "Flu",
  "Dengue",   "RVF",    "Zika",   "YFV",
  "Nipah",
  "Mpox",     "Measles",
  "Cholera"
)

# Build legend breaks and labels filtered to pathogens present in the map data.
# The first active pathogen in each group gets a two-line label:
#   "Group name\n  Pathogen"
# so group names appear inline without needing fake header entries (which
# ggplot2 silently drops from the guide key data, causing override.aes errors).
.build_legend <- function(active_pathogens) {
  breaks       <- .legend_breaks_full[.legend_breaks_full %in% active_pathogens]
  groups_seen  <- character(0)
  labels <- vapply(breaks, function(b) {
    grp      <- unname(.pathogen_group_map[b])
    is_first <- !(grp %in% groups_seen)
    if (is_first) groups_seen <<- c(groups_seen, grp)
    if (is_first) paste0(grp, "\n  ", b) else paste0("  ", b)
  }, character(1L))
  list(breaks = breaks, labels = labels)
}


# --- Internal helpers --------------------------------------------------------

.extract_n <- function(entry) {
  if (!is.null(entry$n)) {
    val <- suppressWarnings(as.numeric(entry$n))
    if (!is.na(val)) return(val)
  }
  if (!is.null(entry$freq_count)) return(sum(entry$freq_count))
  NA_real_
}

# Normalise country name strings to match rnaturalearth name_long field
.normalise_country <- function(x) {
  dplyr::case_when(
    x == "UK"            ~ "United Kingdom",
    x == "South Korea"   ~ "Republic of Korea",
    x == "Brunei"        ~ "Brunei Darussalam",
    x == "USA"           ~ "United States of America",
    x == "DRC"           ~ "Democratic Republic of the Congo",
    x == "Faroe Islands" ~ "Faeroe Islands",   # rnaturalearth spelling
    TRUE                 ~ x
  )
}

# Strings that represent multi-country or non-geocodable entries;
# these are shown in the inset box rather than plotted on the map.
.multi_country_strings <- c("Mixed", "Middle East")

#' @noRd
get_country_centroids <- function(country_names) {
  if (!requireNamespace("rnaturalearth", quietly = TRUE) ||
      !requireNamespace("rnaturalearthdata", quietly = TRUE) ||
      !requireNamespace("sf", quietly = TRUE)) {
    cli::cli_abort(c(
      "Geographic packages are required for mapping.",
      "i" = "Install with: {.code install.packages(c('rnaturalearth', 'rnaturalearthdata', 'sf'))}"
    ))
  }

  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  suppressWarnings(
    centroids <- sf::st_centroid(world)
  )
  coords <- sf::st_coordinates(centroids)

  world_df <- data.frame(
    name_long = world$name_long,
    name      = world$name,
    lon       = coords[, 1],
    lat       = coords[, 2],
    stringsAsFactors = FALSE
  )

  norm_names <- .normalise_country(country_names)

  idx_long <- match(norm_names, world_df$name_long)
  idx_alt  <- match(norm_names, world_df$name)
  idx      <- ifelse(!is.na(idx_long), idx_long, idx_alt)

  unmatched <- unique(norm_names[is.na(idx) & !is.na(country_names) & country_names != "Mixed"])
  if (length(unmatched) > 0) {
    cli::cli_warn("Could not geocode countr{?y/ies}: {.val {unmatched}}")
  }

  data.frame(
    country = country_names,
    lon     = world_df$lon[idx],
    lat     = world_df$lat[idx],
    stringsAsFactors = FALSE
  )
}


# --- Exported functions ------------------------------------------------------

#' Extract a tidy summary of all built-in datasets
#'
#' Loops through all per-pathogen dataset lists and returns one row per dataset
#' entry. Sample size `n` is taken from the `n` field if present, otherwise
#' computed as `sum(freq_count)`.
#'
#' @return A [tibble::tibble()] with columns `pathogen`, `pathogen_group`,
#'   `dataset_id`, `country`, `n`, `n_log` (log1p of n), `subgroup`, and `source`.
#' @export
extract_dataset_summary <- function() {
  pathogen_datasets <- list(
    "Nipah"    = datasets_Nipah,
    "MVD"      = datasets_MVD,
    "EVD"      = datasets_EVD,
    "Lassa"    = datasets_Lassa,
    "SARS"     = datasets_SARS,
    "MERS"     = datasets_MERS,
    "Zika"     = datasets_Zika,
    "Measles"  = datasets_Measles,
    "Mpox"     = datasets_Mpox,
    "Cholera"  = datasets_Cholera,
    "RVF"      = datasets_RVF,
    "CCHF"     = datasets_CCHF,
    "COVID-19" = datasets_COVID_19,
    "Dengue"   = datasets_Dengue,
    "YFV"      = datasets_YFV
  )

  # Flu may not yet be defined
  if (exists("datasets_Flu", envir = asNamespace("ddsynth"), inherits = FALSE)) {
    pathogen_datasets[["Flu"]] <- get("datasets_Flu", envir = asNamespace("ddsynth"))
  }

  rows <- vector("list", 300L)
  row_idx <- 0L

  for (pathogen in names(pathogen_datasets)) {
    dset <- pathogen_datasets[[pathogen]]
    if (length(dset) == 0L) next

    entry_names <- names(dset)

    for (i in seq_along(dset)) {
      entry <- dset[[i]]

      # Handle unnamed entries (when <- is used inside list() instead of =)
      nm <- if (!is.null(entry_names) && nchar(entry_names[i]) > 0) entry_names[i]
            else paste0("d", i)

      n_val    <- .extract_n(entry)
      n_log    <- if (!is.na(n_val) && n_val > 0) log1p(n_val) else NA_real_
      raw_country <- if (!is.null(entry$country)) entry$country else NA_character_
      country  <- if (!is.na(raw_country) && nchar(trimws(raw_country)) > 0) trimws(raw_country) else NA_character_
      subgroup <- if (!is.null(entry$subgroup)) entry$subgroup else NA_character_
      source   <- if (!is.null(entry$source))   entry$source   else NA_character_

      row_idx <- row_idx + 1L
      rows[[row_idx]] <- data.frame(
        pathogen       = pathogen,
        pathogen_group = unname(.pathogen_group_map[pathogen]),
        dataset_id     = nm,
        country        = country,
        n              = n_val,
        n_log          = n_log,
        subgroup       = subgroup,
        source         = source,
        stringsAsFactors = FALSE
      )
    }
  }

  tibble::as_tibble(do.call(rbind, rows[seq_len(row_idx)]))
}


#' Create a global map of dataset geographic distribution
#'
#' Plots each dataset entry as a dot on a world map. Dot size is proportional
#' to `log(n+1)` and dot colour indicates pathogen. Datasets with
#' `country = "Mixed"` are excluded from the map and listed in an inset box.
#'
#' For reproducible jitter, call [set.seed()] before this function.
#'
#' @param data A data frame from [extract_dataset_summary()]. If `NULL`,
#'   [extract_dataset_summary()] is called internally.
#' @param exclude_no_n Logical. If `TRUE`, entries with missing sample size are
#'   dropped. If `FALSE` (default), they appear at size `default_n_log`.
#' @param default_n_log Numeric. Log-scale dot size used when `n` is missing.
#' @param alpha Numeric (0–1). Point transparency.
#' @param jitter_amount Numeric. Maximum spatial jitter in degrees lat/lon,
#'   applied to separate overlapping points from the same country.
#' @param title Character. Plot title.
#'
#' @return A [ggplot2::ggplot()] object.
#' @export
create_data_map <- function(
    data           = NULL,
    exclude_no_n   = FALSE,
    default_n_log  = 1.5,
    alpha          = 0.6,
    jitter_amount  = 3.0,
    title          = "Geographic distribution of delay distribution datasets"
) {
  if (!requireNamespace("rnaturalearth",     quietly = TRUE) ||
      !requireNamespace("rnaturalearthdata", quietly = TRUE) ||
      !requireNamespace("sf",                quietly = TRUE)) {
    cli::cli_abort(c(
      "Geographic packages are required for mapping.",
      "i" = "Install with: {.code install.packages(c('rnaturalearth', 'rnaturalearthdata', 'sf'))}"
    ))
  }

  if (is.null(data)) data <- extract_dataset_summary()

  # --- Handle missing n -------------------------------------------------------
  n_missing <- sum(is.na(data$n_log))
  if (exclude_no_n) {
    data <- data[!is.na(data$n_log), ]
  } else if (n_missing > 0) {
    data$n_log[is.na(data$n_log)] <- default_n_log
    cli::cli_inform(
      "{n_missing} entr{?y/ies} with missing n shown at default size (log(n+1) = {default_n_log})."
    )
  }

  # --- Separate Mixed/multi-country from geocodable entries ------------------
  is_mixed <- !is.na(data$country) & data$country %in% .multi_country_strings
  is_na    <- is.na(data$country)

  mixed_df <- data[is_mixed, ]
  map_df   <- data[!is_mixed & !is_na, ]

  n_na <- sum(is_na)
  if (n_na > 0) {
    cli::cli_inform("{n_na} entr{?y/ies} with missing country excluded from map.")
  }

  # --- Geocode ----------------------------------------------------------------
  unique_countries <- unique(map_df$country)
  centroids        <- get_country_centroids(unique_countries)
  map_df           <- merge(map_df, centroids, by = "country", all.x = TRUE)
  map_df           <- map_df[!is.na(map_df$lon), ]

  # --- Spatial jitter (call set.seed() before create_data_map() to fix) -------
  map_df$lon_j <- map_df$lon + stats::runif(nrow(map_df), -jitter_amount, jitter_amount)
  map_df$lat_j <- map_df$lat + stats::runif(nrow(map_df), -jitter_amount, jitter_amount)
  map_df$lat_j <- pmax(pmin(map_df$lat_j, 70), -43)  # clip to visible range

  # --- World basemap ----------------------------------------------------------
  world_sf <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

  # --- Legend structure (built from pathogens actually in the map data) -------
  legend <- .build_legend(unique(map_df$pathogen))

  # --- Mixed inset label text -------------------------------------------------
  mixed_label <- NULL
  if (nrow(mixed_df) > 0) {
    mixed_lines <- vapply(seq_len(nrow(mixed_df)), function(i) {
      row   <- mixed_df[i, ]
      n_str <- if (!is.na(row$n)) paste0("n=", row$n) else "n=?"
      sub_str <- if (!is.na(row$subgroup)) paste0(" [", row$subgroup, "]") else ""
      paste0(row$pathogen, sub_str, ": ", n_str)
    }, character(1L))
    mixed_label <- paste(c("Mixed/multi-country:", mixed_lines), collapse = "\n")
  }

  # --- Build plot -------------------------------------------------------------
  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(
      data      = world_sf,
      fill      = "#eeeeee",
      colour    = "white",
      linewidth = 0.2
    ) +
    ggplot2::geom_point(
      data    = map_df,
      mapping = ggplot2::aes(
        x      = .data$lon_j,
        y      = .data$lat_j,
        colour = .data$pathogen,
        size   = .data$n_log
      ),
      alpha = alpha
    ) +
    ggplot2::scale_colour_manual(
      name     = NULL,
      values   = .pathogen_colours,
      breaks   = legend$breaks,
      labels   = legend$labels,
      na.value = "grey60"
    ) +
    ggplot2::scale_size_continuous(
      range = c(0.5, 4),
      name  = "log(n+1)"
    ) +
    ggplot2::coord_sf(
      crs  = sf::st_crs(4326),
      xlim = c(-125, 155),
      ylim = c(-45, 72)
    ) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        ncol         = 1,
        override.aes = list(size = rep(3, length(legend$breaks)))
      )
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text        = ggplot2::element_blank(),
      axis.ticks       = ggplot2::element_blank(),
      panel.grid       = ggplot2::element_line(colour = "white"),
      legend.position  = "right",
      legend.key.size  = ggplot2::unit(0.4, "lines"),
      legend.text      = ggplot2::element_text(size = 8)
    )

  # Place Mixed inset in the Indian Ocean (open ocean within cropped extent)
  if (!is.null(mixed_label)) {
    p <- p + ggplot2::annotate(
      geom       = "label",
      x          = -110,
      y          = -30,
      label      = mixed_label,
      hjust      = 0,
      vjust      = 1,
      size       = 2.3,
      fill       = "white",
      colour     = "grey30",
      label.size = 0.3,
      lineheight = 1.1
    )
  }

  p
}
