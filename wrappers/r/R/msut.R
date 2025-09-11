# R/msut.R

#' Parse mzML into compact bin (raw bytes)
#' @param data raw vector with mzML file content
#' @param slim logical; pass TRUE to slim parsing
#' @return
#' @export
parse_mzml <- function(data, slim = FALSE) {
  stopifnot(is.raw(data))
  .Call("C_parse_mzml", data, as.logical(slim), PACKAGE = "msut")
}

#' Convert bin blob back to JSON text
#' @param bin
#' @param pretty logical pretty-print
#' @return character(1) JSON
#' @export
bin_to_json <- function(bin, pretty = FALSE) {
  stopifnot(is.raw(bin))
  .Call("C_bin_to_json", bin, as.logical(pretty), PACKAGE = "msut")
}

bin_to_df <- function(bin, pretty = FALSE) {
  stopifnot(is.raw(bin))
  txt <- bin_to_json(bin, pretty = pretty)

  x <- jsonlite::fromJSON(
    txt,
    simplifyVector   = TRUE,
    simplifyDataFrame= TRUE,
    simplifyMatrix   = TRUE
  )

  if (is.null(x$Ok)) return(x)
  ok  <- x$Ok
  run <- ok$run

  ensure_list_col <- function(df, col) {
    if (!nrow(df)) { df[[col]] <- I(vector("list", 0)); return(df) }
    v <- df[[col]]
    if (is.null(v)) v <- vector("list", nrow(df))
    if (!is.list(v)) v <- lapply(seq_len(nrow(df)), function(i) v[[i]])
    df[[col]] <- I(v)
    df
  }

  chrom <- run$chromatograms
  if (is.null(chrom)) {
    chrom <- data.frame(
      index = integer(), array_length = integer(),
      time_array = I(list()), intensity_array = I(list()),
      id = character(), stringsAsFactors = FALSE
    )
  } else if (!is.data.frame(chrom)) {
    chrom <- as.data.frame(chrom, stringsAsFactors = FALSE)
  }
  chrom <- ensure_list_col(chrom, "time_array")
  chrom <- ensure_list_col(chrom, "intensity_array")

  spec <- run$spectra
  if (is.null(spec)) {
    spec <- data.frame(
      index = integer(), array_length = integer(),
      ms_level = integer(), polarity = integer(), spectrum_type = integer(),
      retention_time = numeric(),
      scan_window_lower_limit = numeric(),
      scan_window_upper_limit = numeric(),
      total_ion_current = numeric(),
      base_peak_intensity = numeric(),
      base_peak_mz = numeric(),
      mz_array = I(list()), intensity_array = I(list()),
      stringsAsFactors = FALSE
    )
  } else if (!is.data.frame(spec)) {
    spec <- as.data.frame(spec, stringsAsFactors = FALSE)
  }
  spec <- ensure_list_col(spec, "mz_array")
  spec <- ensure_list_col(spec, "intensity_array")

  run$chromatograms <- chrom
  run$spectra       <- spec
  ok$run            <- run
  x$Ok              <- ok
  x
}


#' Parse mzML and return {json metadata, arrays blob}
#' @param data raw mzML bytes
#' @param slim logical
#' @return list(json=character(1), blob=raw)
#' @export
parse_mzml_to_json <- function(data, slim = FALSE) {
  stopifnot(is.raw(data))
  .Call("C_parse_mzml_to_json", data, as.logical(slim), PACKAGE = "msut")
}

#' Get peak from an (x, y) chromatogram
#'
#' @param x numeric vector of retention times
#' @param y numeric vector of intensities (will be converted to f32 internally)
#' @param rt numeric target retention time (center)
#' @param range numeric RT range scale used to window/score the peak
#' @param options list of thresholds (optional): \cr
#' @return named list with numeric/integer fields:
#'   \code{from, to, rt, integral, intensity, ratio, np}
#' @export
get_peak <- function(x, y, rt, range, options = NULL) {
  stopifnot(is.numeric(x), is.numeric(y))
  if (length(x) != length(y) || length(x) < 3) {
    stop("x and y must have the same length (>= 3)")
  }
  out_json <- .Call(
    "C_get_peak",
    as.numeric(x), as.numeric(y),
    as.numeric(rt), as.numeric(range),
    options,
    PACKAGE = "msut"
  )
  jsonlite::fromJSON(out_json, simplifyVector = TRUE)
}

#' Get peaks for many (rt, mz, range) items
#'
#' @param bin raw BIN1
#' @param rts,mzs,ranges numeric vectors (same length)
#' @param from_left,to_right numeric window around each `rt`
#' @param options list, same structure as in `get_peak()`
#' @return data.frame with columns: from,to,rt,integral,intensity,ratio,np (rows aligned to inputs)
#' @export
get_peaks_from_eic <- function(
  bin,
  df,
  from_left = 0.5,
  to_right  = 0.5,
  options = NULL,
  cores = 1L
) {
  stopifnot(is.raw(bin))
  if (!is.data.frame(df)) stop("`df` must be a data.frame")
  req <- c("id", "rt", "mz", "ranges")
  miss <- setdiff(req, names(df))
  if (length(miss)) stop("missing columns: ", paste(miss, collapse = ", "))

  id     <- as.character(df$id)
  rts    <- suppressWarnings(as.numeric(df$rt))
  mzs    <- suppressWarnings(as.numeric(df$mz))
  ranges <- suppressWarnings(as.numeric(df$ranges))

  n <- length(id)
  if (!(length(rts) == n && length(mzs) == n && length(ranges) == n)) {
    stop("id, rt, mz, ranges must have the same length")
  }
  id[is.na(id)] <- ""

  out_json <- .Call(
    "C_get_peaks_from_eic",
    bin,
    as.numeric(rts),
    as.numeric(mzs),
    as.numeric(ranges),
    as.character(id),
    as.numeric(from_left),
    as.numeric(to_right),
    options,
    as.integer(cores),
    PACKAGE = "msut"
  )

  res <- jsonlite::fromJSON(out_json, simplifyVector = TRUE)
  if (!is.data.frame(res)) res <- as.data.frame(res)
  want <- c("id", "mz", "ort", "rt", "from", "to", "intensity", "integral")
  present <- intersect(want, names(res))
  extras  <- setdiff(names(res), present)
  res <- res[, c(present, extras), drop = FALSE]
  rownames(res) <- NULL
  res
}



#' Get peaks from chromatograms (BIN1 + items)
#' @param bin raw BIN1 (output of parse_mzml)
#' @param items list/data.frame with numeric columns rt, range (optional mz)
#' @param options list of thresholds (optional)
#' @return data.frame with columns index, id, rt, from, to, intensity, integral
#' @export
get_peaks_from_chrom <- function(bin, items, options = NULL, cores = 1L) {
  stopifnot(is.raw(bin))
  if (is.null(items) || !(is.list(items) || is.data.frame(items))) stop("items must be a list/data.frame")
  idxs    <- suppressWarnings(as.integer(items$idx %||% items$index))
  rts     <- suppressWarnings(as.numeric(items$rt))
  windows <- suppressWarnings(as.numeric(items$window %||% items$range))
  if (length(idxs) != length(rts) || length(windows) != length(rts)) stop("idx, rt, range length mismatch")

  out_json <- .Call("C_get_peaks_from_chrom", bin, idxs, rts, windows, options, as.integer(cores), PACKAGE = "msut")

  df <- jsonlite::fromJSON(out_json, simplifyVector = TRUE)
  if (!is.data.frame(df)) df <- as.data.frame(df)
  if (!"ort" %in% names(df)) stop("internal error: 'ort' missing from result")

  desired <- c("index", "id", "ort", "rt", "from", "to", "intensity", "integral")
  present <- desired[desired %in% names(df)]
  extras  <- setdiff(names(df), present)
  df <- df[, c(present, extras), drop = FALSE]
  if ("index" %in% names(df)) df <- df[order(df$index), , drop = FALSE]
  rownames(df) <- NULL
  df
}


#' Extracted Ion Chromatogram (EIC)
#' @param bin raw BIN1
#' @param targets character(1) comma/space-separated masses
#' @param from,to numeric RT window
#' @param ppm_tolerance,mz_tolerance numeric tolerances
#' @return list(x=numeric, y=numeric)
#' @export
bin_to_eic <- function(
  bin, targets, from, to, ppm_tolerance = 20, mz_tolerance = 0.005
) {
  stopifnot(is.raw(bin), is.character(targets), length(targets) == 1)
  .Call(
    "C_bin_to_eic",
    bin, targets, as.numeric(from), as.numeric(to),
    as.numeric(ppm_tolerance), as.numeric(mz_tolerance),
    PACKAGE = "msut"
  )
}

#' Find peaks from an (x,y) chromatogram
#'
#' @param x numeric vector (RT, f64)
#' @param y numeric vector (intensity; will be converted to f32 internally)
#' @param options list, same structure as in `get_peak()`
#' @return data.frame with columns: from,to,rt,integral,intensity,ratio,np
#' @export
find_peaks <- function(x, y, options = NULL) {
  stopifnot(is.numeric(x), is.numeric(y))
  if (length(x) != length(y) || length(x) < 3) {
    stop("x and y must have the same length >= 3")
  }
  out_json <- .Call("C_find_peaks_json", as.numeric(x), as.numeric(y), options, PACKAGE = "msut")
  jsonlite::fromJSON(out_json, simplifyVector = TRUE)
}
