# R/msut.R

.pack_opts <- function(options) {
  if (is.null(options)) return(NULL)
  if (is.raw(options)) {
    if (length(options) != 48L) stop("options raw blob must be length 48")
    return(options)
  }
  if (!is.list(options) || is.null(names(options))) stop("options must be a named list in snake_case")

  allow <- c(
    "integral_threshold","intensity_threshold","width_threshold","noise",
    "auto_noise","allow_overlap","window_size","sn_ratio"
  )
  bad <- setdiff(names(options), allow)
  if (length(bad)) stop(paste0("unrecognized option(s): ", paste(bad, collapse = ", ")))

  num <- function(k, d) { v <- options[[k]]; if (is.null(v) || is.na(v[1])) d else as.numeric(v[1]) }
  int <- function(k, d) { v <- options[[k]]; if (is.null(v) || is.na(v[1])) d else as.integer(v[1]) }
  lgc <- function(k, d=FALSE) {
    v <- options[[k]]
    if (is.null(v)) return(as.logical(d))
    if (!is.logical(v) || length(v) != 1 || is.na(v)) stop(paste0("option '", k, "' must be logical TRUE/FALSE"))
    v
  }

  raw48 <- raw(48)
  con <- rawConnection(raw48, "r+"); on.exit(close(con), add=TRUE)
  w <- function(val, off, what=c("double","int")) {
    seek(con, off)
    if (what[1] == "double") writeBin(as.double(val), con, size=8, endian="little")
    else writeBin(as.integer(val), con, size=4, endian="little")
  }

  w(num("integral_threshold",  NaN), 0,  "double")
  w(num("intensity_threshold", NaN), 8,  "double")
  w(int("width_threshold",     0L),  16, "int")
  w(num("noise",               NaN), 24, "double")
  w(if (lgc("auto_noise",   FALSE)) 1L else 0L, 32, "int")
  w(if (lgc("allow_overlap",FALSE)) 1L else 0L, 36, "int")
  w(int("window_size",         0L),  40, "int")
  w(int("sn_ratio",            0L),  44, "int")

  rawConnectionValue(con)
}

#' Parse mzML into compact bin (raw bytes)
#' @param data raw vector with mzML file content
#' @return raw vector (BIN1 bytes)
#' @export
parse_mzml <- function(data) {
  stopifnot(is.raw(data))
  .Call("C_parse_mzml", data, PACKAGE="msut")
}

#' Convert bin blob back to JSON text
#' @param bin raw BIN1 blob
#' @return character(1) JSON
#' @export
bin_to_json <- function(bin) {
  stopifnot(is.raw(bin))
  .Call("C_bin_to_json", bin, PACKAGE="msut")
}

#' Convert bin blob into R-friendly list/data.frames
#' @param bin raw BIN1 blob
#' @return list with $Ok$run$chromatograms and $Ok$run$spectra as data.frames
#' @export
bin_to_df <- function(bin) {
  stopifnot(is.raw(bin))
  txt <- bin_to_json(bin)
  x <- jsonlite::fromJSON(txt, simplifyVector=TRUE, simplifyDataFrame=TRUE, simplifyMatrix=TRUE)
  if (is.null(x$Ok)) return(x)

  ensure_list_col <- function(df, col) {
    if (nrow(df) == 0) { df[[col]] <- I(vector("list", 0)); return(df) }
    v <- df[[col]]; if (is.null(v)) v <- vector("list", nrow(df))
    if (!is.list(v)) v <- lapply(seq_len(nrow(df)), function(i) v[[i]])
    df[[col]] <- I(v); df
  }

  ok <- x$Ok; run <- ok$run

  chrom <- run$chromatograms
  if (is.null(chrom)) {
    chrom <- data.frame(index=integer(), array_length=integer(),
                        time_array=I(list()), intensity_array=I(list()),
                        id=character(), stringsAsFactors=FALSE)
  } else if (!is.data.frame(chrom)) chrom <- as.data.frame(chrom, stringsAsFactors=FALSE)
  chrom <- ensure_list_col(chrom, "time_array")
  chrom <- ensure_list_col(chrom, "intensity_array")

  spec <- run$spectra
  if (is.null(spec)) {
    spec <- data.frame(index=integer(), array_length=integer(),
                       ms_level=integer(), polarity=integer(), spectrum_type=integer(),
                       retention_time=numeric(), scan_window_lower_limit=numeric(),
                       scan_window_upper_limit=numeric(), total_ion_current=numeric(),
                       base_peak_intensity=numeric(), base_peak_mz=numeric(),
                       mz_array=I(list()), intensity_array=I(list()),
                       stringsAsFactors=FALSE)
  } else if (!is.data.frame(spec)) spec <- as.data.frame(spec, stringsAsFactors=FALSE)
  spec <- ensure_list_col(spec, "mz_array")
  spec <- ensure_list_col(spec, "intensity_array")

  run$chromatograms <- chrom
  run$spectra <- spec
  ok$run <- run
  x$Ok <- ok
  x
}

#' Get peak from a chromatogram
#' @param x numeric vector of retention times
#' @param y numeric vector of intensities
#' @param rt numeric target retention time
#' @param range numeric RT window size
#' @param integral_threshold numeric integrated area threshold (or NaN)
#' @param intensity_threshold numeric height threshold (or NaN)
#' @param width_threshold integer minimum width (points)
#' @param noise numeric noise level (or NaN)
#' @param auto_noise logical auto noise TRUE/FALSE
#' @param allow_overlap logical allow overlaps TRUE/FALSE
#' @param window_size integer window size (odd recommended)
#' @param sn_ratio integer minimum S/N
#' @return list(from,to,rt,integral,intensity,ratio,np)
#' @export
get_peak <- function(
  x, y, rt, range,
  integral_threshold=NaN, intensity_threshold=NaN, width_threshold=0L,
  noise=NaN, auto_noise=FALSE, allow_overlap=FALSE, window_size=0L, sn_ratio=0L
) {
  stopifnot(is.numeric(x), is.numeric(y))
  if (length(x) != length(y) || length(x) < 3) stop("x and y must have the same length (>= 3)")
  if (!is.logical(auto_noise) || length(auto_noise) != 1 || is.na(auto_noise)) stop("auto_noise must be logical TRUE/FALSE")
  if (!is.logical(allow_overlap) || length(allow_overlap) != 1 || is.na(allow_overlap)) stop("allow_overlap must be logical TRUE/FALSE")

  opt <- .pack_opts(list(
    integral_threshold=integral_threshold, intensity_threshold=intensity_threshold,
    width_threshold=width_threshold, noise=noise,
    auto_noise=auto_noise, allow_overlap=allow_overlap,
    window_size=window_size, sn_ratio=sn_ratio
  ))

  out_json <- .Call("C_get_peak", as.numeric(x), as.numeric(y), as.numeric(rt), as.numeric(range), opt, PACKAGE="msut")
  jsonlite::fromJSON(out_json, simplifyVector=TRUE)
}

#' Get peaks for many (rt,mz,range) items from BIN1
#' @param bin raw BIN1 blob
#' @param df data.frame with columns: id, rt, mz, ranges
#' @param from_left numeric left RT window
#' @param to_right numeric right RT window
#' @param cores integer number of threads
#' @param integral_threshold numeric integrated area threshold (or NaN)
#' @param intensity_threshold numeric height threshold (or NaN)
#' @param width_threshold integer minimum width (points)
#' @param noise numeric noise level (or NaN)
#' @param auto_noise logical auto noise TRUE/FALSE
#' @param allow_overlap logical allow overlaps TRUE/FALSE
#' @param window_size integer window size (odd recommended)
#' @param sn_ratio integer minimum S/N
#' @return data.frame(id,mz,ort,rt,from,to,intensity,integral, ...)
#' @export
get_peaks_from_eic <- function(
  bin, df, from_left=0.5, to_right=0.5, cores=1L,
  integral_threshold=NaN, intensity_threshold=NaN, width_threshold=0L,
  noise=NaN, auto_noise=FALSE, allow_overlap=FALSE, window_size=0L, sn_ratio=0L
) {
  stopifnot(is.raw(bin))
  if (!is.data.frame(df)) stop("`df` must be a data.frame")
  req <- c("id","rt","mz","ranges"); miss <- setdiff(req, names(df))
  if (length(miss)) stop(paste0("missing columns: ", paste(miss, collapse=", ")))

  id <- as.character(df$id)
  rts <- suppressWarnings(as.numeric(df$rt))
  mzs <- suppressWarnings(as.numeric(df$mz))
  ranges <- suppressWarnings(as.numeric(df$ranges))
  n <- length(id)
  if (!(length(rts)==n && length(mzs)==n && length(ranges)==n)) stop("id, rt, mz, ranges must have the same length")
  id[is.na(id)] <- ""

  if (!is.logical(auto_noise) || length(auto_noise) != 1 || is.na(auto_noise)) stop("auto_noise must be logical TRUE/FALSE")
  if (!is.logical(allow_overlap) || length(allow_overlap) != 1 || is.na(allow_overlap)) stop("allow_overlap must be logical TRUE/FALSE")

  opt <- .pack_opts(list(
    integral_threshold=integral_threshold, intensity_threshold=intensity_threshold,
    width_threshold=width_threshold, noise=noise,
    auto_noise=auto_noise, allow_overlap=allow_overlap,
    window_size=window_size, sn_ratio=sn_ratio
  ))

  out_json <- .Call("C_get_peaks_from_eic",
    bin, as.numeric(rts), as.numeric(mzs), as.numeric(ranges), as.character(id),
    as.numeric(from_left), as.numeric(to_right), opt, as.integer(cores),
    PACKAGE="msut"
  )

  res <- jsonlite::fromJSON(out_json, simplifyVector=TRUE)
  if (!is.data.frame(res)) res <- as.data.frame(res)
  want <- c("id","mz","ort","rt","from","to","intensity","integral")
  present <- intersect(want, names(res)); extras <- setdiff(names(res), present)
  res <- res[, c(present, extras), drop=FALSE]
  rownames(res) <- NULL
  res
}

#' Get peaks from chromatograms in BIN1
#' @param bin raw BIN1 blob
#' @param items list/data.frame with columns: idx or index, rt, window or range
#' @param cores integer number of threads
#' @param integral_threshold numeric integrated area threshold (or NaN)
#' @param intensity_threshold numeric height threshold (or NaN)
#' @param width_threshold integer minimum width (points)
#' @param noise numeric noise level (or NaN)
#' @param auto_noise logical auto noise TRUE/FALSE
#' @param allow_overlap logical allow overlaps TRUE/FALSE
#' @param window_size integer window size (odd recommended)
#' @param sn_ratio integer minimum S/N
#' @return data.frame(index,id,ort,rt,from,to,intensity,integral, ...)
#' @export
get_peaks_from_chrom <- function(
  bin, items, cores=1L,
  integral_threshold=NaN, intensity_threshold=NaN, width_threshold=0L,
  noise=NaN, auto_noise=FALSE, allow_overlap=FALSE, window_size=0L, sn_ratio=0L
) {
  stopifnot(is.raw(bin))
  if (is.null(items) || !(is.list(items) || is.data.frame(items))) stop("items must be a list/data.frame")

  idxs <- suppressWarnings(as.integer(if (!is.null(items$idx)) items$idx else items$index))
  rts  <- suppressWarnings(as.numeric(items$rt))
  wins <- suppressWarnings(as.numeric(if (!is.null(items$window)) items$window else items$range))
  if (length(idxs) != length(rts) || length(wins) != length(rts)) stop("idx, rt, range length mismatch")

  if (!is.logical(auto_noise) || length(auto_noise) != 1 || is.na(auto_noise)) stop("auto_noise must be logical TRUE/FALSE")
  if (!is.logical(allow_overlap) || length(allow_overlap) != 1 || is.na(allow_overlap)) stop("allow_overlap must be logical TRUE/FALSE")

  opt <- .pack_opts(list(
    integral_threshold=integral_threshold, intensity_threshold=intensity_threshold,
    width_threshold=width_threshold, noise=noise,
    auto_noise=auto_noise, allow_overlap=allow_overlap,
    window_size=window_size, sn_ratio=sn_ratio
  ))

  out_json <- .Call("C_get_peaks_from_chrom",
    bin, idxs, rts, wins, opt, as.integer(cores), PACKAGE="msut"
  )

  df <- jsonlite::fromJSON(out_json, simplifyVector=TRUE)
  if (!is.data.frame(df)) df <- as.data.frame(df)
  if (!"ort" %in% names(df)) stop("internal error: 'ort' missing from result")
  want <- c("index","id","ort","rt","from","to","intensity","integral")
  present <- want[want %in% names(df)]; extras <- setdiff(names(df), present)
  df <- df[, c(present, extras), drop=FALSE]
  if ("index" %in% names(df)) df <- df[order(df$index), , drop=FALSE]
  rownames(df) <- NULL
  df
}

#' Extracted Ion Chromatogram (EIC)
#' @param bin raw BIN1 blob
#' @param targets character(1) of comma/space-separated m/z targets
#' @param from numeric start RT
#' @param to numeric end RT
#' @param ppm_tolerance numeric ppm tolerance
#' @param mz_tolerance numeric absolute m/z tolerance
#' @return list(x=numeric, y=numeric)
#' @export
calculate_eic <- function(bin, targets, from, to, ppm_tolerance=20, mz_tolerance=0.005) {
  stopifnot(is.raw(bin), is.character(targets), length(targets) == 1)
  .Call("C_calculate_eic",
    bin, targets, as.numeric(from), as.numeric(to),
    as.numeric(ppm_tolerance), as.numeric(mz_tolerance),
    PACKAGE="msut"
  )
}

#' Find peaks from a chromatogram
#' @param x numeric vector of retention times
#' @param y numeric vector of intensities
#' @param integral_threshold numeric integrated area threshold (or NaN)
#' @param intensity_threshold numeric height threshold (or NaN)
#' @param width_threshold integer minimum width (points)
#' @param noise numeric noise level (or NaN)
#' @param auto_noise logical auto noise TRUE/FALSE
#' @param allow_overlap logical allow overlaps TRUE/FALSE
#' @param window_size integer window size (odd recommended)
#' @param sn_ratio integer minimum S/N
#' @return data.frame(from,to,rt,integral,intensity,ratio,np)
#' @export
find_peaks <- function(
  x, y,
  integral_threshold=NaN, intensity_threshold=NaN, width_threshold=0L,
  noise=NaN, auto_noise=FALSE, allow_overlap=FALSE, window_size=0L, sn_ratio=0L
) {
  stopifnot(is.numeric(x), is.numeric(y))
  if (length(x) != length(y) || length(x) < 3) stop("x and y must have the same length >= 3")
  if (!is.logical(auto_noise) || length(auto_noise) != 1 || is.na(auto_noise)) stop("auto_noise must be logical TRUE/FALSE")
  if (!is.logical(allow_overlap) || length(allow_overlap) != 1 || is.na(allow_overlap)) stop("allow_overlap must be logical TRUE/FALSE")

  opt <- .pack_opts(list(
    integral_threshold=integral_threshold, intensity_threshold=intensity_threshold,
    width_threshold=width_threshold, noise=noise,
    auto_noise=auto_noise, allow_overlap=allow_overlap,
    window_size=window_size, sn_ratio=sn_ratio
  ))

  out_json <- .Call("C_find_peaks_json", as.numeric(x), as.numeric(y), opt, PACKAGE="msut")
  jsonlite::fromJSON(out_json, simplifyVector=TRUE)
}
