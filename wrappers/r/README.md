# msut

## Install

You can install the R wrapper directly from GitHub:

Using **remotes**:

```r
remotes::install_github("josoriom/msut", subdir = "wrappers/r", force = TRUE, upgrade = "never")
```

## Parse an mzML (parse_mzml)

```r
path <- "/path/to/file.mzML"
bin <- readBin(path, "raw", file.info(path)\$size)
mzfile <- msut::parse_mzml(bin)
```

## Optional: view bin as data frames

```r
df <- msut::bin_to_df(mzfile)
chroms <- df\$Ok\$run\$chromatograms # columns: index, id, array_length, time_array, intensity_array
```

## Build the list you want to quantify

get_peaks_from_chrom() expects a data frame with columns:

- idx — chromatogram index (0-based; first chrom = 0)
- rt — target retention time (center)
- window — RT window half-width used to score the peak
- id — optional, only for readability in your records

### Minimal manual example

```r
items <- data.frame(
  idx    = c(0L, 5L, 18L),             # 0-based indices
  id     = chroms\$id\[c(1, 6, 19)],     # optional
  rt     = c(1.20, 2.85, 3.50),
  window = c(0.40, 0.50, 0.60),
  stringsAsFactors = FALSE
)
```

## Run peak picking from Chromatogram

```r
library(parallel)
cores <- max(1L, detectCores(logical = FALSE) - 1L)
opts <- list(
  integral_threshold  = 50,
  intensity_threshold = 1e3,
  width_threshold     = 5L,
  noise               = 100,
  auto_noise          = FALSE,
  allow_overlap       = FALSE,
  window_size         = 11L,
  sn_ratio            = 3L
)

peaks <- msut::get_peaks_from_chrom(mzml, items, options = opts, cores = cores)
```

## Run peak picking from EIC

```r
path <- "path/to/your/file.mzML"

targets <- data.frame(
  id = c(
    "1-methylhistidine",
    "Leucine",
    "Glutamine",
    "Alanine",
    "Arginine",
    "Asparagine"
  ),
  rt = c(2.61, 4.15, 1.91, 2.39, 2.59, 1.75),
  mz = c(340.1404, 302.1499, 317.1244, 260.1030, 345.1670, 303.1088),
  ranges = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2),
  stringsAsFactors = FALSE
)

bin <- readBin(path, "raw", file.info(path)$size)
mzml <- msut::parse_mzml(bin)

peaks <- msut::get_peaks_from_eic(mzml, targets, cores = 2)

```
