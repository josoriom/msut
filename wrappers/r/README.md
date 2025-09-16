# msut

## Install

You can install the R wrapper directly from GitHub:

Using **remotes**:

```r
remotes::install_github("josoriom/msut", subdir = "wrappers/r")
```

## Parse an mzML (parse_mzml)

```r
path <- "/path/to/file.mzML"
bin <- readBin(path, "raw", file.info(path)\$size)
file <- msut::parse_mzml(bin)
```

## Optional: view bin as data frames

```r
df <- msut::bin_to_df(file)
chroms <- df\$Ok\$run\$chromatograms
```

## Run peak picking from Chromatogram

```r
library(parallel)
transitions <- data.frame(
  idx    = c(0L, 5L, 18L),
  id     = chroms\$id\[c(1, 6, 19)],
  rt     = c(1.20, 2.85, 3.50),
  window = c(0.40, 0.50, 0.60),
  stringsAsFactors = FALSE
)

cores <- max(1L, detectCores(logical = FALSE) - 1L)

peaks <- msut::get_peaks_from_chrom(file, transitions, cores = cores, auto_noise = FALSE)
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
file <- msut::parse_mzml(bin)
peaks <- msut::get_peaks_from_eic(file, targets, cores = 2)
```
