#!/usr/bin/env Rscript
# Small unit tests for counts/coldata sample-name alignment logic

align_sample_names <- function(cn, rn) {
  # cn: character vector of column names from counts
  # rn: character vector of rownames from coldata
  if (length(cn) == 0) return(character(0))

  if (all(cn %in% rn)) return(cn)

  # try stripping leading 'X'
  cn_stripX <- sub('^X', '', cn)
  if (all(cn_stripX %in% rn)) return(cn_stripX)

  # if numeric columns, try prefixing with 'Sample'
  cn_pref_sample <- ifelse(grepl('^[0-9]+$', cn_stripX), paste0('Sample', cn_stripX), cn_stripX)
  if (all(cn_pref_sample %in% rn)) return(cn_pref_sample)

  # relaxed matching: remove non-alphanumeric and lowercase
  clean <- function(x) tolower(gsub('[^0-9A-Za-z]', '', x))
  rn_clean_map <- setNames(rn, clean(rn))
  cn_clean <- clean(cn)
  if (all(cn_clean %in% names(rn_clean_map))) {
    return(as.character(rn_clean_map[cn_clean]))
  }

  stop('Could not align counts column names to coldata rownames')
}

cat('Running tests for align_sample_names...\n')

# Test 1: exact match
cn <- c('Sample1', 'Sample2')
rn <- c('Sample1', 'Sample2')
res <- align_sample_names(cn, rn)
stopifnot(identical(res, cn))

# Test 2: numeric columns -> SampleN mapping
cn <- c('1', '2')
rn <- c('Sample1', 'Sample2')
res <- align_sample_names(cn, rn)
stopifnot(identical(res, c('Sample1', 'Sample2')))

# Test 3: leading X
cn <- c('XSample1', 'XSample2')
rn <- c('Sample1', 'Sample2')
res <- align_sample_names(cn, rn)
stopifnot(identical(res, c('Sample1', 'Sample2')))

# Test 4: cleaned mapping (punctuation / case)
cn <- c('sample-1', 'SAMPLE2')
rn <- c('Sample1', 'Sample2')
res <- align_sample_names(cn, rn)
stopifnot(identical(res, c('Sample1', 'Sample2')))

# Test 5: expected failure (no possible mapping)
cn <- c('foo', 'bar')
rn <- c('Sample1', 'Sample2')
failed <- FALSE
tryCatch({
  align_sample_names(cn, rn)
}, error = function(e) {
  failed <<- TRUE
})
stopifnot(isTRUE(failed))

cat('All align_sample_names tests passed.\n')

invisible(NULL)
