library(stringr)
library(purrr)

mk_mapper <- function(patch) {
  m <- str_split(patch, "\n")[[1]] |>
    str_match("^@@ *-([0-9]+),([0-9]+) *\\+([0-9]+),([0-9]+) *@@") |>
    na.omit() |>
    _[, -1, drop = FALSE]
  storage.mode(m) <- "integer"

  ends <- m[, 1L] + m[, 2L]
  shifts <- cumsum(c(0L, m[, 4L] - m[, 2L]))

  function(i) {
    pos <- cut(i, c(0, ends, Inf), labels = FALSE, right = FALSE)
    i + shifts[pos]
  }
}

map_lines <- function(warnings, mapper) {
  str_replace_all(warnings, "line=[0-9]+", function(m) {
    str_c("line=", mapper(as.integer(str_sub(m, 6))))
  })
}

lintr_filter <- function(files) {
  changed_files <- files |>
    keep(\(f) nchar(f$patch) > 0) |>
    map_chr("filename")

  base_files <- list.files("base", recursive = TRUE)
  base_exclude <- setdiff(base_files, changed_files)
  base <- lintr::lint_package("base", exclusions = base_exclude)

  head_files <- list.files("head", recursive = TRUE)
  head_exclude <- setdiff(head_files, changed_files)
  head <- lintr::lint_package("head", exclusions = head_exclude)

  mappers <- map(files, \(f) mk_mapper(f$patch)) |>
    setNames(map_chr(files, "filename"))

  base <- structure(
    map(base, function(w) {
      w$line_number <- mappers[[w$filename]](w$line_number)
      w
    }),
    class = class(base)
  )

  setdiff(capture.output(head), capture.output(base))
}
