library(stringr)
library(purrr)

mk_mapper <- function(patch) {
  m <- str_split(patch, "\n")[[1]] |>
    str_match("^([@+ -])(@ *-([0-9]+))?") |>
    _[, -c(1L, 3L), drop = FALSE] |>
    as.data.frame() |>
    setNames(c("t", "s"))

  starts <- as.integer(m$s[m$t == "@"])
  lines <- split(m$t, cumsum(m$t == "@")) |> map(`[`, -1L)
  lens <- lengths(lines)
  ends <- starts + lens

  lmap <- numeric(max(ends))

  pos1 <- 1L
  pos2 <- 1L
  for (ch in seq_along(starts)) {
    gap <- starts[ch] - pos1 - 1L
    lmap[pos1 + seq_len(gap)] <- pos2 + seq_len(gap)

    pos1 <- pos1 + gap
    pos2 <- pos2 + gap
    for (d in lines[[ch]]) {
      if (d %in% c(" ", "-")) pos1 <- pos1 + 1L
      if (d %in% c(" ", "+")) pos2 <- pos2 + 1L
      if (d == " ") lmap[pos1] <- pos2
    }
  }

  gap <- length(lmap) - pos1
  lmap[pos1 + seq_len(gap)] <- pos2 + seq_len(gap)

  shift <- pos2 - pos1

  function(i) {
    ifelse(i <= length(lmap), lmap[i], i + shift)
  }
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
