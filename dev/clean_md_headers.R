clean_md_headers <- function(input, output = input) {
  x <- readLines(input, warn = FALSE)

  # 1) Strip tabset-specific attributes from headers
  x <- gsub("\\s*\\{\\.tabset[^}]*\\}", "", x)

  # 2) Strip any remaining Pandoc-style header attributes: ### Title {#id .class}
  x <- gsub(
    pattern = "^(#{1,6}[^\\{\\n]*?)\\s*\\{[^}]*\\}\\s*$",
    replacement = "\\1",
    x
  )

  writeLines(x, output, useBytes = TRUE)
  invisible(output)
}
