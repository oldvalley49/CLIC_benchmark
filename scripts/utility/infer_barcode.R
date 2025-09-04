library(stringr)

infer_barcode <- function(s) {
  ifelse(str_detect(s, "_"),
         str_replace(s, "^[^_]+_", ""),
         s)
}