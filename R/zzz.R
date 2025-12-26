.onLoad <- function(libname, pkgname) {
  reticulate::py_require(c("pandas", "numpy"))
}