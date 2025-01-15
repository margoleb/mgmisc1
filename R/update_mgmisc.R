#' @export update_mgmisc

update_mgmisc <- function(path="/home/mgoleb/software/mgmisc1_0.1.0.30092024_R_x86_64-pc-linux-gnu.tar.gz"){
  unloadNamespace("mgmisc1")
  remove.packages("mgmisc1")
  install.packages(path)
  library(mgmisc1)
}
