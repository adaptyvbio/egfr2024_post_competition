packages <- c(
  "ggplot2", "data.table", "showtext", "scales", "ggtext",
  "dplyr", "logger", "optparse", "jsonlite"
)

install_if_missing <- function(pkgs) {
  to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(to_install)) {
    install.packages(to_install, repos = "https://cloud.r-project.org")
  }
}

install_if_missing(packages)

if (!requireNamespace("sysfonts", quietly = TRUE)) {
  install.packages("sysfonts", repos = "https://cloud.r-project.org")
}

if (!requireNamespace("showtextdb", quietly = TRUE)) {
  install.packages("showtextdb", repos = "https://cloud.r-project.org")
}

message("R plotting dependencies installed.")

packages <- c(
  "ggplot2", "data.table", "showtext", "scales", "ggtext",
  "dplyr", "logger", "optparse", "jsonlite"
)

install_if_missing <- function(pkgs) {
  to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(to_install)) {
    install.packages(to_install, repos = "https://cloud.r-project.org")
  }
}

install_if_missing(packages)

if (!requireNamespace("sysfonts", quietly = TRUE)) {
  install.packages("sysfonts", repos = "https://cloud.r-project.org")
}

if (!requireNamespace("showtextdb", quietly = TRUE)) {
  install.packages("showtextdb", repos = "https://cloud.r-project.org")
}

message("R plotting dependencies installed.")


