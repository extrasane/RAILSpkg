#' @keywords internal
"_PACKAGE"

## Package-level imports. External calls are namespace-qualified in each
## function (Matrix::, survey::, stats::, utils::); only the pipe, the .data
## pronoun, the dplyr verbs used inside pipes, and the `weights` generic are
## imported here.
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom dplyr mutate filter bind_rows
#' @importFrom stats weights
## survey::calibrate(calfun = "raking") calls MASS::ginv, but survey only lists
## MASS in Suggests. Raking is the whole method here, so MASS is a hard
## requirement for us and is imported to make that explicit.
#' @importFrom MASS ginv
NULL
