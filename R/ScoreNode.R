#' Compute conditional entropy
#' @description
#' For discrete variables, the computation is based on the
#' empirical frequency plus a complexity cost (computed as the
#' Normalized Maximum Likelihood).
#' @param X [a vector]
#' A vector that contains the observational data of the variable.
#' @param df_conditioning [a data frame]
#' The data frame of the observations of the set of conditioning variables
#' @return A list that contains :
#' \itemize{
#' \item entro: The estimation of conditional entropy without the
#' complexity cost.
#' \item entrok: The estimation of conditional entropy with the
#' complexity cost (\eqn{entrok = entro + cplx}).}
#' @export
#' @useDynLib miic
#' @importFrom stats density sd
#'
#' @examples
#' library(miic)
#' N <- 1000
#' X <- sample(x = c(1, 2, 3, 4, 5), size = N, replace = TRUE)
#' Y <- sample(x = c(1, 2, 3), size = N, replace = TRUE)
#' res <- computeEntropy(as.factor(X), df_conditioning = as.factor(Y))
#' message("H(X|Y) = ", res$entro)
#' message("Hk(X|Y) = ", res$entrok)
#'
computeEntropy <- function(X, df_conditioning = NULL) {
  input_data <- data.frame(X, df_conditioning)

  n_samples <- nrow(input_data)
  n_nodes <- ncol(input_data)

  is_continuous <- sapply(input_data, is.numeric)


  # Numeric factor matrix, level starts from 0
  input_factor <- as.matrix(apply(
    input_data, 2,
    function(x) (as.numeric(factor(x, levels = unique(x))) - 1)
  ))
  max_level_list <- as.numeric(apply(input_factor, 2, max)) + 1
  # Data list, numeric for continuous columns, -1 for discrete columns
  input_double <- matrix(nrow = n_samples, ncol = n_nodes)
  # Order list, order(column) for continuous columns (index starting from 0),
  # -1 for discrete columns
  input_order <- matrix(nrow = n_samples, ncol = n_nodes)
  for (i in c(1:n_nodes)) {
    input_double[, i] <- rep_len(-1, n_samples)
    input_order[, i] <- rep_len(-1, n_samples)
  }

  arg_list <- list(
    "is_continuous" = is_continuous,
    "levels" = max_level_list,
    "n_eff" = -1,
    "n_nodes" = n_nodes,
    "n_samples" = n_samples
  )

  cpp_input <- list(
    "factor" = as.vector(input_factor),
    "double" = as.vector(input_double),
    "order" = as.vector(input_order)
  )
  # Call cpp code
  rescpp <- mydiscretizeEntropy(cpp_input, arg_list)
  result <- list()

  result$entro <- rescpp$entro
  result$entrok <- rescpp$entrok
  return(result)
}