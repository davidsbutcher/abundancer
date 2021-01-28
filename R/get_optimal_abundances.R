#' get_optimal_abundances
#'
#' @param score_matrix
#'
#' @return
#' @export
#'
#' @examples

get_optimal_abundances <-
   function(
      score_matrix
   ){

      score_matrix_melt <-
         reshape2::melt(score_matrix) %>%
         tibble::as_tibble()

      names(score_matrix_melt) <-
         c("14N Abundance", "12C Abundance", "Cosine\nSimilarity")

      optimal_12c <-
         score_matrix_melt %>%
         dplyr::arrange(desc(`Cosine\nSimilarity`)) %>%
         dplyr::pull(`12C Abundance`) %>%
         dplyr::first()

      optimal_14n <-
         score_matrix_melt %>%
         dplyr::arrange(desc(`Cosine\nSimilarity`)) %>%
         dplyr::pull(`14N Abundance`) %>%
         dplyr::first()

      return(c(optimal_14n, optimal_12c))

   }
