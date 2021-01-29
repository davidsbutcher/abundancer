#' generate_score_heatmap
#'
#' @param score_matrix
#' @param fillOption
#' @param scaleFillValues
#' @param scaleFillLimits
#' @param compFunc
#'
#' @return
#' @export
#'
#' @examples

generate_score_heatmap <-
   function(
      score_matrix,
      fillOption = "C",
      scaleFillValues = c(0,1),
      scaleFillLimits = c(0,signif(max(score_matrix), digits = 1)),
      compFunc = "dotproduct"
   ) {

      canonical_12c <- 0.9893
      canonical_14n <- 0.99636

      if (compFunc == "dotproduct") {

         compFunc_translation <-
            "Cosine\nSimilarity"

      } else if (compFunc == "scoremfa") {

         compFunc_translation <-
            "ScoreMFA"

      } else {

         compFunc_translation <-
            "ERROR, CHECK COMPFUNC"

      }

      score_matrix_melt <-
         reshape2::melt(score_matrix) %>%
         tibble::as_tibble()

      names(score_matrix_melt) <-
         c("14N Abundance", "12C Abundance", "COMP_SCORE")

      optimal_12c <-
         score_matrix_melt %>%
         dplyr::arrange(desc(COMP_SCORE)) %>%
         dplyr::pull(`12C Abundance`) %>%
         dplyr::first()

      optimal_14n <-
         score_matrix_melt %>%
         dplyr::arrange(desc(COMP_SCORE)) %>%
         dplyr::pull(`14N Abundance`) %>%
         dplyr::first()

      ggplot2::ggplot(score_matrix_melt) +
         ggplot2::geom_raster(
            ggplot2::aes(
               x = `14N Abundance`,
               y = `12C Abundance`,
               fill = COMP_SCORE
            )
         ) +
         ggplot2::annotate(
            "point",
            x = optimal_14n,
            y = optimal_12c,
            color = "red",
            shape = 10,
            size = 8
         ) +
         ggplot2::annotate(
            "label",
            x = min(score_matrix_melt$`14N Abundance`),
            y = min(score_matrix_melt$`12C Abundance`),
            color = "red",
            label = paste("14N: ", optimal_14n, "\n12C: ", optimal_12c),
            size = 4,
            vjust="inward",
            hjust="inward"
         ) +
         ggplot2::theme_minimal() +
         ggplot2::scale_fill_viridis_c(
            option = fillOption,
            values = scaleFillValues,
            limits = scaleFillLimits,
            guide =
               "colourbar",
               ggplot2::guide_legend(
                  title = compFunc_translation
               )
         ) +
         ggplot2::theme(
            text = ggplot2::element_text(size = 16)
         )


   }
