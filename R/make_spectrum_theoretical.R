#' make_spectrum_theoretical
#'
#' @param score_matrix
#'
#' @return
#'
#' @examples

make_spectrum_theoretical <-
   function(
      score_matrix,
      sequence = NULL,
      charge = NULL,
      isoAbundCutoff = 1,
      hClust_height = 0.005,
      xrange = NULL
   ) {

      # Get isotope indices -----------------------------------------------------

      data(isotopes, package = "enviPat")

      index_12C <-
         which(isotopes$isotope == "12C") %>%
         .[[1]]

      index_13C <-
         which(isotopes$isotope == "13C") %>%
         .[[1]]

      index_14N <-
         which(isotopes$isotope == "14N") %>%
         .[[1]]

      index_15N <-
         which(isotopes$isotope == "15N") %>%
         .[[1]]

      # Get optimal abundances --------------------------------------------------

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

      # Replace standard abundances with optimal values --------------------------

      isotopes$abundance[index_12C] <- optimal_12c

      isotopes$abundance[index_13C] <- 1 - optimal_12c

      isotopes$abundance[index_14N] <- optimal_14n

      isotopes$abundance[index_15N] <- 1 - optimal_14n


      # Make chemical formula for sequence --------------------------------------

      chemform <-
         OrgMassSpecR::ConvertPeptide(sequence) %>%
         unlist() %>%
         purrr::imap(~paste0(.y, .x, collapse = '')) %>%
         paste0(collapse = '') %>%
         enviPat::mergeform(
            paste0("H", charge)
         )

      # Remove C0|H0|N0|O0|P0|S0 from formula to prevent errors

      if (
         stringr::str_detect(chemform, "C0|H0|N0|O0|P0|S0") == TRUE
      ) {
         chemform <-
            stringr::str_remove_all(chemform, "C0|H0|N0|O0|P0|S0")
      }


      # Generate theoretical isotopic distribution ------------------------------

      isopat_temp <-
         enviPat::isopattern(
            isotopes,
            chemform,
            charge = charge,
            verbose = F
         ) %>%
         .[[1]] %>%
         tibble::as_tibble() %>%
         dplyr::filter(abundance > isoAbundCutoff)


      # Cluster theoretical peaks -----------------------------------------------

      isopat_cluster <-
         dplyr::mutate(
            isopat_temp,
            cluster =
               cutree(
                  hclust(
                     dist(
                        `m/z`, method = "maximum"), method = "centroid"
                  ),
                  h = hClust_height
               )
         ) %>%
         dplyr::group_by(cluster) %>%
         dplyr::summarise(
            `m/z` = mean(`m/z`),
            abundance = sum(abundance), # All clustered peak abundances are summed
            charge = mean(charge)
         ) %>%
         dplyr::filter(charge %% 1 == 0) %>% # Remove all incorrectly clustered peaks
         dplyr::ungroup()


      # Generate spectrum -------------------------------------------------------

      ggplot2::ggplot(
         isopat_cluster,
         ggplot2::aes(`m/z`, abundance)
      ) +
         ggplot2::geom_segment(
            ggplot2::aes(x = `m/z`, xend = `m/z`, y = 0, yend = abundance),
            color = 'red',
            size = 1,
            alpha = 1
         ) +
         ggplot2::theme_minimal() +
         ggplot2::theme(
            text = ggplot2::element_text(size = 16)
         ) +
         ggplot2::labs(
            x = "m/z",
            y = "Relative Abundance"
         ) +
         ggplot2::xlim(xrange)


   }
