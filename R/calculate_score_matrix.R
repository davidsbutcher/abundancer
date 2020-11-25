#' calculate_score_matrix
#'
#' @param mz
#' @param intensity
#' @param sequence
#' @param charge
#' @param start12C
#' @param start14N
#' @param abundStep
#' @param isoAbundCutoff
#' @param binSize
#'
#' @return
#' @export
#' @importFrom magrittr %>%
#' @import MSnbase
#'
#' @examples

calculate_score_matrix <-
   function(
      mz,
      intensity,
      sequence = NULL,
      charge = 1,
      start12C = 0.98,
      start14N = 0.98,
      abundStep = 0.001,
      isoAbundCutoff = 5,
      binSize = 0.05
   ) {

      # Initialize progress bar

      p <- progressr::progressor(along = seq(start12C, 1, by = abundStep))

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


      # Peak picking ------------------------------------------------------------

      spectrum <-
         new(
            "Spectrum1",
            mz = mz,
            intensity = intensity,
            centroided = FALSE
         )

      peaks <-
         MSnbase::pickPeaks(
            spectrum,
            SNR = 10,
            refineMz = "kNeighbors",
            k = 2
         )


      # Make chemical formula ---------------------------------------------------

      chemform <-
         OrgMassSpecR::ConvertPeptide(sequence) %>%
         unlist() %>%
         purrr::imap(~paste0(.y, .x, collapse = '')) %>%
         paste0(collapse = '') %>%
         enviPat::mergeform(
            paste0("H", charge)
         )


      # Calculate matrix --------------------------------------------------------

      iso_matrix <-
         matrix(
            ncol = length(seq(start12C, by = abundStep)),
            nrow = length(seq(start14N, by = abundStep))
         )

      colnames(iso_matrix) <-
         seq(start12C, by = abundStep)

      rownames(iso_matrix) <-
         seq(start14N, by = abundStep)


      for (i in seq_along(seq(start12C, 1, by = abundStep))) {

         for (j in seq_along(seq(start14N, 1, by = abundStep))) {

            k <- seq(start14N, 1, by = abundStep)[i]

            l <- seq(start14N, 1, by = abundStep)[j]

            isotopes$abundance[index_12C] <- k

            isotopes$abundance[index_13C] <- 1 - k

            isotopes$abundance[index_14N] <- l

            isotopes$abundance[index_15N] <- 1 - l

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

            peaks_IsoPat <-
               new(
                  "Spectrum1",
                  mz = isopat_temp$`m/z`,
                  intensity = isopat_temp$abundance,
                  centroided = TRUE
               )

            compSpec_temp <-
               MSnbase::compareSpectra(
                  peaks,
                  peaks_IsoPat,
                  fun = "dotproduct",
                  binSize = binSize
               )

            iso_matrix[j,i] <- compSpec_temp

         }

         p(message = paste0("12C abund: ", k))

      }

      return(iso_matrix)

   }
