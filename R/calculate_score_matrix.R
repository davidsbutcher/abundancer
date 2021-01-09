#' calculate_score_matrix
#'
#' @param MSspectrum
#' @param mz
#' @param intensity
#' @param sequence
#' @param charge
#' @param start12C
#' @param start14N
#' @param abundStep
#' @param isoAbundCutoff
#' @param SNR
#' @param resolvingPower
#' @param compFunc
#' @param binSize
#' @param hClust_height
#'
#' @return
#' @export
#' @importFrom magrittr %>%
#' @import MSnbase
#'
#' @examples

calculate_score_matrix <-
   function(
      MSspectrum = NULL,
      mz = NULL,
      intensity = NULL,
      sequence = NULL,
      charge = 1,
      start12C = 0.98,
      start14N = 0.98,
      abundStep = 0.001,
      isoAbundCutoff = 5,
      SNR = 10,
      resolvingPower = 300000,
      compFunc = "dotproduct",
      binSize = 0.05,
      hClust_height = 0.005
   ) {

      # Assertions ---------

      if (!is.null(MSspectrum)) {

         assertthat::assert_that(
            class(MSspectrum)[[1]] == "Spectrum1",
            msg = "MSspectrum not an MSnbase Spectrum1 object"
         )

      }

      if (!is.null(mz)) {

         assertthat::assert_that(
            is.numeric(mz),
            msg = "mz is not a numeric vector"
         )

      }

      if (!is.null(intensity)) {

         assertthat::assert_that(
            is.numeric(intensity),
            msg = "intensity is not a numeric vector"
         )

      }

      assertthat::assert_that(
         assertthat::is.string(sequence),
         msg = "Sequence is not a string"
      )

      assertthat::assert_that(
         assertthat::is.count(charge),
         msg = "Charge is not a positive integer"
      )

      assertthat::assert_that(
         assertthat::is.number(start12C),
         msg = "start12C is not a length 1 numeric vector"
      )

      assertthat::assert_that(
         start12C >= 0 && start12C <= 1,
         msg = "start12C is not in the range 0 to 1"
      )

      assertthat::assert_that(
         assertthat::is.number(start14N),
         msg = "start14N is not a length 1 numeric vector"
      )

      assertthat::assert_that(
         start14N >= 0 && start14N <= 1,
         msg = "start14N is not in the range 0 to 1"
      )

      assertthat::assert_that(
         assertthat::is.number(abundStep),
         msg = "abundStep is not a length 1 numeric vector"
      )

      assertthat::assert_that(
         assertthat::is.number(isoAbundCutoff),
         msg = "isoAbundCutoff is not a length 1 numeric vector"
      )

      assertthat::assert_that(
         assertthat::is.number(SNR),
         msg = "SNR is not a length 1 numeric vector"
      )

      assertthat::assert_that(
         compFunc == "dotproduct" | compFunc == "cor" | compFunc == "scoremfa"
      )

      assertthat::assert_that(
         assertthat::is.number(binSize),
         msg = "binSize is not a length 1 numeric vector"
      )



      # Calculate matrix --------------------------------------------------------

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


      # Extract spectrum from raw file ------------------------------------------

      if (
         !is.null(mz) &
         !is.null(intensity)
      ) {

         spectrum <-
            new(
               "Spectrum1",
               mz = mz,
               intensity = intensity,
               centroided = FALSE
            )

      } else if (!is.null(MSspectrum)) {

         spectrum <-
            MSspectrum

      }

      # Peak picking ------------------------------------------------------------

      peaks <-
         MSnbase::pickPeaks(
            spectrum,
            SNR = SNR,
            method = "MAD",
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

      # Remove C0|H0|N0|O0|P0|S0 from formula to prevent errors

      if (
         stringr::str_detect(chemform, "C0|H0|N0|O0|P0|S0") == TRUE
      ) {
         chemform <-
            stringr::str_remove_all(chemform, "C0|H0|N0|O0|P0|S0")
      }


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

      if (compFunc == "dotproduct") {

         for (i in seq_along(seq(start12C, 1, by = abundStep))) {

            for (j in seq_along(seq(start14N, 1, by = abundStep))) {

               k <- seq(start12C, 1, by = abundStep)[i]

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

               # Cluster temporary isotopic distribution

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
                  # dplyr::arrange(desc(abundance)) %>%
                  dplyr::ungroup()

               peaks_IsoPat <-
                  new(
                     "Spectrum1",
                     mz = isopat_cluster$`m/z`,
                     intensity = isopat_cluster$abundance,
                     centroided = TRUE
                  )

               compSpec_temp <-
                  MSnbase::compareSpectra(
                     peaks,
                     peaks_IsoPat,
                     fun = compFunc,
                     binSize = binSize
                  )

               iso_matrix[j,i] <- compSpec_temp

            }

            p(message = paste0("Calculating matrix\n", k))

         }

      } else if (compFunc == "scoremfa") {

         for (i in seq_along(seq(start12C, 1, by = abundStep))) {

            for (j in seq_along(seq(start14N, 1, by = abundStep))) {

               k <- seq(start12C, 1, by = abundStep)[i]

               l <- seq(start14N, 1, by = abundStep)[j]

               isotopes$abundance[index_12C] <- k

               isotopes$abundance[index_13C] <- 1 - k

               isotopes$abundance[index_14N] <- l

               isotopes$abundance[index_15N] <- 1 - l

               # Create temporary isotopic distribution

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

               # Cluster temporary isotopic distribution

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
                  # dplyr::arrange(desc(abundance)) %>%
                  dplyr::ungroup()

               # Create experimental vectors for scoreMFA input

               vexp_tbl <-
                  tibble::tibble(
                     `m/z` = MSnbase::mz(peaks),
                     intensity = MSnbase::intensity(peaks)
                  )
                  # dplyr::arrange(desc(intensity))

               vexp_mz <-
                  dplyr::pull(vexp_tbl, `m/z`)

               vexp_sn <-
                  dplyr::pull(vexp_tbl, intensity)


               # Create theoretical vectors for scoreMFA input, match
               # lengths to experimental vectors

               vtheo_mz <-
                  isopat_cluster %>%
                  dplyr::pull(`m/z`) %>%
                  fix_vector_length(length(vexp_mz))

               vtheo_sn <-
                  isopat_cluster %>%
                  dplyr::pull(abundance) %>%
                  fix_vector_length(length(vexp_mz))

               # Compare spectra with scoreMFA

               compSpec_temp <-
                  ScoreMFA(
                     vexp_mz = vexp_mz,
                     vtheo_mz = vtheo_mz,
                     vrp = rep(resolvingPower, length(vexp_mz)),
                     vexp_sn = vexp_sn,
                     vtheo_sn = vtheo_sn
                  )

               iso_matrix[j,i] <- compSpec_temp

            }

            # Increment the progress bar

            p(message = paste0("Calculating matrix\n", k))

         }

      }

      return(iso_matrix)

   }
