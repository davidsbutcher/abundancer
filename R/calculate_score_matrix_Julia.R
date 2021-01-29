#' calculate_score_matrix_Julia
#'
#' @param MSspectrum
#' @param mz
#' @param intensity
#' @param sequence
#' @param charge
#' @param start12C
#' @param start14N
#' @param abundStep
#' @param SNR
#' @param resolvingPower
#' @param compFunc
#' @param binSize
#'
#' @return
#' @export
#' @importFrom magrittr %>%
#' @import MSnbase
#'
#' @examples

calculate_score_matrix_Julia <-
   function(
      MSspectrum = NULL,
      mz = NULL,
      intensity = NULL,
      sequence = NULL,
      charge = 1,
      start12C = 0.989,
      start14N = 0.995,
      abundStep = 0.001,
      SNR = 10,
      method = "MAD",
      refineMz = "kNeighbors",
      k = 2,
      resolvingPower = 300000,
      compFunc = "dotproduct",
      binSize = 0.05
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

      # assertthat::assert_that(
      #    assertthat::is.number(isoAbundCutoff),
      #    msg = "isoAbundCutoff is not a length 1 numeric vector"
      # )

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

      # Need to determine these so they can be replaced later

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

      # This is only run if mz and intensity vectors are not provided

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

      # Peak picking for user-supplied/experimental spectrum

      peaks_exp_picked_binned <-
         MSnbase::pickPeaks(
            spectrum,
            SNR = SNR,
            method = method,
            refineMz = refineMz,
            k = k
         ) %>%
         MSnbase::bin(
            binSize = binSize
         )

      peaks_exp_picked_binned_inten <-
         intensity(peaks_exp_picked_binned)[intensity(peaks_exp_picked_binned) != 0]

      # Make chemical formula ---------------------------------------------------

      # For use in calculating theoretical isotopic distributions

      chemform <-
         OrgMassSpecR::ConvertPeptide(sequence) %>%
         unlist() %>%
         purrr::imap(~paste0(.y, .x, collapse = '')) %>%
         paste0(collapse = '') %>%
         enviPat::mergeform(
            paste0("H", charge)
         )

      # Remove C0|H0|N0|O0|P0|S0 from formula to prevent errors

      chemform <-
         stringr::str_remove_all(chemform, "C0|H0|N0|O0|P0|S0")


      # Prepare TIDs for score matrix calc --------------------------------------

      isopat_list <-
         vector(
            "list",
            length = length(seq_along(seq(start12C, 1, by = abundStep)))
         ) %>%
         purrr::set_names(
            seq(start12C, 1, by = abundStep)
         ) %>%
         purrr::map(
            ~vector(
               "list",
               length = length(seq_along(seq(start14N, 1, by = abundStep)))
            ) %>%
               purrr::set_names(
                  seq(start14N, 1, by = abundStep)
               )
         )


      for (i in seq_along(seq(start12C, 1, by = abundStep))) {

         for (j in seq_along(seq(start14N, 1, by = abundStep))) {

            l <- seq(start14N, 1, by = abundStep)[j]

            k <- seq(start12C, 1, by = abundStep)[i]

            isotopes$abundance[index_14N] <- l

            isotopes$abundance[index_15N] <- 1 - l

            isotopes$abundance[index_12C] <- k

            isotopes$abundance[index_13C] <- 1 - k

            # Calculate theoretical isotopic distribution for the current set
            # of isotopic abundances (i and j)

            isopat_temp <-
               enviPat::isopattern(
                  isotopes,
                  chemform,
                  charge = charge,
                  verbose = F
               )

            isopat_cluster <-
               enviPat::envelope(
                  isopat_temp,
                  dmz  = "get",
                  resolution = resolvingPower,
                  verbose = F
               ) %>%
               .[[1]] %>%
               tibble::as_tibble() %>%
               dplyr::filter(abundance > 0)

            peaks_IsoPat <-
               new(
                  "Spectrum1",
                  mz = isopat_cluster$`m/z`,
                  intensity = isopat_cluster$abundance,
                  centroided = FALSE
               )

            # use peak picking on the theoretical isotopic distribution

            peaks_IsoPat_picked_binned <-
               MSnbase::pickPeaks(
                  peaks_IsoPat,
                  SNR = SNR,
                  method = method,
                  refineMz = refineMz,
                  k = k
               ) %>%
               MSnbase::bin(
                  binSize = binSize
               )

            peaks_IsoPat_picked_binned_inten <-
               intensity(peaks_IsoPat_picked_binned)[intensity(peaks_IsoPat_picked_binned) != 0]


            isopat_list[[i]][[j]] <- peaks_IsoPat_picked_binned_inten

         }

      }

      # Calculate score matrix --------------------------------------------------------

      # Initialize matrix with appropriate dimensions and name rows and cols

      score_matrix <-
         matrix(
            ncol = length(seq(start12C, by = abundStep)),
            nrow = length(seq(start14N, by = abundStep))
         )

      colnames(score_matrix) <-
         seq(start12C, by = abundStep)

      rownames(score_matrix) <-
         seq(start14N, by = abundStep)


      # Main loop which calculates scores. i iterates 12C abundance,
      # j iterates 14N abundance. THIS SHOULD BE IMPLEMENTED ENTIRELY IN JULIA

      for (i in seq_along(seq(start12C, 1, by = abundStep))) {

         for (j in seq_along(seq(start14N, 1, by = abundStep))) {

            if (compFunc == "dotproduct") {

               compSpec_temp <-
                  dotproduct(
                     isopat_list[[i]][[j]],
                     peaks_exp_picked_binned_inten
                  )

               score_matrix[j,i] <- compSpec_temp

            } else if (compFunc == "scoremfa") {

               stop("This isn't implemented yet. Stop it")

            }

         }

         p(message = paste0("Calculating matrix\n", k))

      }

      return(score_matrix)

   }
