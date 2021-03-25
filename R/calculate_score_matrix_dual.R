#' calculate_score_matrix_dual
#'
#' @param MSspectrum An MSnbase Spectrum1 object.
#' @param mz A numeric vector containing m/z values for a spectrum.
#' @param intensity A numeric vector containing intensity values for a spectrum.
#' @param sequence Sequence of the protein represented by the m/z and intensity vectors.
#' @param PTMformula Chemical formula of PTMs of the proteoform represented by the m/z and intensity vectors. Formulas for all PTMs should be combined.
#' @param charge Charge state of the proteoform represented by the m/z and intensity vectors.
#' @param start12C Initial abundance of 12C to use for calculating the score matrix. Should be lower than the expected value.
#' @param start14N Initial abundance of 14N to use for calculating the score matrix. Should be lower than the expected value.
#' @param abundStepCoarse Abundance step to use for the coarse score matrix. Must be larger than abundStepFine.
#' @param abundStepFine Abundance step to use for the fine score matrix. Must be smaller than abundStepCoarse.
#' @param SNR Signal-to-noise cutoff to use for peak picking. See ?MSnbase::pickPeaks.
#' @param method Method to use for peak picking. See ?MSnbase::pickPeaks.
#' @param refineMz Method for m/z refinement for peak picking. See ?MSnbase::pickPeaks.
#' @param k Number of neighboring signals to use for m/z refinement if refineMz = "kNeighbors". See ?MSnbase::pickPeaks.
#' @param binSize Bin size to use for peak binning prior to comparing spectra. See ?MSnbase::bin.
#' @param resolvingPower Resolving power to be used for generating the initial theoretical spectrum. This parameter does not need to match the resolving power of the experimental spectrum.
#' @param compFunc Function to use for comparison of experimental and theoretical spectra. Acceptable values are "dotproduct", "scoremfa", and "scoremfacpp".
#'
#' @return
#' @export
#' @importFrom magrittr %>%
#' @import MSnbase
#' @useDynLib abundancer
#'

calculate_score_matrix_dual <-
   function(
      MSspectrum = NULL,
      mz = NULL,
      intensity = NULL,
      sequence = NULL,
      PTMformula = "C0",
      charge = 1,
      start12C = 0.987,
      start14N = 0.994,
      abundStepCoarse = 0.001,
      abundStepFine = 0.0001,
      SNR = 25,
      method = "MAD",
      refineMz = "kNeighbors",
      k = 2,
      resolvingPower = 300000,
      compFunc = "dotproduct",
      binSize = 0.05
   ) {

      # Assertions --------------------------------------

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
         assertthat::is.number(SNR),
         msg = "SNR is not a length 1 numeric vector"
      )

      assertthat::assert_that(
         compFunc == "dotproduct" | compFunc == "cor" | compFunc == "scoremfa" | compFunc == "scoremfacpp"
      )

      assertthat::assert_that(
         assertthat::is.number(binSize),
         msg = "binSize is not a length 1 numeric vector"
      )

      # Initialize parameters ---------------------------------------------------

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

      peaks_exp_picked <-
         MSnbase::pickPeaks(
            spectrum,
            SNR = SNR,
            method = method,
            refineMz = refineMz,
            k = k
         )

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


      # Add PTM chem form

      chemform <-
         enviPat::mergeform(chemform, PTMformula)


      # Remove C0|H0|N0|O0|P0|S0 from formula to prevent errors

      chemform <-
         stringr::str_remove_all(chemform, "C0|H0|N0|O0|P0|S0")


      # Calculate coarse matrix --------------------------------------------------------

      # Initialize progress bar

      p <- progressr::progressor(steps = 100)

      # Initialize matrix with appropriate dimensions and name rows and cols

      iso_matrix_coarse <-
         matrix(
            ncol = length(seq(start12C, by = abundStepCoarse)),
            nrow = length(seq(start14N, by = abundStepCoarse))
         )

      colnames(iso_matrix_coarse) <-
         seq(start12C, by = abundStepCoarse)

      rownames(iso_matrix_coarse) <-
         seq(start14N, by = abundStepCoarse)

      # Main loop which calculates scores for COARSE matrix.
      # i iterates 12C abundance, j iterates 14N abundance

      for (i in seq_along(seq(start12C, 1, by = abundStepCoarse))) {

         for (j in seq_along(seq(start14N, 1, by = abundStepCoarse))) {

            k <- seq(start12C, 1, by = abundStepCoarse)[i]

            l <- seq(start14N, 1, by = abundStepCoarse)[j]

            isotopes$abundance[index_12C] <- k

            isotopes$abundance[index_13C] <- 1 - k

            isotopes$abundance[index_14N] <- l

            isotopes$abundance[index_15N] <- 1 - l


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

            peaks_IsoPat_picked <-
               MSnbase::pickPeaks(
                  peaks_IsoPat,
                  SNR = SNR,
                  method = method,
                  refineMz = refineMz,
                  k = k
               )

            if (compFunc == "dotproduct" | compFunc == "cor") {

               compSpec_temp <-
                  MSnbase::compareSpectra(
                     peaks_exp_picked,
                     peaks_IsoPat_picked,
                     fun = compFunc,
                     binSize = binSize
                  )

               iso_matrix_coarse[j,i] <- compSpec_temp

            } else if (compFunc == "scoremfa") {

               compSpec_pared <-
                  pare_spectra(
                     peaks_exp_picked,
                     peaks_IsoPat_picked
                  )

               scaling_factor <-
                  max(
                     MSnbase::intensity(compSpec_pared[[1]])
                  )/
                  max(
                     MSnbase::intensity(compSpec_pared[[2]])
                  )

               compSpec_pared[[2]] <-
                  new(
                     "Spectrum1",
                     mz = MSnbase::mz(compSpec_pared[[2]]),
                     intensity = MSnbase::intensity(compSpec_pared[[2]]) * scaling_factor,
                     centroided = TRUE
                  )

               compSpec_temp <-
                  ScoreMFA(
                     vexp_mz = MSnbase::mz(compSpec_pared[[1]]),
                     vtheo_mz = MSnbase::mz(compSpec_pared[[2]]),
                     vrp = rep(resolvingPower, length(MSnbase::mz(compSpec_pared[[1]]))),
                     vexp_sn = MSnbase::intensity(compSpec_pared[[1]]),
                     vtheo_sn = MSnbase::intensity(compSpec_pared[[2]])
                  )

               iso_matrix_coarse[j,i] <- compSpec_temp

            } else if (compFunc == "scoremfacpp") {

               compSpec_pared <-
                  pare_spectra(
                     peaks_exp_picked,
                     peaks_IsoPat_picked
                  )

               scaling_factor <-
                  max(
                     MSnbase::intensity(compSpec_pared[[1]])
                  )/
                  max(
                     MSnbase::intensity(compSpec_pared[[2]])
                  )

               compSpec_pared[[2]] <-
                  new(
                     "Spectrum1",
                     mz = MSnbase::mz(compSpec_pared[[2]]),
                     intensity = MSnbase::intensity(compSpec_pared[[2]]) * scaling_factor,
                     centroided = TRUE
                  )

               compSpec_temp <-
                  ScoreMFA_cpp(
                     vexp_mz = MSnbase::mz(compSpec_pared[[1]]),
                     vtheo_mz = MSnbase::mz(compSpec_pared[[2]]),
                     vrp = rep(resolvingPower, length(MSnbase::mz(compSpec_pared[[1]]))),
                     vexp_sn = MSnbase::intensity(compSpec_pared[[1]]),
                     vtheo_sn = MSnbase::intensity(compSpec_pared[[2]])
                  )

               iso_matrix_coarse[j,i] <- compSpec_temp

            }

         }

         p(
            50/length(seq(start12C, 1, by = abundStepCoarse)),
            message = paste0("Calculating coarse matrix")
         )

      }


      # Get/set parameters for fine matrix ------------------------------------------

      optimal_abund <-
         get_optimal_abundances(iso_matrix_coarse)

      start14Nfine <-
         optimal_abund[[1]] - abundStepCoarse/2

      # if (start14Nfine == 1) start14Nfine <- 1 - abundStepCoarse

      end14Nfine <-
         optimal_abund[[1]] + abundStepCoarse/2

      if (end14Nfine > 1) end14Nfine <- 1

      start12Cfine <-
         optimal_abund[[2]] - abundStepCoarse/2

      # if (start12Cfine == 1) start12Cfine <- 1 - abundStepCoarse

      end12Cfine <-
         optimal_abund[[2]] + abundStepCoarse/2

      if (end12Cfine > 1) end12Cfine <- 1

      # Calculate fine matrix --------------------------------------------------------

      # q <- progressr::progressor(along = seq(start12Cfine, end12Cfine, by = abundStepFine))

      # Initialize matrix with appropriate dimensions and name rows and cols

      iso_matrix_fine <-
         matrix(
            ncol = length(seq(start12Cfine, end12Cfine, by = abundStepFine)),
            nrow = length(seq(start14Nfine, end14Nfine, by = abundStepFine))
         )

      colnames(iso_matrix_fine) <-
         seq(start12Cfine, end12Cfine, by = abundStepFine)

      rownames(iso_matrix_fine) <-
         seq(start14Nfine, end14Nfine, by = abundStepFine)

      # Main loop which calculates scores for COARSE matrix.
      # i iterates 12C abundance, j iterates 14N abundance

      for (i in seq_along(seq(start12Cfine, end12Cfine, by = abundStepFine))) {

         for (j in seq_along(seq(start14Nfine, end14Nfine, by = abundStepFine))) {

            k <- seq(start12Cfine, end12Cfine, by = abundStepFine)[i]

            l <- seq(start14Nfine, end14Nfine, by = abundStepFine)[j]

            isotopes$abundance[index_12C] <- k

            isotopes$abundance[index_13C] <- 1 - k

            isotopes$abundance[index_14N] <- l

            isotopes$abundance[index_15N] <- 1 - l


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

            peaks_IsoPat_picked <-
               MSnbase::pickPeaks(
                  peaks_IsoPat,
                  SNR = SNR,
                  method = method,
                  refineMz = refineMz,
                  k = k
               )

            if (compFunc == "dotproduct" | compFunc == "cor") {

               compSpec_temp <-
                  MSnbase::compareSpectra(
                     peaks_exp_picked,
                     peaks_IsoPat_picked,
                     fun = compFunc,
                     binSize = binSize
                  )

               iso_matrix_fine[j,i] <- compSpec_temp

            } else if (compFunc == "scoremfa") {

               compSpec_pared <-
                  pare_spectra(
                     peaks_exp_picked,
                     peaks_IsoPat_picked
                  )

               scaling_factor <-
                  max(
                     MSnbase::intensity(compSpec_pared[[1]])
                  )/
                  max(
                     MSnbase::intensity(compSpec_pared[[2]])
                  )

               compSpec_pared[[2]] <-
                  new(
                     "Spectrum1",
                     mz = MSnbase::mz(compSpec_pared[[2]]),
                     intensity = MSnbase::intensity(compSpec_pared[[2]]) * scaling_factor,
                     centroided = TRUE
                  )

               compSpec_temp <-
                  ScoreMFA(
                     vexp_mz = MSnbase::mz(compSpec_pared[[1]]),
                     vtheo_mz = MSnbase::mz(compSpec_pared[[2]]),
                     vrp = rep(resolvingPower, length(MSnbase::mz(compSpec_pared[[1]]))),
                     vexp_sn = MSnbase::intensity(compSpec_pared[[1]]),
                     vtheo_sn = MSnbase::intensity(compSpec_pared[[2]])
                  )

               iso_matrix_fine[j,i] <- compSpec_temp

            } else if (compFunc == "scoremfacpp") {

               compSpec_pared <-
                  pare_spectra(
                     peaks_exp_picked,
                     peaks_IsoPat_picked
                  )

               scaling_factor <-
                  max(
                     MSnbase::intensity(compSpec_pared[[1]])
                  )/
                  max(
                     MSnbase::intensity(compSpec_pared[[2]])
                  )

               compSpec_pared[[2]] <-
                  new(
                     "Spectrum1",
                     mz = MSnbase::mz(compSpec_pared[[2]]),
                     intensity = MSnbase::intensity(compSpec_pared[[2]]) * scaling_factor,
                     centroided = TRUE
                  )

               compSpec_temp <-
                  ScoreMFA_cpp(
                     vexp_mz = MSnbase::mz(compSpec_pared[[1]]),
                     vtheo_mz = MSnbase::mz(compSpec_pared[[2]]),
                     vrp = rep(resolvingPower, length(MSnbase::mz(compSpec_pared[[1]]))),
                     vexp_sn = MSnbase::intensity(compSpec_pared[[1]]),
                     vtheo_sn = MSnbase::intensity(compSpec_pared[[2]])
                  )

               iso_matrix_fine[j,i] <- compSpec_temp

            }


         }

         p(
            50/length(seq(start12Cfine, end12Cfine, by = abundStepFine)),
            message = paste0("Calculating fine matrix")
         )

      }

      return(
         list(
            iso_matrix_coarse,
            iso_matrix_fine
         )
      )

   }
