#' calculate_score_array_DEV
#'
#' @param MSspectrum An MSnbase Spectrum1 object.
#' @param mz A numeric vector containing m/z values for a spectrum.
#' @param intensity A numeric vector containing intensity values for a spectrum.
#' @param sequence Sequence of the protein represented by the m/z and intensity vectors.
#' @param PTMformula Chemical formula of PTMs of the proteoform represented by the m/z and intensity vectors. Formulas for all PTMs should be combined.
#' @param charge Charge state of the proteoform represented by the m/z and intensity vectors.
#' @param range12C
#' @param range14N
#' @param abundStepCoarse Abundance step to use for the coarse score array Must be larger than abundStepFine.
#' @param abundStepFine Abundance step to use for the fine score array Must be smaller than abundStepCoarse.
#' @param SNR Signal-to-noise cutoff to use for peak picking. See ?MSnbase::pickPeaks.
#' @param method Method to use for peak picking. See ?MSnbase::pickPeaks.
#' @param refineMz Method for m/z refinement for peak picking. See ?MSnbase::pickPeaks.
#' @param k Number of neighboring signals to use for m/z refinement if refineMz = "kNeighbors". See ?MSnbase::pickPeaks.
#' @param binSize Bin size to use for peak binning prior to comparing spectra. See ?MSnbase::bin.
#' @param resolvingPower Resolving power to be usedy for generating the initial theoretical spectrum. This parameter does not need to match the resolving power of the experimental spectrum.
#' @param compFunc Function to use for comparison of experimental and theoretical spectra. Acceptable values are "dotproduct" and "scoremfa".
#'
#' @return
#' @export
#' @importFrom magrittr %>%
#' @import MSnbase
#'

calculate_score_array_DEV <-
   function(
      MSspectrum = NULL,
      mz = NULL,
      intensity = NULL,
      sequence = NULL,
      PTMformula = "C0",
      charge = 1,
      range12C = c(0.986,1),
      range14N = c(0.994,1),
      range16O = c(0.996,1),
      abund17Ofixed = c(0.00038),
      range32S = c(0.944,0.953),
      abund33Sfixed = c(0.00749),
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

      # Initialize parameters ---------------------------------------------------

      # Get isotope indices -----------------------------------------------------

      # Need to determine these so they can be replaced later

      data(isotopes, package = "enviPat")

      indices <-
         list(
            "12C" = which(isotopes$isotope == "12C") %>% .[[1]],
            "13C" = which(isotopes$isotope == "13C") %>% .[[1]],
            "14N" = which(isotopes$isotope == "14N") %>% .[[1]],
            "15N" = which(isotopes$isotope == "15N") %>% .[[1]],
            "16O" =which(isotopes$isotope == "16O") %>% .[[1]],
            "17O" = which(isotopes$isotope == "17O") %>% .[[1]],
            "18O" = which(isotopes$isotope == "18O") %>% .[[1]],
            "32S" = which(isotopes$isotope == "32S") %>% .[[1]],
            "33S" = which(isotopes$isotope == "33S") %>% .[[1]],
            "34S" = which(isotopes$isotope == "34S") %>% .[[1]],
            "36S" = which(isotopes$isotope == "36S") %>% .[[1]]
         )


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


      # Calculate coarse array for C and N ---------------------------------------------------

      # Initialize progress bar

      p <- progressr::progressor(steps = 100)

      # Generate sequences for isotope abundances

      seq_coarse <- list()

      seq_coarse[["12C"]] <- seq(from = range12C[1], to = range12C[2], by = abundStepCoarse)
      seq_coarse[["14N"]] <- seq(from = range14N[1], to = range14N[2], by = abundStepCoarse)
      seq_coarse[["16O"]] <- seq(from = range16O[1], to = range16O[2], by = abundStepCoarse)
      seq_coarse[["32S"]] <- seq(from = range32S[1], to = range32S[2], by = abundStepCoarse)


      # Initialize array with appropriate dimensions and name rows and cols

      iso_array_coarse_4D <-
         array(
            dim = c(
               length(seq_coarse[["12C"]]),
               length(seq_coarse[["14N"]]),
               length(seq_coarse[["16O"]]),
               length(seq_coarse[["32S"]])
            ),
            dimnames =
               list(
                  seq_coarse[["12C"]],
                  seq_coarse[["14N"]],
                  seq_coarse[["16O"]],
                  seq_coarse[["32S"]]
               )
         )

      # Main loop which calculates scores for first coarse array for carbon/nitrogen
      # i iterates 12C abundance, j iterates 14N abundance

      for (i in seq_along(seq_coarse[["12C"]])) {

         for (j in seq_along(seq_coarse[["14N"]])) {

            for (k in seq_along(seq_coarse[["16O"]])) {

               for (l in seq_along(seq_coarse[["32S"]])) {

                  # Set isotope abundances for this iteration -------------------------------

                  isotopes$abundance[indices[["12C"]]] <- seq_coarse[["12C"]][i]
                  isotopes$abundance[indices[["13C"]]] <- 1 - seq_coarse[["12C"]][i]
                  isotopes$abundance[indices[["14N"]]] <- seq_coarse[["14N"]][j]
                  isotopes$abundance[indices[["15N"]]] <- 1 - seq_coarse[["14N"]][j]
                  isotopes$abundance[indices[["16O"]]] <- seq_coarse[["16O"]][k]
                  isotopes$abundance[indices[["17O"]]] <- abund17Ofixed
                  isotopes$abundance[indices[["18O"]]] <- 1 - seq_coarse[["16O"]][k] - abund17Ofixed
                  isotopes$abundance[indices[["32S"]]] <- seq_coarse[["32S"]][l]
                  isotopes$abundance[indices[["33S"]]] <- abund33Sfixed
                  isotopes$abundance[indices[["34S"]]] <- 1 - seq_coarse[["32S"]][l] - abund33Sfixed
                  isotopes$abundance[indices[["36S"]]] <- 0

                  # Check for errors with abundances

                  if (isotopes$abundance[indices[["18O"]]] < 0) {

                     isotopes$abundance[indices[["18O"]]] <- 0

                  }

                  if (isotopes$abundance[indices[["34S"]]] < 0) {

                     isotopes$abundance[indices[["34S"]]] <-  0

                  }

                  # Check for sums of abundances != 0

                  # if (
                  #    isotopes$abundance[indices[["16O"]]] +
                  #    isotopes$abundance[indices[["16O"]]] +
                  #    isotopes$abundance[indices[["16O"]]] != 1
                  # ) {
                  #
                  #    iso_array_coarse_4D[i,j,k,l] <- NA
                  #    next()
                  #
                  # }
                  #
                  # if (
                  #    isotopes$abundance[indices[["32S"]]] +
                  #    isotopes$abundance[indices[["33S"]]] +
                  #    isotopes$abundance[indices[["34S"]]] != 1
                  # ) {
                  #
                  #    iso_array_coarse_4D[i,j,k,l] <- NA
                  #    next()
                  #
                  # }

                  # Generate and compare theo. and obs. spectra, get score ------------------

                  iso_array_coarse_4D[i,j,k,l] <-
                     generate_and_compare_spectra(
                        peaks_exp_picked,
                        isotopes = isotopes,
                        chemform = chemform,
                        charge = charge,
                        resolvingPower = resolvingPower,
                        compFunc = compFunc,
                        SNR = SNR,
                        method = method,
                        refineMz = refineMz,
                        k = k,
                        binSize = binSize
                     )

               }

            }

         }

         p(
            100/length(seq_coarse[["12C"]]),
            message = paste0("Calculating coarse array")
         )

      }



      # # Get/set parameters for fine array ---------------------------------
      #
      # # Create ranges for isotope abundances
      #
      # coarse_optvals <-
      #    c(
      #       get_array_max_names(iso_array_coarse_CN)
      #    ) %>%
      #    purrr::set_names("12C", "14N")
      #
      # range_fine <-
      #    purrr::map(
      #       as.list(coarse_optvals),
      #       ~c(
      #          .x - abundStepCoarse, # Adjust to abundStepCoarse/2 if too expensive
      #          .x + abundStepCoarse
      #       ) %>%
      #          {if (.[[2]]>1) c(.[[1]], 1) else .}
      #    )
      #
      # # If range_fine for any isotope is wider than the range specified in the
      # # function argument, replace with function argument
      #
      # range_argument <-
      #    list(
      #       range12C,
      #       range14N
      #    ) %>%
      #    purrr::set_names("12C", "14N")
      #
      # range_fine <-
      #    purrr::map2(
      #       range_fine,
      #       range_argument,
      #       ~{if (abs(.x[1] - .x[2]) > abs(.y[1] - .y[2])) .y else .x}
      #    )
      #
      # # Create ranges for 34S, 33S, 16O, 18O
      #
      #
      #
      #
      #
      #
      #
      #
      #
      # # Create sequences to use for making array using isotope ranges
      #
      # seq_fine <-
      #    purrr::map(
      #       range_fine,
      #       ~seq(.x[[1]], .x[[2]], by = abundStepFine)
      #    )
      #
      # # Initialize array with appropriate dimensions and name rows and cols
      #
      # iso_array_fine <-
      #    array(
      #       dim = c(
      #          length(seq_fine[["12C"]]),
      #          length(seq_fine[["14N"]]),
      #          length(seq_fine[["16O"]]),
      #          length(seq_fine[["32S"]]),
      #          length(seq_fine[["34S"]]),
      #          length(seq_fine[["33S"]])
      #       ),
      #       dimnames =
      #          list(
      #             seq_fine[["12C"]],
      #             seq_fine[["14N"]],
      #             seq_fine[["16O"]],
      #             seq_fine[["32S"]],
      #             seq_fine[["34S"]],
      #             seq_fine[["33S"]]
      #          )
      #    )
      #
      #
      # # Calculate fine array --------------------------------------------------------
      #
      # for (i in seq_along(seq_fine[["12C"]])) {
      #
      #    for (j in seq_along(seq_fine[["14N"]])) {
      #
      #       for (k in seq_along(seq_fine[["16O"]])) {
      #
      #          for (l in seq_along(seq_fine[["32S"]])) {
      #
      #             for (m in seq_along(seq_fine[["34S"]])) {
      #
      #                for (n in seq_along(seq_fine[["33S"]])) {
      #
      #                   if (seq_fine[["32S"]][[l]] + seq_fine[["33S"]][[n]] + seq_fine[["34S"]][[m]] != 1) next()
      #
      #                   if (seq_fine[["16O"]][[k]] + seq_fine[["18O"]][[k]] + abund17O != 1) next()
      #
      #                   isotopes$abundance[indices[["12C"]]] <- seq_fine[["12C"]][[i]]
      #
      #                   isotopes$abundance[indices[["13C"]]] <- 1 - seq_fine[["12C"]][[i]]
      #
      #                   isotopes$abundance[indices[["14N"]]] <- seq_fine[["14N"]][[j]]
      #
      #                   isotopes$abundance[indices[["15N"]]] <- 1 - seq_fine[["14N"]][[j]]
      #
      #                   isotopes$abundance[indices[["16O"]]] <- seq_fine[["16O"]][[k]]
      #
      #                   isotopes$abundance[indices[["17O"]]] <- abund17O
      #
      #                   isotopes$abundance[indices[["18O"]]] <- seq_fine[["18O"]][[k]]
      #
      #                   isotopes$abundance[indices[["32S"]]] <- seq_fine[["32S"]][[l]]
      #
      #                   isotopes$abundance[indices[["33S"]]] <- seq_fine[["33S"]][[n]]
      #
      #                   isotopes$abundance[indices[["34S"]]] <- seq_fine[["34S"]][[m]]
      #
      #                   # Calculate theoretical isotopic distribution for the current set
      #                   # of isotopic abundances (i and j)
      #
      #                   iso_array_fine[i,j,k,l,m,n] <-
      #                      generate_and_compare_spectra(
      #                         peaks_exp_picked,
      #                         isotopes = isotopes,
      #                         chemform = chemform,
      #                         charge = charge,
      #                         resolvingPower = resolvingPower,
      #                         compFunc = compFunc,
      #                         SNR = SNR,
      #                         method = method,
      #                         refineMz = refineMz,
      #                         k = k,
      #                         binSize = binSize
      #                      )
      #
      #                }
      #
      #             }
      #
      #          }
      #
      #       }
      #
      #    }
      #
      #    p(
      #       50/length(seq_fine[['12C']]),
      #       message = paste0("Calculating fine array")
      #    )
      #
      # }


      return(
         list(
            iso_array_coarse_4D
         )
      )

   }
