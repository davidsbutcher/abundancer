#' make_spectrum_theoretical
#'
#' @description
#' Used to create a theoretical spectrum from a protein sequence.
#'
#' @param sequence Sequence of the protein represented by the m/z and intensity vectors.
#' @param PTMformula Chemical formula of PTMs of the proteoform represented by the m/z and intensity vectors. Formulas for all PTMs should be combined.
#' @param charge Charge state of the proteoform represented by the m/z and intensity vectors.
#' @param SNR Signal-to-noise cutoff to use for peak picking. See ?MSnbase::pickPeaks.
#' @param method Method to use for peak picking. See `?MSnbase::pickPeaks`.
#' @param refineMz Method for m/z refinement for peak picking. See ?MSnbase::pickPeaks.
#' @param k Number of neighboring signals to use for m/z refinement if refineMz = "kNeighbors". See ?MSnbase::pickPeaks.
#' @param resolvingPower Resolving power to be used for generating the initial theoretical spectrum. This parameter does not need to match the resolving power of the experimental spectrum.
#' @param pctAbundCutoff Cutoff of percent of maximum abundance to be used for culling low abundance peaks.
#' @param isoAbund Two-element vector containing abundances of 12C and 14N to use for theoretical modeling.
#' @param mzRange Range to use for X axis of theoretical spectrum.
#'
#' @return
#' @export
#' @importFrom magrittr %>%
#' @import MSnbase
#'

make_spectrum_theoretical <-
   function(
      sequence = NULL,
      PTMformula = "C0",
      charge = NULL,
      SNR = 10,
      method = "MAD",
      refineMz = "kNeighbors",
      k = 2,
      resolvingPower = 300000,
      pctAbundCutoff = 1,
      isoAbund = c(0.9893,0.99636),
      mzRange = NULL
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

      # Replace standard abundances with optimal values --------------------------

      isotopes$abundance[index_12C] <- isoAbund[[1]]

      isotopes$abundance[index_13C] <- 1 - isoAbund[[1]]

      isotopes$abundance[index_14N] <- isoAbund[[2]]

      isotopes$abundance[index_15N] <- 1 - isoAbund[[2]]

      # Make chemical formula for sequence --------------------------------------

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
         )

      # Cluster theoretical peaks -----------------------------------------------

      isopat_cluster <-
         enviPat::envelope(
            isopat_temp,
            dmz  = "get",
            resolution = resolvingPower,
            verbose = F
         ) %>%
         .[[1]] %>%
         tibble::as_tibble()

      peaks_IsoPat <-
         new(
            "Spectrum1",
            mz = isopat_cluster$`m/z`,
            intensity = isopat_cluster$abundance,
            centroided = FALSE
         )

      peaks_IsoPat_picked <-
         MSnbase::pickPeaks(
            peaks_IsoPat,
            SNR = SNR,
            method = method,
            refineMz = refineMz,
            k = k
         )


      # Prepare data for plotting -----------------------------------------------

      peaks_IsoPat_tibble <-
         tibble::tibble(
            mz = MSnbase::mz(peaks_IsoPat_picked),
            intensity = MSnbase::intensity(peaks_IsoPat_picked),
            pctMaxAbund = (intensity/max(MSnbase::intensity(peaks_IsoPat_picked)))*100
         ) %>%
         dplyr::filter(pctMaxAbund >= pctAbundCutoff)

      # Make spectrum -----------------------------------------------------------

      plot <-
      ggplot2::ggplot() +
         ggplot2::geom_segment(
            data = peaks_IsoPat_tibble,
            ggplot2::aes(
               x = mz,
               xend = mz,
               y = 0,
               yend = intensity
            )
         ) +
         ggplot2::geom_point(
            data = peaks_IsoPat_tibble,
            ggplot2::aes(x = mz, y = intensity),
            color = 'blue',
            size = 3,
            shape = 18,
            alpha = 0.35
         ) +
         # ggplot2::geom_point(
         #    ggplot2::aes(x = mz, y = intensity),
         #    color = 'red',
         #    size = 3,
         #    alpha = 0.5
         # ) +
         ggplot2::theme_minimal() +
         ggplot2::theme(
            text = ggplot2::element_text(size = 16)
         ) +
         ggplot2::labs(
            x = "m/z",
            y = "Relative Abundance"
         )

      if (!is.null(mzRange)) {

         plot <-
            plot +
            ggplot2::lims(
               x = mzRange
            )


      }

      return(plot)

   }
