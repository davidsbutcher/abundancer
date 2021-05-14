#' make_spectrum_array
#'
#' @description
#' Used to create a combined spectrum showing an experimental spectrum and peaks
#' picked from a theoretical spectrum generated using a score matrix, sequence,
#' optional PTM formula, charge, and resolving power.
#'
#' @param mz A numeric vector containing m/z values for an experimental spectrum.
#' @param intensity A numeric vector containing intensity values for an experimental spectrum.
#' @param SNR Signal-to-noise cutoff to use for peak picking. See ?MSnbase::pickPeaks.
#' @param method Method to use for peak picking. See `?MSnbase::pickPeaks`.
#' @param refineMz Method for m/z refinement for peak picking. See ?MSnbase::pickPeaks.
#' @param k Number of neighboring signals to use for m/z refinement if refineMz = "kNeighbors". See ?MSnbase::pickPeaks.
#' @param score_matrix A matrix containing 12C and 14N abundances.
#' @param sequence Sequence of the protein represented by the m/z and intensity vectors.
#' @param PTMformula Chemical formula of PTMs of the proteoform represented by the m/z and intensity vectors. Formulas for all PTMs should be combined.
#' @param charge Charge state of the proteoform represented by the m/z and intensity vectors.
#' @param resolvingPower Resolving power to be used for generating the initial theoretical spectrum. This parameter does not need to match the resolving power of the experimental spectrum.
#'
#' @return
#' @export
#' @importFrom magrittr %>%
#' @import MSnbase
#'

make_spectrum_array <-
   function(
      mz = NULL,
      intensity = NULL,
      SNR = 10,
      method = "MAD",
      refineMz = "kNeighbors",
      k = 2,
      isoAbund = list(),
      sequence = NULL,
      PTMformula = "C0",
      charge = NULL,
      resolvingPower = 300000
   ) {

      # Set isotope abundances -----------------------------------------------------

      data(isotopes, package = "enviPat")

      if (!is.null(names(isoAbund))) {

         if (!any(names(isoAbund) %in% isotopes$isotope)) {
            stop("Isotope specified in isoAbund not found in enviPat isotopes table")
         }

         for (i in seq_along(isoAbund)) {

            isotopes$abundance[which(isotopes$isotope == names(isoAbund)[[i]] )] <- isoAbund[[i]]

         }

      }


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

      chemform <-
         stringr::str_remove_all(chemform, "C0|H0|N0|O0|P0|S0")

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


      # Extract spectrum from raw file ------------------------------------------

      spectrum <-
         new(
            "Spectrum1",
            mz = mz,
            intensity = intensity,
            centroided = FALSE
         )

      # Peak picking ------------------------------------------------------------

      peaks <-
         MSnbase::pickPeaks(
            spectrum,
            SNR = SNR,
            method = method,
            refineMz = refineMz,
            k = k
         )

      # Pare spectra to only common peaks ---------------------------------------

      peaks_pared <-
         pare_spectra_closest_match(
            peaks_IsoPat_picked,
            peaks
         )

      peaks_IsoPat_picked <-
         peaks_pared[[1]]

      peaks <-
         peaks_pared[[2]]

      # Prepare data for plotting -----------------------------------------------

      peaks_tibble <-
         tibble::tibble(
            mz = MSnbase::mz(peaks),
            intensity = MSnbase::intensity(peaks)
         )

      peaks_IsoPat_tibble <-
         tibble::tibble(
            mz = MSnbase::mz(peaks_IsoPat_picked),
            intensity = MSnbase::intensity(peaks_IsoPat_picked)
         )

      scaling_factor <-
         max(peaks_tibble$intensity)/max(peaks_IsoPat_tibble$intensity)

      peaks_IsoPat_tibble <-
         peaks_IsoPat_tibble %>%
         dplyr::mutate(scaled_intensity = intensity*scaling_factor)

      # Make spectrum -----------------------------------------------------------

      ggplot2::ggplot(
         peaks_tibble,
         ggplot2::aes(mz, intensity)
      ) +
         ggplot2::geom_line(
            data =
               tibble::tibble(
                  mz = spectrum@mz,
                  intensity = spectrum@intensity
               ),
            ggplot2::aes(mz, intensity)
         ) +
         ggplot2::geom_point(
            data = peaks_IsoPat_tibble,
            ggplot2::aes(x = mz, y = scaled_intensity),
            color = 'blue',
            size = 3,
            shape = 4,
            alpha = 0.75
         ) +
         ggplot2::geom_point(
            ggplot2::aes(x = mz, y = intensity),
            color = 'red',
            size = 3,
            alpha = 0.5
         ) +
         ggplot2::theme_minimal() +
         ggplot2::theme(
            text = ggplot2::element_text(size = 16)
         ) +
         ggplot2::labs(
            x = "m/z",
            y = "Intensity"
         )


   }
