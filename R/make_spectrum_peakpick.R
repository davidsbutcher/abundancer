#' make_spectrum_peakpick
#'
#' @param MSspectrum An MSnbase Spectrum1 object.
#' @param mz A numeric vector containing m/z values for a spectrum.
#' @param intensity A numeric vector containing intensity values for a spectrum.
#' @param SNR Signal-to-noise cutoff to use for peak picking. See ?MSnbase::pickPeaks.
#' @param method Method to use for peak picking. See ?MSnbase::pickPeaks.
#' @param refineMz Method for m/z refinement for peak picking. See ?MSnbase::pickPeaks.
#' @param k Number of neighboring signals to use for m/z refinement if refineMz = "kNeighbors". See ?MSnbase::pickPeaks.
#' @param binSize Bin size to use for peak binning prior to comparing spectra. See ?MSnbase::bin.
#'
#' @return
#' @export
#' @importFrom magrittr %>%
#' @import MSnbase
#'

make_spectrum_peakpick <-
   function(
      MSspectrum = NULL,
      mz = NULL,
      intensity = NULL,
      SNR = 10,
      method = "MAD",
      refineMz = "kNeighbors",
      k = 2,
      binSize = 0.05
   ) {

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
            method = method,
            refineMz = refineMz,
            k = k
         )

      # Make spectrum -----------------------------------------------------------

      peaks_tibble <-
         tibble::tibble(
            mz = MSnbase::mz(peaks),
            intensity = MSnbase::intensity(peaks)
         )

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
         # ggplot2::geom_segment(
         #    ggplot2::aes(x = mz, xend = mz, y = 0, yend = intensity),
         #    color = 'red',
         #    size = 1,
         #    alpha = 0.75
         # ) +
         ggplot2::geom_point(
            ggplot2::aes(x = mz, y = intensity),
            color = 'red',
            size = 3,
            alpha = 0.75
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
