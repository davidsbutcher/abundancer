#' make_spectrum_peakpick
#'
#' @param mz
#' @param intensity
#' @param SNR
#' @param binSize
#'
#' @return
#'
#' @examples

make_spectrum_peakpick <-
   function(
      mz,
      intensity,
      SNR = 10,
      method = "MAD",
      refineMz = "kNeighbors",
      k = 2,
      binSize = 0.05
   ) {

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
