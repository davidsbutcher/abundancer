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
#' @param SNR
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
      mz = NULL,
      intensity = NULL,
      rawFileDir = NULL,
      rawFileName = NULL,
      scanNum = NULL,
      mzRange = NULL,
      sequence = NULL,
      charge = 1,
      start12C = 0.98,
      start14N = 0.98,
      abundStep = 0.001,
      isoAbundCutoff = 5,
      SNR = 10,
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


      # Extract spectrum from raw file ------------------------------------------

      if (
         !is.null(rawFileDir) &
         !is.null(rawFileName) &
         !is.null(scanNum) &
         !is.null(mzRange)
      ) {

         # Get raw file path -------------------------------------------------------


         rawFilesInDir <-
            fs::dir_ls(
               rawFileDir,
               recurse = TRUE,
               type = "file",
               regexp = c("[.]raw$")
            )

         if (length(rawFilesInDir) == 0) {
            stop("No .raw files found in raw file directory")
         }

         if (rawFilesInDir %>%
             stringr::str_detect(rawFileName) %>%
             any() == FALSE) {
            stop("Input raw file not found in raw file directory")
         }

         rawFile <-
            rawFilesInDir %>%
            stringr::str_subset(rawFileName)

         library(rawR)

         spectrum <-
            rawR::readSpectrum(
               rawfile = rawFile,
               scan = scanNum
            )

         spectrum <-
            new(
               "Spectrum1",
               mz = spectrum[[1]]$mZ[spectrum[[1]]$mZ >= mzRange[[1]] & spectrum[[1]]$mZ <= mzRange[[2]]],
               intensity = spectrum[[1]]$intensity[spectrum[[1]]$mZ >= mzRange[[1]] & spectrum[[1]]$mZ <= mzRange[[2]]],
               centroided = FALSE
            )

      } else if (
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
