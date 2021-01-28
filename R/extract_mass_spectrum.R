#' extract_mass_spectrum
#'
#' @param rawFileDir
#' @param rawFileName
#' @param scanNum
#' @param mzRange
#'
#' @return
#' @export
#'

extract_mass_spectrum <-
   function(
      rawFileDir = NULL,
      rawFileName = NULL,
      scanNum = NULL,
      mzRange = NULL
   ) {

      # Assertions --------------------------------------------------------------

      assertthat::assert_that(
         assertthat::is.dir(rawFileDir),
         msg = "rawFileDir is not a recognized path"
      )

      assertthat::assert_that(
         fs::is_file(
            fs::path(
               rawFileDir, rawFileName
            )
         ),
         msg = "Raw file not found"
      )

      assertthat::assert_that(
         assertthat::is.readable(
            fs::path(
               rawFileDir, rawFileName
            )
         ),
         msg = "Raw file not readable"
      )

      assertthat::assert_that(
         assertthat::is.count(scanNum),
         msg = "scanNum not a positive integer"
      )

      assertthat::assert_that(
         length(mzRange) == 2,
         msg = "mzRange doesn't have a length of 2"
      )

      assertthat::assert_that(
         is.numeric(mzRange),
         msg = "mzRange not a numeric vector"
      )

      # Extract mass spec -------------------------------------------------------

      rawFilesInDir <-
         fs::dir_ls(
            rawFileDir,
            recurse = TRUE,
            type = "file",
            regexp = c("[.]raw$")
         )

      rawFile <-
         rawFilesInDir %>%
         stringr::str_subset(rawFileName)

      library(rawrr)

      spectrum <-
         rawrr::readSpectrum(
            rawfile = rawFile,
            scan = scanNum
         )

      new(
         "Spectrum1",
         mz = spectrum[[1]]$mZ[spectrum[[1]]$mZ >= mzRange[[1]] & spectrum[[1]]$mZ <= mzRange[[2]]],
         intensity = spectrum[[1]]$intensity[spectrum[[1]]$mZ >= mzRange[[1]] & spectrum[[1]]$mZ <= mzRange[[2]]],
         centroided = FALSE
      )


   }
