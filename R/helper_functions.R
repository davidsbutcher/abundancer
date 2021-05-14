ovln <- function(theo_mz, exp_mz, rp, rp_mult = 6){
   sd <- exp_mz / (rp_mult * rp)
   Zn <- -abs(theo_mz - exp_mz) / (2 * sd)
   return(2*pnorm(Zn, mean=0, sd=1))
}

ovlhn <- function(theo_mz, exp_mz, theo_sn, exp_sn, rp, rp_mult = 6){

   # Added by DSB to handle bug with missing values in exp_sn

   if (length(theo_sn) != length(exp_sn)) return(0)
   if (is.na(theo_sn) | is.na(exp_sn)) return(0)
   if (is.null(theo_sn) | is.null(exp_sn)) return(0)

   if(theo_sn==exp_sn){
      return(ovln(theo_mz,exp_mz, rp, rp_mult = rp_mult))
   } else {
      sd <- exp_mz / (rp_mult * rp)

      mu1t <- theo_mz - trunc(exp_mz)
      mu2t <- exp_mz - trunc(exp_mz)

      Z <- (mu1t + mu2t) / 2 + ((sd * sd) * log(theo_sn / exp_sn)) / (mu2t - mu1t)
      #Note log=ln

      muy <- min(mu1t, mu2t)
      mux <- max(mu1t, mu2t)

      if (muy == mu1t) {
         sn_muy <- theo_sn
         sn_mux <- exp_sn
      } else{
         sn_muy <- exp_sn
         sn_mux <- theo_sn
      }

      #Overlap of Normals without h
      ovln_sem_sn = pnorm(Z, mux, sd) + 1 - pnorm(Z, muy, sd)

      Iarea <- sn_mux*pnorm(Z, mux, sd) + sn_muy*(1 - pnorm(Z, muy, sd))

      Tarea <- sn_mux*(1 - pnorm(Z, mux, sd)) +  sn_muy*pnorm(Z, muy, sd)

      return(Iarea/Tarea)
   }

}

ScoreMFA <- function(vexp_mz, vtheo_mz, vrp, vexp_sn, vtheo_sn, rp_mult = 6) {

   if (length(vexp_sn) != length(vtheo_sn)) return(0)

   vsd <- vexp_mz / (rp_mult*vrp)
   ntheo <- length(vtheo_mz)
   nexp <- length(vexp_mz)
   if(nexp < ntheo){
      vtheo_mz <- vtheo_mz[1:nexp]
      vtheo_sn <- vtheo_sn[1:nexp]
   }
   vovlhn <- numeric(nexp)
   for(i in 1:ntheo){
      vovlhn[i] <- ovlhn(theo_mz=vtheo_mz[i], exp_mz=vexp_mz[i], theo_sn=vtheo_sn[i], exp_sn=vexp_sn[i], rp=vrp[i], rp_mult = rp_mult)
   }
   Score <- sum(vtheo_sn*vovlhn)/sum(vtheo_sn)
   return(Score)

}

fix_vector_length <-
   function(
      vector, desired_length, fill = c(0)
   ) {

      if (length(vector) < desired_length) {

         while (length(vector) < desired_length) {
            vector[[length(vector)+1]] <- fill
         }

      } else if (length(vector) > desired_length) {

         while (length(vector) > desired_length) {
            vector[[length(vector)]] <- NULL
         }

      }

      return(vector)

   }


pare_spectra <-
   function(
      spectrum1,
      spectrum2
   ) {

      # Get info about spec1

      maxint_index_spec1 <-
         which(
            MSnbase::intensity(spectrum1) == max(MSnbase::intensity(spectrum1))
         )

      mz_maxint <-
         MSnbase::mz(spectrum1)[maxint_index_spec1]


      peaks_belowmax_spec1 <-
         length(which(MSnbase::mz(spectrum1) < mz_maxint))

      peaks_abovemax_spec1 <-
         length(which(MSnbase::mz(spectrum1) > mz_maxint))


      # Get info about spec2

      maxint_index_spec2 <-
         which(
            MSnbase::intensity(spectrum2) == max(MSnbase::intensity(spectrum2))
         )

      mz_maxint <-
         MSnbase::mz(spectrum2)[maxint_index_spec2]

      peaks_belowmax_spec2 <-
         length(which(MSnbase::mz(spectrum2) < mz_maxint))

      peaks_abovemax_spec2 <-
         length(which(MSnbase::mz(spectrum2) > mz_maxint))

      # Choose the lower number of peaks below max as the allowed number

      if (peaks_belowmax_spec1 > peaks_belowmax_spec2) {

         peaks_allowed_belowmax <- peaks_belowmax_spec2

      } else if (peaks_belowmax_spec1 < peaks_belowmax_spec2) {

         peaks_allowed_belowmax <- peaks_belowmax_spec1

      } else if (peaks_belowmax_spec1 == peaks_belowmax_spec2) {

         peaks_allowed_belowmax <- peaks_belowmax_spec1

      }

      # Choose the lower number of peaks above max as the allowed number


      if (peaks_abovemax_spec1 > peaks_abovemax_spec2) {

         peaks_allowed_abovemax <- peaks_abovemax_spec2

      } else if (peaks_abovemax_spec1 < peaks_abovemax_spec2) {

         peaks_allowed_abovemax <- peaks_abovemax_spec1

      } else if (peaks_abovemax_spec1 == peaks_abovemax_spec2) {

         peaks_allowed_abovemax <- peaks_abovemax_spec1

      }

      # Get indices for allowed values for mz and intensity for both spectra

      start_index_spec1 <-
         maxint_index_spec1 - peaks_allowed_belowmax

      end_index_spec1 <-
         maxint_index_spec1 + peaks_allowed_abovemax

      indices_spec1 <-
         seq(start_index_spec1, end_index_spec1)


      start_index_spec2 <-
         maxint_index_spec2 - peaks_allowed_belowmax

      end_index_spec2 <-
         maxint_index_spec2 + peaks_allowed_abovemax

      indices_spec2 <-
         seq(start_index_spec2, end_index_spec2)

      # Pare both spectra to the allowed numbers of peaks above and below max

      new_spectrum1 <-
         new(
            "Spectrum1",
            mz =
               MSnbase::mz(spectrum1)[indices_spec1],
            intensity =
               MSnbase::intensity(spectrum1)[indices_spec1],
            centroided = TRUE
         )

      new_spectrum2 <-
         new(
            "Spectrum1",
            mz =
               MSnbase::mz(spectrum2)[indices_spec2],
            intensity =
               MSnbase::intensity(spectrum2)[indices_spec2],
            centroided = TRUE
         )

      # return new spectra

      return(
         list(
            new_spectrum1,
            new_spectrum2
         )
      )

   }

pare_spectra_closest_match <-
   function(
      spectrum_theo,
      spectrum_obs,
      resPowerMS1 = 300000,
      isoWinMultiplier = 1
   ) {

      mz_at_max_theo <-
         which(MSnbase::intensity(spectrum_theo) == max(MSnbase::intensity(spectrum_theo)))

      mz_at_max_obs <-
         which(MSnbase::intensity(spectrum_obs) == max(MSnbase::intensity(spectrum_obs)))

      indi_obs <-
         MALDIquant::match.closest(
            MSnbase::mz(spectrum_theo),
            MSnbase::mz(spectrum_obs),
            tolerance = (MSnbase::mz(spectrum_theo)[mz_at_max_theo]/resPowerMS1)*isoWinMultiplier
         ) %>%
         unique() %>%
         na.exclude()

      indi_theo <-
         MALDIquant::match.closest(
            MSnbase::mz(spectrum_obs),
            MSnbase::mz(spectrum_theo),
            tolerance = (MSnbase::mz(spectrum_obs)[mz_at_max_obs]/resPowerMS1)*isoWinMultiplier
         ) %>%
         unique() %>%
         na.exclude()


      new_spectrum1 <-
         new(
            "Spectrum1",
            mz =
               MSnbase::mz(spectrum_theo)[indi_theo],
            intensity =
               MSnbase::intensity(spectrum_theo)[indi_theo],
            centroided = TRUE
         )

      new_spectrum2 <-
         new(
            "Spectrum1",
            mz =
               MSnbase::mz(spectrum_obs)[indi_obs],
            intensity =
               MSnbase::intensity(spectrum_obs)[indi_obs],
            centroided = TRUE
         )

      return(
         list(
            new_spectrum1,
            new_spectrum2
         )
      )

   }

dotproduct <- function(x, y) {
   #Taken from MSnbase package, 20210127
   maxlength <-
      min(
         c(length(x), length(y))
      )

   x <- x[1:maxlength]
   y <- y[1:maxlength]

   as.vector(x %*% y) / (sqrt(sum(x*x)) * sqrt(sum(y*y)))
}

ceiling_dec <- function(x, digits=1) round(x + 5*10^(-digits-1), digits)


#' get_array_max_names
#'
#' @param array
#'
#' @return
#' @export
#'
#' @examples
#'

get_array_max_names <-
   function(
      array
   ) {

      # Replace all NA values with 0

      array[is.na(array)] <- 0

      dim_names <-
         dimnames(array)

      # WARNING - IF MULTIPLE VALUES HAVE THE SAME "MAX" SCORE,
      # the top row is arbitrarily chosen. This is a problem!

      max_indices <-
         tibble::as_tibble(which(array == max(array), arr.ind = T)) %>%
         dplyr::slice(1) %>%
         as.list()

      # max_indices <-
      #    as.list(which(array == max(array), arr.ind = T)) %>%
      #    purrr::map(
      #       magrittr::extract,
      #       1
      #    )

      out <-
         purrr::map2(
            dim_names,
            max_indices[1:length(dim_names)], # To account for weird edge cases where length(max_indices) >> length(dim_names)
            ~as.double(.x[[.y]])
         )

      unlist(out)

   }


calculate_mma_ppm <-
   function(
      obs_mass,
      theo_mass
   ) {

      ((obs_mass - theo_mass)/theo_mass) * 1E6

   }
