ovln <- function(theo_mz, exp_mz, rp){
   sd <- exp_mz / (6 * rp)
   Zn <- -abs(theo_mz - exp_mz) / (2 * sd)
   return(2*pnorm(Zn, mean=0, sd=1))
}

ovlhn <- function(theo_mz, exp_mz, theo_sn, exp_sn, rp){

   if (theo_sn == exp_sn){
      return(ovln(theo_mz,exp_mz,rp))
   } else {
      sd <- exp_mz / (6 * rp)

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

ScoreMFA <-
   function(vexp_mz, vtheo_mz, vrp, vexp_sn, vtheo_sn) {

      vsd <- vexp_mz / (6*vrp)
      ntheo <- length(vtheo_mz)
      nexp <- length(vtheo_mz)
      if(nexp < ntheo){
         vtheo_mz <- vtheo_mz[1:nexp]
         vtheo_sn <- vtheo_sn[1:nexp]
      }
      vovlhn <- numeric(nexp)
      for(i in 1:ntheo){
         vovlhn[i] <- ovlhn(theo_mz=vtheo_mz[i], exp_mz=vexp_mz[i], theo_sn=vtheo_sn[i], exp_sn=vexp_sn[i], rp=vrp[i])
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
