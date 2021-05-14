generate_and_compare_spectra <-
   function(
      peaks_exp_picked,
      isotopes = NULL,
      chemform = NULL,
      charge = NULL,
      resolvingPower = NULL,
      compFunc = NULL,
      SNR = NULL,
      method = NULL,
      refineMz = NULL,
      k = NULL,
      binSize = NULL
   ) {

      # Calculate theoretical isotopic distribution for the current set
      # of isotopic abundances (i and j)

      isopat_cluster <-
         enviPat::isopattern(
            isotopes,
            chemform,
            charge = charge,
            verbose = F
         ) %>%
         enviPat::envelope(
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

      pared_spec <-
         pare_spectra_closest_match(
            peaks_IsoPat_picked,
            peaks_exp_picked
         )

      if (compFunc == "dotproduct") {

         compSpec_temp <-
            MSnbase::compareSpectra(
               pared_spec[[1]],
               pared_spec[[2]],
               fun = compFunc,
               binSize = binSize
            )

      } else if (compFunc == "scoremfa") {

         scaling_factor <-
            max(
               MSnbase::intensity(pared_spec[[2]])
            )/
            max(
               MSnbase::intensity(pared_spec[[1]])
            )

         pared_spec[[1]] <-
            new(
               "Spectrum1",
               mz = MSnbase::mz(pared_spec[[1]]),
               intensity = MSnbase::intensity(pared_spec[[1]]) * scaling_factor,
               centroided = TRUE
            )

         compSpec_temp <-
            ScoreMFA(
               vexp_mz = MSnbase::mz(pared_spec[[2]]),
               vtheo_mz = MSnbase::mz(pared_spec[[1]]),
               vrp = rep(resolvingPower, length(MSnbase::mz(pared_spec[[1]]))),
               vexp_sn = MSnbase::intensity(pared_spec[[2]]),
               vtheo_sn = MSnbase::intensity(pared_spec[[1]]),
               rp_mult = 6
            )

      }

      return(compSpec_temp)

   }
