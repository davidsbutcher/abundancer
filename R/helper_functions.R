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
