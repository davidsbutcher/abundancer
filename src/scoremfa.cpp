#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// from Dirk Eddelbuettel
// https://gallery.rcpp.org/articles/vector-minimum/index.html
// [[Rcpp::export]]
double vecmin(NumericVector x) {
   // Rcpp supports STL-style iterators
   NumericVector::iterator it = std::min_element(x.begin(), x.end());
   // we want the value so dereference 
   return *it;
}

// [[Rcpp::export]]
double vecmax(NumericVector x) {
   // Rcpp supports STL-style iterators
   NumericVector::iterator it = std::max_element(x.begin(), x.end());
   // we want the value so dereference 
   return *it;
}

// [[Rcpp::export]]
bool all_sug(LogicalVector x) {
   // Note the use of is_true to return a bool type.
   return is_true(all(x == TRUE));
}

// [[Rcpp::export]]
NumericVector ovln_cpp_old(
      NumericVector theo_mz,
      NumericVector exp_mz,
      NumericVector rp
) {
   NumericVector sd = exp_mz / (6 * rp);
   NumericVector Zn = -abs(theo_mz - exp_mz) / (2 * sd);
   
   int size = Zn.size();
   NumericVector out(size);
   
   for (int i = 0; i < size; i++) {
      out[i] = 2 * R::pnorm5(Zn[i], 0.0, 1.0, 1, 0);
   }
   
   return out;
}

// This is the good one
// [[Rcpp::export]]
NumericVector ovln_cpp(
      NumericVector theo_mz,
      NumericVector exp_mz,
      NumericVector rp
) {
   
   NumericVector Zn = -Rcpp::abs(theo_mz - exp_mz) / (2 * (exp_mz / (6 * rp)));
   
   NumericVector out = 2 * pnorm(Zn, 0.0, 1.0, 1, 0);
   
   return out;
}

// [[Rcpp::export]]
double ovlhn_cpp_old(
      NumericVector theo_mz_ovlhn,
      NumericVector exp_mz_ovlhn,
      NumericVector theo_sn_ovlhn,
      NumericVector exp_sn_ovlhn,
      NumericVector rp_ovlhn
) {
   
   // if (all_sug(theo_sn_ovlhn == exp_sn_ovlhn) == TRUE) {
   //    return ovln_cpp(theo_mz_ovlhn, exp_mz_ovlhn, rp_ovlhn);
   // } 
   
   NumericVector sd = exp_mz_ovlhn / (6 * rp_ovlhn);
   double ave_sd = mean(sd);
   
   NumericVector mu1t = theo_mz_ovlhn - Rcpp::trunc(exp_mz_ovlhn);
   NumericVector mu2t = exp_mz_ovlhn - Rcpp::trunc(exp_mz_ovlhn);
   
   NumericVector Z = (mu1t + mu2t) / 2 + ((sd * sd) * Rcpp::log(theo_sn_ovlhn / exp_sn_ovlhn)) / (mu2t - mu1t);
   
   double mu1t_min = vecmin(mu1t);
   double mu1t_max = vecmax(mu1t);
   
   double mu2t_min = vecmin(mu2t);
   double mu2t_max = vecmax(mu2t);
   
   NumericVector muy(1);
   NumericVector mux(1);
   
   if (mu1t_min <= mu2t_min) {
      
      muy[0] = mu1t_min;
      
   } else {
      
      muy[0] = mu2t_min;
      
   }
   
   if (mu1t_max >= mu2t_max) {
      
      mux[0] = mu1t_max;
      
   } else {
      
      mux[0] = mu2t_max;
      
   }
   
   NumericVector sn_muy;
   NumericVector sn_mux;
   
   if (all_sug(muy == mu1t)) {
      sn_muy = theo_sn_ovlhn;
      sn_mux = exp_sn_ovlhn;
   } else {
      sn_muy = exp_sn_ovlhn;
      sn_mux = theo_sn_ovlhn;
   }
   
   // the use of ave_sd in these pnorm calls instead of sd is slightly different
   // than Mel's implementation in R and leads to slightly different results
   
   NumericVector ovln_sem_sn(Z.size());
   ovln_sem_sn = (pnorm(Z, mux[0], ave_sd, 1, 0) + 1) - pnorm(Z, muy[0], ave_sd, 1, 0);
   
   NumericVector Iarea(Z.size());
   Iarea = sn_mux * pnorm(Z, mux[0], ave_sd, 1, 0) + sn_muy*(1 - pnorm(Z, muy[0], ave_sd, 1, 0));
   
   NumericVector  Tarea(Z.size());
   Tarea = sn_mux * (1 - pnorm(Z, mux[0], ave_sd, 1, 0)) +  sn_muy*pnorm(Z, muy[0], ave_sd, 1, 0);
   
   double out = (Iarea/Tarea)[0];
   
   return(out);
   
}

// [[Rcpp::export]]
double ovlhn_cpp(
      double theo_mz_ovlhn,
      double exp_mz_ovlhn,
      double theo_sn_ovlhn,
      double exp_sn_ovlhn,
      double rp_ovlhn
) {
   
   // if (all_sug(theo_sn_ovlhn == exp_sn_ovlhn) == TRUE) {
   //    return ovln_cpp(theo_mz_ovlhn, exp_mz_ovlhn, rp_ovlhn);
   // } 
   
   double sd = exp_mz_ovlhn / (6 * rp_ovlhn);

   double mu1t = theo_mz_ovlhn - trunc(exp_mz_ovlhn);
   double mu2t = exp_mz_ovlhn - trunc(exp_mz_ovlhn);
   
   double Z;
   Z = (mu1t + mu2t) / 2 + ((sd * sd) * log(theo_sn_ovlhn / exp_sn_ovlhn)) / (mu2t - mu1t);
   
   double muy;
   double mux;
   
   if (mu1t < mu2t) {
      muy = mu1t;
   } else {
      muy = mu2t;
   }
   
   if (mu1t > mu2t) {
      mux = mu1t;
   } else {
      mux = mu2t;
   }
   
   double sn_muy;
   double sn_mux;
   
   if (muy == mu1t) {
      sn_muy = theo_sn_ovlhn;
      sn_mux = exp_sn_ovlhn;
   } else{
      sn_muy = exp_sn_ovlhn;
      sn_mux = theo_sn_ovlhn;
   }
   

   double ovln_sem_sn;
   ovln_sem_sn = (R::pnorm5(Z, mux, sd, 1, 0) + 1) - R::pnorm5(Z, muy, sd, 1, 0);
   
   double Iarea;
   Iarea = sn_mux * R::pnorm5(Z, mux, sd, 1, 0) + sn_muy*(1 - R::pnorm5(Z, muy, sd, 1, 0));
   
   double  Tarea;
   Tarea = sn_mux * (1 - R::pnorm5(Z, mux, sd, 1, 0)) +  sn_muy*R::pnorm5(Z, muy, sd, 1, 0);
   
   double out = (Iarea/Tarea);
   
   return(out);
   
}

// [[Rcpp::export]]
NumericVector ScoreMFA_cpp(
      NumericVector vexp_mz,
      NumericVector vtheo_mz,
      NumericVector vrp,
      NumericVector vexp_sn,
      NumericVector vtheo_sn
) {
   
   NumericVector vsd = vexp_mz / (6*vrp);
   int ntheo = vtheo_mz.size();
   int nexp = vexp_mz.size();
   
   if(nexp < ntheo) {
      vtheo_mz = vtheo_mz[Range(0, nexp-1)];
      vtheo_sn = vtheo_sn[Range(0, nexp-1)];
   }
   
   NumericVector vovlhn(ntheo);
   
   for (int i = 0; i < ntheo; ++i) {
      vovlhn[i] = ovlhn_cpp(vtheo_mz[i], vexp_mz[i], vtheo_sn[i], vexp_sn[i], vrp[i]);
   }

   NumericVector score(ntheo);
   score = sum(vtheo_sn*vovlhn)/sum(vtheo_sn);
   
   return(score);
   
}

