#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector gaussianLocalDensity(NumericVector distance, int nrow, double dc) {
   int size = distance.size();
   NumericVector half(size);
   for (int i = 0; i < size; i++) {
      double combOver = distance[i] / dc;
      double negSq = pow(combOver, 2) * -1;
      half[i] = exp(negSq);
   }     
   int ncol = nrow;
   
   NumericVector result(nrow);
   
   int i = 0;
   for (int col = 0; col < ncol; col++) {
      for (int row = col + 1; row < nrow; row++) {
         double temp = half[i];
         result[row] += temp;
         result[col] += temp;
         i++;
      }
   }
   
   return result;
}

// [[Rcpp::export]]
NumericVector nonGaussianLocalDensity(NumericVector distance, int nrow, double dc) {
   int ncol = nrow;
   NumericVector result(nrow);
   int i = 0;
   for (int col = 0; col < ncol; col++) {
      for (int row = col + 1; row < nrow; row++) {
         if (distance[i] < dc) {
            result[row] += 1;
            result[col] += 1;
         } else {
            // do nothing
         }
         i++;
      }
   }
   return result;
}