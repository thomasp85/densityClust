#include <Rcpp.h>
using namespace Rcpp;

NumericMatrix distToMatrix(NumericVector x, int size) {
   NumericMatrix df(size, size); 
   int nrow = df.nrow();
   int ncol = df.ncol();
   
   int i = 0;
   for (int col = 0; col < ncol; col++) {
      for (int row = col + 1; row < nrow; row++) {
         df(row, col) = x[i];
         df(col, row) = x[i];
         i++;
      }
   }
   
   for (int cell = 0; cell < ncol; cell++) {
      df(cell, cell) = 0;
   }
   
   return df;
}

// [[Rcpp::export]]
NumericVector gaussianLocalDensity(NumericVector distance, int nrow, double dc) {
   // The computations computed on each element of `distance` are expensive 
   // enough that it's faster to compute them for each element of distance and
   // then convert that to matrix form, rather than the other way around.
   int size = distance.size();
   NumericVector half(size);
   for (int i = 0; i < size; i++) {
      double combOver = distance[i] / dc;
      double negSq = pow(combOver, 2) * -1;
      half[i] = exp(negSq);
   }     
   NumericMatrix comb = distToMatrix(half, nrow);
   nrow = comb.nrow();
   int ncol = comb.ncol();
   
   NumericVector result(nrow);
   for (int row = 0; row < nrow; row++) {
      double sum = 0;
      for (int col = 0; col < ncol; col++) {
         double temp = comb(row, col);
         sum += temp;
      }
      result[row] = sum;
   }
   return result;
}

// [[Rcpp::export]]
NumericVector nonGaussianLocalDensity(NumericVector distance, int size, double dc) {
   // The computations computed in the non-Gaussian case are trivial and first
   // converting to a matrix is equally fast
   NumericMatrix comb = distToMatrix(distance, size);
   int nrow = comb.nrow();
   int ncol = comb.ncol();
   NumericVector result(nrow);
   for (int row = 0; row < nrow; row++) {
      double sum = 0;
      for (int col = 0; col < ncol; col++) {
         double temp = comb(row, col);
         if (temp < dc) {
            sum += 1;
         }
      }
      result[row] = sum - 1;
   }
   return result;
}