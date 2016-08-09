#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector findDistValueByRowColInd(NumericVector distance, size_t num_row, NumericVector row_inds,  NumericVector col_inds) {
  size_t row_inds_len = row_inds.size();
  size_t col_inds_len = col_inds.size();
  NumericVector res(row_inds_len * col_inds_len);
  
  size_t i = 0;
  
  for (size_t col = 0; col < col_inds_len; col++) {
    size_t col_ind = col_inds[col];
    for (size_t row = col + 1; row < row_inds_len; row++) {
      size_t row_ind = row_inds[row];

      if(row_ind == col_ind){
        res[i] = 0;
      }
      else{
        size_t row_ind_new; 
        size_t col_ind_new; 

        if(col_ind > row_ind) {
          size_t row_ind_tmp = row_ind; 
          size_t col_ind_tmp = col_ind;
          row_ind_new = col_ind_tmp;
          col_ind_new = row_ind_tmp;
        }
        else{
          row_ind_new = row_ind;
          col_ind_new = col_ind;          
        }
        size_t dist_ind = num_row * (col_ind_new - 1) + row_ind_new - 1/2 * (1 + col_ind_new) * col_ind_new;
        res[i] = distance[dist_ind]; 
      }
      i++; 
    }
  }
  
  return res;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
findDistValueByRowColInd(test, attr(test, 'Size'), 1:100, 1:100)
*/
