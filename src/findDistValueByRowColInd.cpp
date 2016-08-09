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

  // Rcout << "distance is " << distance  << "col_inds index is " << col_inds << std::endl; 
  // Rcout << "row_inds is " << row_inds  << "col_inds index is " << col_inds << std::endl;
  // Rcout << "col_inds_len is " << col_inds_len  << "col_inds_len index is " << col_inds_len  << "res length is " << res.size() << std::endl;
  
  size_t i = 0;
  size_t dist_ind; 

  for (size_t row = 0; row < row_inds_len; row++) {
    size_t row_ind = row_inds[row];
    for (size_t col = 0; col < col_inds_len; col++) {
      size_t col_ind = col_inds[col];
      // Rcout << "col is " << col  << "row is " << row << std::endl;

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
         dist_ind = num_row * (col_ind_new - 1) + row_ind_new - 0.5 * (1 + col_ind_new) * col_ind_new - 1;
         // if(row_ind == 3 && col_ind == 2){
          // Rcout << "num_row * (col_ind_new - 1) is " << num_row * (col_ind_new - 1)  << " 1/2 * (1 + col_ind_new) * col_ind_new is " << 0.5 * (1 + 1) * col_ind_new  << " num_row is " << num_row << " dist_ind is " << num_row * (col_ind_new - 1) + row_ind_new - 1/2 * (1 + col_ind_new) * col_ind_new - 1 << std::endl;

          // Rcout << "col_ind_new is " << col_ind_new  << " row_ind_new index is " << row_ind_new  << " num_row is " << num_row << " dist_ind is " << dist_ind << " distance under the current index " << distance[dist_ind] << std::endl;
         // }
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
