#include <cpp11/doubles.hpp>
#include <cpp11/integers.hpp>
#include <cpp11/list.hpp>
#include <vector>

using namespace cpp11::literals;

[[cpp11::register]]
cpp11::writable::doubles findDistValueByRowColInd(cpp11::doubles distance, int num_row, cpp11::integers row_inds,  cpp11::integers col_inds) {
  int row_inds_len = row_inds.size();
  int col_inds_len = col_inds.size();
  cpp11::writable::doubles res(row_inds_len * col_inds_len);
  
  int i = 0;
  int dist_ind; 
  
  for (int row = 0; row < row_inds_len; row++) {
    int row_ind = row_inds[row];
    for (int col = 0; col < col_inds_len; col++) {
      int col_ind = col_inds[col];
      
      if(row_ind == col_ind){
        res[i] = 0;
      }
      else{
        int row_ind_new; 
        int col_ind_new; 
        
        if(col_ind > row_ind) {
          int row_ind_tmp = row_ind; 
          int col_ind_tmp = col_ind;
          row_ind_new = col_ind_tmp;
          col_ind_new = row_ind_tmp;
        }
        else{
          row_ind_new = row_ind;
          col_ind_new = col_ind;          
        }
        dist_ind = ((unsigned long long) num_row) * (col_ind_new - 1) + row_ind_new - 0.5 * (1 + col_ind_new) * ((unsigned long long) col_ind_new) - 1;
        res[i] = distance[dist_ind]; 
      }
      i++; 
    }
  }
  
  return res;
}

std::vector<double> all_finite(cpp11::doubles x) {
  std::vector<double> res;
  for (int i = 0; i < x.size(); i++) {
    if (x[i] < R_PosInf) {
      res.push_back(x[i]);
    }
  }
  return res;
}

[[cpp11::register]]
cpp11::writable::list smallest_dist_rho_order_coords(cpp11::doubles ordered_rho, cpp11::doubles ordered_coords) {
  int sample_size = ordered_rho.size();
  int dim_num = ordered_coords.size() / sample_size;
  cpp11::writable::doubles smallest_dist(sample_size);
  cpp11::writable::doubles nearest_higher_density_sample(sample_size);
  
  double current_dist;
  
  for (int cell_ind = 0; cell_ind < sample_size; cell_ind ++) {
    smallest_dist[cell_ind] = R_PosInf;
    nearest_higher_density_sample[cell_ind] = cell_ind;
    
    if(cell_ind == sample_size - 1) { // assign the last distance to the highest density peak cell
      std::vector<double> all_finite_vals = all_finite(smallest_dist);
      auto maximal =  std::max_element(all_finite_vals.begin(), all_finite_vals.end());
      smallest_dist[cell_ind] = all_finite_vals[maximal - all_finite_vals.begin()]; 
      nearest_higher_density_sample[cell_ind] = maximal - all_finite_vals.begin();
    }
    
    for (int higher_local_density_cell_ind = cell_ind + 1; higher_local_density_cell_ind < sample_size; higher_local_density_cell_ind ++) {
      cpp11::writable::doubles source_coord(dim_num);
      cpp11::writable::doubles target_coord(dim_num);
      
      current_dist = 0;
      double tmp;
      for(int dim_num_tmp = 0; dim_num_tmp < dim_num; dim_num_tmp ++) {
        source_coord[dim_num_tmp] = ordered_coords[cell_ind + dim_num_tmp * sample_size];
        target_coord[dim_num_tmp] = ordered_coords[higher_local_density_cell_ind + dim_num_tmp  * sample_size];
        
        tmp = source_coord[dim_num_tmp] - target_coord[dim_num_tmp];
        current_dist += tmp * tmp;
      }
      
      current_dist = sqrt(current_dist);
      if(smallest_dist[cell_ind] > current_dist) {
        smallest_dist[cell_ind] = current_dist;
        nearest_higher_density_sample[cell_ind] = higher_local_density_cell_ind;
      }
    }
  }
  return {"smallest_dist"_nm = smallest_dist,
          "nearest_higher_density_sample"_nm = nearest_higher_density_sample};
}
