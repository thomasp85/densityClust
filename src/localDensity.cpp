#include <cpp11/doubles.hpp>

[[cpp11::register]]
cpp11::writable::doubles gaussianLocalDensity(cpp11::doubles distance, int nrow, double dc) {
  int size = distance.size();
  cpp11::writable::doubles half(size);
  for (int i = 0; i < size; i++) {
    double combOver = distance[i] / dc;
    double negSq = pow(combOver, 2) * -1;
    half[i] = exp(negSq);
  }     
  int ncol = nrow;
  
  cpp11::writable::doubles result(nrow);
  std::fill(result.begin(), result.end(), 0.0);
  
  int i = 0;
  for (int col = 0; col < ncol; col++) {
    for (int row = col + 1; row < nrow; row++) {
      if(i > distance.size()){
        break;
      }
      double temp = half[i];
      result[row] += temp;
      result[col] += temp;
      i++;
    }
  }
  
  return result;
}

[[cpp11::register]]
cpp11::writable::doubles nonGaussianLocalDensity(cpp11::doubles distance, int nrow, double dc) {
  int ncol = nrow;
  cpp11::writable::doubles result(nrow);
  std::fill(result.begin(), result.end(), 0.0);
  int i = 0;
  for (int col = 0; col < ncol; col++) {
    for (int row = col + 1; row < nrow; row++) {
      if(i > distance.size()){
        break;
      }
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
