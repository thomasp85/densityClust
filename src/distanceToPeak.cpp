#include <cpp11/doubles.hpp>

[[cpp11::register]]
cpp11::writable::doubles distanceToPeakCpp(cpp11::doubles distance, cpp11::doubles rho) {
  int size = rho.size();
  cpp11::writable::doubles peaks(size);
  std::fill(peaks.begin(), peaks.end(), 0.0);
  cpp11::writable::doubles maximum(size);
  std::fill(maximum.begin(), maximum.end(), 0.0);
  
  int i = 0;
  for (int col = 0; col < size; col++) {
    for (int row = col + 1; row < size; row++) {
      double newValue = distance[i];
      double rhoRow = rho[row];
      double rhoCol = rho[col];
      
      if (rhoRow > rhoCol) {
        double peaksCol = peaks[col];
        if (newValue < peaksCol || peaksCol == 0) {
          peaks[col] = newValue;
        }
      } else if (newValue > maximum[col]) {
        maximum[col] = newValue;
      }
      
      if (rhoCol > rhoRow) {
        double peaksRow = peaks[row];
        if (newValue < peaksRow || peaksRow == 0) {
          peaks[row] = newValue;
        }
      } else if (newValue > maximum[row]) {
        maximum[row] = newValue;
      }
      i++;
    }
  }
  
  for (int j = 0; j < size; j++) {
    if (peaks[j] == 0) {
      peaks[j] = double(maximum[j]);
    } else {
      // do nothing, peaks is already min
    }
  } 
  
  return peaks;
}
