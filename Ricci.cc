/*
 * Copyright 2014 Barry Wardell
 *
 * This file is distributed under the University of Illinois/NCSA
 * Open Source License. See LICENSE file for details.
 *
 */

#include <vector>
#include <complex>
#include <iostream>
#include <sstream>
#include <boost/multi_array.hpp>
#ifdef _OPENMP
#  include <omp.h>
#endif
#include "h1.h"
#include "R2.h"
#include "utils.h"
#include "h5wrapper.h"

using namespace std;
using boost::multi_array;
typedef boost::multi_array_types::extent_range range;

const double M = 1.0;

int main(int argc, char* argv[])
{
  if(argc != 2)
  {
    cout << "Usage: " << argv[0] << " <dir>" << endl;
    exit(1);
  }

#ifdef _OPENMP
  int num_threads = omp_get_max_threads();
  cout << "Running with " << num_threads << " OpenMP threads" << endl;
#endif

  /* Read in first-order fields */
  string dir(argv[1]);
  vector<double> r, f, fp;
  double r0;
  multi_array<complex<double>, 4> h, dh, ddh;
  read_h1(dir, r0, r, f, fp, h, dh, ddh);
  const int lmax = h[1].size()-1;
  const int N = r.size();

  /* Compute the source */
  multi_array<complex<double>,4> src(boost::extents[range(1,11)][lmax+1][lmax+1][N]);
  fill(src.data(), src.data() + src.num_elements(), 0.0);

  /* Loop over i3, l3, m3 */
  vector<ilm_mode> modes;
  for(int l3=0; l3<=lmax; ++l3) {
    for(int m3=0; m3<=l3; ++m3) {
      for(int i3=1; i3<=10; ++i3) {
        if(((i3<=7) && isOdd(l3+m3)) || ((i3>7) && isEven(l3+m3)))
          continue;
        modes.push_back(ilm_mode({i3, l3, m3}));
      }
    }
  }

  int status = 0;
  const int status_frequency = 4*num_threads;
  const int num_modes = modes.size();
#pragma omp parallel for
  for(vector<ilm_mode>::iterator mode = modes.begin(); mode < modes.end(); mode++) {
    int i3 = mode->i;
    int l3 = mode->l;
    int m3 = mode->m;

    /* Sum over l1, l2, m1, m2 = m3-m1 */
    vector<complex<double>> tmp(r.size(), 0.0);
    for(int l1=0; l1<=lmax; ++l1) {
      for(int l2=max(0,abs(l3-l1)); l2<=min(l3+l1, lmax); ++l2) {
        for(int m1=-l1; m1<=l1; ++m1) {
          if(abs(m3 - m1) > l2)
            continue;
          switch(i3) {
            case 1:
              tmp = R2_1(M, r0, r, f, fp, h, dh, ddh, l3, m3, l1, m1, l2, m3-m1);
              break;
            case 2:
              tmp = R2_2(M, r0, r, f, fp, h, dh, ddh, l3, m3, l1, m1, l2, m3-m1);
              break;
            case 3:
              tmp = R2_3(M, r0, r, f, fp, h, dh, ddh, l3, m3, l1, m1, l2, m3-m1);
              break;
            case 4:
              tmp = R2_4(M, r0, r, f, fp, h, dh, ddh, l3, m3, l1, m1, l2, m3-m1);
              break;
            case 5:
              tmp = R2_5(M, r0, r, f, fp, h, dh, ddh, l3, m3, l1, m1, l2, m3-m1);
              break;
            case 6:
              tmp = R2_6(M, r0, r, f, fp, h, dh, ddh, l3, m3, l1, m1, l2, m3-m1);
              break;
            case 7:
              tmp = R2_7(M, r0, r, f, fp, h, dh, ddh, l3, m3, l1, m1, l2, m3-m1);
              break;
            case 8:
              tmp = R2_8(M, r0, r, f, fp, h, dh, ddh, l3, m3, l1, m1, l2, m3-m1);
              break;
            case 9:
              tmp = R2_9(M, r0, r, f, fp, h, dh, ddh, l3, m3, l1, m1, l2, m3-m1);
              break;
            case 10:
              tmp = R2_10(M, r0, r, f, fp, h, dh, ddh, l3, m3, l1, m1, l2, m3-m1);
              break;
          }
          for(size_t j = 0; j < tmp.size(); ++j)
            src[i3][l3][m3][j] += tmp[j];
        }
      }
    }

#pragma omp critical
    {
      if(++status % status_frequency == 0)
        cout << "Finished " << status << " of " << num_modes << " modes" << endl;
    }
  }
  if(num_modes % status_frequency != 0)
    cout << "Finished " << status << " of " << num_modes << " modes" << endl;
  
  /* Output result */
  cout << "Saving results" << endl;
  H5F src_h5("src.h5", H5F_ACC_TRUNC);
  for(vector<ilm_mode>::iterator mode = modes.begin(); mode < modes.end(); mode++) {
    int i = mode->i;
    int l = mode->l;
    int m = mode->m;

    stringstream dsName;
    dsName << "src i=" << i << " l=" << l << " m=" << m;

    multi_array<double,2> data(boost::extents[N][2]);
    for(int j=0; j<N; ++j) {
      data[j][0] = real(src[i][l][m][j]);
      data[j][1] = imag(src[i][l][m][j]);
    }
    
    src_h5.write_dataset(dsName.str(), data);
  }
  cout << "Done" << endl;
  return 0;
}
