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
  cout << "Version: " << __GIT_VERSION << endl;

  string dirA, dirB;
  switch(argc) {
    case 2:
      dirA = dirB = argv[1];
      break;
    case 5:
      dirA = argv[1];
      dirB = argv[2];
      break;
    default:
      cout << "Usage: " << argv[0] << " <dirA> [<dirB>] [<delta lmax>] [<lmax>]" << endl;
      exit(1);
  }

#ifdef _OPENMP
  int num_threads = omp_get_max_threads();
  cout << "Running with " << num_threads << " OpenMP threads" << endl;
#endif

  /* Read in first-order fields */
  vector<double> r, f, fp;
  double r0;
  multi_array<complex<double>, 4> hA, dhA, ddhA, hB, dhB, ddhB;
  read_h1(dirA, r0, r, f, fp, hA, dhA, ddhA, std::numeric_limits<int>::max());
  read_h1(dirB, r0, r, f, fp, hB, dhB, ddhB, std::numeric_limits<int>::max());
  int h1lmax = min(hA[1].size()-1, hB[1].size()-1);
  const int N = r.size();

  int lmax, dlmax;
  switch(argc) {
    case 2:
      dlmax = lmax = h1lmax;
      break;
    case 5:
      dlmax = atoi(argv[3]);

      if(dlmax<0) {
        h1lmax += dlmax;
        dlmax = h1lmax;
      } else {
        dlmax = min(dlmax, h1lmax);
      }
      lmax = min(atoi(argv[4]), h1lmax);
      break;
    default:
      cout << "Error, inconsistent command line arguments." << endl;
      exit(1);
  }

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
        if(((i3==4)||(i3==5)||(i3==8)||(i3==9)) && l3<1)
          continue;
        if(((i3==7)||(i3==10)) && l3<2)
          continue;
        modes.push_back(ilm_mode({i3, l3, m3}));
      }
    }
  }

  cout << "Computing source: " << modes.size() << " (i,l,m) modes, l_max = " << lmax
       << ", delta l_max = " << dlmax << endl;

  int status = 0;
  const int status_frequency = 4*num_threads;
  const int num_modes = modes.size();
#pragma omp parallel for
  for(vector<ilm_mode>::iterator mode = modes.begin(); mode < modes.end(); mode++) {
    int i3 = mode->i;
    int l3 = mode->l;
    int m3 = mode->m;

    /* Sum over l1, l2, m1, m2 = m3-m1.
     * The sum is done in the range l3-dlmax <= (l1,l2) <= l3+dlmax, within
     * physical (l1,l2 >=0) and practical (l1,l2 <= lmax for the first order
     * fields) constraints.*/
    vector<complex<double>> tmp(r.size(), 0.0);
    for(int l1=max(0,l3-dlmax); l1<=min(h1lmax,l3+dlmax); ++l1) {
      for(int l2=max(abs(l3-l1),l3-dlmax); l2<=min(min(l3+l1, h1lmax),l3+dlmax); ++l2) {
        for(int m1=-l1; m1<=l1; ++m1) {
          if(abs(m3 - m1) > l2)
            continue;
          switch(i3) {
            case 1:
              tmp = R2_1(M, r0, r, f, fp, hA, dhA, ddhA, hB, dhB, ddhB, l3, m3, l1, m1, l2, m3-m1);
              break;
            case 2:
              tmp = R2_2(M, r0, r, f, fp, hA, dhA, ddhA, hB, dhB, ddhB, l3, m3, l1, m1, l2, m3-m1);
              break;
            case 3:
              tmp = R2_3(M, r0, r, f, fp, hA, dhA, ddhA, hB, dhB, ddhB, l3, m3, l1, m1, l2, m3-m1);
              break;
            case 4:
              tmp = R2_4(M, r0, r, f, fp, hA, dhA, ddhA, hB, dhB, ddhB, l3, m3, l1, m1, l2, m3-m1);
              break;
            case 5:
              tmp = R2_5(M, r0, r, f, fp, hA, dhA, ddhA, hB, dhB, ddhB, l3, m3, l1, m1, l2, m3-m1);
              break;
            case 6:
              tmp = R2_6(M, r0, r, f, fp, hA, dhA, ddhA, hB, dhB, ddhB, l3, m3, l1, m1, l2, m3-m1);
              break;
            case 7:
              tmp = R2_7(M, r0, r, f, fp, hA, dhA, ddhA, hB, dhB, ddhB, l3, m3, l1, m1, l2, m3-m1);
              break;
            case 8:
              tmp = R2_8(M, r0, r, f, fp, hA, dhA, ddhA, hB, dhB, ddhB, l3, m3, l1, m1, l2, m3-m1);
              break;
            case 9:
              tmp = R2_9(M, r0, r, f, fp, hA, dhA, ddhA, hB, dhB, ddhB, l3, m3, l1, m1, l2, m3-m1);
              break;
            case 10:
              tmp = R2_10(M, r0, r, f, fp, hA, dhA, ddhA, hB, dhB, ddhB, l3, m3, l1, m1, l2, m3-m1);
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

    multi_array<double,2> data(boost::extents[N][2]);

    /* Second order source */
    stringstream srcds;
    srcds << "src i=" << i << " l=" << l << " m=" << m;
    for(int j=0; j<N; ++j) {
      data[j][0] = real(src[i][l][m][j]);
      data[j][1] = imag(src[i][l][m][j]);
    }
    src_h5.write_dataset(srcds.str(), data);
  }
  /* Grid coordinates */
  stringstream rds;
  rds << "r";
  multi_array<double,1> rdata(boost::extents[N]);
  for(int j=0; j<N; ++j) {
    rdata[j] = r[j];
  }
  src_h5.write_dataset(rds.str(), rdata);
  cout << "Done" << endl;
  return 0;
}
