/*
 * Copyright 2014 Barry Wardell
 *
 * This file is distributed under the University of Illinois/NCSA
 * Open Source License. See LICENSE file for details.
 *
 */

#include <vector>
#include <complex>
#include <boost/multi_array.hpp>
#include "h1.h"
#include "R2.h"
#include "utils.h"

using namespace std;
using boost::multi_array;
typedef boost::multi_array_types::extent_range range;

int main(int argc, char* argv[])
{
  /* Read in first-order fields */
  vector<double> r, f, fp;
  multi_array<complex<double>, 4> h, dh, ddh;
  read_h1(r, f, fp, h, dh, ddh);
  const int lmax = h[1].size()-1;
  const int N = r.size();

  /* Compute the source */
  multi_array<complex<double>,4> src(boost::extents[range(1,11)][lmax+1][range(-lmax,lmax+1)][N]);
  fill(src.data(), src.data() + src.num_elements(), 0.0);

  /* Loop over i3, l3, m3 */
  for(int l3=0; l3<=3; ++l3)
  for(int m3=-l3; m3<=l3; ++m3)
  for(int i3=1; i3<=10; ++i3)
  {
    if(((i3<=7) && isOdd(l3+m3)) || ((i3>7) && isEven(l3+m3)))
      continue;

    /* Sum over l1, l2, m1, m2 = m3-m1 */
    vector<complex<double>> tmp(r.size(), 0.0);
    for(int l1=0; l1<=lmax; ++l1)
      for(int l2=max(0,abs(l3-l1)); l2<=min(l3+l1, lmax); ++l2)
        for(int m1=-l1; m1<=l1; ++m1) {
          if(abs(m3 - m1) > l2)
            continue;
          switch(i3) {
            case 1:
              tmp = R2_1(r, f, fp, h, dh, l3, m3, l1, m1, l2, m3-m1);
              break;
            case 2:
              tmp = R2_2(r, f, fp, h, dh, l3, m3, l1, m1, l2, m3-m1);
              break;
            case 3:
              tmp = R2_3(r, f, fp, h, dh, l3, m3, l1, m1, l2, m3-m1);
              break;
            case 4:
              tmp = R2_4(r, f, fp, h, dh, l3, m3, l1, m1, l2, m3-m1);
              break;
            case 5:
              tmp = R2_5(r, f, fp, h, dh, l3, m3, l1, m1, l2, m3-m1);
              break;
            case 6:
              tmp = R2_6(r, f, fp, h, dh, l3, m3, l1, m1, l2, m3-m1);
              break;
            case 7:
              tmp = R2_7(r, f, fp, h, dh, l3, m3, l1, m1, l2, m3-m1);
              break;
            case 8:
              tmp = R2_8(r, f, fp, h, dh, l3, m3, l1, m1, l2, m3-m1);
              break;
            case 9:
              tmp = R2_9(r, f, fp, h, dh, l3, m3, l1, m1, l2, m3-m1);
              break;
            case 10:
              tmp = R2_10(r, f, fp, h, dh, l3, m3, l1, m1, l2, m3-m1);
              break;
          }
          for(int j=0; j< tmp.size(); ++j)
            src[i3][l3][m3][j] += tmp[j];
        }
  }

  return 0;
}
