/*
 * Copyright 2014 Barry Wardell
 *
 * This file is distributed under the University of Illinois/NCSA
 * Open Source License. See LICENSE file for details.
 *
 */

#include "h1_M.h"
#include "h5wrapper.h"
#include <complex>
#include <iostream>

using namespace std;
using boost::multi_array;
typedef boost::multi_array_types::extent_range range;
typedef boost::multi_array_types::index_range irange;
typedef boost::multi_array<complex<double>,4> field_type;

const double M = 1.0;

/* Compute first order mass perturbation */
void h1_S(const double &r0, const vector<double> &r, const vector<double> &f, const vector<double> &fp,
             field_type &h, field_type &dh, field_type &ddh)
{
  cout << "Computing first order angular momentum perturbation fields" << endl;

  const int l = 1;
  const int m = 0;
  int l_max = 1;

  size_t N = r.size();  /* Number of grid points */

  cout << "Grid size: [" << r.front() << ", " << r.back() << "] (" << N << " points)" << endl;
  cout << "Worldline: r_0 = " << r0  << endl;

  /* Evaluate the perturbation on the grid */
  h.resize(boost::extents[range(1,11)][l_max+1][range(-l_max,l_max+1)][N]);
  dh.resize(boost::extents[range(1,11)][l_max+1][range(-l_max,l_max+1)][N]);
  ddh.resize(boost::extents[range(1,11)][l_max+1][range(-l_max,l_max+1)][N]);

  fill(h.data(), h.data() + h.num_elements(), 0.0);
  fill(dh.data(), dh.data() + dh.num_elements(), 0.0);
  fill(ddh.data(), ddh.data() + ddh.num_elements(), 0.0);

  for(size_t j=0; j<N; ++j) {
    // i=8
    h[8][l][m][j] = -8*pow(M_PI/3.,0.5)*pow(r[j],-2);
    dh[8][l][m][j] = 16*pow(M_PI/3.,0.5)*pow(r[j],-3);
    ddh[8][l][m][j] = -48*pow(M_PI/3.,0.5)*pow(r[j],-4);

    // i=9
    h[9][l][m][j] = -16*pow(M_PI/3.,0.5)*pow(r[j],-3);
    dh[9][l][m][j] = 48*pow(M_PI/3.,0.5)*pow(r[j],-4);
    ddh[9][l][m][j] = -192*pow(M_PI/3.,0.5)*pow(r[j],-5);
  }
}
