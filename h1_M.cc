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
void h1_M(const double &r0, const vector<double> &r, const vector<double> &f, const vector<double> &fp,
             field_type &h, field_type &dh, field_type &ddh)
{
  cout << "Computing first order mass perturbation fields" << endl;

  const int l = 0;
  const int m = 0;
  int l_max = 0;

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
    // i=1
    h[1][l][m][j] = (2*pow(M,-1)*pow(2*M_PI,0.5)*pow(r[j],-4)*(16*log((r[j]*pow(M,-1))/2.)*(M - r[j])*pow(M,3) - 
       2*(4*pow(M,2) - pow(r[j],2))*(-6*M*r[j] + 2*f[j]*(2*M - r[j])*r[j] + 5*pow(M,2) + 2*pow(r[j],2)) + 
       log(f[j])*(-16*r[j]*pow(M,3) + 24*pow(M,4) + M*f[j]*r[j]*(2*M*r[j] + 4*pow(M,2) + pow(r[j],2)) + 
          3*M*pow(r[j],3) - pow(r[j],4))))/3.;
    dh[1][l][m][j] = (-2*pow(2*M_PI,0.5)*pow(f[j],-1)*pow(r[j],-6)*(32*r[j]*pow(M,3) - 48*pow(M,4) - 6*M*pow(r[j],3) + 
       f[j]*r[j]*(8*f[j]*r[j]*(-2*M + r[j])*(6*M + r[j]) + 16*log(r[j]*pow(M,-1))*(4*M - 3*r[j])*pow(M,2) + 
          4*(31 + 6*log(4))*r[j]*pow(M,2) - 8*(15 + log(256))*pow(M,3) - 30*M*pow(r[j],2) - 4*pow(r[j],3) + 
          log(f[j])*(-52*r[j]*pow(M,2) + 88*pow(M,3) - 2*M*pow(r[j],2) + 
             f[j]*r[j]*(4*M*r[j] + 12*pow(M,2) + pow(r[j],2)) + 3*pow(r[j],3))) + 2*pow(r[j],4)))/3.;
    ddh[1][l][m][j] = (4*pow(2*M_PI,0.5)*pow(f[j],-2)*pow(r[j],-8)*(2*M*(16*r[j]*pow(M,3) - 24*pow(M,4) - 3*M*pow(r[j],3) + pow(r[j],4)) + 
       f[j]*r[j]*(f[j]*r[j]*(32*log((r[j]*pow(M,-1))/2.)*(5*M - 3*r[j])*pow(M,2) + 236*r[j]*pow(M,2) - 
             248*pow(M,3) + log(f[j])*(2*M - r[j])*(-2*M*r[j] + 104*pow(M,2) - 3*pow(r[j],2)) - 54*M*pow(r[j],2) + 
             f[j]*r[j]*(8*(6*M*r[j] - 24*pow(M,2) + pow(r[j],2)) + 
                log(f[j])*(6*M*r[j] + 24*pow(M,2) + pow(r[j],2))) - 4*pow(r[j],3)) + 
          2*(66*r[j]*pow(M,3) - 116*pow(M,4) + pow(M,2)*pow(r[j],2) - 6*M*pow(r[j],3) + pow(r[j],4)))))/3.;

    // i=3
    h[3][l][m][j] = (2*pow(M,-1)*pow(2*M_PI,0.5)*pow(f[j],-1)*pow(r[j],-4)*
     (16*log((r[j]*pow(M,-1))/2.)*(-2*M + r[j])*pow(M,3) + 
       f[j]*r[j]*(64*pow(M,3) + M*log(f[j])*(2*M*r[j] + 4*pow(M,2) + pow(r[j],2)) - 2*pow(r[j],3)) + 
       log(f[j])*(16*r[j]*pow(M,3) - 24*pow(M,4) - 3*M*pow(r[j],3) + pow(r[j],4)) + 
       2*(3*M*r[j] + 10*pow(M,2) + pow(r[j],2))*pow(-2*M + r[j],2)))/3.;
    dh[3][l][m][j] = (-2*pow(2*M_PI,0.5)*pow(f[j],-2)*pow(r[j],-6)*(pow(f[j],2)*pow(r[j],2)*
        (192*pow(M,2) + log(f[j])*(4*M*r[j] + 12*pow(M,2) + pow(r[j],2))) + 
       2*(2*M - r[j])*(-10*r[j]*pow(M,2) + (52 + 8*log(4))*pow(M,3) - 16*log(r[j]*pow(M,-1))*pow(M,3) - 
          3*M*pow(r[j],2) - pow(r[j],3)) + f[j]*r[j]*
        (-4*(47 + 6*log(4))*r[j]*pow(M,2) + 16*log(r[j]*pow(M,-1))*(-8*M + 3*r[j])*pow(M,2) + 
          8*(43 + 8*log(4))*pow(M,3) + 6*M*pow(r[j],2) - 2*pow(r[j],3) - 
          3*log(f[j])*(-16*r[j]*pow(M,2) + 32*pow(M,3) + pow(r[j],3))) + 
       2*log(f[j])*(16*r[j]*pow(M,3) - 24*pow(M,4) - 3*M*pow(r[j],3) + pow(r[j],4))))/3.;
    ddh[3][l][m][j] = (4*pow(2*M_PI,0.5)*pow(f[j],-3)*pow(r[j],-8)*(2*M*(2*M - r[j])*
        (-22*r[j]*pow(M,2) + 4*(29 + log(256))*pow(M,3) - 32*log(r[j]*pow(M,-1))*pow(M,3) - 7*M*pow(r[j],2) - 
          pow(r[j],3)) + pow(f[j],3)*(384*pow(M,2) + log(f[j])*(6*M*r[j] + 24*pow(M,2) + pow(r[j],2)))*
        pow(r[j],3) + pow(f[j],2)*pow(r[j],2)*(-4*(101 + 12*log(4))*r[j]*pow(M,2) + 
          32*log(r[j]*pow(M,-1))*(-10*M + 3*r[j])*pow(M,2) + 16*(57 + 10*log(4))*pow(M,3) + 8*M*pow(r[j],2) - 
          2*pow(r[j],3) - 3*log(f[j])*(-32*r[j]*pow(M,2) + 80*pow(M,3) + pow(r[j],3))) + 
       4*M*log(f[j])*(16*r[j]*pow(M,3) - 24*pow(M,4) - 3*M*pow(r[j],3) + pow(r[j],4)) + 
       2*f[j]*r[j]*(-2*(153 + 16*log(4))*r[j]*pow(M,3) + 32*log(r[j]*pow(M,-1))*(-5*M + 2*r[j])*pow(M,3) + 
          (548 + 80*log(4))*pow(M,4) + 11*pow(M,2)*pow(r[j],2) + 2*M*pow(r[j],3) + pow(r[j],4) + 
          log(f[j])*(64*r[j]*pow(M,3) - 120*pow(M,4) - 6*M*pow(r[j],3) + pow(r[j],4)))))/3.;

    // i=6
    h[6][l][m][j] = (-4*pow(M,-1)*pow(2*M_PI,0.5)*pow(r[j],-3)*(-4*r[j]*pow(M,2) + 4*(5 + log(4))*pow(M,3) - 
       8*log(r[j]*pow(M,-1))*pow(M,3) - M*pow(r[j],2) + 4*f[j]*r[j]*(2*M*r[j] + 4*pow(M,2) + pow(r[j],2)) - 
       4*pow(r[j],3) + log(f[j])*(-8*pow(M,3) + pow(r[j],3))))/3.;
    dh[6][l][m][j] = (4*pow(2*M_PI,0.5)*pow(f[j],-1)*pow(r[j],-5)*(16*pow(M,3) + 
       f[j]*r[j]*(8*f[j]*r[j]*(4*M + r[j]) + 3*
           (-8*M*r[j] + 4*(3 + log(4))*pow(M,2) - 8*(log(f[j]) + log(r[j]*pow(M,-1)))*pow(M,2) - 3*pow(r[j],2))) - 
       2*pow(r[j],3)))/3.;
    ddh[6][l][m][j] = (-8*pow(2*M_PI,0.5)*pow(f[j],-2)*pow(r[j],-7)*(16*pow(M,4) + 
       f[j]*r[j]*(64*pow(M,3) + f[j]*r[j]*(-44*M*r[j] + 8*f[j]*r[j]*(6*M + r[j]) + 4*(13 + 6*log(4))*pow(M,2) - 
             48*(log(f[j]) + log(r[j]*pow(M,-1)))*pow(M,2) - 9*pow(r[j],2)) - 2*pow(r[j],3)) - 2*M*pow(r[j],3)))/3.;
  }
}
