/*
 * Copyright 2014 Barry Wardell
 *
 * This file is distributed under the University of Illinois/NCSA
 * Open Source License. See LICENSE file for details.
 *
 */

#include <boost/multi_array.hpp>
#include <complex>
#include <math.h>
#include <iostream>

#include "Coupling.h"
#include "hh.h"

using namespace std;

vector<complex<double>> hh_9(const double M, const double r0,
                     const vector<double> &r, const vector<double> &f, const vector<double> &fp,
                     const field_type &h1, const field_type &h2,
                     int l3, int m3, int l1, int m1, int l2, int m2) {
  vector<complex<double>> hh(r.size());
  const complex<double> I(0.0, 1.0);

  const double c10 = Coupling(l3,m3,-1,l1,m1,-1,l2,m2,0);
  const double c11 = Coupling(l3,m3,1,l1,m1,1,l2,m2,0);
  const double c12 = Coupling(l3,m3,-1,l1,m1,1,l2,m2,-2);
  const double c13 = Coupling(l3,m3,1,l1,m1,-1,l2,m2,2);

  if(l1>=2&&l2>=2) {
     for(vector<complex<double>>::size_type j=0; j!=hh.size(); ++j) {
       hh[j] = (0.35355339059327373*c10*(I*(h1[4][l1][m1][j] + I*h1[8][l1][m1][j])*h2[2][l2][m2][j] + (-1.*I*h1[5][l1][m1][j] + h1[9][l1][m1][j])*(h2[1][l2][m2][j] + f[j]*(-h2[3][l2][m2][j] + h2[6][l2][m2][j]))))/f[j] + (0.35355339059327373*c11*((-1.*I*h1[4][l1][m1][j] - h1[8][l1][m1][j])*h2[2][l2][m2][j] + (I*h1[5][l1][m1][j] + h1[9][l1][m1][j])*(h2[1][l2][m2][j] + f[j]*(-h2[3][l2][m2][j] + h2[6][l2][m2][j]))))/f[j] + 0.35355339059327373*c13*(-1.*I*h1[5][l1][m1][j] + h1[9][l1][m1][j])*(h2[7][l2][m2][j] - I*h2[10][l2][m2][j]) + 0.35355339059327373*c12*(I*h1[5][l1][m1][j] + h1[9][l1][m1][j])*(h2[7][l2][m2][j] + I*h2[10][l2][m2][j]);
     }
  } else if(l1==0&&l2==0) {
    fill(hh.begin(), hh.end(), 0.0);
  } else if(l1==1&&l2==0) {
    for(vector<complex<double>>::size_type j=0; j!=hh.size(); ++j) {
      hh[j] = c10*((0.35355339059327373*h1[9][1][m1][j]*h2[1][0][0][j])/f[j] + (0.35355339059327373*I*h1[4][1][m1][j]*h2[2][0][0][j])/f[j] - (0.35355339059327373*h1[8][1][m1][j]*h2[2][0][0][j])/f[j] + h1[5][1][m1][j]*((-0.35355339059327373*I*h2[1][0][0][j])/f[j] + 0.35355339059327373*I*h2[3][0][0][j] - 0.35355339059327373*I*h2[6][0][0][j]) + h1[9][1][m1][j]*(-0.35355339059327373*h2[3][0][0][j] + 0.35355339059327373*h2[6][0][0][j])) + c11*((0.35355339059327373*h1[9][1][m1][j]*h2[1][0][0][j])/f[j] - (0.35355339059327373*I*h1[4][1][m1][j]*h2[2][0][0][j])/f[j] - (0.35355339059327373*h1[8][1][m1][j]*h2[2][0][0][j])/f[j] + h1[5][1][m1][j]*((0.35355339059327373*I*h2[1][0][0][j])/f[j] - 0.35355339059327373*I*h2[3][0][0][j] + 0.35355339059327373*I*h2[6][0][0][j]) + h1[9][1][m1][j]*(-0.35355339059327373*h2[3][0][0][j] + 0.35355339059327373*h2[6][0][0][j]));
    }
  } else if(l1==0&&l2==1) {
    fill(hh.begin(), hh.end(), 0.0);
  } else if(l1==1&&l2==1) {
    for(vector<complex<double>>::size_type j=0; j!=hh.size(); ++j) {
      hh[j] = c10*((0.35355339059327373*h1[9][1][m1][j]*h2[1][1][m2][j])/f[j] + (0.35355339059327373*I*h1[4][1][m1][j]*h2[2][1][m2][j])/f[j] - (0.35355339059327373*h1[8][1][m1][j]*h2[2][1][m2][j])/f[j] + h1[5][1][m1][j]*((-0.35355339059327373*I*h2[1][1][m2][j])/f[j] + 0.35355339059327373*I*h2[3][1][m2][j] - 0.35355339059327373*I*h2[6][1][m2][j]) + h1[9][1][m1][j]*(-0.35355339059327373*h2[3][1][m2][j] + 0.35355339059327373*h2[6][1][m2][j])) + c11*((0.35355339059327373*h1[9][1][m1][j]*h2[1][1][m2][j])/f[j] - (0.35355339059327373*I*h1[4][1][m1][j]*h2[2][1][m2][j])/f[j] - (0.35355339059327373*h1[8][1][m1][j]*h2[2][1][m2][j])/f[j] + h1[5][1][m1][j]*((0.35355339059327373*I*h2[1][1][m2][j])/f[j] - 0.35355339059327373*I*h2[3][1][m2][j] + 0.35355339059327373*I*h2[6][1][m2][j]) + h1[9][1][m1][j]*(-0.35355339059327373*h2[3][1][m2][j] + 0.35355339059327373*h2[6][1][m2][j]));
    }
  } else if(l1==0) {
    fill(hh.begin(), hh.end(), 0.0);
  } else if(l2==0) {
    for(vector<complex<double>>::size_type j=0; j!=hh.size(); ++j) {
      hh[j] = c10*((0.35355339059327373*h1[9][l1][m1][j]*h2[1][0][0][j])/f[j] + (0.35355339059327373*I*h1[4][l1][m1][j]*h2[2][0][0][j])/f[j] - (0.35355339059327373*h1[8][l1][m1][j]*h2[2][0][0][j])/f[j] + h1[5][l1][m1][j]*((-0.35355339059327373*I*h2[1][0][0][j])/f[j] + 0.35355339059327373*I*h2[3][0][0][j] - 0.35355339059327373*I*h2[6][0][0][j]) + h1[9][l1][m1][j]*(-0.35355339059327373*h2[3][0][0][j] + 0.35355339059327373*h2[6][0][0][j])) + c11*((0.35355339059327373*h1[9][l1][m1][j]*h2[1][0][0][j])/f[j] - (0.35355339059327373*I*h1[4][l1][m1][j]*h2[2][0][0][j])/f[j] - (0.35355339059327373*h1[8][l1][m1][j]*h2[2][0][0][j])/f[j] + h1[5][l1][m1][j]*((0.35355339059327373*I*h2[1][0][0][j])/f[j] - 0.35355339059327373*I*h2[3][0][0][j] + 0.35355339059327373*I*h2[6][0][0][j]) + h1[9][l1][m1][j]*(-0.35355339059327373*h2[3][0][0][j] + 0.35355339059327373*h2[6][0][0][j]));
    }
  } else if(l1==1) {
    for(vector<complex<double>>::size_type j=0; j!=hh.size(); ++j) {
      hh[j] = c10*((0.35355339059327373*h1[9][1][m1][j]*h2[1][l2][m2][j])/f[j] + (0.35355339059327373*I*h1[4][1][m1][j]*h2[2][l2][m2][j])/f[j] - (0.35355339059327373*h1[8][1][m1][j]*h2[2][l2][m2][j])/f[j] + h1[5][1][m1][j]*((-0.35355339059327373*I*h2[1][l2][m2][j])/f[j] + 0.35355339059327373*I*h2[3][l2][m2][j] - 0.35355339059327373*I*h2[6][l2][m2][j]) + h1[9][1][m1][j]*(-0.35355339059327373*h2[3][l2][m2][j] + 0.35355339059327373*h2[6][l2][m2][j])) + c11*((0.35355339059327373*h1[9][1][m1][j]*h2[1][l2][m2][j])/f[j] - (0.35355339059327373*I*h1[4][1][m1][j]*h2[2][l2][m2][j])/f[j] - (0.35355339059327373*h1[8][1][m1][j]*h2[2][l2][m2][j])/f[j] + h1[5][1][m1][j]*((0.35355339059327373*I*h2[1][l2][m2][j])/f[j] - 0.35355339059327373*I*h2[3][l2][m2][j] + 0.35355339059327373*I*h2[6][l2][m2][j]) + h1[9][1][m1][j]*(-0.35355339059327373*h2[3][l2][m2][j] + 0.35355339059327373*h2[6][l2][m2][j])) + c13*(h1[5][1][m1][j]*(-0.35355339059327373*I*h2[7][l2][m2][j] - 0.35355339059327373*h2[10][l2][m2][j]) + h1[9][1][m1][j]*(0.35355339059327373*h2[7][l2][m2][j] - 0.35355339059327373*I*h2[10][l2][m2][j])) + c12*(h1[5][1][m1][j]*(0.35355339059327373*I*h2[7][l2][m2][j] - 0.35355339059327373*h2[10][l2][m2][j]) + h1[9][1][m1][j]*(0.35355339059327373*h2[7][l2][m2][j] + 0.35355339059327373*I*h2[10][l2][m2][j]));
    }
  } else if(l2==1) {
    for(vector<complex<double>>::size_type j=0; j!=hh.size(); ++j) {
      hh[j] = c10*((0.35355339059327373*h1[9][l1][m1][j]*h2[1][1][m2][j])/f[j] + (0.35355339059327373*I*h1[4][l1][m1][j]*h2[2][1][m2][j])/f[j] - (0.35355339059327373*h1[8][l1][m1][j]*h2[2][1][m2][j])/f[j] + h1[5][l1][m1][j]*((-0.35355339059327373*I*h2[1][1][m2][j])/f[j] + 0.35355339059327373*I*h2[3][1][m2][j] - 0.35355339059327373*I*h2[6][1][m2][j]) + h1[9][l1][m1][j]*(-0.35355339059327373*h2[3][1][m2][j] + 0.35355339059327373*h2[6][1][m2][j])) + c11*((0.35355339059327373*h1[9][l1][m1][j]*h2[1][1][m2][j])/f[j] - (0.35355339059327373*I*h1[4][l1][m1][j]*h2[2][1][m2][j])/f[j] - (0.35355339059327373*h1[8][l1][m1][j]*h2[2][1][m2][j])/f[j] + h1[5][l1][m1][j]*((0.35355339059327373*I*h2[1][1][m2][j])/f[j] - 0.35355339059327373*I*h2[3][1][m2][j] + 0.35355339059327373*I*h2[6][1][m2][j]) + h1[9][l1][m1][j]*(-0.35355339059327373*h2[3][1][m2][j] + 0.35355339059327373*h2[6][1][m2][j]));
    }
  } else {
    /* We should never reach here */
    cerr << "Invalid value for l1 (" << l1 << ") or l2 (" << l2 << ")." << endl;
    exit(EXIT_FAILURE);
  }

  return hh;
}
