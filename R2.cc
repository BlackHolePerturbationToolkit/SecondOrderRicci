/*
 * Copyright 2014 Barry Wardell
 *
 * This file is distributed under the University of Illinois/NCSA
 * Open Source License. See LICENSE file for details.
 *
 */

#include "Coupling.h"
#include <boost/multi_array.hpp>
#include <complex>
#include <math.h>

using std::complex;
using std::vector;
typedef boost::multi_array<complex<double>,4> field_type;

vector<complex<double>> R2_1(const vector<double> &r, const vector<double> &f, const vector<double> &fp,
                     const field_type &h, const field_type &dh,
                     int l3, int m3, int l1, int m1, int l2, int m2) {
  vector<complex<double>> ricci(r.size());
  complex<double> I(0.0, 1.0);
  double Omega = sqrt(1./1000.);
  double l1d = l1, m1d = m1, l2d = l2, m2d = m2;

  for(vector<complex<double>>::size_type j=0; j!=ricci.size(); ++j) {
    ricci[j] = (Coupling(l3,m3,0,l1,m1,0,l2,m2,0)*(-(pow(f[j],3)*pow(r[j],2)*(dh[1][l2][m2][j]*dh[3][l1][m1][j] +
      dh[1][l1][m1][j]*dh[3][l2][m2][j])) - 2.0*(1.0*I*m1d*Omega*h[1][l2][m2][j]*h[2][l1][m1][j] -
      1.0*I*m2d*Omega*h[1][l2][m2][j]*h[2][l1][m1][j] - 1.0*I*m1d*Omega*h[1][l1][m1][j]*h[2][l2][m2][j] +
      1.0*I*m2d*Omega*h[1][l1][m1][j]*h[2][l2][m2][j] - m1d*m2d*pow(Omega,2)*pow(r[j],2)*h[2][l1][m1][j]*h[2][l2][m2][j]) +
      f[j]*(-2.0*dh[1][l2][m2][j]*h[1][l1][m1][j] - 1.0*I*m1d*Omega*pow(r[j],2)*dh[2][l2][m2][j]*h[1][l1][m1][j] -
      2.0*dh[1][l1][m1][j]*h[1][l2][m2][j] - 1.0*I*m2d*Omega*pow(r[j],2)*dh[2][l1][m1][j]*h[1][l2][m2][j] -
      1.0*I*m1d*Omega*pow(r[j],2)*dh[1][l2][m2][j]*h[2][l1][m1][j] + 2.0*dh[2][l2][m2][j]*h[2][l1][m1][j] -
      1.0*I*m2d*Omega*pow(r[j],2)*dh[1][l1][m1][j]*h[2][l2][m2][j] + 2.0*dh[2][l1][m1][j]*h[2][l2][m2][j] -
      m1d*m2d*pow(Omega,2)*pow(r[j],2)*h[1][l2][m2][j]*h[3][l1][m1][j] -
      m1d*m2d*pow(Omega,2)*pow(r[j],2)*h[1][l1][m1][j]*h[3][l2][m2][j]) + pow(f[j],4)*r[j]*(dh[6][l2][m2][j]*(r[j]*dh[6][l1][m1][j] +
      2.0*h[6][l1][m1][j]) - 2.0*dh[6][l1][m1][j]*h[6][l2][m2][j]) + pow(f[j],2)*(pow(r[j],2)*(-2.0*dh[2][l1][m1][j]*dh[2][l2][m2][j]
      - 1.0*I*m1d*Omega*dh[3][l2][m2][j]*h[2][l1][m1][j] - 1.0*I*m2d*Omega*dh[3][l1][m1][j]*h[2][l2][m2][j] +
      1.0*I*m1d*Omega*dh[2][l2][m2][j]*h[3][l1][m1][j] + 1.0*I*m2d*Omega*dh[2][l1][m1][j]*h[3][l2][m2][j] -
      m1d*m2d*pow(Omega,2)*h[6][l1][m1][j]*h[6][l2][m2][j]) + sqrt(l1d*(1 + l1d)*l2d*(1 + l2d))*h[8][l1][m1][j]*h[8][l2][m2][j] +
      sqrt(l1d*(1 + l1d)*l2d*(1 + l2d))*h[9][l1][m1][j]*h[9][l2][m2][j])))/(2.*sqrt(2)*pow(f[j],2)*pow(r[j],2)); +
      (Coupling(l3,m3,0,l1,m1,-1,l2,m2,1)*(r[j]*(1.0*I*sqrt(l1d*(1 + l1d))*m2d*Omega*r[j]*h[1][l1][m1][j]*h[4][l2][m2][j] -
      m1d*m2d*pow(Omega,2)*pow(r[j],2)*h[4][l1][m1][j]*h[4][l2][m2][j] + 1.0*I*m2d*Omega*h[4][l2][m2][j]*h[5][l1][m1][j] +
      1.0*I*sqrt(l1d*(1 + l1d))*m2d*Omega*r[j]*h[2][l1][m1][j]*h[5][l2][m2][j] -
      1.0*I*m1d*m2d*pow(Omega,2)*pow(r[j],2)*h[4][l2][m2][j]*h[8][l1][m1][j] - h[5][l2][m2][j]*(-1.0*I*m1d*Omega*h[4][l1][m1][j] +
      m1d*Omega*h[8][l1][m1][j]) + sqrt(l1d*(1 + l1d))*m2d*Omega*r[j]*h[1][l1][m1][j]*h[8][l2][m2][j] +
      1.0*I*m1d*m2d*pow(Omega,2)*pow(r[j],2)*h[4][l1][m1][j]*h[8][l2][m2][j] + m2d*Omega*h[5][l1][m1][j]*h[8][l2][m2][j] -
      m1d*m2d*pow(Omega,2)*pow(r[j],2)*h[8][l1][m1][j]*h[8][l2][m2][j] - m2d*Omega*h[4][l2][m2][j]*h[9][l1][m1][j] +
      1.0*I*m2d*Omega*h[8][l2][m2][j]*h[9][l1][m1][j] + sqrt(l1d*(1 + l1d))*m2d*Omega*r[j]*h[2][l1][m1][j]*h[9][l2][m2][j] +
      1.0*I*(-1.0*I*m1d*Omega*h[4][l1][m1][j] + m1d*Omega*h[8][l1][m1][j])*h[9][l2][m2][j]) -
      pow(f[j],2)*r[j]*(pow(r[j],2)*dh[5][l1][m1][j]*dh[5][l2][m2][j] + 1.0*I*pow(r[j],2)*dh[5][l2][m2][j]*dh[9][l1][m1][j] -
      1.0*I*pow(r[j],2)*dh[5][l1][m1][j]*dh[9][l2][m2][j] + pow(r[j],2)*dh[9][l1][m1][j]*dh[9][l2][m2][j] +
      r[j]*dh[4][l1][m1][j]*h[4][l2][m2][j] + 1.0*I*r[j]*dh[8][l1][m1][j]*h[4][l2][m2][j] + r[j]*dh[5][l1][m1][j]*h[5][l2][m2][j] +
      1.0*I*r[j]*dh[9][l1][m1][j]*h[5][l2][m2][j] - 1.0*I*r[j]*dh[4][l1][m1][j]*h[8][l2][m2][j] +
      r[j]*dh[8][l1][m1][j]*h[8][l2][m2][j] + 1.0*I*r[j]*dh[5][l2][m2][j]*h[9][l1][m1][j] + r[j]*dh[9][l2][m2][j]*h[9][l1][m1][j] +
      2.0*I*h[5][l2][m2][j]*h[9][l1][m1][j] + sqrt(l1d*(1 + l1d))*h[3][l1][m1][j]*(r[j]*dh[5][l2][m2][j] - 1.0*I*r[j]*dh[9][l2][m2][j]
      + 2.0*h[5][l2][m2][j] - 2.0*I*h[9][l2][m2][j]) + h[5][l1][m1][j]*(r[j]*dh[5][l2][m2][j] - 1.0*I*r[j]*dh[9][l2][m2][j] +
      2.0*h[5][l2][m2][j] - 2.0*I*h[9][l2][m2][j]) - 1.0*I*r[j]*dh[5][l1][m1][j]*h[9][l2][m2][j] +
      r[j]*dh[9][l1][m1][j]*h[9][l2][m2][j] + 2.0*h[9][l1][m1][j]*h[9][l2][m2][j]) + f[j]*(sqrt(l1d*(1 +
      l1d))*pow(r[j],2)*dh[5][l2][m2][j]*h[1][l1][m1][j] - 1.0*I*sqrt(l1d*(1 + l1d))*pow(r[j],2)*dh[9][l2][m2][j]*h[1][l1][m1][j] +
      sqrt(l1d*(1 + l1d))*pow(r[j],2)*dh[4][l2][m2][j]*h[2][l1][m1][j] - 1.0*I*sqrt(l1d*(1 +
      l1d))*pow(r[j],2)*dh[8][l2][m2][j]*h[2][l1][m1][j] + sqrt(l1d*(1 + l1d)*l2d*(1 + l2d))*r[j]*h[1][l1][m1][j]*h[3][l2][m2][j] -
      1.0*I*m1d*Omega*pow(r[j],2)*h[4][l2][m2][j]*h[5][l1][m1][j] + r[j]*dh[5][l1][m1][j]*h[5][l2][m2][j] +
      1.0*I*r[j]*dh[9][l1][m1][j]*h[5][l2][m2][j] + 2.0*sqrt(l1d*(1 + l1d))*r[j]*h[1][l1][m1][j]*h[5][l2][m2][j] -
      m1d*Omega*pow(r[j],2)*h[5][l1][m1][j]*h[8][l2][m2][j] + 1.0*I*r[j]*dh[5][l2][m2][j]*h[9][l1][m1][j] +
      r[j]*dh[9][l2][m2][j]*h[9][l1][m1][j] + m1d*Omega*pow(r[j],2)*h[4][l2][m2][j]*h[9][l1][m1][j] +
      2.0*I*h[5][l2][m2][j]*h[9][l1][m1][j] - 1.0*I*m1d*Omega*pow(r[j],2)*h[8][l2][m2][j]*h[9][l1][m1][j] +
      h[5][l1][m1][j]*(r[j]*(dh[5][l2][m2][j] - 1.0*I*dh[9][l2][m2][j]) + 2.0*h[5][l2][m2][j] - 2.0*I*h[9][l2][m2][j]) -
      1.0*I*r[j]*dh[5][l1][m1][j]*h[9][l2][m2][j] + r[j]*dh[9][l1][m1][j]*h[9][l2][m2][j] - 2.0*I*sqrt(l1d*(1 +
      l1d))*r[j]*h[1][l1][m1][j]*h[9][l2][m2][j] + 2.0*h[9][l1][m1][j]*h[9][l2][m2][j] + h[3][l1][m1][j]*(2.0*sqrt(l1d*(1 +
      l1d))*h[5][l2][m2][j] + r[j]*(sqrt(l1d*(1 + l1d)*l2d*(1 + l2d))*h[1][l2][m2][j] - sqrt(l1d*(1 +
      l1d))*r[j]*(-1.0*I*m2d*Omega*h[4][l2][m2][j] - m2d*Omega*h[8][l2][m2][j])) - 2.0*I*sqrt(l1d*(1 +
      l1d))*h[9][l2][m2][j]))))/(2.*sqrt(2)*f[j]*pow(r[j],3)) + (Coupling(l3,m3,0,l1,m1,1,l2,m2,-1)*(r[j]*(1.0*I*sqrt(l1d*(1 +
      l1d))*m2d*Omega*r[j]*h[1][l1][m1][j]*h[4][l2][m2][j] - m1d*m2d*pow(Omega,2)*pow(r[j],2)*h[4][l1][m1][j]*h[4][l2][m2][j] +
      1.0*I*m2d*Omega*h[4][l2][m2][j]*h[5][l1][m1][j] + 1.0*I*sqrt(l1d*(1 + l1d))*m2d*Omega*r[j]*h[2][l1][m1][j]*h[5][l2][m2][j] +
      1.0*I*m1d*m2d*pow(Omega,2)*pow(r[j],2)*h[4][l2][m2][j]*h[8][l1][m1][j] - h[5][l2][m2][j]*(-1.0*I*m1d*Omega*h[4][l1][m1][j] -
      m1d*Omega*h[8][l1][m1][j]) - sqrt(l1d*(1 + l1d))*m2d*Omega*r[j]*h[1][l1][m1][j]*h[8][l2][m2][j] -
      1.0*I*m1d*m2d*pow(Omega,2)*pow(r[j],2)*h[4][l1][m1][j]*h[8][l2][m2][j] - m2d*Omega*h[5][l1][m1][j]*h[8][l2][m2][j] -
      m1d*m2d*pow(Omega,2)*pow(r[j],2)*h[8][l1][m1][j]*h[8][l2][m2][j] + m2d*Omega*h[4][l2][m2][j]*h[9][l1][m1][j] +
      1.0*I*m2d*Omega*h[8][l2][m2][j]*h[9][l1][m1][j] - sqrt(l1d*(1 + l1d))*m2d*Omega*r[j]*h[2][l1][m1][j]*h[9][l2][m2][j] +
      (-(m1d*Omega*h[4][l1][m1][j]) + 1.0*I*m1d*Omega*h[8][l1][m1][j])*h[9][l2][m2][j]) -
      pow(f[j],2)*r[j]*(pow(r[j],2)*dh[5][l1][m1][j]*dh[5][l2][m2][j] - 1.0*I*pow(r[j],2)*dh[5][l2][m2][j]*dh[9][l1][m1][j] +
      1.0*I*pow(r[j],2)*dh[5][l1][m1][j]*dh[9][l2][m2][j] + pow(r[j],2)*dh[9][l1][m1][j]*dh[9][l2][m2][j] +
      r[j]*dh[4][l1][m1][j]*h[4][l2][m2][j] - 1.0*I*r[j]*dh[8][l1][m1][j]*h[4][l2][m2][j] + r[j]*dh[5][l1][m1][j]*h[5][l2][m2][j] -
      1.0*I*r[j]*dh[9][l1][m1][j]*h[5][l2][m2][j] + 1.0*I*r[j]*dh[4][l1][m1][j]*h[8][l2][m2][j] +
      r[j]*dh[8][l1][m1][j]*h[8][l2][m2][j] - 1.0*I*r[j]*dh[5][l2][m2][j]*h[9][l1][m1][j] + r[j]*dh[9][l2][m2][j]*h[9][l1][m1][j] -
      2.0*I*h[5][l2][m2][j]*h[9][l1][m1][j] + sqrt(l1d*(1 + l1d))*h[3][l1][m1][j]*(r[j]*dh[5][l2][m2][j] + 1.0*I*r[j]*dh[9][l2][m2][j]
      + 2.0*h[5][l2][m2][j] + 2.0*I*h[9][l2][m2][j]) + h[5][l1][m1][j]*(r[j]*dh[5][l2][m2][j] + 1.0*I*r[j]*dh[9][l2][m2][j] +
      2.0*h[5][l2][m2][j] + 2.0*I*h[9][l2][m2][j]) + 1.0*I*r[j]*dh[5][l1][m1][j]*h[9][l2][m2][j] +
      r[j]*dh[9][l1][m1][j]*h[9][l2][m2][j] + 2.0*h[9][l1][m1][j]*h[9][l2][m2][j]) + f[j]*(sqrt(l1d*(1 +
      l1d))*pow(r[j],2)*dh[5][l2][m2][j]*h[1][l1][m1][j] + 1.0*I*sqrt(l1d*(1 + l1d))*pow(r[j],2)*dh[9][l2][m2][j]*h[1][l1][m1][j] +
      sqrt(l1d*(1 + l1d))*pow(r[j],2)*dh[4][l2][m2][j]*h[2][l1][m1][j] + 1.0*I*sqrt(l1d*(1 +
      l1d))*pow(r[j],2)*dh[8][l2][m2][j]*h[2][l1][m1][j] + sqrt(l1d*(1 + l1d)*l2d*(1 + l2d))*r[j]*h[1][l1][m1][j]*h[3][l2][m2][j] -
      1.0*I*m1d*Omega*pow(r[j],2)*h[4][l2][m2][j]*h[5][l1][m1][j] + r[j]*dh[5][l1][m1][j]*h[5][l2][m2][j] -
      1.0*I*r[j]*dh[9][l1][m1][j]*h[5][l2][m2][j] + 2.0*sqrt(l1d*(1 + l1d))*r[j]*h[1][l1][m1][j]*h[5][l2][m2][j] +
      m1d*Omega*pow(r[j],2)*h[5][l1][m1][j]*h[8][l2][m2][j] - 1.0*I*r[j]*dh[5][l2][m2][j]*h[9][l1][m1][j] +
      r[j]*dh[9][l2][m2][j]*h[9][l1][m1][j] - m1d*Omega*pow(r[j],2)*h[4][l2][m2][j]*h[9][l1][m1][j] -
      2.0*I*h[5][l2][m2][j]*h[9][l1][m1][j] - 1.0*I*m1d*Omega*pow(r[j],2)*h[8][l2][m2][j]*h[9][l1][m1][j] +
      h[5][l1][m1][j]*(r[j]*(dh[5][l2][m2][j] + 1.0*I*dh[9][l2][m2][j]) + 2.0*h[5][l2][m2][j] + 2.0*I*h[9][l2][m2][j]) +
      1.0*I*r[j]*dh[5][l1][m1][j]*h[9][l2][m2][j] + r[j]*dh[9][l1][m1][j]*h[9][l2][m2][j] + 2.0*I*sqrt(l1d*(1 +
      l1d))*r[j]*h[1][l1][m1][j]*h[9][l2][m2][j] + 2.0*h[9][l1][m1][j]*h[9][l2][m2][j] + h[3][l1][m1][j]*(2.0*sqrt(l1d*(1 +
      l1d))*h[5][l2][m2][j] + r[j]*(sqrt(l1d*(1 + l1d)*l2d*(1 + l2d))*h[1][l2][m2][j] - sqrt(l1d*(1 +
      l1d))*r[j]*(-1.0*I*m2d*Omega*h[4][l2][m2][j] + m2d*Omega*h[8][l2][m2][j])) + 2.0*I*sqrt(l1d*(1 +
      l1d))*h[9][l2][m2][j]))))/(2.*sqrt(2)*f[j]*pow(r[j],3)) +
      (Coupling(l3,m3,0,l1,m1,2,l2,m2,-2)*(r[j]*(-1.0*I*m1d*Omega*h[7][l1][m1][j] -
      m1d*Omega*h[10][l1][m1][j])*(-1.0*I*m2d*Omega*h[7][l2][m2][j] + m2d*Omega*h[10][l2][m2][j]) +
      pow(f[j],2)*(-2.0*(dh[7][l1][m1][j] - 1.0*I*dh[10][l1][m1][j])*h[7][l2][m2][j] + (dh[7][l2][m2][j] +
      1.0*I*dh[10][l2][m2][j])*(r[j]*dh[7][l1][m1][j] - 1.0*I*r[j]*dh[10][l1][m1][j] + 2.0*h[7][l1][m1][j] - 2.0*I*h[10][l1][m1][j]) +
      (-2.0*I*dh[7][l1][m1][j] - 2.0*dh[10][l1][m1][j])*h[10][l2][m2][j])))/(4.*sqrt(2)*r[j]) +
      (Coupling(l3,m3,0,l1,m1,-2,l2,m2,2)*(r[j]*(-1.0*I*m1d*Omega*h[7][l1][m1][j] +
      m1d*Omega*h[10][l1][m1][j])*(-1.0*I*m2d*Omega*h[7][l2][m2][j] - m2d*Omega*h[10][l2][m2][j]) +
      pow(f[j],2)*(-2.0*(dh[7][l1][m1][j] + 1.0*I*dh[10][l1][m1][j])*h[7][l2][m2][j] + (dh[7][l2][m2][j] -
      1.0*I*dh[10][l2][m2][j])*(r[j]*dh[7][l1][m1][j] + 1.0*I*r[j]*dh[10][l1][m1][j] + 2.0*h[7][l1][m1][j] + 2.0*I*h[10][l1][m1][j]) +
      2.0*I*(dh[7][l1][m1][j] + 1.0*I*dh[10][l1][m1][j])*h[10][l2][m2][j])))/(4.*sqrt(2)*r[j]);
  }

  return ricci;
}
