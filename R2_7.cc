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

#include "Coupling.h"
#include "R2.h"

using std::complex;
using std::vector;

vector<complex<double>> R2_7(const double M, const double r0,
                     const vector<double> &r, const vector<double> &f, const vector<double> &fp,
                     const field_type &h, const field_type &dh, const field_type &ddh,
                     int l3, int m3, int l1, int m1, int l2, int m2) {
  vector<complex<double>> ricci(r.size());
  const complex<double> I(0.0, 1.0);
  const double Omega = sqrt(M/(r0*r0*r0));
  const double Omega_2 = Omega*Omega;
  const double l1d = l1, m1d = m1, l2d = l2, m2d = m2;

  const double sqrtl1 = sqrt(l1d*(1. + l1d));
  const double sqrtl2 = sqrt(l1d*(1. + l1d)*l2d*(1. + l2d));
  const double sqrtl3 = sqrt(-2. + l1d + l1d*l1d);
  const double sqrtl4 = sqrt(1. + l1d);
  const double sqrtl5 = sqrt(l1d);
  const double sqrtl6 = sqrt(1. + l2d);
  const double sqrtl7 = sqrt((-2. + l1d + l1d*l1d)*(1. + l2d));
  const double sqrtl8 = sqrt(l2d);
  const double sqrtl9 = sqrt(-6. + l1d + l1d*l1d);
  const double sqrtl10 = sqrt(l1d*(1. + l1d)*(-2. + l2d + l2d*l2d));
  const double sqrtl11 = sqrt((-2. + l1d + l1d*l1d)*l2d);
  const double sqrtl12 = 1/sqrt(1. + l2d);
  const double sqrtl13 = sqrt(-6. + l2d + l2d*l2d);
  const double sqrtl14 = sqrt(-2. + l2d + l2d*l2d);
  const double sqrtl15 = 1/sqrt(1. + l1d);
  const double sqrtl16 = sqrt((-2. + l1d + l1d*l1d)*l2d*(1. + l2d));
  const double sqrtl17 = sqrt(l1d*(-2. + l2d + l2d*l2d));
  const double sqrtl18 = sqrt(l1d*(-2. + l1d*l1d + l1d));
  const double sqrtl19 = sqrt(-2. + 2.*l1d*l1d + l1d*l1d*l1d - l1d);
  const double c1 = Coupling(l3,m3,2,l1,m1,2,l2,m2,0);
  const double c2 = Coupling(l3,m3,-2,l1,m1,-2,l2,m2,0);
  const double c3 = Coupling(l3,m3,2,l1,m1,-1,l2,m2,3);
  const double c4 = Coupling(l3,m3,-2,l1,m1,1,l2,m2,-3);
  const double c5 = Coupling(l3,m3,2,l1,m1,3,l2,m2,-1);
  const double c6 = Coupling(l3,m3,-2,l1,m1,-3,l2,m2,1);
  const double c7 = Coupling(l3,m3,-2,l1,m1,-1,l2,m2,-1);
  const double c8 = Coupling(l3,m3,2,l1,m1,1,l2,m2,1);
  const double c9 = Coupling(l3,m3,-2,l1,m1,0,l2,m2,-2);
  const double c10 = Coupling(l3,m3,2,l1,m1,0,l2,m2,2);

  for(vector<complex<double>>::size_type j=0; j!=ricci.size(); ++j) {
    double r_2 = r[j]*r[j];
    double r_3 = r_2*r[j];
    double f_2 = f[j]*f[j];
    double f_3 = f_2*f[j];

    ricci[j] = (0.17677669529663687*sqrtl12*c1*(2.*f_2*r_2*sqrtl6*dh[6][l2][m2][j]*dh[7][l1][m1][j] - 2.0*I*f_2*r_2*sqrtl6*dh[6][l2][m2][j]*dh[10][l1][m1][j] - 4.*sqrtl6*f[j]*r[j]*dh[7][l1][m1][j]*h[1][l2][m2][j] + 4.0*I*sqrtl6*f[j]*r[j]*dh[10][l1][m1][j]*h[1][l2][m2][j] + 4.*f_2*sqrtl6*r[j]*dh[7][l1][m1][j]*h[3][l2][m2][j] - 4.0*I*f_2*sqrtl6*r[j]*dh[10][l1][m1][j]*h[3][l2][m2][j] - 2.*sqrtl7*f[j]*r[j]*dh[6][l2][m2][j]*h[5][l1][m1][j] + 4.*sqrtl7*h[1][l2][m2][j]*h[5][l1][m1][j] - 4.*sqrtl7*f[j]*h[3][l2][m2][j]*h[5][l1][m1][j] - sqrtl11*h[5][l1][m1][j]*h[5][l2][m2][j] - l2d*sqrtl3*sqrtl8*h[5][l1][m1][j]*h[5][l2][m2][j] - 4.*sqrtl7*f[j]*h[5][l1][m1][j]*h[6][l2][m2][j] + h[4][l1][m1][j]*((sqrtl11 + l2d*sqrtl3*sqrtl8)*h[4][l2][m2][j] - 2.0*I*m2d*sqrtl7*Omega*r[j]*h[6][l2][m2][j]) + 4.*f_2*sqrtl6*r[j]*dh[6][l2][m2][j]*h[7][l1][m1][j] - 8.*sqrtl6*f[j]*h[1][l2][m2][j]*h[7][l1][m1][j] + 8.*f_2*sqrtl6*h[3][l2][m2][j]*h[7][l1][m1][j] + 4.*f_2*sqrtl6*h[6][l2][m2][j]*h[7][l1][m1][j] + 2.*m1d*m2d*r_2*sqrtl6*Omega_2*h[6][l2][m2][j]*h[7][l1][m1][j] - I*sqrtl11*h[4][l2][m2][j]*h[8][l1][m1][j] - I*l2d*sqrtl3*sqrtl8*h[4][l2][m2][j]*h[8][l1][m1][j] - 2.*m2d*sqrtl7*Omega*r[j]*h[6][l2][m2][j]*h[8][l1][m1][j] + 2.0*I*sqrtl7*f[j]*r[j]*dh[6][l2][m2][j]*h[9][l1][m1][j] - 4.0*I*sqrtl7*h[1][l2][m2][j]*h[9][l1][m1][j] + 4.0*I*sqrtl7*f[j]*h[3][l2][m2][j]*h[9][l1][m1][j] + I*sqrtl11*h[5][l2][m2][j]*h[9][l1][m1][j] + I*l2d*sqrtl3*sqrtl8*h[5][l2][m2][j]*h[9][l1][m1][j] + 4.0*I*sqrtl7*f[j]*h[6][l2][m2][j]*h[9][l1][m1][j] - 4.0*I*f_2*sqrtl6*r[j]*dh[6][l2][m2][j]*h[10][l1][m1][j] + 8.0*I*sqrtl6*f[j]*h[1][l2][m2][j]*h[10][l1][m1][j] - 8.0*I*f_2*sqrtl6*h[3][l2][m2][j]*h[10][l1][m1][j] - 4.0*I*f_2*sqrtl6*h[6][l2][m2][j]*h[10][l1][m1][j] - 2.0*I*m1d*m2d*r_2*sqrtl6*Omega_2*h[6][l2][m2][j]*h[10][l1][m1][j] + h[2][l2][m2][j]*(-4.*sqrtl7*h[4][l1][m1][j] + 4.0*I*sqrtl7*h[8][l1][m1][j] + 4.*sqrtl6*r[j]*(-1.0*I*m1d*Omega*h[7][l1][m1][j] - m1d*Omega*h[10][l1][m1][j]))))/(r_2*f[j]) + (0.17677669529663687*sqrtl12*c2*(2.*f_2*r_2*sqrtl6*dh[6][l2][m2][j]*dh[7][l1][m1][j] + 2.0*I*f_2*r_2*sqrtl6*dh[6][l2][m2][j]*dh[10][l1][m1][j] - 4.*sqrtl6*f[j]*r[j]*dh[7][l1][m1][j]*h[1][l2][m2][j] - 4.0*I*sqrtl6*f[j]*r[j]*dh[10][l1][m1][j]*h[1][l2][m2][j] + 4.*f_2*sqrtl6*r[j]*dh[7][l1][m1][j]*h[3][l2][m2][j] + 4.0*I*f_2*sqrtl6*r[j]*dh[10][l1][m1][j]*h[3][l2][m2][j] - 2.*sqrtl7*f[j]*r[j]*dh[6][l2][m2][j]*h[5][l1][m1][j] + 4.*sqrtl7*h[1][l2][m2][j]*h[5][l1][m1][j] - 4.*sqrtl7*f[j]*h[3][l2][m2][j]*h[5][l1][m1][j] - sqrtl11*h[5][l1][m1][j]*h[5][l2][m2][j] - l2d*sqrtl3*sqrtl8*h[5][l1][m1][j]*h[5][l2][m2][j] - 4.*sqrtl7*f[j]*h[5][l1][m1][j]*h[6][l2][m2][j] + h[4][l1][m1][j]*((sqrtl11 + l2d*sqrtl3*sqrtl8)*h[4][l2][m2][j] - 2.0*I*m2d*sqrtl7*Omega*r[j]*h[6][l2][m2][j]) + 4.*f_2*sqrtl6*r[j]*dh[6][l2][m2][j]*h[7][l1][m1][j] - 8.*sqrtl6*f[j]*h[1][l2][m2][j]*h[7][l1][m1][j] + 8.*f_2*sqrtl6*h[3][l2][m2][j]*h[7][l1][m1][j] + 4.*f_2*sqrtl6*h[6][l2][m2][j]*h[7][l1][m1][j] + 2.*m1d*m2d*r_2*sqrtl6*Omega_2*h[6][l2][m2][j]*h[7][l1][m1][j] + I*sqrtl11*h[4][l2][m2][j]*h[8][l1][m1][j] + I*l2d*sqrtl3*sqrtl8*h[4][l2][m2][j]*h[8][l1][m1][j] + 2.*m2d*sqrtl7*Omega*r[j]*h[6][l2][m2][j]*h[8][l1][m1][j] - 2.0*I*sqrtl7*f[j]*r[j]*dh[6][l2][m2][j]*h[9][l1][m1][j] + 4.0*I*sqrtl7*h[1][l2][m2][j]*h[9][l1][m1][j] - 4.0*I*sqrtl7*f[j]*h[3][l2][m2][j]*h[9][l1][m1][j] - I*sqrtl11*h[5][l2][m2][j]*h[9][l1][m1][j] - I*l2d*sqrtl3*sqrtl8*h[5][l2][m2][j]*h[9][l1][m1][j] - 4.0*I*sqrtl7*f[j]*h[6][l2][m2][j]*h[9][l1][m1][j] + 4.0*I*f_2*sqrtl6*r[j]*dh[6][l2][m2][j]*h[10][l1][m1][j] - 8.0*I*sqrtl6*f[j]*h[1][l2][m2][j]*h[10][l1][m1][j] + 8.0*I*f_2*sqrtl6*h[3][l2][m2][j]*h[10][l1][m1][j] + 4.0*I*f_2*sqrtl6*h[6][l2][m2][j]*h[10][l1][m1][j] + 2.0*I*m1d*m2d*r_2*sqrtl6*Omega_2*h[6][l2][m2][j]*h[10][l1][m1][j] + h[2][l2][m2][j]*(-4.*sqrtl7*h[4][l1][m1][j] - 4.0*I*sqrtl7*h[8][l1][m1][j] + 4.*sqrtl6*r[j]*(-1.0*I*m1d*Omega*h[7][l1][m1][j] + m1d*Omega*h[10][l1][m1][j]))))/(r_2*f[j]) + (0.08838834764831843*sqrtl13*sqrtl3*c3*(h[7][l1][m1][j] + I*h[10][l1][m1][j])*(h[7][l2][m2][j] - I*h[10][l2][m2][j]))/r_2 + (0.08838834764831843*sqrtl13*sqrtl3*c4*(h[7][l1][m1][j] - I*h[10][l1][m1][j])*(h[7][l2][m2][j] + I*h[10][l2][m2][j]))/r_2 - (0.08838834764831843*sqrtl9*c5*(h[7][l1][m1][j] - I*h[10][l1][m1][j])*(4.*h[5][l2][m2][j] - sqrtl14*h[7][l2][m2][j] + 4.0*I*h[9][l2][m2][j] - I*sqrtl14*h[10][l2][m2][j]))/r_2 - (0.08838834764831843*sqrtl9*c6*(h[7][l1][m1][j] + I*h[10][l1][m1][j])*(4.*h[5][l2][m2][j] - sqrtl14*h[7][l2][m2][j] - 4.0*I*h[9][l2][m2][j] + I*sqrtl14*h[10][l2][m2][j]))/r_2 + (0.17677669529663687*c7*(-(f_2*r_2*dh[4][l1][m1][j]*dh[4][l2][m2][j]) - I*f_2*r_2*dh[4][l2][m2][j]*dh[8][l1][m1][j] - I*f_2*r_2*dh[4][l1][m1][j]*dh[8][l2][m2][j] + f_2*r_2*dh[8][l1][m1][j]*dh[8][l2][m2][j] + f_2*sqrtl2*h[3][l1][m1][j]*h[3][l2][m2][j] - f_2*r[j]*dh[4][l2][m2][j]*h[4][l1][m1][j] - I*f_2*r[j]*dh[8][l2][m2][j]*h[4][l1][m1][j] + f_2*r[j]*dh[4][l1][m1][j]*h[4][l2][m2][j] + I*f_2*r[j]*dh[8][l1][m1][j]*h[4][l2][m2][j] - f_2*h[4][l1][m1][j]*h[4][l2][m2][j] - I*m1d*r_2*Omega*f[j]*dh[4][l2][m2][j]*h[5][l1][m1][j] + m1d*r_2*Omega*f[j]*dh[8][l2][m2][j]*h[5][l1][m1][j] + I*m1d*Omega*f[j]*r[j]*h[4][l2][m2][j]*h[5][l1][m1][j] - I*m2d*r_2*Omega*f[j]*dh[4][l1][m1][j]*h[5][l2][m2][j] + m2d*r_2*Omega*f[j]*dh[8][l1][m1][j]*h[5][l2][m2][j] + 2.*f_2*sqrtl1*h[3][l1][m1][j]*h[5][l2][m2][j] - I*m2d*Omega*f[j]*r[j]*h[4][l1][m1][j]*h[5][l2][m2][j] + 2.*f_2*h[5][l1][m1][j]*h[5][l2][m2][j] + m1d*m2d*r_2*Omega_2*h[5][l1][m1][j]*h[5][l2][m2][j] + 4.*f_2*sqrtl1*h[5][l2][m2][j]*h[6][l1][m1][j] + 2.*f_2*sqrtl3*h[5][l2][m2][j]*h[7][l1][m1][j] - f_2*sqrtl16*h[6][l2][m2][j]*h[7][l1][m1][j] - f_2*sqrtl10*h[6][l1][m1][j]*h[7][l2][m2][j] - f_2*sqrtl14*sqrtl3*h[7][l1][m1][j]*h[7][l2][m2][j] - I*f_2*r[j]*dh[4][l2][m2][j]*h[8][l1][m1][j] + f_2*r[j]*dh[8][l2][m2][j]*h[8][l1][m1][j] - I*f_2*h[4][l2][m2][j]*h[8][l1][m1][j] + m2d*Omega*f[j]*r[j]*h[5][l2][m2][j]*h[8][l1][m1][j] + h[2][l1][m1][j]*(-(sqrtl2*h[2][l2][m2][j]) + 2.*sqrtl1*f[j]*(h[4][l2][m2][j] + I*h[8][l2][m2][j])) + I*f_2*r[j]*dh[4][l1][m1][j]*h[8][l2][m2][j] - f_2*r[j]*dh[8][l1][m1][j]*h[8][l2][m2][j] - I*f_2*h[4][l1][m1][j]*h[8][l2][m2][j] - m1d*Omega*f[j]*r[j]*h[5][l1][m1][j]*h[8][l2][m2][j] + f_2*h[8][l1][m1][j]*h[8][l2][m2][j] + m1d*r_2*Omega*f[j]*dh[4][l2][m2][j]*h[9][l1][m1][j] + I*m1d*r_2*Omega*f[j]*dh[8][l2][m2][j]*h[9][l1][m1][j] - m1d*Omega*f[j]*r[j]*h[4][l2][m2][j]*h[9][l1][m1][j] + 2.0*I*f_2*h[5][l2][m2][j]*h[9][l1][m1][j] + I*m1d*m2d*r_2*Omega_2*h[5][l2][m2][j]*h[9][l1][m1][j] - I*m1d*Omega*f[j]*r[j]*h[8][l2][m2][j]*h[9][l1][m1][j] + h[1][l1][m1][j]*(sqrtl2*h[1][l2][m2][j] - 2.*sqrtl1*f[j]*(h[5][l2][m2][j] + I*h[9][l2][m2][j])) + m2d*r_2*Omega*f[j]*dh[4][l1][m1][j]*h[9][l2][m2][j] + I*m2d*r_2*Omega*f[j]*dh[8][l1][m1][j]*h[9][l2][m2][j] + 2.0*I*f_2*sqrtl1*h[3][l1][m1][j]*h[9][l2][m2][j] + m2d*Omega*f[j]*r[j]*h[4][l1][m1][j]*h[9][l2][m2][j] + 2.0*I*f_2*h[5][l1][m1][j]*h[9][l2][m2][j] + I*m1d*m2d*r_2*Omega_2*h[5][l1][m1][j]*h[9][l2][m2][j] + 4.0*I*f_2*sqrtl1*h[6][l1][m1][j]*h[9][l2][m2][j] + 2.0*I*f_2*sqrtl3*h[7][l1][m1][j]*h[9][l2][m2][j] + I*m2d*Omega*f[j]*r[j]*h[8][l1][m1][j]*h[9][l2][m2][j] - 2.*f_2*h[9][l1][m1][j]*h[9][l2][m2][j] - m1d*m2d*r_2*Omega_2*h[9][l1][m1][j]*h[9][l2][m2][j] + 2.0*I*f_2*sqrtl3*h[5][l2][m2][j]*h[10][l1][m1][j] - I*f_2*sqrtl16*h[6][l2][m2][j]*h[10][l1][m1][j] - I*f_2*sqrtl14*sqrtl3*h[7][l2][m2][j]*h[10][l1][m1][j] - 2.*f_2*sqrtl3*h[9][l2][m2][j]*h[10][l1][m1][j] - I*f_2*sqrtl10*h[6][l1][m1][j]*h[10][l2][m2][j] - I*f_2*sqrtl14*sqrtl3*h[7][l1][m1][j]*h[10][l2][m2][j] + f_2*sqrtl14*sqrtl3*h[10][l1][m1][j]*h[10][l2][m2][j]))/(f_2*r_2) + (0.17677669529663687*c8*(-(f_2*r_2*dh[4][l1][m1][j]*dh[4][l2][m2][j]) + I*f_2*r_2*dh[4][l2][m2][j]*dh[8][l1][m1][j] + I*f_2*r_2*dh[4][l1][m1][j]*dh[8][l2][m2][j] + f_2*r_2*dh[8][l1][m1][j]*dh[8][l2][m2][j] + f_2*sqrtl2*h[3][l1][m1][j]*h[3][l2][m2][j] - f_2*r[j]*dh[4][l2][m2][j]*h[4][l1][m1][j] + I*f_2*r[j]*dh[8][l2][m2][j]*h[4][l1][m1][j] + f_2*r[j]*dh[4][l1][m1][j]*h[4][l2][m2][j] - I*f_2*r[j]*dh[8][l1][m1][j]*h[4][l2][m2][j] - f_2*h[4][l1][m1][j]*h[4][l2][m2][j] - I*m1d*r_2*Omega*f[j]*dh[4][l2][m2][j]*h[5][l1][m1][j] - m1d*r_2*Omega*f[j]*dh[8][l2][m2][j]*h[5][l1][m1][j] + I*m1d*Omega*f[j]*r[j]*h[4][l2][m2][j]*h[5][l1][m1][j] - I*m2d*r_2*Omega*f[j]*dh[4][l1][m1][j]*h[5][l2][m2][j] - m2d*r_2*Omega*f[j]*dh[8][l1][m1][j]*h[5][l2][m2][j] + 2.*f_2*sqrtl1*h[3][l1][m1][j]*h[5][l2][m2][j] - I*m2d*Omega*f[j]*r[j]*h[4][l1][m1][j]*h[5][l2][m2][j] + 2.*f_2*h[5][l1][m1][j]*h[5][l2][m2][j] + m1d*m2d*r_2*Omega_2*h[5][l1][m1][j]*h[5][l2][m2][j] + 4.*f_2*sqrtl1*h[5][l2][m2][j]*h[6][l1][m1][j] + 2.*f_2*sqrtl3*h[5][l2][m2][j]*h[7][l1][m1][j] - f_2*sqrtl16*h[6][l2][m2][j]*h[7][l1][m1][j] - f_2*sqrtl10*h[6][l1][m1][j]*h[7][l2][m2][j] - f_2*sqrtl14*sqrtl3*h[7][l1][m1][j]*h[7][l2][m2][j] + I*f_2*r[j]*dh[4][l2][m2][j]*h[8][l1][m1][j] + f_2*r[j]*dh[8][l2][m2][j]*h[8][l1][m1][j] + I*f_2*h[4][l2][m2][j]*h[8][l1][m1][j] - m2d*Omega*f[j]*r[j]*h[5][l2][m2][j]*h[8][l1][m1][j] - h[2][l1][m1][j]*(sqrtl2*h[2][l2][m2][j] - 2.*sqrtl1*f[j]*(h[4][l2][m2][j] - I*h[8][l2][m2][j])) - I*f_2*r[j]*dh[4][l1][m1][j]*h[8][l2][m2][j] - f_2*r[j]*dh[8][l1][m1][j]*h[8][l2][m2][j] + I*f_2*h[4][l1][m1][j]*h[8][l2][m2][j] + m1d*Omega*f[j]*r[j]*h[5][l1][m1][j]*h[8][l2][m2][j] + f_2*h[8][l1][m1][j]*h[8][l2][m2][j] - m1d*r_2*Omega*f[j]*dh[4][l2][m2][j]*h[9][l1][m1][j] + I*m1d*r_2*Omega*f[j]*dh[8][l2][m2][j]*h[9][l1][m1][j] + m1d*Omega*f[j]*r[j]*h[4][l2][m2][j]*h[9][l1][m1][j] - 2.0*I*f_2*h[5][l2][m2][j]*h[9][l1][m1][j] - I*m1d*m2d*r_2*Omega_2*h[5][l2][m2][j]*h[9][l1][m1][j] - I*m1d*Omega*f[j]*r[j]*h[8][l2][m2][j]*h[9][l1][m1][j] + h[1][l1][m1][j]*(sqrtl2*h[1][l2][m2][j] - 2.*sqrtl1*f[j]*(h[5][l2][m2][j] - I*h[9][l2][m2][j])) - m2d*r_2*Omega*f[j]*dh[4][l1][m1][j]*h[9][l2][m2][j] + I*m2d*r_2*Omega*f[j]*dh[8][l1][m1][j]*h[9][l2][m2][j] - 2.0*I*f_2*sqrtl1*h[3][l1][m1][j]*h[9][l2][m2][j] - m2d*Omega*f[j]*r[j]*h[4][l1][m1][j]*h[9][l2][m2][j] - 2.0*I*f_2*h[5][l1][m1][j]*h[9][l2][m2][j] - I*m1d*m2d*r_2*Omega_2*h[5][l1][m1][j]*h[9][l2][m2][j] - 4.0*I*f_2*sqrtl1*h[6][l1][m1][j]*h[9][l2][m2][j] - 2.0*I*f_2*sqrtl3*h[7][l1][m1][j]*h[9][l2][m2][j] + I*m2d*Omega*f[j]*r[j]*h[8][l1][m1][j]*h[9][l2][m2][j] - 2.*f_2*h[9][l1][m1][j]*h[9][l2][m2][j] - m1d*m2d*r_2*Omega_2*h[9][l1][m1][j]*h[9][l2][m2][j] - 2.0*I*f_2*sqrtl3*h[5][l2][m2][j]*h[10][l1][m1][j] + I*f_2*sqrtl16*h[6][l2][m2][j]*h[10][l1][m1][j] + I*f_2*sqrtl14*sqrtl3*h[7][l2][m2][j]*h[10][l1][m1][j] - 2.*f_2*sqrtl3*h[9][l2][m2][j]*h[10][l1][m1][j] + I*f_2*sqrtl10*h[6][l1][m1][j]*h[10][l2][m2][j] + I*f_2*sqrtl14*sqrtl3*h[7][l1][m1][j]*h[10][l2][m2][j] + f_2*sqrtl14*sqrtl3*h[10][l1][m1][j]*h[10][l2][m2][j]))/(f_2*r_2) + (0.17677669529663687*sqrtl15*c9*(-(h[5][l1][m1][j]*((sqrtl17 + l1d*sqrtl14*sqrtl5)*h[5][l2][m2][j] + I*(sqrtl17 + l1d*sqrtl14*sqrtl5)*h[9][l2][m2][j] - 2.*(1 + l1d)*sqrtl5*f[j]*(r[j]*dh[7][l2][m2][j] + I*r[j]*dh[10][l2][m2][j] + 2.*h[7][l2][m2][j] + 2.0*I*h[10][l2][m2][j]))) + 2.*(I*(1 + l1d)*sqrtl5*f[j]*r[j]*(dh[7][l2][m2][j] + I*dh[10][l2][m2][j])*h[9][l1][m1][j] + f_2*sqrtl4*(r_2*dh[6][l1][m1][j]*(dh[7][l2][m2][j] + I*dh[10][l2][m2][j]) + 2.*h[6][l1][m1][j]*(r[j]*dh[7][l2][m2][j] + I*r[j]*dh[10][l2][m2][j] + h[7][l2][m2][j] + I*h[10][l2][m2][j])) + r[j]*(-(m1d*sqrtl4*Omega*r[j]*h[6][l1][m1][j]) + (1 + l1d)*sqrtl5*h[8][l1][m1][j])*(-(m2d*Omega*h[7][l2][m2][j]) - I*m2d*Omega*h[10][l2][m2][j])) + h[4][l1][m1][j]*((sqrtl17 + l1d*sqrtl14*sqrtl5)*h[4][l2][m2][j] + I*((sqrtl17 + l1d*sqrtl14*sqrtl5)*h[8][l2][m2][j] + 2.0*I*(1 + l1d)*sqrtl5*r[j]*(-1.0*I*m2d*Omega*h[7][l2][m2][j] + m2d*Omega*h[10][l2][m2][j])))))/(r_2*f[j]) + (0.17677669529663687*sqrtl15*c10*(-(h[5][l1][m1][j]*((sqrtl17 + l1d*sqrtl14*sqrtl5)*h[5][l2][m2][j] - I*(sqrtl17 + l1d*sqrtl14*sqrtl5)*h[9][l2][m2][j] - 2.*(1 + l1d)*sqrtl5*f[j]*(r[j]*dh[7][l2][m2][j] - I*r[j]*dh[10][l2][m2][j] + 2.*h[7][l2][m2][j] - 2.0*I*h[10][l2][m2][j]))) + h[4][l1][m1][j]*((sqrtl17 + l1d*sqrtl14*sqrtl5)*h[4][l2][m2][j] - I*((sqrtl17 + l1d*sqrtl14*sqrtl5)*h[8][l2][m2][j] - 2.0*I*(1 + l1d)*sqrtl5*r[j]*(-1.0*I*m2d*Omega*h[7][l2][m2][j] - m2d*Omega*h[10][l2][m2][j]))) + 2.*(-1.0*I*(1 + l1d)*sqrtl5*f[j]*r[j]*(dh[7][l2][m2][j] - I*dh[10][l2][m2][j])*h[9][l1][m1][j] + r[j]*(m1d*sqrtl4*Omega*r[j]*h[6][l1][m1][j] + (1 + l1d)*sqrtl5*h[8][l1][m1][j])*(m2d*Omega*h[7][l2][m2][j] - I*m2d*Omega*h[10][l2][m2][j]) + f_2*sqrtl4*(r_2*dh[6][l1][m1][j]*(dh[7][l2][m2][j] - I*dh[10][l2][m2][j]) + 2.*h[6][l1][m1][j]*(h[7][l2][m2][j] - I*(I*r[j]*dh[7][l2][m2][j] + r[j]*dh[10][l2][m2][j] + h[10][l2][m2][j]))))))/(r_2*f[j]);
    ricci[j]+= (-0.35355339059327373*sqrtl15*c7*(-(sqrtl4*f[j]*r[j]*dh[4][l1][m1][j]*h[4][l2][m2][j]) - I*sqrtl4*f[j]*r[j]*dh[8][l1][m1][j]*h[4][l2][m2][j] - I*m1d*sqrtl4*Omega*r[j]*h[4][l2][m2][j]*h[5][l1][m1][j] + sqrtl5*f[j]*r[j]*dh[6][l1][m1][j]*h[5][l2][m2][j] + l1d*sqrtl5*f[j]*r[j]*dh[6][l1][m1][j]*h[5][l2][m2][j] + sqrtl19*f[j]*r[j]*dh[7][l1][m1][j]*h[5][l2][m2][j] + I*sqrtl19*f[j]*r[j]*dh[10][l1][m1][j]*h[5][l2][m2][j] - sqrtl5*h[1][l1][m1][j]*h[5][l2][m2][j] - l1d*sqrtl5*h[1][l1][m1][j]*h[5][l2][m2][j] + sqrtl5*f[j]*h[3][l1][m1][j]*h[5][l2][m2][j] + l1d*sqrtl5*f[j]*h[3][l1][m1][j]*h[5][l2][m2][j] + sqrtl4*h[5][l1][m1][j]*h[5][l2][m2][j] + sqrtl4*f[j]*h[5][l1][m1][j]*h[5][l2][m2][j] + I*m1d*sqrtl5*Omega*r[j]*h[4][l2][m2][j]*h[6][l1][m1][j] + I*l1d*m1d*sqrtl5*Omega*r[j]*h[4][l2][m2][j]*h[6][l1][m1][j] + I*m1d*sqrtl19*Omega*r[j]*h[4][l2][m2][j]*h[7][l1][m1][j] - I*sqrtl4*h[4][l2][m2][j]*h[8][l1][m1][j] + I*l1d*l1d*sqrtl4*h[4][l2][m2][j]*h[8][l1][m1][j] + I*l1d*sqrtl4*h[4][l2][m2][j]*h[8][l1][m1][j] + (1 + l1d)*sqrtl5*h[2][l1][m1][j]*(h[4][l2][m2][j] + I*h[8][l2][m2][j]) - sqrtl4*h[4][l1][m1][j]*(h[4][l2][m2][j] + I*h[8][l2][m2][j]) - I*sqrtl4*f[j]*r[j]*dh[4][l1][m1][j]*h[8][l2][m2][j] + sqrtl4*f[j]*r[j]*dh[8][l1][m1][j]*h[8][l2][m2][j] + m1d*sqrtl4*Omega*r[j]*h[5][l1][m1][j]*h[8][l2][m2][j] - m1d*sqrtl5*Omega*r[j]*h[6][l1][m1][j]*h[8][l2][m2][j] - l1d*m1d*sqrtl5*Omega*r[j]*h[6][l1][m1][j]*h[8][l2][m2][j] - m1d*sqrtl19*Omega*r[j]*h[7][l1][m1][j]*h[8][l2][m2][j] + sqrtl4*h[8][l1][m1][j]*h[8][l2][m2][j] - l1d*l1d*sqrtl4*h[8][l1][m1][j]*h[8][l2][m2][j] - l1d*sqrtl4*h[8][l1][m1][j]*h[8][l2][m2][j] + m1d*sqrtl4*Omega*r[j]*h[4][l2][m2][j]*h[9][l1][m1][j] + I*sqrtl4*h[5][l2][m2][j]*h[9][l1][m1][j] - I*l1d*l1d*sqrtl4*h[5][l2][m2][j]*h[9][l1][m1][j] - I*l1d*sqrtl4*h[5][l2][m2][j]*h[9][l1][m1][j] + I*sqrtl4*f[j]*h[5][l2][m2][j]*h[9][l1][m1][j] + I*m1d*sqrtl4*Omega*r[j]*h[8][l2][m2][j]*h[9][l1][m1][j] + I*sqrtl5*f[j]*r[j]*dh[6][l1][m1][j]*h[9][l2][m2][j] + I*l1d*sqrtl5*f[j]*r[j]*dh[6][l1][m1][j]*h[9][l2][m2][j] + I*sqrtl19*f[j]*r[j]*dh[7][l1][m1][j]*h[9][l2][m2][j] - sqrtl19*f[j]*r[j]*dh[10][l1][m1][j]*h[9][l2][m2][j] - I*sqrtl5*h[1][l1][m1][j]*h[9][l2][m2][j] - I*l1d*sqrtl5*h[1][l1][m1][j]*h[9][l2][m2][j] + I*sqrtl5*f[j]*h[3][l1][m1][j]*h[9][l2][m2][j] + I*l1d*sqrtl5*f[j]*h[3][l1][m1][j]*h[9][l2][m2][j] + I*sqrtl4*h[5][l1][m1][j]*h[9][l2][m2][j] + I*sqrtl4*f[j]*h[5][l1][m1][j]*h[9][l2][m2][j] - sqrtl4*h[9][l1][m1][j]*h[9][l2][m2][j] + l1d*l1d*sqrtl4*h[9][l1][m1][j]*h[9][l2][m2][j] + l1d*sqrtl4*h[9][l1][m1][j]*h[9][l2][m2][j] - sqrtl4*f[j]*h[9][l1][m1][j]*h[9][l2][m2][j] - m1d*sqrtl19*Omega*r[j]*h[4][l2][m2][j]*h[10][l1][m1][j] - I*m1d*sqrtl19*Omega*r[j]*h[8][l2][m2][j]*h[10][l1][m1][j]))/(r_2*f[j]) + (0.35355339059327373*sqrtl15*c8*(sqrtl4*f[j]*r[j]*dh[4][l1][m1][j]*h[4][l2][m2][j] - I*sqrtl4*f[j]*r[j]*dh[8][l1][m1][j]*h[4][l2][m2][j] + I*m1d*sqrtl4*Omega*r[j]*h[4][l2][m2][j]*h[5][l1][m1][j] - sqrtl5*f[j]*r[j]*dh[6][l1][m1][j]*h[5][l2][m2][j] - l1d*sqrtl5*f[j]*r[j]*dh[6][l1][m1][j]*h[5][l2][m2][j] - sqrtl19*f[j]*r[j]*dh[7][l1][m1][j]*h[5][l2][m2][j] + I*sqrtl19*f[j]*r[j]*dh[10][l1][m1][j]*h[5][l2][m2][j] + sqrtl5*h[1][l1][m1][j]*h[5][l2][m2][j] + l1d*sqrtl5*h[1][l1][m1][j]*h[5][l2][m2][j] - sqrtl5*f[j]*h[3][l1][m1][j]*h[5][l2][m2][j] - l1d*sqrtl5*f[j]*h[3][l1][m1][j]*h[5][l2][m2][j] - sqrtl4*h[5][l1][m1][j]*h[5][l2][m2][j] - sqrtl4*f[j]*h[5][l1][m1][j]*h[5][l2][m2][j] - I*m1d*sqrtl5*Omega*r[j]*h[4][l2][m2][j]*h[6][l1][m1][j] - I*l1d*m1d*sqrtl5*Omega*r[j]*h[4][l2][m2][j]*h[6][l1][m1][j] - I*m1d*sqrtl19*Omega*r[j]*h[4][l2][m2][j]*h[7][l1][m1][j] - I*sqrtl4*h[4][l2][m2][j]*h[8][l1][m1][j] + I*l1d*l1d*sqrtl4*h[4][l2][m2][j]*h[8][l1][m1][j] + I*l1d*sqrtl4*h[4][l2][m2][j]*h[8][l1][m1][j] - (1 + l1d)*sqrtl5*h[2][l1][m1][j]*(h[4][l2][m2][j] - I*h[8][l2][m2][j]) + sqrtl4*h[4][l1][m1][j]*(h[4][l2][m2][j] - I*h[8][l2][m2][j]) - I*sqrtl4*f[j]*r[j]*dh[4][l1][m1][j]*h[8][l2][m2][j] - sqrtl4*f[j]*r[j]*dh[8][l1][m1][j]*h[8][l2][m2][j] + m1d*sqrtl4*Omega*r[j]*h[5][l1][m1][j]*h[8][l2][m2][j] - m1d*sqrtl5*Omega*r[j]*h[6][l1][m1][j]*h[8][l2][m2][j] - l1d*m1d*sqrtl5*Omega*r[j]*h[6][l1][m1][j]*h[8][l2][m2][j] - m1d*sqrtl19*Omega*r[j]*h[7][l1][m1][j]*h[8][l2][m2][j] - sqrtl4*h[8][l1][m1][j]*h[8][l2][m2][j] + l1d*l1d*sqrtl4*h[8][l1][m1][j]*h[8][l2][m2][j] + l1d*sqrtl4*h[8][l1][m1][j]*h[8][l2][m2][j] + m1d*sqrtl4*Omega*r[j]*h[4][l2][m2][j]*h[9][l1][m1][j] + I*sqrtl4*h[5][l2][m2][j]*h[9][l1][m1][j] - I*l1d*l1d*sqrtl4*h[5][l2][m2][j]*h[9][l1][m1][j] - I*l1d*sqrtl4*h[5][l2][m2][j]*h[9][l1][m1][j] + I*sqrtl4*f[j]*h[5][l2][m2][j]*h[9][l1][m1][j] - I*m1d*sqrtl4*Omega*r[j]*h[8][l2][m2][j]*h[9][l1][m1][j] + I*sqrtl5*f[j]*r[j]*dh[6][l1][m1][j]*h[9][l2][m2][j] + I*l1d*sqrtl5*f[j]*r[j]*dh[6][l1][m1][j]*h[9][l2][m2][j] + I*sqrtl19*f[j]*r[j]*dh[7][l1][m1][j]*h[9][l2][m2][j] + sqrtl19*f[j]*r[j]*dh[10][l1][m1][j]*h[9][l2][m2][j] - I*sqrtl5*h[1][l1][m1][j]*h[9][l2][m2][j] - I*l1d*sqrtl5*h[1][l1][m1][j]*h[9][l2][m2][j] + I*sqrtl5*f[j]*h[3][l1][m1][j]*h[9][l2][m2][j] + I*l1d*sqrtl5*f[j]*h[3][l1][m1][j]*h[9][l2][m2][j] + I*sqrtl4*h[5][l1][m1][j]*h[9][l2][m2][j] + I*sqrtl4*f[j]*h[5][l1][m1][j]*h[9][l2][m2][j] + sqrtl4*h[9][l1][m1][j]*h[9][l2][m2][j] - l1d*l1d*sqrtl4*h[9][l1][m1][j]*h[9][l2][m2][j] - l1d*sqrtl4*h[9][l1][m1][j]*h[9][l2][m2][j] + sqrtl4*f[j]*h[9][l1][m1][j]*h[9][l2][m2][j] - m1d*sqrtl19*Omega*r[j]*h[4][l2][m2][j]*h[10][l1][m1][j] + I*m1d*sqrtl19*Omega*r[j]*h[8][l2][m2][j]*h[10][l1][m1][j]))/(r_2*f[j]) + (0.35355339059327373*sqrtl15*c2*(f_2*(sqrtl4*r[j]*(-2.*h[6][l2][m2][j]*(h[7][l1][m1][j] + I*h[10][l1][m1][j]) + h[1][l2][m2][j]*(r[j]*(r[j]*ddh[7][l1][m1][j] + I*r[j]*ddh[10][l1][m1][j] + 2.*dh[7][l1][m1][j] + 2.0*I*dh[10][l1][m1][j]) - 2.*h[7][l1][m1][j] - 2.0*I*h[10][l1][m1][j])) + h[3][l2][m2][j]*(-4.*sqrtl4*h[7][l1][m1][j] + r[j]*(sqrtl19*r[j]*dh[5][l1][m1][j] - 2.*sqrtl4*dh[7][l1][m1][j] + I*sqrtl19*r[j]*dh[9][l1][m1][j] - 2.0*I*sqrtl4*dh[10][l1][m1][j] + (sqrtl18 + l1d*sqrtl3*sqrtl5)*h[3][l1][m1][j] + sqrtl19*h[5][l1][m1][j] + I*sqrtl19*h[9][l1][m1][j]) - 4.0*I*sqrtl4*h[10][l1][m1][j])) + r[j]*h[1][l2][m2][j]*((sqrtl18 + l1d*sqrtl3*sqrtl5)*h[1][l1][m1][j] + r[j]*(I*m1d*sqrtl19*Omega*h[4][l1][m1][j] + sqrtl19*fp[j]*h[5][l1][m1][j] - m1d*m1d*sqrtl4*Omega_2*r[j]*h[7][l1][m1][j] - m1d*sqrtl19*Omega*h[8][l1][m1][j] + I*sqrtl19*fp[j]*h[9][l1][m1][j] - I*m1d*m1d*sqrtl4*Omega_2*r[j]*h[10][l1][m1][j])) + h[2][l2][m2][j]*(sqrtl19*(-2. + f[j]*r[j])*h[4][l1][m1][j] - I*sqrtl19*(2. - f[j]*r[j])*h[8][l1][m1][j] + r[j]*(sqrtl19*f[j]*r[j]*dh[4][l1][m1][j] + 2.0*I*m1d*r_2*sqrtl4*Omega*f[j]*dh[7][l1][m1][j] + I*sqrtl19*f[j]*r[j]*dh[8][l1][m1][j] - 2.*m1d*r_2*sqrtl4*Omega*f[j]*dh[10][l1][m1][j] - (sqrtl18 + l1d*sqrtl3*sqrtl5)*h[2][l1][m1][j] - I*m1d*sqrtl19*Omega*r[j]*h[5][l1][m1][j] - 2.0*I*m1d*sqrtl4*Omega*(1 - f[j]*r[j])*h[7][l1][m1][j] + m1d*sqrtl19*Omega*r[j]*h[9][l1][m1][j] + 2.*m1d*sqrtl4*Omega*h[10][l1][m1][j] - 2.*m1d*sqrtl4*Omega*f[j]*r[j]*h[10][l1][m1][j])) + f_3*sqrtl4*r[j]*(2.*h[6][l2][m2][j]*(h[7][l1][m1][j] + I*h[10][l1][m1][j]) + h[3][l2][m2][j]*(2.*h[7][l1][m1][j] + I*(I*r[j]*(r[j]*ddh[7][l1][m1][j] + I*r[j]*ddh[10][l1][m1][j] + 2.*dh[7][l1][m1][j] + 2.0*I*dh[10][l1][m1][j]) + 2.*h[10][l1][m1][j]))) + f[j]*(-(sqrtl19*r[j]*h[1][l2][m2][j]*(r[j]*(dh[5][l1][m1][j] + I*dh[9][l1][m1][j]) + h[5][l1][m1][j] + I*h[9][l1][m1][j])) + h[3][l2][m2][j]*(sqrtl19*(2. - r_2*fp[j])*h[5][l1][m1][j] + I*sqrtl19*(2. - r_2*fp[j])*h[9][l1][m1][j] + r_2*(I*m1d*sqrtl19*Omega*h[4][l1][m1][j] - m1d*sqrtl19*Omega*h[8][l1][m1][j] + sqrtl4*r[j]*(-(m1d*m1d*Omega_2*h[7][l1][m1][j]) - I*m1d*m1d*Omega_2*h[10][l1][m1][j]))))))/(f_2*r_3) + (0.35355339059327373*sqrtl15*c1*(f_3*sqrtl4*r[j]*(2.*h[6][l2][m2][j]*(h[7][l1][m1][j] - I*h[10][l1][m1][j]) + h[3][l2][m2][j]*(r[j]*(-(r[j]*ddh[7][l1][m1][j]) + I*r[j]*ddh[10][l1][m1][j] - 2.*dh[7][l1][m1][j] + 2.0*I*dh[10][l1][m1][j]) + 2.*h[7][l1][m1][j] - 2.0*I*h[10][l1][m1][j])) + f_2*(sqrtl4*r[j]*(-2.*h[6][l2][m2][j]*(h[7][l1][m1][j] - I*h[10][l1][m1][j]) + h[1][l2][m2][j]*(r[j]*(r[j]*ddh[7][l1][m1][j] - I*r[j]*ddh[10][l1][m1][j] + 2.*dh[7][l1][m1][j] - 2.0*I*dh[10][l1][m1][j]) - 2.*h[7][l1][m1][j] + 2.0*I*h[10][l1][m1][j])) + h[3][l2][m2][j]*(-4.*sqrtl4*h[7][l1][m1][j] + r[j]*(sqrtl19*r[j]*dh[5][l1][m1][j] - 2.*sqrtl4*dh[7][l1][m1][j] - I*sqrtl19*r[j]*dh[9][l1][m1][j] + 2.0*I*sqrtl4*dh[10][l1][m1][j] + (sqrtl18 + l1d*sqrtl3*sqrtl5)*h[3][l1][m1][j] + sqrtl19*h[5][l1][m1][j] - I*sqrtl19*h[9][l1][m1][j]) + 4.0*I*sqrtl4*h[10][l1][m1][j])) + r[j]*h[1][l2][m2][j]*((sqrtl18 + l1d*sqrtl3*sqrtl5)*h[1][l1][m1][j] + r[j]*(I*m1d*sqrtl19*Omega*h[4][l1][m1][j] + sqrtl19*fp[j]*h[5][l1][m1][j] - m1d*m1d*sqrtl4*Omega_2*r[j]*h[7][l1][m1][j] + m1d*sqrtl19*Omega*h[8][l1][m1][j] - I*sqrtl19*fp[j]*h[9][l1][m1][j] + I*m1d*m1d*sqrtl4*Omega_2*r[j]*h[10][l1][m1][j])) + h[2][l2][m2][j]*(-(sqrtl19*(2. - f[j]*r[j])*h[4][l1][m1][j]) + I*sqrtl19*(2. - f[j]*r[j])*h[8][l1][m1][j] - r[j]*(-(sqrtl19*f[j]*r[j]*dh[4][l1][m1][j]) - 2.0*I*m1d*r_2*sqrtl4*Omega*f[j]*dh[7][l1][m1][j] + I*sqrtl19*f[j]*r[j]*dh[8][l1][m1][j] - 2.*m1d*r_2*sqrtl4*Omega*f[j]*dh[10][l1][m1][j] + (sqrtl18 + l1d*sqrtl3*sqrtl5)*h[2][l1][m1][j] + I*m1d*sqrtl19*Omega*r[j]*h[5][l1][m1][j] + 2.0*I*m1d*sqrtl4*Omega*(1 - f[j]*r[j])*h[7][l1][m1][j] + m1d*sqrtl19*Omega*r[j]*h[9][l1][m1][j] + 2.*m1d*sqrtl4*Omega*h[10][l1][m1][j] - 2.*m1d*sqrtl4*Omega*f[j]*r[j]*h[10][l1][m1][j])) + f[j]*(-(sqrtl19*r[j]*h[1][l2][m2][j]*(h[5][l1][m1][j] - I*(r[j]*(I*dh[5][l1][m1][j] + dh[9][l1][m1][j]) + h[9][l1][m1][j]))) + h[3][l2][m2][j]*(sqrtl19*(2. - r_2*fp[j])*h[5][l1][m1][j] - I*sqrtl19*(2. - r_2*fp[j])*h[9][l1][m1][j] + r_2*(I*m1d*sqrtl19*Omega*h[4][l1][m1][j] + m1d*sqrtl19*Omega*h[8][l1][m1][j] + sqrtl4*r[j]*(-(m1d*m1d*Omega_2*h[7][l1][m1][j]) + I*m1d*m1d*Omega_2*h[10][l1][m1][j]))))))/(f_2*r_3) + (0.35355339059327373*sqrtl15*c10*(4.*sqrtl4*h[1][l1][m1][j] - 2.*sqrtl4*f[j]*(r[j]*dh[6][l1][m1][j] + 2.*h[3][l1][m1][j] + h[6][l1][m1][j]) + (1 + l1d)*sqrtl5*(-2.*h[5][l1][m1][j] + sqrtl1*h[6][l1][m1][j]))*(h[7][l2][m2][j] - I*h[10][l2][m2][j]))/r_2 + (0.35355339059327373*sqrtl15*c9*(4.*sqrtl4*h[1][l1][m1][j] - 2.*sqrtl4*f[j]*(r[j]*dh[6][l1][m1][j] + 2.*h[3][l1][m1][j] + h[6][l1][m1][j]) + (1 + l1d)*sqrtl5*(-2.*h[5][l1][m1][j] + sqrtl1*h[6][l1][m1][j]))*(h[7][l2][m2][j] + I*h[10][l2][m2][j]))/r_2;
  }

  return ricci;
}
