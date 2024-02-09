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
#include "R2.h"

using namespace std;

vector<complex<double>> R2_5(const double M, const double r0,
                     const vector<double> &r_grid, const vector<double> &f_grid, const vector<double> &fp_grid,
                     const field_type &h1_grid, const field_type &dh1_grid, const field_type &ddh1_grid,
                     const field_type &h2_grid, const field_type &dh2_grid, const field_type &ddh2_grid,
                     int l3, int m3, int l1, int m1, int l2, int m2) {
  vector<complex<double>> ricci(r_grid.size());
  const complex<double> I(0.0, 1.0);
  const double Omega = sqrt(M/(r0*r0*r0));
  const double omega1 = m1*Omega;
  const double omega2 = m2*Omega;

  const double cv1 = Coupling(l3,m3,1,l2,m2,0,l1,m1,1);
  const double cv2 = Coupling(l3,m3,-1,l2,m2,-1,l1,m1,0);
  const double cv3 = Coupling(l3,m3,-1,l2,m2,-2,l1,m1,1);
  const double cv4 = Coupling(l3,m3,1,l2,m2,-1,l1,m1,2);
  const double cv5 = Coupling(l3,m3,1,l2,m2,-2,l1,m1,3);
  const double cv6 = Coupling(l3,m3,-1,l2,m2,-3,l1,m1,2);

  const double mu1l1 = sqrt(l1*(l1+1));
  const double mu1l2 = sqrt(l2*(l2+1));
  const double mu2l1 = sqrt((l1-1)*(l1+2));
  const double mu2l2 = sqrt((l2-1)*(l2+2));
  const double mu3l1 = sqrt((l1-2)*(l1+3));
  const double mu3l2 = sqrt((l2-2)*(l2+3));

  const double omega1_2 = omega1*omega1;
  const double omega2_2 = omega2*omega2;
  const double mu1l1_2 = mu1l1*mu1l1;
  const double mu1l2_2 = mu1l2*mu1l2;
  const double mu2l1_2 = mu2l1*mu2l1;
  const double mu2l2_2 = mu2l2*mu2l2;

  for(vector<complex<double>>::size_type j=0; j!=ricci.size(); ++j) {
    const double r = r_grid[j];
    const double f = f_grid[j];

    const double r_2 = r*r;
    const double r_3 = r_2*r;
    const double r_4 = r_2*r_2;
    const double r_5 = r_3*r_2;
    const double f_2 = f*f;
    const double f_3 = f_2*f;
    const double f_4 = f_2*f_2;

    const complex<double> h1_1  = h1_grid[1][l1][m1][j];
    const complex<double> h1_2  = h1_grid[2][l1][m1][j];
    const complex<double> h1_3  = h1_grid[3][l1][m1][j];
    const complex<double> h1_4  = h1_grid[4][l1][m1][j];
    const complex<double> h1_5  = h1_grid[5][l1][m1][j];
    const complex<double> h1_6  = h1_grid[6][l1][m1][j];
    const complex<double> h1_7  = h1_grid[7][l1][m1][j];
    const complex<double> h1_8  = h1_grid[8][l1][m1][j];
    const complex<double> h1_9  = h1_grid[9][l1][m1][j];
    const complex<double> h1_10 = h1_grid[10][l1][m1][j];

    const complex<double> h2_1  = h2_grid[1][l2][m2][j];
    const complex<double> h2_2  = h2_grid[2][l2][m2][j];
    const complex<double> h2_3  = h2_grid[3][l2][m2][j];
    const complex<double> h2_4  = h2_grid[4][l2][m2][j];
    const complex<double> h2_5  = h2_grid[5][l2][m2][j];
    const complex<double> h2_6  = h2_grid[6][l2][m2][j];
    const complex<double> h2_7  = h2_grid[7][l2][m2][j];
    const complex<double> h2_8  = h2_grid[8][l2][m2][j];
    const complex<double> h2_9  = h2_grid[9][l2][m2][j];
    const complex<double> h2_10 = h2_grid[10][l2][m2][j];

    const complex<double> dh1_1  = dh1_grid[1][l1][m1][j];
    const complex<double> dh1_2  = dh1_grid[2][l1][m1][j];
    const complex<double> dh1_3  = dh1_grid[3][l1][m1][j];
    const complex<double> dh1_4  = dh1_grid[4][l1][m1][j];
    const complex<double> dh1_5  = dh1_grid[5][l1][m1][j];
    const complex<double> dh1_6  = dh1_grid[6][l1][m1][j];
    const complex<double> dh1_7  = dh1_grid[7][l1][m1][j];
    const complex<double> dh1_8  = dh1_grid[8][l1][m1][j];
    const complex<double> dh1_9  = dh1_grid[9][l1][m1][j];
    const complex<double> dh1_10 = dh1_grid[10][l1][m1][j];

    const complex<double> dh2_1  = dh2_grid[1][l2][m2][j];
    const complex<double> dh2_2  = dh2_grid[2][l2][m2][j];
    const complex<double> dh2_3  = dh2_grid[3][l2][m2][j];
    const complex<double> dh2_4  = dh2_grid[4][l2][m2][j];
    const complex<double> dh2_5  = dh2_grid[5][l2][m2][j];
    const complex<double> dh2_6  = dh2_grid[6][l2][m2][j];
    const complex<double> dh2_7  = dh2_grid[7][l2][m2][j];
    const complex<double> dh2_8  = dh2_grid[8][l2][m2][j];
    const complex<double> dh2_9  = dh2_grid[9][l2][m2][j];
    const complex<double> dh2_10 = dh2_grid[10][l2][m2][j];

    const complex<double> ddh1_1  = ddh1_grid[1][l1][m1][j];
    const complex<double> ddh1_2  = ddh1_grid[2][l1][m1][j];
    const complex<double> ddh1_3  = ddh1_grid[3][l1][m1][j];
    const complex<double> ddh1_4  = ddh1_grid[4][l1][m1][j];
    const complex<double> ddh1_5  = ddh1_grid[5][l1][m1][j];
    const complex<double> ddh1_6  = ddh1_grid[6][l1][m1][j];
    const complex<double> ddh1_7  = ddh1_grid[7][l1][m1][j];
    const complex<double> ddh1_8  = ddh1_grid[8][l1][m1][j];
    const complex<double> ddh1_9  = ddh1_grid[9][l1][m1][j];
    const complex<double> ddh1_10 = ddh1_grid[10][l1][m1][j];

    const complex<double> ddh2_1  = ddh2_grid[1][l2][m2][j];
    const complex<double> ddh2_2  = ddh2_grid[2][l2][m2][j];
    const complex<double> ddh2_3  = ddh2_grid[3][l2][m2][j];
    const complex<double> ddh2_4  = ddh2_grid[4][l2][m2][j];
    const complex<double> ddh2_5  = ddh2_grid[5][l2][m2][j];
    const complex<double> ddh2_6  = ddh2_grid[6][l2][m2][j];
    const complex<double> ddh2_7  = ddh2_grid[7][l2][m2][j];
    const complex<double> ddh2_8  = ddh2_grid[8][l2][m2][j];
    const complex<double> ddh2_9  = ddh2_grid[9][l2][m2][j];
    const complex<double> ddh2_10 = ddh2_grid[10][l2][m2][j];

    ricci[j] = 0;
    if(l1>=1) ricci[j] += (0.17677669529663687*cv1*(-2.*(omega1*r*(h1_5 - I*h1_9)*((omega1 - omega2)*r_2*h2_1 + 2.*I*M*h2_2) - mu1l1_2*h1_2*(I*omega1*r_2*h2_1 + 2.*M*h2_2) + mu1l1_2*h1_1*(2.*M*h2_1 + I*omega1*r_2*h2_2)) + 2.*f_3*r*(r_2*dh1_5*dh2_6 - I*r_2*dh1_9*dh2_6 + mu1l1_2*r*dh2_3*h1_3 + mu1l1_2*r*dh2_6*h1_3 + mu1l1_2*r*dh2_6*h1_6 + I*r_2*ddh2_6*h1_9 - 2.*I*r*dh2_6*h1_9 + (mu1l1_2*(r*dh1_3 - 2.*h1_3) + 2.*h1_5 - 2.*I*h1_9)*h2_3 - r*dh1_5*h2_6 + mu1l1_2*r*dh1_6*h2_6 + r*dh1_7*h2_6 + I*r*dh1_9*h2_6 - I*r*dh1_10*h2_6 + mu1l1_2*h1_3*h2_6 - 2.*mu1l1_2*h1_6*h2_6 - h1_7*h2_6 - 2.*I*h1_9*h2_6 + I*h1_10*h2_6 + h1_5*(r*(-(r*ddh2_6) + 2.*dh2_6) + 2.*h2_6)) + f_2*(-2.*r_3*dh1_4*dh2_2 + 2.*I*r_3*dh1_8*dh2_2 - 2.*mu1l1_2*r_2*dh2_6*h1_1 + 4.*r_2*dh2_2*h1_4 + r_2*dh2_4*h1_4 - 2.*I*omega2*r_3*dh2_6*h1_4 + 3.*I*r_2*dh2_8*h1_4 - 2.*r_2*dh2_5*h1_5 - 4.*M*r*dh2_6*h1_5 - 4.*I*r_2*dh2_2*h1_8 - I*r_2*dh2_4*h1_8 - 2.*omega2*r_3*dh2_6*h1_8 + 3.*r_2*dh2_8*h1_8 + 2.*I*r_2*dh2_5*h1_9 + 4.*I*M*r*dh2_6*h1_9 + 2.*mu1l1_2*r_2*dh1_3*h2_1 - 6.*mu1l1_2*r*h1_3*h2_1 - 4.*r*h1_5*h2_1 + 4.*I*r*h1_9*h2_1 + 2.*r_3*ddh1_4*h2_2 - 2.*I*r_3*ddh1_8*h2_2 - 2.*r_2*dh1_4*h2_2 + 2.*I*r_2*dh1_8*h2_2 + (-((8.*M + mu1l2_2*r)*h1_5) + 2.*r*(-3.*mu1l1_2*h1_1 + r*(mu1l1_2*dh1_1 + I*omega1*r*dh1_4 + I*omega2*r*dh1_4 + omega1*r*dh1_8 + omega2*r*dh1_8 - 2.*I*(omega1 + omega2)*h1_4 - 2.*(omega1 + omega2)*h1_8)) + I*(8.*M + mu1l2_2*r)*h1_9)*h2_3 + r_2*dh1_4*h2_4 - I*r_2*dh1_8*h2_4 - 2.*r*h1_4*h2_4 + 2.*I*r*h1_8*h2_4 + mu1l1_2*r*h1_3*h2_5 + 2.*r*h1_5*h2_5 - 2.*I*r*h1_9*h2_5 - 2.*mu1l1_2*r*h1_1*h2_6 + 2.*I*omega2*r_2*h1_4*h2_6 + 4.*M*h1_5*h2_6 + 2.*omega2*r_2*h1_8*h2_6 - 4.*I*M*h1_9*h2_6 + 2.*I*mu1l1_2*r*h1_9*h2_6 + I*r_2*dh1_4*h2_8 + r_2*dh1_8*h2_8 - 6.*I*r*h1_4*h2_8 - 6.*r*h1_8*h2_8 - 2.*I*r_2*dh1_5*h2_9 - 2.*r_2*dh1_9*h2_9 - I*mu1l1_2*r*h1_3*h2_9 - 2.*I*r*h1_5*h2_9 + 2.*I*mu1l1_2*r*h1_6*h2_9 + 2.*I*r*h1_7*h2_9 - 2.*r*h1_9*h2_9 + 2.*r*h1_10*h2_9) + f*(2.*mu1l1_2*r_2*dh2_1*h1_1 - 2.*mu1l1_2*r_2*dh2_2*h1_2 - 2.*omega1*r_3*dh2_2*h1_9 + 2.*mu1l1_2*r_2*dh1_1*h2_1 + 2.*I*omega1*r_3*dh1_4*h2_1 - 2.*I*omega2*r_3*dh1_4*h2_1 + 2.*omega1*r_3*dh1_8*h2_1 - 2.*omega2*r_3*dh1_8*h2_1 - 4.*mu1l1_2*r*h1_1*h2_1 - 4.*omega1*r_2*h1_8*h2_1 + 4.*omega2*r_2*h1_8*h2_1 - I*mu1l2_2*r*h1_9*h2_1 - 2.*mu1l1_2*r_2*dh1_2*h2_2 + 2.*I*omega1*r_3*dh1_5*h2_2 + 2.*omega1*r_3*dh1_9*h2_2 + 4.*mu1l1_2*r*h1_2*h2_2 + 2.*I*mu1l1_2*omega1*r_2*h1_3*h2_2 + I*mu1l2_2*r*h1_8*h2_2 + 2.*omega1*r_2*h1_9*h2_2 + 2.*I*mu1l1_2*omega1*r_2*h1_2*h2_3 + 2.*I*omega1_2*r_3*h1_9*h2_3 + 2.*I*omega1*omega2*r_3*h1_9*h2_3 + mu1l1_2*r*h1_2*h2_4 + 2.*I*M*h1_8*h2_4 + omega1*r_2*h1_9*h2_4 - mu1l1_2*r*h1_1*h2_5 - omega2*r_2*h1_8*h2_5 - 2.*I*M*h1_9*h2_5 - 2.*I*mu1l1_2*omega2*r_2*h1_2*h2_6 - 2.*I*omega1*omega2*r_3*h1_9*h2_6 - I*mu1l1_2*r*h1_2*h2_8 - 2.*M*h1_8*h2_8 - I*omega1*r_2*h1_9*h2_8 + I*mu1l1_2*r*h1_1*h2_9 + 3.*I*omega2*r_2*h1_8*h2_9 + 2.*M*h1_9*h2_9 + h1_5*(2.*M*h2_5 + r*(mu1l2_2*h2_1 + omega1*r*(-2.*I*r*dh2_2 + 2.*I*h2_2 - 2.*omega1*r*h2_3 - 2.*omega2*r*h2_3 + I*h2_4 + 2.*omega2*r*h2_6 + h2_8)) + 2.*I*M*h2_9) - h1_4*(2.*M*h2_4 + 2.*I*M*h2_8 + r*(4.*I*(omega1 - omega2)*r*h2_1 + mu1l2_2*h2_2 + omega2*r*(I*h2_5 + 3.*h2_9))))))/(f_2*mu1l1*r_5);
    if(l2>=1) ricci[j] += (0.17677669529663687*cv2*(h1_1*(-4.*M*mu1l2_2*h2_1 + 2.*omega2*r_2*(I*mu1l2_2*h2_2 + (omega1 - omega2)*r*(h2_5 + I*h2_9))) + 2.*h1_2*(2.*M*mu1l2_2*h2_2 + omega2*r*(-1.*I*mu1l2_2*r*h2_1 - 2.*I*M*h2_5 + 2.*M*h2_9)) - f*(-2.*mu1l2_2*r_2*dh2_1*h1_1 + 2.*I*omega1*r_3*dh2_4*h1_1 - 2.*I*omega2*r_3*dh2_4*h1_1 - 2.*omega1*r_3*dh2_8*h1_1 + 2.*omega2*r_3*dh2_8*h1_1 + 2.*mu1l2_2*r_2*dh2_2*h1_2 - 2.*I*omega2*r_3*dh2_5*h1_2 + 2.*omega2*r_3*dh2_9*h1_2 - 2.*mu1l2_2*r_2*dh1_1*h2_1 + 4.*mu1l2_2*r*h1_1*h2_1 + I*mu1l2_2*r*h1_9*h2_1 + 2.*mu1l2_2*r_2*dh1_2*h2_2 - 4.*mu1l2_2*r*h1_2*h2_2 - 2.*I*mu1l2_2*omega2*r_2*h1_3*h2_2 + 2.*I*mu1l2_2*omega1*r_2*h1_6*h2_2 - I*mu1l2_2*r*h1_8*h2_2 - 2.*I*mu1l2_2*omega2*r_2*h1_2*h2_3 - 4.*I*omega1*r_2*h1_1*h2_4 + 4.*I*omega2*r_2*h1_1*h2_4 + mu1l1_2*r*h1_2*h2_4 - 2.*I*M*h1_8*h2_4 - 3.*omega1*r_2*h1_9*h2_4 + 2.*I*omega2*r_3*dh1_2*h2_5 - mu1l1_2*r*h1_1*h2_5 - 2.*I*omega2*r_2*h1_2*h2_5 + 2.*omega1*omega2*r_3*h1_3*h2_5 + 2.*omega2_2*r_3*h1_3*h2_5 - 2.*omega1*omega2*r_3*h1_6*h2_5 + omega2*r_2*h1_8*h2_5 + 2.*I*M*h1_9*h2_5 + 4.*omega1*r_2*h1_1*h2_8 - 4.*omega2*r_2*h1_1*h2_8 + I*mu1l1_2*r*h1_2*h2_8 + 2.*M*h1_8*h2_8 - 3.*I*omega1*r_2*h1_9*h2_8 - 2.*omega2*r_3*dh1_2*h2_9 - I*mu1l1_2*r*h1_1*h2_9 + 2.*omega2*r_2*h1_2*h2_9 + 2.*I*omega1*omega2*r_3*h1_3*h2_9 + 2.*I*omega2_2*r_3*h1_3*h2_9 - 2.*I*omega1*omega2*r_3*h1_6*h2_9 + I*omega2*r_2*h1_8*h2_9 - 2.*M*h1_9*h2_9 + h1_5*(-2.*M*h2_5 + r*(mu1l2_2*h2_1 + I*omega1*r*h2_4 - omega1*r*h2_8) - 2.*I*M*h2_9) + h1_4*(2.*M*h2_4 + 2.*I*M*h2_8 + r*(-(mu1l2_2*h2_2) + omega2*r*(-1.*I*h2_5 + h2_9)))) + 2.*f_3*r*(r_2*dh1_6*dh2_5 + I*r_2*dh1_6*dh2_9 - r*dh2_5*h1_6 + mu1l2_2*r*dh2_6*h1_6 + r*dh2_7*h1_6 - I*r*dh2_9*h1_6 + I*r*dh2_10*h1_6 + mu1l2_2*r*dh1_3*h2_3 + mu1l2_2*r*dh1_6*h2_3 + mu1l2_2*h1_6*h2_3 + (r*(-(r*ddh1_6) + 2.*dh1_6) + h1_6)*h2_5 + mu1l2_2*r*dh1_6*h2_6 - 2.*mu1l2_2*h1_6*h2_6 - h1_6*h2_7 + h1_3*(mu1l2_2*(r*dh2_3 - 2.*h2_3) + 2.*h2_5 + 2.*I*h2_9) - I*r_2*ddh1_6*h2_9 + 2.*I*r*dh1_6*h2_9 + I*h1_6*h2_9 - I*h1_6*h2_10) - f_2*(h1_3*((8.*M + mu1l1_2*r)*h2_5 + 2.*r*(3.*mu1l2_2*h2_1 + r*(-(mu1l2_2*dh2_1) - I*omega1*r*dh2_4 - I*omega2*r*dh2_4 + omega1*r*dh2_8 + omega2*r*dh2_8 + 2.*I*(omega1 + omega2)*h2_4 - 2.*(omega1 + omega2)*h2_8)) + I*(8.*M + mu1l1_2*r)*h2_9) + r*(2.*r_2*dh1_2*dh2_4 + 2.*I*r_2*dh1_2*dh2_8 - 2.*r_2*ddh2_4*h1_2 - 2.*I*r_2*ddh2_8*h1_2 + 2.*r*dh2_4*h1_2 + 2.*I*r*dh2_8*h1_2 + I*r*dh2_4*h1_8 - r*dh2_8*h1_8 - 2.*I*r*dh2_5*h1_9 + 2.*r*dh2_9*h1_9 + 2.*mu1l2_2*r*dh1_6*h2_1 + 2.*mu1l2_2*h1_6*h2_1 - mu1l2_2*h1_5*h2_3 - I*mu1l2_2*h1_9*h2_3 - 4.*r*dh1_2*h2_4 - r*dh1_4*h2_4 + 2.*I*omega1*r_2*dh1_6*h2_4 + 3.*I*r*dh1_8*h2_4 - 2.*I*omega1*r*h1_6*h2_4 - 6.*I*h1_8*h2_4 + 2.*r*dh1_5*h2_5 + 4.*M*dh1_6*h2_5 - 2.*h1_5*h2_5 - mu1l2_2*h1_6*h2_5 + mu2l2_2*h1_6*h2_5 - 2.*I*h1_9*h2_5 + 2.*I*mu1l2_2*h1_9*h2_6 + 2.*I*h1_9*h2_7 + h1_4*(-(r*(dh2_4 + I*dh2_8)) + 2.*h2_4 + 2.*I*h2_8) - 4.*I*r*dh1_2*h2_8 - I*r*dh1_4*h2_8 - 2.*omega1*r_2*dh1_6*h2_8 - 3.*r*dh1_8*h2_8 + 2.*omega1*r*h1_6*h2_8 + 6.*h1_8*h2_8 + h1_1*(2.*mu1l2_2*(-(r*dh2_3) + 3.*h2_3) + 4.*h2_5 + 4.*I*h2_9) + 2.*I*r*dh1_5*h2_9 + 4.*I*M*dh1_6*h2_9 - 2.*I*h1_5*h2_9 + I*mu1l2_2*h1_6*h2_9 + I*mu2l2_2*h1_6*h2_9 + 2.*h1_9*h2_9 - 2.*h1_9*h2_10))))/(f_2*mu1l2*r_5);
    if(l1>=1 && l2>=2) ricci[j] += (0.17677669529663687*cv3*(mu1l2_2*mu2l2_2*r*h1_5*h2_1 - I*mu1l2_2*mu2l2_2*r*h1_9*h2_1 + I*mu1l2_2*mu2l2_2*r*h1_8*h2_2 + mu1l1_2*mu2l2_2*r*h1_2*h2_4 + I*mu2l2_2*omega1*r_2*h1_5*h2_4 + 2.*I*M*mu2l2_2*h1_8*h2_4 + mu2l2_2*omega1*r_2*h1_9*h2_4 - mu1l1_2*mu2l2_2*r*h1_1*h2_5 + 2.*M*mu2l2_2*h1_5*h2_5 - mu2l2_2*omega2*r_2*h1_8*h2_5 - 2.*I*M*mu2l2_2*h1_9*h2_5 + 2.*I*mu1l1_2*omega2*r_2*h1_2*h2_7 - 2.*omega1*omega2*r_3*h1_5*h2_7 + 2.*I*omega1*omega2*r_3*h1_9*h2_7 + I*mu1l1_2*mu2l2_2*r*h1_2*h2_8 - mu2l2_2*omega1*r_2*h1_5*h2_8 - 2.*M*mu2l2_2*h1_8*h2_8 + I*mu2l2_2*omega1*r_2*h1_9*h2_8 - I*mu1l1_2*mu2l2_2*r*h1_1*h2_9 + 2.*I*M*mu2l2_2*h1_5*h2_9 - I*mu2l2_2*omega2*r_2*h1_8*h2_9 + 2.*M*mu2l2_2*h1_9*h2_9 - mu2l2_2*h1_4*(2.*M*h2_4 + 2.*I*M*h2_8 + r*(mu1l2_2*h2_2 + I*omega2*r*h2_5 - omega2*r*h2_9)) - 2.*mu1l1_2*omega2*r_2*h1_2*h2_10 - 2.*I*omega1*omega2*r_3*h1_5*h2_10 - 2.*omega1*omega2*r_3*h1_9*h2_10 + f_2*r*(-2.*r_2*dh1_5*dh2_7 + 2.*I*r_2*dh1_9*dh2_7 - 2.*I*r_2*dh1_5*dh2_10 - 2.*r_2*dh1_9*dh2_10 - 2.*mu1l1_2*r*dh2_7*h1_3 - 2.*I*mu1l1_2*r*dh2_10*h1_3 - 2.*I*r_2*ddh2_7*h1_9 + 2.*r_2*ddh2_10*h1_9 + 4.*I*r*dh2_7*h1_9 - 4.*r*dh2_10*h1_9 - I*r*dh2_7*h1_10 + r*dh2_10*h1_10 + 2.*r*dh1_5*h2_7 + 2.*mu1l1_2*r*dh1_6*h2_7 + 2.*r*dh1_7*h2_7 - 2.*I*r*dh1_9*h2_7 - 2.*I*r*dh1_10*h2_7 + 2.*mu1l1_2*h1_3*h2_7 - 2.*mu1l1_2*h1_6*h2_7 - 4.*I*h1_9*h2_7 + 3.*I*h1_10*h2_7 + 2.*h1_5*(r*(r*ddh2_7 + I*r*ddh2_10 - 2.*dh2_7 - 2.*I*dh2_10) + 2.*h2_7 + 2.*I*h2_10) + h1_7*(r*(dh2_7 + I*dh2_10) - 3.*h2_7 - 3.*I*h2_10) + 2.*I*r*dh1_5*h2_10 + 2.*I*mu1l1_2*r*dh1_6*h2_10 + 2.*I*r*dh1_7*h2_10 + 2.*r*dh1_9*h2_10 + 2.*r*dh1_10*h2_10 + 2.*I*mu1l1_2*h1_3*h2_10 - 2.*I*mu1l1_2*h1_6*h2_10 + 4.*h1_9*h2_10 - 3.*h1_10*h2_10) + f*(2.*mu1l1_2*r_2*dh2_7*h1_1 + 2.*I*mu1l1_2*r_2*dh2_10*h1_1 + mu2l2_2*r_2*dh2_4*h1_4 + 2.*I*omega2*r_3*dh2_7*h1_4 + I*mu2l2_2*r_2*dh2_8*h1_4 - 2.*omega2*r_3*dh2_10*h1_4 - I*mu2l2_2*r_2*dh2_4*h1_8 + 2.*omega2*r_3*dh2_7*h1_8 + mu2l2_2*r_2*dh2_8*h1_8 + 2.*I*omega2*r_3*dh2_10*h1_8 + 2.*I*mu2l2_2*r_2*dh2_5*h1_9 - 4.*I*M*r*dh2_7*h1_9 - 2.*mu2l2_2*r_2*dh2_9*h1_9 + 4.*M*r*dh2_10*h1_9 + I*mu1l2_2*mu2l2_2*r*h1_9*h2_3 + mu2l2_2*r_2*dh1_4*h2_4 - I*mu2l2_2*r_2*dh1_8*h2_4 - 2.*mu2l2_2*r*h1_4*h2_4 + 2.*I*mu2l2_2*r*h1_8*h2_4 + mu1l1_2*mu2l2_2*r*h1_3*h2_5 - 2.*I*mu2l2_2*r*h1_9*h2_5 + 2.*I*(I*r*(mu1l1_2*h1_1 + omega2*r*(I*h1_4 + h1_8)) + (2.*M + mu1l1_2*r)*h1_9)*h2_7 + I*mu2l2_2*r_2*dh1_4*h2_8 + mu2l2_2*r_2*dh1_8*h2_8 - 2.*I*mu2l2_2*r*h1_4*h2_8 - 2.*mu2l2_2*r*h1_8*h2_8 + I*mu1l1_2*mu2l2_2*r*h1_3*h2_9 + 2.*mu2l2_2*r*h1_9*h2_9 - 2.*I*mu1l1_2*r*h1_1*h2_10 + 2.*omega2*r_2*h1_4*h2_10 - 2.*I*omega2*r_2*h1_8*h2_10 - 4.*M*h1_9*h2_10 - 2.*mu1l1_2*r*h1_9*h2_10 - h1_5*(4.*M*h2_7 + r*(2.*mu2l2_2*r*dh2_5 - 4.*M*dh2_7 + 2.*I*mu2l2_2*r*dh2_9 - 4.*I*M*dh2_10 + mu1l2_2*mu2l2_2*h2_3 - 2.*mu2l2_2*h2_5 - 2.*I*mu2l2_2*h2_9) + 4.*I*M*h2_10))))/(f*mu1l1*mu1l2*mu2l2*r_5);
    if(l1>=2 && l2>=1) ricci[j] += (0.17677669529663687*cv4*(-(mu1l2_2*mu2l1_2*r*h1_5*h2_1) + I*mu1l2_2*mu2l1_2*r*h1_9*h2_1 + 2.*I*mu1l2_2*omega1*r_2*h1_7*h2_2 - I*mu1l2_2*mu2l1_2*r*h1_8*h2_2 + 2.*mu1l2_2*omega1*r_2*h1_10*h2_2 - mu1l1_2*mu2l1_2*r*h1_2*h2_4 - I*mu2l1_2*omega1*r_2*h1_5*h2_4 + 2.*I*M*mu2l1_2*h1_8*h2_4 - mu2l1_2*omega1*r_2*h1_9*h2_4 + mu1l1_2*mu2l1_2*r*h1_1*h2_5 + 2.*M*mu2l1_2*h1_5*h2_5 - 2.*omega1*omega2*r_3*h1_7*h2_5 + mu2l1_2*omega2*r_2*h1_8*h2_5 - 2.*I*M*mu2l1_2*h1_9*h2_5 + 2.*I*omega1*omega2*r_3*h1_10*h2_5 - I*mu1l1_2*mu2l1_2*r*h1_2*h2_8 + mu2l1_2*omega1*r_2*h1_5*h2_8 - 2.*M*mu2l1_2*h1_8*h2_8 - I*mu2l1_2*omega1*r_2*h1_9*h2_8 + I*mu1l1_2*mu2l1_2*r*h1_1*h2_9 + 2.*I*M*mu2l1_2*h1_5*h2_9 - 2.*I*omega1*omega2*r_3*h1_7*h2_9 + I*mu2l1_2*omega2*r_2*h1_8*h2_9 + 2.*M*mu2l1_2*h1_9*h2_9 - 2.*omega1*omega2*r_3*h1_10*h2_9 + mu2l1_2*h1_4*(-2.*M*h2_4 - 2.*I*M*h2_8 + r*(mu1l2_2*h2_2 + I*omega2*r*h2_5 - omega2*r*h2_9)) - f*(((8.*M - mu1l2_2*r + mu2l2_2*r)*h1_7 + r*(2.*mu2l1_2*r*dh1_5 - 4.*M*dh1_7 - 2.*I*mu2l1_2*r*dh1_9 + 4.*I*M*dh1_10 + mu1l1_2*mu2l1_2*h1_3 - 2.*mu2l1_2*h1_5 + 2.*I*mu2l1_2*h1_9) - I*(8.*M - mu1l2_2*r + mu2l2_2*r)*h1_10)*h2_5 + I*r*(mu2l1_2*r*dh2_4*h1_8 + I*mu2l1_2*r*dh2_8*h1_8 + 2.*I*mu1l2_2*r*dh1_7*h2_1 + 2.*mu1l2_2*r*dh1_10*h2_1 - 2.*mu1l2_2*h1_10*h2_1 + I*mu1l2_2*mu2l1_2*h1_5*h2_3 + mu1l2_2*mu2l1_2*h1_9*h2_3 + (-2.*mu2l1_2*h1_8 + r*(mu2l1_2*(I*dh1_4 + dh1_8) - 2.*omega1*r*(dh1_7 - I*dh1_10) - 2.*I*omega1*h1_10))*h2_4 - mu2l1_2*r*dh1_4*h2_8 - 2.*I*omega1*r_2*dh1_7*h2_8 + I*mu2l1_2*r*dh1_8*h2_8 - 2.*omega1*r_2*dh1_10*h2_8 - 2.*I*mu2l1_2*h1_8*h2_8 + 2.*omega1*r*h1_10*h2_8 + mu2l1_2*h1_4*(I*r*dh2_4 - r*dh2_8 - 2.*I*h2_4 + 2.*h2_8)) + (r*(2.*I*mu2l1_2*r*dh1_5 - 4.*I*M*dh1_7 + 2.*mu2l1_2*r*dh1_9 - 4.*M*dh1_10 + I*mu1l1_2*mu2l1_2*h1_3 - 2.*I*mu2l1_2*h1_5 - 2.*mu2l1_2*h1_9) + (8.*M + mu1l2_2*r + mu2l2_2*r)*h1_10)*h2_9 + I*h1_7*(2.*r*(-1.*I*mu1l2_2*h2_1 + omega1*r*(h2_4 + I*h2_8)) + (8.*M + mu1l2_2*r + mu2l2_2*r)*h2_9)) + f_2*r*(-2.*r_2*dh1_7*dh2_5 + 2.*I*r_2*dh1_10*dh2_5 - 2.*I*r_2*dh1_7*dh2_9 - 2.*r_2*dh1_10*dh2_9 - 2.*I*r*dh2_5*h1_10 - 2.*I*mu1l2_2*r*dh2_6*h1_10 - 2.*I*r*dh2_7*h1_10 + 2.*r*dh2_9*h1_10 + 2.*r*dh2_10*h1_10 - 2.*mu1l2_2*r*dh1_7*h2_3 + 2.*I*mu1l2_2*r*dh1_10*h2_3 - 2.*I*mu1l2_2*h1_10*h2_3 + 2.*(r*(r*ddh1_7 - I*r*ddh1_10 - 2.*dh1_7 + 2.*I*dh1_10) + h1_7 - I*h1_10)*h2_5 + 2.*I*mu1l2_2*h1_10*h2_6 + r*dh1_7*h2_7 - I*r*dh1_10*h2_7 + 3.*I*h1_10*h2_7 + 2.*I*r_2*ddh1_7*h2_9 + 2.*r_2*ddh1_10*h2_9 - 4.*I*r*dh1_7*h2_9 - 4.*r*dh1_10*h2_9 + 2.*h1_10*h2_9 + h1_7*(2.*r*dh2_5 + 2.*mu1l2_2*r*dh2_6 + 2.*r*dh2_7 + 2.*I*r*dh2_9 + 2.*I*r*dh2_10 + 2.*mu1l2_2*h2_3 - 2.*mu1l2_2*h2_6 - 3.*h2_7 + 2.*I*h2_9 - 3.*I*h2_10) + I*r*dh1_7*h2_10 + r*dh1_10*h2_10 - 3.*h1_10*h2_10)))/(f*mu1l1*mu1l2*mu2l1*r_5);
    if(l1>=2 && l2>=2) ricci[j] += (0.17677669529663687*cv5*f*mu3l1*(h1_7 - I*h1_10)*(r*(dh2_7 + I*dh2_10) - h2_7 - I*h2_10))/(mu1l1*mu1l2*mu2l1*mu2l2*r_4);
    if(l1>=2 && l2>=2) ricci[j] += (0.17677669529663687*cv6*f*mu3l2*(r*(dh1_7 - I*dh1_10) - h1_7 + I*h1_10)*(h2_7 + I*h2_10))/(mu1l1*mu1l2*mu2l1*mu2l2*r_4);
  }

  return ricci;
}
