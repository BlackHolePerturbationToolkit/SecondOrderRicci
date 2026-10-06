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
typedef Field field_type;

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

  /* Evaluate the perturbation on the grid (packed storage: only the l = 1, m = 0 mode) */
  h.init(l_max, N);
  dh.init(l_max, N);
  ddh.init(l_max, N);

  for(size_t j=0; j<N; ++j) {
    // i=8
    h.slot(8, l, m)[j] = -8*pow(M_PI/3.,0.5)*pow(r[j],-2);
    dh.slot(8, l, m)[j] = 16*pow(M_PI/3.,0.5)*pow(r[j],-3);
    ddh.slot(8, l, m)[j] = -48*pow(M_PI/3.,0.5)*pow(r[j],-4);

    // i=9
    h.slot(9, l, m)[j] = -16*pow(M_PI/3.,0.5)*pow(r[j],-3);
    dh.slot(9, l, m)[j] = 48*pow(M_PI/3.,0.5)*pow(r[j],-4);
    ddh.slot(9, l, m)[j] = -192*pow(M_PI/3.,0.5)*pow(r[j],-5);
  }

  /* The expressions above give the field in the convention of the original coupling code
   * (Barack-Lousto a_il and 1/r factors included). The kernels generated from
   * PerturbationEquations use the field without these factors: h -> h r / a_il. */
  for (int i = 1; i <= 10; ++i) {
    if (!Field::present(i, l, m)) continue;
    double a = pow(2, -0.5);
    if (i == 4 || i == 5 || i == 8 || i == 9) a *= pow(l*(l+1), -0.5);
    if (i == 7 || i == 10) a *= pow((l-1)*l*(l+1)*(l+2), -0.5);
    complex<double> *ph = h.slot(i, l, m), *pdh = dh.slot(i, l, m), *pddh = ddh.slot(i, l, m);
    for (size_t j = 0; j < N; ++j) {
      const complex<double> h0 = ph[j], d0 = pdh[j], dd0 = pddh[j];
      ph[j] = h0*r[j]/a;
      pdh[j] = (h0 + r[j]*d0)/a;
      pddh[j] = (2.0*d0 + r[j]*dd0)/a;
    }
  }
}
