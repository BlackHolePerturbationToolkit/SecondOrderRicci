/*
 * Copyright 2014 Barry Wardell
 *
 * This file is distributed under the University of Illinois/NCSA
 * Open Source License. See LICENSE file for details.
 *
 */

#include "utils.h"
#include <math.h>
extern "C"
{
#include <gsl/gsl_sf_coupling.h>
}

double Coupling(int lc, int mc, int sc, int la, int ma, int sa, int lb, int mb, int sb) {
  double sign = isOdd(mc+sc) ? -1.0 : 1.0;
  return sign*sqrt((2.*la+1.)*(2.*lb+1.)*(2.*lc+1.)/(4.0*M_PI))*
    gsl_sf_coupling_3j(2*lc, 2*la, 2*lb, 2*sc, -2*sa, -2*sb)*
    gsl_sf_coupling_3j(2*lc, 2*la, 2*lb, -2*mc, 2*ma, 2*mb);
}