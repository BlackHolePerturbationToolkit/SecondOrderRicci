/*
 * Copyright 2014 Barry Wardell
 *
 * This file is distributed under the University of Illinois/NCSA
 * Open Source License. See LICENSE file for details.
 *
 */

#ifndef WIGNERDMATRIX_H
#define WIGNERDMATRIX_H

#include <array>
#include <complex>

/** WignerDMatrix class.
 *  Compute the Wigner-D matrix element for a rotation by pi, pi/2, pi/2.
 */
class WignerDMatrix
{
public:
  /** Constructor.
   */
  WignerDMatrix();

  /** operator().
   *  Compute the matrix element for a specific l, m, m'.
   *  @param l Integer value for l
   *  @param m Integer value for m
   *  @param mp Integer value for m'
   */
  std::complex<double> operator()(int l, int m, int mp) const;

private:
  const double * const D0;
  const double * const D1;
  const double * const D2;
  const double * const D3;
  const double * const D4;
};

#endif // WIGNERDMATRIX_H
