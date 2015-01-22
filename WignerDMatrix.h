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
  std::array<double, 5151> D0;
  std::array<double, 5150> D1;
  std::array<double, 5148> D2;
  std::array<double, 5145> D3;
  std::array<double, 5141> D4;
};

#endif // WIGNERDMATRIX_H
