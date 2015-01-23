/*
 * Copyright 2014 Barry Wardell
 *
 * This file is distributed under the University of Illinois/NCSA
 * Open Source License. See LICENSE file for details.
 *
 */

#include <complex>
#include <assert.h>

#include "WignerDMatrix.h"
#include "WignerD0.h"
#include "WignerD1.h"
#include "WignerD2.h"
#include "WignerD3.h"
#include "WignerD4.h"

static bool isOdd(int n)
{
  return n % 2;
}

WignerDMatrix::WignerDMatrix()
  : D0((const double * const)WignerD0_bin),
    D1((const double * const)WignerD1_bin),
    D2((const double * const)WignerD2_bin),
    D3((const double * const)WignerD3_bin),
    D4((const double * const)WignerD4_bin)
{}

std::complex<double> WignerDMatrix::operator()(int l, int mp, int m) const
{
  const std::complex<double> I(0.0, 1.0);
  int am = abs(m);
  int amp = abs(mp);

  assert(am<=l && amp<=l && amp <= 4 && l <= 100);

  int index = ((l + amp)*(l - amp + 1))/2 + am - amp;
  std::complex<double> val;

  switch(amp) {
    case 0:
      val = D0[index];
      break;

    case 1:
      val = I*D1[index];
      break;

    case 2:
      val = D2[index];
      break;

    case 3:
      val = I*D3[index];
      break;

    case 4:
      val = D4[index];
      break;
  }

  if(m < 0 && isOdd(l+mp))
    val *= -1;

  if(mp < 0 && isOdd(l + m + mp))
    val *= -1;

  return val;
}

