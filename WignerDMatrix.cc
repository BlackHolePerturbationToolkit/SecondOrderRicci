/*
 * Copyright 2014 Barry Wardell
 *
 * This file is distributed under the University of Illinois/NCSA
 * Open Source License. See LICENSE file for details.
 *
 */

#include <complex>
#include <assert.h>

#include "h5wrapper.h"
#include "WignerDMatrix.h"
#include "utils.h"

WignerDMatrix::WignerDMatrix()
{
  H5F WignerD_h5("WignerD.h5", H5F_ACC_RDONLY);

  H5D d0(WignerD_h5, "m'=0");
  H5Dread(d0.getId(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, D0.data());

  H5D d1(WignerD_h5, "m'=1");
  H5Dread(d1.getId(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, D1.data());

  H5D d2(WignerD_h5, "m'=2");
  H5Dread(d2.getId(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, D2.data());

  H5D d3(WignerD_h5, "m'=3");
  H5Dread(d3.getId(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, D3.data());

  H5D d4(WignerD_h5, "m'=4");
  H5Dread(d4.getId(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, D4.data());
}

std::complex<double> WignerDMatrix::operator()(int l, int mp, int m) const
{
  const std::complex<double> I(0.0, 1.0);
  int am = abs(m);
  int amp = abs(mp);

  assert(am<=l && amp<=l && amp <= 4);

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

