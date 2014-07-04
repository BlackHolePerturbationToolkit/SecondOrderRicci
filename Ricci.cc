/*
 * Copyright 2014 Barry Wardell
 *
 * This file is distributed under the University of Illinois/NCSA
 * Open Source License. See LICENSE file for details.
 *
 */

#include <vector>
#include <complex>
#include <boost/multi_array.hpp>
#include "h1.h"
#include "R2.h"

using namespace std;
using boost::multi_array;

int main(int argc, char* argv[])
{
  /* Read in first-order fields */
  vector<double> r, f, fp;
  multi_array<complex<double>, 4> h, dh, ddh;
  read_h1(r, f, fp, h, dh, ddh);
  vector<complex<double>> r2_1 = R2_1(r, f, fp, h, dh, 2, 2, 2, 0, 2, 2);
  for(auto r2: r2_1)
    cout << r2 << endl;

  return 0;
}
