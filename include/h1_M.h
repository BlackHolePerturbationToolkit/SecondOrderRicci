/*
 * Copyright 2020 Barry Wardell
 *
 * This file is distributed under the University of Illinois/NCSA
 * Open Source License. See LICENSE file for details.
 *
 */

#include <complex>
#include <vector>
#include "field.h"

void h1_M(const double &r0,
             const std::vector<double> &r, const std::vector<double> &f, const std::vector<double> &fp,
             Field &h,
             Field &dh,
             Field &ddh);
