/*
 * Copyright 2020 Barry Wardell
 *
 * This file is distributed under the University of Illinois/NCSA
 * Open Source License. See LICENSE file for details.
 *
 */

#include <boost/multi_array.hpp>

void h1_M(const double &r0,
             const std::vector<double> &r, const std::vector<double> &f, const std::vector<double> &fp,
             boost::multi_array<std::complex<double>,4> &h,
             boost::multi_array<std::complex<double>,4> &dh,
             boost::multi_array<std::complex<double>,4> &ddh);
