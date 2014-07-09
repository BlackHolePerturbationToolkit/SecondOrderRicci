/*
 * Copyright 2014 Barry Wardell
 *
 * This file is distributed under the University of Illinois/NCSA
 * Open Source License. See LICENSE file for details.
 *
 */

#include <boost/multi_array.hpp>

void read_h1(const std::string dir, double &r0,
             std::vector<double> &r, std::vector<double> &f, std::vector<double> &fp,
             boost::multi_array<std::complex<double>,4> &h,
             boost::multi_array<std::complex<double>,4> &dh,
             boost::multi_array<std::complex<double>,4> &ddh);

void compute_h1P(const double r0, const std::vector<double> &r, const int l_max,
             boost::multi_array<std::complex<double>,4> &hP,
             boost::multi_array<std::complex<double>,4> &dhP,
             boost::multi_array<std::complex<double>,4> &ddhP);
