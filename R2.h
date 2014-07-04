/*
 * Copyright 2014 Barry Wardell
 *
 * This file is distributed under the University of Illinois/NCSA
 * Open Source License. See LICENSE file for details.
 *
 */

std::vector<std::complex<double>> R2_1(const std::vector<double> &r, const std::vector<double> &f, const std::vector<double> &fp,
  const boost::multi_array<std::complex<double>,4> &h, const boost::multi_array<std::complex<double>,4> &dh,
  int l3, int m3, int l1, int m1, int l2, int m2);
