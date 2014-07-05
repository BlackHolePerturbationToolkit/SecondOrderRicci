/*
 * Copyright 2014 Barry Wardell
 *
 * This file is distributed under the University of Illinois/NCSA
 * Open Source License. See LICENSE file for details.
 *
 */

typedef boost::multi_array<std::complex<double>,4> field_type;

std::vector<std::complex<double>> R2_1(const std::vector<double> &r, const std::vector<double> &f, const std::vector<double> &fp,
  const field_type &h, const field_type &dh, const field_type &ddh, int l3, int m3, int l1, int m1, int l2, int m2);
std::vector<std::complex<double>> R2_2(const std::vector<double> &r, const std::vector<double> &f, const std::vector<double> &fp,
  const field_type &h, const field_type &dh, const field_type &ddh, int l3, int m3, int l1, int m1, int l2, int m2);
std::vector<std::complex<double>> R2_3(const std::vector<double> &r, const std::vector<double> &f, const std::vector<double> &fp,
  const field_type &h, const field_type &dh, const field_type &ddh, int l3, int m3, int l1, int m1, int l2, int m2);
std::vector<std::complex<double>> R2_4(const std::vector<double> &r, const std::vector<double> &f, const std::vector<double> &fp,
  const field_type &h, const field_type &dh, const field_type &ddh, int l3, int m3, int l1, int m1, int l2, int m2);
std::vector<std::complex<double>> R2_5(const std::vector<double> &r, const std::vector<double> &f, const std::vector<double> &fp,
  const field_type &h, const field_type &dh, const field_type &ddh, int l3, int m3, int l1, int m1, int l2, int m2);
std::vector<std::complex<double>> R2_6(const std::vector<double> &r, const std::vector<double> &f, const std::vector<double> &fp,
  const field_type &h, const field_type &dh, const field_type &ddh, int l3, int m3, int l1, int m1, int l2, int m2);
std::vector<std::complex<double>> R2_7(const std::vector<double> &r, const std::vector<double> &f, const std::vector<double> &fp,
  const field_type &h, const field_type &dh, const field_type &ddh, int l3, int m3, int l1, int m1, int l2, int m2);
std::vector<std::complex<double>> R2_8(const std::vector<double> &r, const std::vector<double> &f, const std::vector<double> &fp,
  const field_type &h, const field_type &dh, const field_type &ddh, int l3, int m3, int l1, int m1, int l2, int m2);
std::vector<std::complex<double>> R2_9(const std::vector<double> &r, const std::vector<double> &f, const std::vector<double> &fp,
  const field_type &h, const field_type &dh, const field_type &ddh, int l3, int m3, int l1, int m1, int l2, int m2);
std::vector<std::complex<double>> R2_10(const std::vector<double> &r, const std::vector<double> &f, const std::vector<double> &fp,
  const field_type &h, const field_type &dh, const field_type &ddh, int l3, int m3, int l1, int m1, int l2, int m2);
