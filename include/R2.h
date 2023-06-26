/*
 * Copyright 2014 Barry Wardell
 *
 * This file is distributed under the University of Illinois/NCSA
 * Open Source License. See LICENSE file for details.
 *
 */

typedef boost::multi_array<std::complex<double>,4> field_type;

std::vector<std::complex<double>> R2_1(const double M, const double r0,
  const std::vector<double> &r, const std::vector<double> &f, const std::vector<double> &fp,
  const field_type &h1, const field_type &dh1, const field_type &ddh1,
  const field_type &h2, const field_type &dh2, const field_type &ddh2, int l3, int m3, int l1, int m1, int l2, int m2);
std::vector<std::complex<double>> R2_2(const double M, const double r0,
  const std::vector<double> &r, const std::vector<double> &f, const std::vector<double> &fp,
  const field_type &h1, const field_type &dh1, const field_type &ddh1,
  const field_type &h2, const field_type &dh2, const field_type &ddh2, int l3, int m3, int l1, int m1, int l2, int m2);
std::vector<std::complex<double>> R2_3(const double M, const double r0,
  const std::vector<double> &r, const std::vector<double> &f, const std::vector<double> &fp,
  const field_type &h1, const field_type &dh1, const field_type &ddh1,
  const field_type &h2, const field_type &dh2, const field_type &ddh2, int l3, int m3, int l1, int m1, int l2, int m2);
std::vector<std::complex<double>> R2_4(const double M, const double r0,
  const std::vector<double> &r, const std::vector<double> &f, const std::vector<double> &fp,
  const field_type &h1, const field_type &dh1, const field_type &ddh1,
  const field_type &h2, const field_type &dh2, const field_type &ddh2, int l3, int m3, int l1, int m1, int l2, int m2);
std::vector<std::complex<double>> R2_5(const double M, const double r0,
  const std::vector<double> &r, const std::vector<double> &f, const std::vector<double> &fp,
  const field_type &h1, const field_type &dh1, const field_type &ddh1,
  const field_type &h2, const field_type &dh2, const field_type &ddh2, int l3, int m3, int l1, int m1, int l2, int m2);
std::vector<std::complex<double>> R2_6(const double M, const double r0,
  const std::vector<double> &r, const std::vector<double> &f, const std::vector<double> &fp,
  const field_type &h1, const field_type &dh1, const field_type &ddh1,
  const field_type &h2, const field_type &dh2, const field_type &ddh2, int l3, int m3, int l1, int m1, int l2, int m2);
std::vector<std::complex<double>> R2_7(const double M, const double r0,
  const std::vector<double> &r, const std::vector<double> &f, const std::vector<double> &fp,
  const field_type &h1, const field_type &dh1, const field_type &ddh1,
  const field_type &h2, const field_type &dh2, const field_type &ddh2, int l3, int m3, int l1, int m1, int l2, int m2);
std::vector<std::complex<double>> R2_8(const double M, const double r0,
  const std::vector<double> &r, const std::vector<double> &f, const std::vector<double> &fp,
  const field_type &h1, const field_type &dh1, const field_type &ddh1,
  const field_type &h2, const field_type &dh2, const field_type &ddh2, int l3, int m3, int l1, int m1, int l2, int m2);
std::vector<std::complex<double>> R2_9(const double M, const double r0,
  const std::vector<double> &r, const std::vector<double> &f, const std::vector<double> &fp,
  const field_type &h1, const field_type &dh1, const field_type &ddh1,
  const field_type &h2, const field_type &dh2, const field_type &ddh2, int l3, int m3, int l1, int m1, int l2, int m2);
std::vector<std::complex<double>> R2_10(const double M, const double r0,
  const std::vector<double> &r, const std::vector<double> &f, const std::vector<double> &fp,
  const field_type &h1, const field_type &dh1, const field_type &ddh1,
  const field_type &h2, const field_type &dh2, const field_type &ddh2, int l3, int m3, int l1, int m1, int l2, int m2);
