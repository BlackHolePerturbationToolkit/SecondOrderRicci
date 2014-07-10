/*
 * Copyright 2014 Barry Wardell
 *
 * This file is distributed under the University of Illinois/NCSA
 * Open Source License. See LICENSE file for details.
 *
 */

typedef boost::multi_array<std::complex<double>,4> field_type;

std::vector<std::complex<double>> hh_1(const double M, const double r0,
  const std::vector<double> &r, const std::vector<double> &f, const std::vector<double> &fp,
  const field_type &h1, const field_type &h2, int l3, int m3, int l1, int m1, int l2, int m2);
std::vector<std::complex<double>> hh_2(const double M, const double r0,
  const std::vector<double> &r, const std::vector<double> &f, const std::vector<double> &fp,
  const field_type &h1, const field_type &h2, int l3, int m3, int l1, int m1, int l2, int m2);
std::vector<std::complex<double>> hh_3(const double M, const double r0,
  const std::vector<double> &r, const std::vector<double> &f, const std::vector<double> &fp,
  const field_type &h1, const field_type &h2, int l3, int m3, int l1, int m1, int l2, int m2);
std::vector<std::complex<double>> hh_4(const double M, const double r0,
  const std::vector<double> &r, const std::vector<double> &f, const std::vector<double> &fp,
  const field_type &h1, const field_type &h2, int l3, int m3, int l1, int m1, int l2, int m2);
std::vector<std::complex<double>> hh_5(const double M, const double r0,
  const std::vector<double> &r, const std::vector<double> &f, const std::vector<double> &fp,
  const field_type &h1, const field_type &h2, int l3, int m3, int l1, int m1, int l2, int m2);
std::vector<std::complex<double>> hh_6(const double M, const double r0,
  const std::vector<double> &r, const std::vector<double> &f, const std::vector<double> &fp,
  const field_type &h1, const field_type &h2, int l3, int m3, int l1, int m1, int l2, int m2);
std::vector<std::complex<double>> hh_7(const double M, const double r0,
  const std::vector<double> &r, const std::vector<double> &f, const std::vector<double> &fp,
  const field_type &h1, const field_type &h2, int l3, int m3, int l1, int m1, int l2, int m2);
std::vector<std::complex<double>> hh_8(const double M, const double r0,
  const std::vector<double> &r, const std::vector<double> &f, const std::vector<double> &fp,
  const field_type &h1, const field_type &h2, int l3, int m3, int l1, int m1, int l2, int m2);
std::vector<std::complex<double>> hh_9(const double M, const double r0,
  const std::vector<double> &r, const std::vector<double> &f, const std::vector<double> &fp,
  const field_type &h1, const field_type &h2, int l3, int m3, int l1, int m1, int l2, int m2);
std::vector<std::complex<double>> hh_10(const double M, const double r0,
  const std::vector<double> &r, const std::vector<double> &f, const std::vector<double> &fp,
  const field_type &h1, const field_type &h2, int l3, int m3, int l1, int m1, int l2, int m2);
