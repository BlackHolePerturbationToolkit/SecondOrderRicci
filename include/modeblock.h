/*
 * One resolved (l, m) mode of a first-order field: real and imaginary planes (N grid values
 * each) of every component i = 1..10 (index 0 unused). Absent components point at zeros;
 * for m < 0 the planes are a copy of (-1)^m conj(h_{l |m|}) held in a ModeBuffer.
 *
 * GBuffer holds the kernel outputs G_k (real and imaginary planes) for the couplings k.
 */
#pragma once

#include <complex>
#include <cstddef>
#include <vector>
#include <algorithm>
#include "field.h"

struct ModeBlock {
  const double *re[11];
  const double *im[11];
};

/* Per-thread scratch space for negative-m modes. */
class ModeBuffer {
public:
  void reserve(std::size_t N) { re_.resize(11*N); im_.resize(11*N); N_ = N; }

  /* Resolve mode (l, m) of `field` into `block`. For m < 0 the values on jb <= j < je are
   * copied (with the sign and conjugation) into this buffer (only that range is valid). */
  void resolve(const Field &field, int l, int m, std::size_t jb, std::size_t je, ModeBlock &block) {
    block.re[0] = block.im[0] = nullptr;
    for (int i = 1; i <= 10; ++i) {
      const Field::Mode md = field.mode(i, l, m);
      if (!md.negate_conj || md.re == field.zeros()) {
        block.re[i] = md.re;
        block.im[i] = md.im;
      } else {
        double *dr = re_.data() + i*N_, *di = im_.data() + i*N_;
        const double s = md.sign;
        for (std::size_t j = jb; j < je; ++j) {
          dr[j] = s*md.re[j];
          di[j] = -s*md.im[j];
        }
        block.re[i] = dr;
        block.im[i] = di;
      }
    }
  }

private:
  std::vector<double> re_, im_;
  std::size_t N_ = 0;
};

/* Kernel outputs G_k, k < K, as real and imaginary planes of N points. */
class GBuffer {
public:
  GBuffer(int Kmax, std::size_t N) : Kmax_(Kmax), N_(N), re_(Kmax*N), im_(Kmax*N), pre_(Kmax), pim_(Kmax) {
    for (int k = 0; k < Kmax; ++k) { pre_[k] = re_.data() + k*N; pim_[k] = im_.data() + k*N; }
  }
  double *const *re() { return pre_.data(); }
  double *const *im() { return pim_.data(); }
  void zero(int K, std::size_t jb, std::size_t je) {
    for (int k = 0; k < K; ++k) {
      std::fill(pre_[k] + jb, pre_[k] + je, 0.0);
      std::fill(pim_[k] + jb, pim_[k] + je, 0.0);
    }
  }
  /* sum_k c_k G_k at grid point j */
  inline std::complex<double> combine(int K, const double *c, std::size_t j) const {
    double vr = 0.0, vi = 0.0;
    for (int k = 0; k < K; ++k) { vr += c[k]*pre_[k][j]; vi += c[k]*pim_[k][j]; }
    return std::complex<double>(vr, vi);
  }

private:
  int Kmax_;
  std::size_t N_;
  std::vector<double> re_, im_;
  std::vector<double*> pre_, pim_;
};
