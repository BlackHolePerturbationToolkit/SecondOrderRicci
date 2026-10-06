/*
 * One resolved (l, m) mode of a first-order field: pointers to the N grid values of each
 * component i = 1..10 (index 0 unused). Absent components point at zeros; for m < 0 the
 * pointers refer to a conjugated copy (h_{l,-m} = (-1)^m conj(h_{l m})) held in a ModeBuffer.
 */
#pragma once

#include <complex>
#include <cstddef>
#include <vector>
#include "field.h"

struct ModeBlock {
  const std::complex<double> *p[11];
  inline const std::complex<double> *operator[](int i) const { return p[i]; }
};

/* Per-thread scratch space for conjugated negative-m modes. */
class ModeBuffer {
public:
  void reserve(std::size_t N) { buf_.resize(11*N); N_ = N; }

  /* Resolve mode (l, m) of `field` into `block`. For m < 0 the values on jb <= j < je are
   * conjugated into this buffer (only that range is valid). */
  void resolve(const Field &field, int l, int m, std::size_t jb, std::size_t je, ModeBlock &block) {
    block.p[0] = nullptr;
    for (int i = 1; i <= 10; ++i) {
      Field::Mode md = field.mode(i, l, m);
      if (!md.negate_conj || md.p == field.zeros()) {
        block.p[i] = md.p;
      } else {
        std::complex<double> *dst = buf_.data() + i*N_;
        for (std::size_t j = jb; j < je; ++j)
          dst[j] = md.sign*std::conj(md.p[j]);
        block.p[i] = dst;
      }
    }
  }

private:
  std::vector<std::complex<double>> buf_;
  std::size_t N_ = 0;
};
