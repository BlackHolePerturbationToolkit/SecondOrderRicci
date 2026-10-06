/*
 * Packed storage for the first-order field modes h_{i l m}(r) (and their radial derivatives).
 *
 * Compared with a dense boost::multi_array indexed [i][l][m][j] with i = 1..10,
 * 0 <= l <= lmax, -lmax <= m <= lmax, this stores only
 *   - m >= 0 (negative m follow from h_{i l -m} = (-1)^m conj(h_{i l m})),
 *   - |m| <= l,
 *   - the components present for the mode's parity (i = 1..7 for l+m even,
 *     i = 8..10 for l+m odd, fewer at l = 0, 1),
 * which is ~8x less memory at lmax = 50.
 *
 * Real and imaginary parts are stored in separate planes, so that each component is a
 * contiguous array of doubles over the grid and the kernels can be vectorised across grid
 * points.
 */
#pragma once

#include <complex>
#include <cstddef>
#include <vector>

class Field {
public:
  /* A resolved (i, l, m) mode: real and imaginary planes of the stored m >= 0 mode (or
   * zeros for an absent component); for m < 0 the values are sign*conj(stored). */
  struct Mode {
    const double *re, *im;
    bool negate_conj;   /* true for m < 0 */
    double sign;        /* (-1)^m for m < 0 */
  };

  Field() : lmax_(-1), N_(0) {}

  /* Set up the layout for modes 0 <= l <= lmax, 0 <= m <= l on N grid points. Storage is
   * allocated for the components listed by present(i, l, m). */
  void init(int lmax, std::size_t N);

  /* Components stored for an (l, m >= 0) mode. */
  static bool present(int i, int l, int m);

  /* Writable planes of a stored (i, l, m >= 0) component. */
  double *slot_re(int i, int l, int m);
  double *slot_im(int i, int l, int m);
  void set(int i, int l, int m, std::size_t j, std::complex<double> v) {
    slot_re(i, l, m)[j] = v.real();
    slot_im(i, l, m)[j] = v.imag();
  }
  std::complex<double> get(int i, int l, int m, std::size_t j) {
    return std::complex<double>(slot_re(i, l, m)[j], slot_im(i, l, m)[j]);
  }

  inline Mode mode(int i, int l, int m) const {
    const int am = m < 0 ? -m : m;
    const double *re = zeros_.data(), *im = zeros_.data();
    if (l >= 0 && l <= lmax_ && am <= l && i >= 1 && i <= 10) {
      const long off = offset_[index(i, l, am)];
      if (off >= 0) { re = re_.data() + off; im = im_.data() + off; }
    }
    if (m < 0)
      return Mode{re, im, true, (am % 2) ? -1.0 : 1.0};
    return Mode{re, im, false, 1.0};
  }

  const double *zeros() const { return zeros_.data(); }
  int lmax() const { return lmax_; }
  std::size_t size() const { return N_; }
  std::size_t bytes() const { return (re_.size() + im_.size() + zeros_.size())*sizeof(double); }

private:
  inline std::size_t index(int i, int l, int m) const {
    return (static_cast<std::size_t>(i - 1)*(lmax_ + 1) + l)*(lmax_ + 1) + m;
  }
  long offset(int i, int l, int m) const;
  int lmax_;
  std::size_t N_;
  std::vector<long> offset_;           /* [i][l][m>=0] -> offset into re_/im_, or -1 */
  std::vector<double> re_, im_;
  std::vector<double> zeros_;
};
