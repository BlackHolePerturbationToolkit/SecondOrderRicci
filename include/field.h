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
 * Access keeps the dense array syntax h[i][l][m][j] used by the generated Ricci
 * expressions: h[i][l][m] resolves (once per mode) to a pointer into the packed data, or
 * to a shared buffer of zeros for absent components, plus the conjugation needed for
 * m < 0; [j] then returns the value at grid point j.
 */
#pragma once

#include <complex>
#include <cstddef>
#include <vector>

class Field {
public:
  typedef std::complex<double> value_type;

  /* A resolved (i, l, m) mode: indexing by grid point j. */
  struct Mode {
    const value_type *p;
    bool negate_conj;   /* true for m < 0 */
    double sign;        /* (-1)^m for m < 0 */
    inline value_type operator[](std::size_t j) const {
      return negate_conj ? sign*std::conj(p[j]) : p[j];
    }
    inline const value_type *data() const { return p; }
  };

  struct LIndex {
    const Field *f; int i; int l;
    inline Mode operator[](int m) const { return f->mode(i, l, m); }
  };
  struct IIndex {
    const Field *f; int i;
    inline LIndex operator[](int l) const { return LIndex{f, i, l}; }
  };

  Field() : lmax_(-1), N_(0) {}

  /* Set up the layout for modes 0 <= l <= lmax, 0 <= m <= l on N grid points. Storage is
   * allocated for the components listed by present(i, l, m). */
  void init(int lmax, std::size_t N);

  /* Components stored for an (l, m >= 0) mode. */
  static bool present(int i, int l, int m);

  /* Writable pointer to the N values of a stored (i, l, m >= 0) component. */
  value_type *slot(int i, int l, int m);

  inline Mode mode(int i, int l, int m) const {
    const int am = m < 0 ? -m : m;
    const value_type *p = zeros_.data();
    if (l >= 0 && l <= lmax_ && am <= l && i >= 1 && i <= 10) {
      const long off = offset_[index(i, l, am)];
      if (off >= 0) p = data_.data() + off;
    }
    if (m < 0)
      return Mode{p, true, (am % 2) ? -1.0 : 1.0};
    return Mode{p, false, 1.0};
  }

  inline IIndex operator[](int i) const { return IIndex{this, i}; }

  const value_type *zeros() const { return zeros_.data(); }
  int lmax() const { return lmax_; }
  std::size_t size() const { return N_; }
  std::size_t bytes() const { return data_.size()*sizeof(value_type) + zeros_.size()*sizeof(value_type); }

private:
  inline std::size_t index(int i, int l, int m) const {
    return (static_cast<std::size_t>(i - 1)*(lmax_ + 1) + l)*(lmax_ + 1) + m;
  }
  int lmax_;
  std::size_t N_;
  std::vector<long> offset_;           /* [i][l][m>=0] -> offset into data_, or -1 */
  std::vector<value_type> data_;
  std::vector<value_type> zeros_;
};
