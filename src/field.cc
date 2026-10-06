#include "field.h"
#include <stdexcept>

bool Field::present(int i, int l, int m)
{
  if (m < 0 || m > l) return false;
  const bool even = ((l + m) % 2) == 0;
  if (even) {
    /* even parity: 1, 2, 3, 6 always; 4, 5 need l >= 1; 7 needs l >= 2 */
    if (i == 1 || i == 2 || i == 3 || i == 6) return true;
    if (i == 4 || i == 5) return l >= 1;
    if (i == 7) return l >= 2;
    return false;
  }
  /* odd parity: 8, 9 need l >= 1; 10 needs l >= 2 */
  if (i == 8 || i == 9) return l >= 1;
  if (i == 10) return l >= 2;
  return false;
}

void Field::init(int lmax, std::size_t N)
{
  lmax_ = lmax;
  N_ = N;
  offset_.assign(static_cast<std::size_t>(10)*(lmax + 1)*(lmax + 1), -1);
  long next = 0;
  for (int i = 1; i <= 10; ++i)
    for (int l = 0; l <= lmax; ++l)
      for (int m = 0; m <= l; ++m)
        if (present(i, l, m)) {
          offset_[index(i, l, m)] = next;
          next += static_cast<long>(N);
        }
  data_.assign(static_cast<std::size_t>(next), value_type(0.0, 0.0));
  zeros_.assign(N, value_type(0.0, 0.0));
}

Field::value_type *Field::slot(int i, int l, int m)
{
  if (l < 0 || l > lmax_ || m < 0 || m > l || i < 1 || i > 10)
    throw std::out_of_range("Field::slot: mode outside the stored range");
  const long off = offset_[index(i, l, m)];
  if (off < 0)
    throw std::out_of_range("Field::slot: component not stored for this parity");
  return data_.data() + off;
}
