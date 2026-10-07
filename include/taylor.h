/*
 * Truncated Taylor series sum_{k < N} c_k x^k with real coefficients, for evaluating the
 * generated kernels as power series about the horizon (x = r - 2M): the kernels only use
 * +, -, * and 1/r, so instantiating them with this type gives the exact Taylor coefficients
 * of the output to order N - 1.
 */
#pragma once

template<int N>
struct Taylor {
  double c[N];
  Taylor(double v = 0.0) { c[0] = v; for (int k = 1; k < N; ++k) c[k] = 0.0; }
  Taylor &operator+=(const Taylor &b) { for (int k = 0; k < N; ++k) c[k] += b.c[k]; return *this; }
  Taylor &operator-=(const Taylor &b) { for (int k = 0; k < N; ++k) c[k] -= b.c[k]; return *this; }
};

template<int N> inline Taylor<N> operator+(Taylor<N> a, const Taylor<N> &b) { return a += b; }
template<int N> inline Taylor<N> operator-(Taylor<N> a, const Taylor<N> &b) { return a -= b; }
template<int N> inline Taylor<N> operator-(const Taylor<N> &a) { Taylor<N> r; for (int k = 0; k < N; ++k) r.c[k] = -a.c[k]; return r; }
template<int N> inline Taylor<N> operator*(const Taylor<N> &a, const Taylor<N> &b) {
  Taylor<N> r;
  for (int k = 0; k < N; ++k) { double s = 0.0; for (int j = 0; j <= k; ++j) s += a.c[j]*b.c[k-j]; r.c[k] = s; }
  return r;
}
template<int N> inline Taylor<N> operator*(double s, Taylor<N> a) { for (int k = 0; k < N; ++k) a.c[k] *= s; return a; }
template<int N> inline Taylor<N> operator*(Taylor<N> a, double s) { return s*a; }
template<int N> inline Taylor<N> operator+(double s, Taylor<N> a) { a.c[0] += s; return a; }
template<int N> inline Taylor<N> operator+(Taylor<N> a, double s) { a.c[0] += s; return a; }
template<int N> inline Taylor<N> operator-(Taylor<N> a, double s) { a.c[0] -= s; return a; }
template<int N> inline Taylor<N> operator-(double s, const Taylor<N> &a) { return s + (-a); }
/* s / a, a.c[0] != 0 */
template<int N> inline Taylor<N> operator/(double s, const Taylor<N> &a) {
  Taylor<N> r;
  for (int k = 0; k < N; ++k) {
    double v = (k == 0) ? s : 0.0;
    for (int j = 1; j <= k; ++j) v -= a.c[j]*r.c[k-j];
    r.c[k] = v/a.c[0];
  }
  return r;
}

typedef Taylor<3> Taylor3;
