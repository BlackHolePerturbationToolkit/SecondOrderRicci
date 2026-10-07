/*
 * Truncated expansions at large r for evaluating the u-slice kernels as series:
 *   sum_{n < YLOG_N} sum_{p <= YLOG_P} c[n][p] y^n (log r)^p,   y = 1/r,
 * with real coefficients. The kernels use only +, -, * and 1/r; r is represented by the tag type
 * InfR, for which 1/r gives y. YLOG_N = 5 gives d2R ~ y^2, y^3, y^4 exactly; products of two
 * first-order fields with (log r)^1 terms give at most (log r)^2.
 */
#pragma once

#define YLOG_N 5
#define YLOG_P 2

struct YLog {
  double c[YLOG_N][YLOG_P + 1];
  YLog(double v = 0.0) { for (int n = 0; n < YLOG_N; ++n) for (int p = 0; p <= YLOG_P; ++p) c[n][p] = 0.0; c[0][0] = v; }
  YLog &operator+=(const YLog &b) { for (int n = 0; n < YLOG_N; ++n) for (int p = 0; p <= YLOG_P; ++p) c[n][p] += b.c[n][p]; return *this; }
  YLog &operator-=(const YLog &b) { for (int n = 0; n < YLOG_N; ++n) for (int p = 0; p <= YLOG_P; ++p) c[n][p] -= b.c[n][p]; return *this; }
};

/* r at large radius: only 1/r is used */
struct InfR {};

inline YLog operator/(double s, const InfR &) { YLog y; y.c[1][0] = s; return y; }
inline YLog operator+(YLog a, const YLog &b) { return a += b; }
inline YLog operator-(YLog a, const YLog &b) { return a -= b; }
inline YLog operator-(const YLog &a) { YLog r; for (int n = 0; n < YLOG_N; ++n) for (int p = 0; p <= YLOG_P; ++p) r.c[n][p] = -a.c[n][p]; return r; }
inline YLog operator*(const YLog &a, const YLog &b) {
  YLog r;
  for (int n = 0; n < YLOG_N; ++n)
    for (int k = 0; k <= n; ++k)
      for (int p = 0; p <= YLOG_P; ++p)
        if (a.c[k][p] != 0.0)
          for (int q = 0; p + q <= YLOG_P; ++q) r.c[n][p+q] += a.c[k][p]*b.c[n-k][q];
  return r;
}
inline YLog operator*(double s, YLog a) { for (int n = 0; n < YLOG_N; ++n) for (int p = 0; p <= YLOG_P; ++p) a.c[n][p] *= s; return a; }
inline YLog operator*(YLog a, double s) { return s*a; }
inline YLog operator+(double s, YLog a) { a.c[0][0] += s; return a; }
inline YLog operator+(YLog a, double s) { a.c[0][0] += s; return a; }
inline YLog operator-(YLog a, double s) { a.c[0][0] -= s; return a; }
inline YLog operator-(double s, const YLog &a) { return s + (-a); }
