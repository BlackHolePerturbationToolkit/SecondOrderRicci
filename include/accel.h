/*
 * Per-point convergence control and Wynn-epsilon acceleration for the shell-ordered
 * mode-coupling sum.
 *
 * For one output mode, each grid point j keeps the last `window` partial sums
 * S_n(j) = sum_{s <= n} c_s(j) over the non-empty shells, and the Wynn epsilon estimate
 * W_n(j) of the limit. A point is converged when, on `confirm` consecutive shells (and
 * after at least `min_shells` non-empty shells),
 *   |W_n - W_{n-1}| <= tol * scale   and   |c_n| <= sqrt(tol) * scale,
 * where scale = max(|W_n(j)|, floor) and the floor is `floor_rel` times the largest |S_n|
 * on the grid (so that points where the mode passes through zero do not demand
 * unattainable relative accuracy). The second condition keeps the acceleration from being
 * trusted while the plain partial sums are still far from converged: acceleration then
 * supplies the last ~half of the digits. Converged points keep W_n(j) and are no longer
 * computed. With tol = 1e-10, min_shells = 6, confirm = 3 the error at r0 = 8 (l <= 50) is
 * <= 1.1e-10 relative, at ~6% of the work of summing all shells at every point.
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>

/* Wynn epsilon algorithm on s[0..n-1]; returns the estimate from the highest even column
 * (its last element, i.e. the one using the most recent partial sums). Stops early if a
 * difference vanishes or the table becomes non-finite. */
inline std::complex<double> wynn_epsilon(const std::complex<double> *s, int n)
{
  typedef std::complex<double> C;
  if (n <= 0) return C(0.0, 0.0);
  C prev[64], cur[64], next[64];
  n = std::min(n, 63);
  for (int i = 0; i <= n; ++i) prev[i] = C(0.0, 0.0);
  for (int i = 0; i < n; ++i) cur[i] = s[i];
  C best = s[n-1];
  int len = n;
  for (int k = 1; k < n; ++k) {
    for (int i = 0; i + 1 < len; ++i) {
      const C d = cur[i+1] - cur[i];
      if (d == C(0.0, 0.0)) return best;
      next[i] = prev[i+1] + 1.0/d;
      if (!std::isfinite(next[i].real()) || !std::isfinite(next[i].imag())) return best;
    }
    len -= 1;
    for (int i = 0; i <= len; ++i) prev[i] = cur[i];
    for (int i = 0; i < len; ++i) cur[i] = next[i];
    if (k % 2 == 0) best = cur[len-1];
  }
  return best;
}

class ShellConvergence {
public:
  typedef std::complex<double> C;

  ShellConvergence(std::size_t N, int window, double tol, double floor_rel, int confirm, bool accelerate,
                   int min_shells = 6)
    : N_(N), window_(window), tol_(tol), floor_rel_(floor_rel), confirm_(confirm), accelerate_(accelerate),
      min_shells_(min_shells), ctol_(std::sqrt(tol)),
      hist_(static_cast<std::size_t>(window)*N), nhist_(0), est_(N), ok_count_(N, 0), active_(N, 1),
      shells_used_(N, 0) {}

  /* Mark points that are not to be computed at all (e.g. inside the worldtube). */
  void deactivate(std::size_t jb, std::size_t je) {
    for (std::size_t j = jb; j < je; ++j) active_[j] = 0;
  }

  bool active(std::size_t j) const { return active_[j] != 0; }
  std::size_t num_active() const { return static_cast<std::size_t>(std::count(active_.begin(), active_.end(), 1)); }

  /* Record the partial sums S (all points; only active ones are used) after one more
   * non-empty shell, update estimates and convergence flags. Returns the number of points
   * that became converged. */
  std::size_t update(const C *S) {
    const std::size_t slot = static_cast<std::size_t>(nhist_ % window_);
    for (std::size_t j = 0; j < N_; ++j)
      if (active_[j]) hist_[slot*N_ + j] = S[j];
    ++nhist_;
    const int n = static_cast<int>(std::min<long>(nhist_, window_));

    double smax = 0.0;
    for (std::size_t j = 0; j < N_; ++j)
      if (active_[j]) smax = std::max(smax, std::abs(S[j]));
    max_seen_ = std::max(max_seen_, smax);
    const double floor = floor_rel_*max_seen_;

    std::size_t newly = 0;
    C seq[64];
    for (std::size_t j = 0; j < N_; ++j) {
      if (!active_[j]) continue;
      shells_used_[j] = static_cast<int>(nhist_);
      /* the last n partial sums in order */
      for (int k = 0; k < n; ++k) {
        const long idx = (nhist_ - n + k) % window_;
        seq[k] = hist_[static_cast<std::size_t>(idx)*N_ + j];
      }
      const double cn = (n >= 2) ? std::abs(seq[n-1] - seq[n-2]) : std::abs(seq[n-1]);
      C w = accelerate_ ? wynn_epsilon(seq, n) : seq[n-1];
      /* Trust the acceleration only for a plausible tail: a geometric tail with ratio
       * q <= 0.9 changes the partial sum by at most 9 |c_n|. Larger corrections come from
       * rounding noise once the contributions are negligible (e.g. near the horizon). */
      if (std::abs(w - seq[n-1]) > 10.0*cn) w = seq[n-1];
      const double scale = std::max(std::abs(w), floor);
      if (nhist_ >= min_shells_ && std::abs(w - est_[j]) <= tol_*scale && cn <= ctol_*scale)
        ++ok_count_[j];
      else
        ok_count_[j] = 0;
      est_[j] = w;
      if (ok_count_[j] >= confirm_ && floor > 0.0) {
        active_[j] = 0;
        ++newly;
      }
    }
    return newly;
  }

  /* Final values: the estimate for every point that was computed (converged or not). */
  const std::vector<C> &estimates() const { return est_; }
  const std::vector<int> &shells_used() const { return shells_used_; }
  /* Points still active when the shells ran out (not converged to tol). */
  std::vector<char> unconverged() const { return active_; }

private:
  std::size_t N_;
  int window_;
  double tol_, floor_rel_;
  int confirm_;
  bool accelerate_;
  int min_shells_;
  double ctol_;
  std::vector<C> hist_;
  long nhist_;
  std::vector<C> est_;
  std::vector<int> ok_count_;
  std::vector<char> active_;
  std::vector<int> shells_used_;
  double max_seen_ = 0.0;
};

/*
 * Levin u-transform of the partial sums S[0..L] of a series with terms c[0..L], using the last
 * k + 1 partial sums and the remainder estimates (n + beta) c_n:
 *   L_k = sum_j w_j S_{n0+j}/omega_{n0+j} / sum_j w_j/omega_{n0+j},   n0 = L - k,
 *   w_j = (-1)^j binom(k, j) ((n0 + j + beta)/(n0 + k + beta))^(k-1).
 * It is exact for remainders S - S_n = omega_n sum_{i<k} a_i/(n + beta)^i, which covers both the
 * algebraic convergence (c_s ~ s^-p) of the worldtube pieces at r = r0 and the mixed geometric x
 * algebraic convergence (c_s ~ rho^s s^-p) away from it. Returns S[L] if a term vanishes.
 */
inline std::complex<double> levin_u(const std::complex<double> *S, const std::complex<double> *c, int L, int k, double beta = 1.0)
{
  typedef std::complex<double> C;
  const int n0 = L - k;
  if (n0 < 0) return S[L];
  C num(0.0, 0.0), den(0.0, 0.0);
  double binom = 1.0;
  for (int j = 0; j <= k; ++j) {
    const int n = n0 + j;
    if (c[n] == C(0.0, 0.0)) return S[L];
    const C om = (n + beta)*c[n];
    const double w = ((j % 2) ? -1.0 : 1.0)*binom*std::pow((n0 + j + beta)/(n0 + k + beta), k - 1);
    num += w*S[n]/om;
    den += w/om;
    binom = binom*(k - j)/(j + 1.0);
  }
  if (den == C(0.0, 0.0) || !std::isfinite(std::abs(num/den))) return S[L];
  return num/den;
}
