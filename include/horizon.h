/*
 * d2R near the horizon, from the expansion of the first-order modes about the horizon,
 *   h^i_lm(r) = exp(-i m Omega r*) sum_k b_k (r - 2M)^k,   r* = r + 2M log(r/2M - 1)
 * (t slice, BLS components, PerturbationEquations numbering), with the kernels R2V_i that are
 * regular at the horizon (scripts/generate_kernels.wl vslice). They take the regular
 * (v, r)-frame fields a_j of each mode,
 *   a1 = h1 + f h3,  a2 = -h3 - (h1 - h2)/f,  a3 = 2 (h1 - h2)/f^2,  a5 = (h5 - h4)/f,
 *   a9 = (h9 - h8)/f,  a_j = h_j otherwise   (phase factored out),
 * and return exp(i m3 Omega r*) d2R_i3 (the normalisation of Ricci's src.h5). Two forms:
 *   - values on grid points near the horizon (where the t-slice sum loses accuracy to the
 *     cancellation of terms that diverge at the horizon);
 *   - the exact Taylor coefficients of order 0, 1, 2 in r - 2M (kernels evaluated in
 *     truncated power-series arithmetic), as written by NearHorizonRicci.m.
 * Piece::RetRet is d2R[h1, h1] summed over all modes l1, l2 <= lmax1; Piece::S1 is
 * d2R[hS, h1] + d2R[h1, hS] with hS the spin perturbation (as Ricci_S),
 *   h^8_10 = -sqrt(pi/3) 16 S1/r,  h^9_10 = -sqrt(pi/3) 32 M S1/r^2   (S1 = 1).
 */
#pragma once

#include <array>
#include <complex>
#include <map>
#include <string>
#include <tuple>
#include <vector>
#include "utils.h"
#include "field.h"

class HorizonData {
public:
  enum class Piece { RetRet, S1 };

  /* Horizon expansion coefficients of the first-order modes with l <= lmax. `source` is a
   * table "l m i Re b0 Im b0 Re b1 Im b1 ..." (i in the h1 file numbering, as written by the
   * boundary_data stage), or a directory (or several, separated by ':', later ones taking
   * precedence) of h1 mode files h1-l<l>m<m>.h5 with a "horizon" dataset (as written by
   * h1Lorenz: row n, Re/Im of the coefficient c_n of each field in the order of inhom_left). */
  void load(const std::string &source, int lmax);
  int lmax() const { return lmax_; }
  int num_coefficients() const { return ncoef_; }

  /* exp(i m3 Omega r*) d2R_ilm on the points r (one vector per mode) */
  std::vector<std::vector<std::complex<double>>> values(Piece piece, double r0, const std::vector<double> &r,
                                                        const std::vector<ilm_mode> &modes) const;
  /* Taylor coefficients of order 0, 1, 2 in r - 2M of exp(i m3 Omega r*) d2R_ilm */
  std::vector<std::array<std::complex<double>, 3>> coefficients(Piece piece, double r0,
                                                                const std::vector<ilm_mode> &modes) const;

private:
  typedef std::vector<std::complex<double>> Poly;
  std::map<std::tuple<int,int,int>, Poly> b_;   /* (l, m, i) -> coefficients, PE numbering */
  int lmax_ = -1, ncoef_ = 0;
  template<class Fill> void build(Field &A, Field &dA, Field &ddA, std::size_t N, bool spin_only, Fill fill) const;
  template<bool Series> void sum(Piece piece, double r0, const std::vector<double> &r, const std::vector<ilm_mode> &modes,
                                 std::vector<std::vector<std::complex<double>>> &out) const;
};

/* r* = r + 2M log(r/2M - 1) */
double rstar(double r);
