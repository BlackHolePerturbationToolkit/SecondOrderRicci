/*
 * d2R at large r on the outgoing (u) slice, from the expansion of the first-order modes at infinity,
 *   h^i_lm(r) = exp(i m Omega r*) sum_k (b_k + c_k log r)/r^k,   r* = r + 2M log(r/2M - 1)
 * (t slice, BLS components, PerturbationEquations numbering; log terms only for static modes),
 * with the kernels R2U_i (scripts/generate_kernels.wl uslice). They take the (u, r)-frame fields
 *   a1 = h1 + f h3,  a2 = h3 + (h1 + h2)/f,  a3 = 2 (h1 + h2)/f^2,  a5 = (h5 + h4)/f,
 *   a9 = (h9 + h8)/f,  a_j = h_j otherwise   (phase factored out),
 * and return exp(-i m3 Omega r*) d2R_i3 (the normalisation of Ricci's src.h5). Every term of
 * these kernels falls off at least as fast as d2R itself (r^-2), so nothing cancels at large r,
 * unlike the t-slice sum, whose O(omega^2) terms cancel down to O(r^-2). Two forms:
 *   - values on grid points at large r;
 *   - the exact coefficients of the expansion at infinity,
 *       exp(-i m3 Omega r*) d2R = A1/r^2 + A2/r^3 + (A3 + A3log log r)/r^4 + ...,
 *     from the kernels evaluated in series arithmetic in 1/r and log r (include/yseries.h), as
 *     written by FarZoneRicci.m.
 * Piece::RetRet is d2R[h1, h1] summed over l1, l2 <= lmax1; Piece::S1 the spin coupling (as Ricci_S).
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

class InfinityData {
public:
  enum class Piece { RetRet, S1 };

  /* Expansion coefficients at infinity of the first-order modes with l <= lmax. `source` is
   * "<coefficients table>[:<log coefficients table>]" (rows "l m i Re b0 Im b0 Re b1 ...", i in
   * the h1 file numbering, as written by the boundary_data stage). */
  void load(const std::string &source, int lmax);
  int lmax() const { return lmax_; }

  /* exp(-i m3 Omega r*) d2R_ilm on the points r */
  std::vector<std::vector<std::complex<double>>> values(Piece piece, double r0, const std::vector<double> &r,
                                                        const std::vector<ilm_mode> &modes) const;
  /* {A1, A2, A3, A3log} and the largest coefficient that should vanish (orders 1/r^0, 1/r^1,
   * log terms other than A3log), for each mode */
  std::vector<std::array<std::complex<double>, 4>> coefficients(Piece piece, double r0, const std::vector<ilm_mode> &modes,
                                                                double *residual = nullptr) const;

private:
  typedef std::vector<std::complex<double>> Poly;
  std::map<std::tuple<int,int,int>, Poly> b_, c_;   /* (l, m, i) -> coefficients of 1/r^k and log r/r^k */
  int lmax_ = -1;
  template<bool Series> void sum(Piece piece, double r0, const std::vector<double> &r, const std::vector<ilm_mode> &modes,
                                 std::vector<std::vector<std::complex<double>>> &out) const;
};
