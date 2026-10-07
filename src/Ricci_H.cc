/*
 * d2R[h1, h1] near the horizon on the ingoing (v) slice.
 *
 * The first-order modes are built from their expansion about the horizon,
 *   h^i_lm(r) = exp(-i m Omega r*) sum_k b_k (r - 2M)^k   (t slice, BLS components),
 * converted to the regular (v, r)-frame fields a_j of scripts/generate_kernels.wl (vslice), and
 * coupled with the kernels R2V_i, whose coefficients are all regular at the horizon. The output
 *   src i= l= m=  :  exp(i m Omega r*) d2R_ilm(r)
 * (original normalisation, as Ricci's src.h5 times the phase) on the grid points r < rmax of an
 * h1 mode file, replaces the symbolic near-horizon expansion (NearHorizonRicci.m).
 *
 * With the optional last argument S1, computes instead the coupling to the spin perturbation
 * (as NearHorizonRicci's RicciSumS1): the odd dipole (l, m) = (1, 0) is replaced by
 *   h^8 = -sqrt(pi/3) 16 S1/r,  h^9 = -sqrt(pi/3) 32 M S1/r^2   (S1 = 1),
 * and only the terms d2R[hS, h] + d2R[h, hS] with hS that dipole are summed. Output src_H_S1.h5.
 *
 * The grid is either that of an h1 mode file (points r < rmax), or, with grid = log:N, N points
 * with r - 2M logarithmically spaced between 1e-6 M and rmax - 2M.
 *
 * Usage: Ricci_H <horizon coefficients> <r0> <h1 mode file | log:N> <rmax> <lmax> <lmax_out> [S1]
 */
#include <vector>
#include <complex>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cassert>
#include <map>
#include <boost/multi_array.hpp>
#ifdef _OPENMP
#  include <omp.h>
#endif
#include <hdf5.h>
#include "field.h"
#include "modeblock.h"
#include "R2V.h"
#include "utils.h"
#include "h5wrapper.h"

using namespace std;
typedef complex<double> cplx;

const double M = 1.0;

typedef void (*R2G_function)(const double, const double, const vector<double>&, const vector<double>&,
                             const vector<double>&, const ModeBlock&, const ModeBlock&, const ModeBlock&,
                             const ModeBlock&, const ModeBlock&, const ModeBlock&,
                             int, int, int, int, double *const*, double *const*, size_t, size_t);
typedef void (*R2C_function)(int, int, int, int, int, int, double*);
static const R2G_function R2V[11] = {nullptr, R2V_1, R2V_2, R2V_3, R2V_4, R2V_5, R2V_6, R2V_7, R2V_8, R2V_9, R2V_10};
static const R2C_function R2C[11] = {nullptr, R2C_1, R2C_2, R2C_3, R2C_4, R2C_5, R2C_6, R2C_7, R2C_8, R2C_9, R2C_10};

/* polynomial in x = r - 2M */
typedef vector<cplx> Poly;

static Poly mul(const Poly &a, const Poly &b) {
  if (a.empty() || b.empty()) return Poly();
  Poly c(a.size() + b.size() - 1, 0.0);
  for (size_t i = 0; i < a.size(); ++i)
    for (size_t j = 0; j < b.size(); ++j)
      c[i+j] += a[i]*b[j];
  return c;
}
static Poly sub(const Poly &a, const Poly &b) {
  Poly c(max(a.size(), b.size()), 0.0);
  for (size_t i = 0; i < a.size(); ++i) c[i] += a[i];
  for (size_t i = 0; i < b.size(); ++i) c[i] -= b[i];
  return c;
}
static Poly scale(const Poly &a, cplx s) { Poly c(a); for (auto &v : c) v *= s; return c; }
/* divide by x^n: the dropped coefficients must vanish (regularity on the future horizon) */
static Poly shift(const Poly &a, size_t n, double tol, const string &what) {
  double amax = 0, low = 0;
  for (size_t i = 0; i < a.size(); ++i) amax = max(amax, abs(a[i]));
  for (size_t i = 0; i < min(n, a.size()); ++i) low = max(low, abs(a[i]));
  if (low > tol*amax && low > 1e-280) {
    cerr << "Warning: " << what << " is not regular at the horizon (" << low << " vs " << amax << ")" << endl;
  }
  return a.size() > n ? Poly(a.begin() + n, a.end()) : Poly();
}
/* value and first two derivatives */
static void eval(const Poly &p, double x, cplx &v, cplx &d, cplx &dd) {
  v = d = dd = 0.0;
  for (size_t k = p.size(); k-- > 0;) {
    dd = dd*x + 2.0*d;
    d = d*x + v;
    v = v*x + p[k];
  }
}

int main(int argc, char* argv[])
{
  cout << "Version: " << __GIT_VERSION << endl;
  if (argc != 7 && !(argc == 8 && string(argv[7]) == "S1")) {
    cout << "Usage: " << argv[0] << " <horizon coefficients> <r0> <h1 mode file | log:N> <rmax> <lmax> <lmax_out> [S1]" << endl;
    exit(1);
  }
  const bool spin = (argc == 8);
  const string coeffFile = argv[1];
  const double r0 = atof(argv[2]);
  const string gridFile = argv[3];
  const double rmax = atof(argv[4]);
  const int lmax = atoi(argv[5]);
  const int lmax_out = atoi(argv[6]);
  const double Omega = sqrt(M/(r0*r0*r0));

  /* grid points r < rmax */
  vector<double> r, f, fp;
  if (gridFile.compare(0, 4, "log:") == 0) {
    const int n = atoi(gridFile.c_str() + 4);
    const double x0 = 1e-6*M, x1 = rmax - 2.0*M;
    assert(n >= 2 && x1 > x0);
    for (int k = 0; k < n; ++k) {
      const double x = x0*pow(x1/x0, k/(n - 1.0));
      r.push_back(2.0*M + x); f.push_back(1.0 - 2.0*M/(2.0*M + x)); fp.push_back(2.0*M/((2.0*M + x)*(2.0*M + x)));
    }
  } else {
    H5F h5_file(gridFile, H5F_ACC_RDONLY);
    H5D dataset(h5_file, "grid");
    H5S dataspace(dataset);
    vector<hsize_t> gridSize = dataspace.getSimpleExtentDims();
    vector<double> grid(gridSize[0]);
    H5S memspace(gridSize);
    H5Dread(dataset.getId(), H5T_NATIVE_DOUBLE, memspace.getId(), dataspace.getId(), H5P_DEFAULT, grid.data());
    for (double x : grid)
      if (x < rmax) { r.push_back(x); f.push_back(1.0 - 2.0*M/x); fp.push_back(2.0*M/(x*x)); }
  }
  const size_t N = r.size();
  cout << "Grid: " << N << " points in [" << r.front() << ", " << r.back() << "]" << endl;

  /* horizon coefficients: l m i, then pairs (Re, Im) of b_0, b_1, ...; i in the file numbering,
   * which has components 3 and 6 swapped relative to PerturbationEquations */
  map<tuple<int,int,int>, Poly> b;
  {
    ifstream in(coeffFile);
    if (!in) { cerr << "Cannot read " << coeffFile << endl; exit(1); }
    string line;
    while (getline(in, line)) {
      istringstream ss(line);
      double l, m, i;
      if (!(ss >> l >> m >> i)) continue;
      Poly c;
      double re, im;
      while (ss >> re >> im) c.push_back(cplx(re, im));
      const int ii = (int) lround(i);
      const int ipe = (ii == 3) ? 6 : ((ii == 6) ? 3 : ii);
      if ((int) lround(l) <= lmax) b[make_tuple((int) lround(l), (int) lround(m), ipe)] = c;
    }
  }

  /* regular (v, r)-frame fields a_j and their r-derivatives */
  Field A, dA, ddA;
  A.init(lmax, N); dA.init(lmax, N); ddA.init(lmax, N);
  const double tol = 1e-12;
  const Poly rpoly = {2.0*M, 1.0};            /* r = 2M + x */
  for (int l = 0; l <= lmax; ++l) {
    for (int m = 0; m <= l; ++m) {
      Poly h[11];
      for (int i = 1; i <= 10; ++i) {
        auto it = b.find(make_tuple(l, m, i));
        if (it != b.end()) h[i] = it->second;
      }
      stringstream lm; lm << "(l, m) = (" << l << ", " << m << ")";
      const Poly d12 = sub(h[1], h[2]);
      Poly a[11];
      /* a1 = h1 + f h3 is evaluated directly below */
      a[2] = sub(scale(h[3], -1.0), mul(rpoly, shift(d12, 1, tol, "h1 - h2 " + lm.str())));
      a[3] = scale(mul(mul(rpoly, rpoly), shift(d12, 2, tol, "h1 - h2 " + lm.str())), 2.0);
      a[4] = h[4];
      a[5] = mul(rpoly, shift(sub(h[5], h[4]), 1, tol, "h5 - h4 " + lm.str()));
      a[6] = h[6]; a[7] = h[7]; a[8] = h[8];
      a[9] = mul(rpoly, shift(sub(h[9], h[8]), 1, tol, "h9 - h8 " + lm.str()));
      a[10] = h[10];
      for (int i = 1; i <= 10; ++i) {
        if (!Field::present(i, l, m)) continue;
        double *vr = A.slot_re(i, l, m), *vi = A.slot_im(i, l, m);
        double *dr = dA.slot_re(i, l, m), *di = dA.slot_im(i, l, m);
        double *ddr = ddA.slot_re(i, l, m), *ddi = ddA.slot_im(i, l, m);
        for (size_t j = 0; j < N; ++j) {
          const double x = r[j] - 2.0*M;
          cplx v, d, dd;
          if (i == 1) {
            cplx h1, dh1, ddh1, h3, dh3, ddh3;
            eval(h[1], x, h1, dh1, ddh1);
            eval(h[3], x, h3, dh3, ddh3);
            const double ff = f[j], f1 = 2.0*M/(r[j]*r[j]), f2 = -4.0*M/(r[j]*r[j]*r[j]);
            v = h1 + ff*h3;
            d = dh1 + f1*h3 + ff*dh3;
            dd = ddh1 + f2*h3 + 2.0*f1*dh3 + ff*ddh3;
          } else {
            eval(a[i], x, v, d, dd);
          }
          vr[j] = v.real(); vi[j] = v.imag();
          dr[j] = d.real(); di[j] = d.imag();
          ddr[j] = dd.real(); ddi[j] = dd.imag();
        }
      }
    }
  }

  /* spin perturbation: (1, 0) odd dipole, regular frame a8 = h8, a9 = (h9 - h8)/f (exact) */
  if (spin) {
    const double c = 16.0*sqrt(M_PI/3.0);
    double *p[3][2][2] = {{{A.slot_re(8,1,0), A.slot_im(8,1,0)}, {A.slot_re(9,1,0), A.slot_im(9,1,0)}},
                          {{dA.slot_re(8,1,0), dA.slot_im(8,1,0)}, {dA.slot_re(9,1,0), dA.slot_im(9,1,0)}},
                          {{ddA.slot_re(8,1,0), ddA.slot_im(8,1,0)}, {ddA.slot_re(9,1,0), ddA.slot_im(9,1,0)}}};
    for (size_t j = 0; j < N; ++j) {
      const double x = r[j];
      const double v[3] = {1.0/x, -1.0/(x*x), 2.0/(x*x*x)};   /* 1/r and its derivatives */
      for (int s = 0; s < 3; ++s) {
        p[s][0][0][j] = -c*v[s]; p[s][0][1][j] = 0.0;     /* a8 = -c/r */
        p[s][1][0][j] =  c*v[s]; p[s][1][1][j] = 0.0;     /* a9 = +c/r */
      }
    }
  }

  /* output modes */
  vector<ilm_mode> modes;
  for (int l3 = 0; l3 <= lmax_out; ++l3)
    for (int m3 = 0; m3 <= l3; ++m3)
      for (int i3 = 1; i3 <= 10; ++i3)
        if (Field::present(i3, l3, m3))
          modes.push_back(ilm_mode({i3, l3, m3}));
  cout << "Computing source: " << modes.size() << " (i,l,m) modes, first-order l <= " << lmax << endl;

  boost::multi_array<cplx, 2> src(boost::extents[modes.size()][N]);
  fill(src.data(), src.data() + src.num_elements(), 0.0);
  const int num_modes = modes.size();
#pragma omp parallel
  {
    ModeBuffer bX[3], bY[3];
    for (int k = 0; k < 3; ++k) { bX[k].reserve(N); bY[k].reserve(N); }
    ModeBlock X[3], Y[3];
    GBuffer G(R2_NC_MAX, N);
    double c[R2_NC_MAX];
#pragma omp for schedule(dynamic)
    for (int im = 0; im < num_modes; ++im) {
      const int i3 = modes[im].i, l3 = modes[im].l, m3 = modes[im].m;
      const int K = R2_NC[i3];
      for (int l1 = 0; l1 <= lmax; ++l1) {
        for (int m1 = -l1; m1 <= l1; ++m1) {
          /* S1: only the pairs containing the spin dipole, in either slot */
          bX[0].resolve(A, l1, m1, 0, N, X[0]); bX[1].resolve(dA, l1, m1, 0, N, X[1]); bX[2].resolve(ddA, l1, m1, 0, N, X[2]);
          for (int l2 = abs(l3 - l1); l2 <= min(lmax, l3 + l1); ++l2) {
            const int m2 = m3 - m1;
            if (abs(m2) > l2) continue;
            if (spin && !((l1 == 1 && m1 == 0) || (l2 == 1 && m2 == 0))) continue;
            bY[0].resolve(A, l2, m2, 0, N, Y[0]); bY[1].resolve(dA, l2, m2, 0, N, Y[1]); bY[2].resolve(ddA, l2, m2, 0, N, Y[2]);
            G.zero(K, 0, N);
            R2V[i3](M, r0, r, f, fp, X[0], X[1], X[2], Y[0], Y[1], Y[2], l1, m1, l2, m2, G.re(), G.im(), 0, N);
            R2C[i3](l3, m3, l1, m1, l2, m2, c);
            /* RicciSumS1 counts the (S, S) pair once in each of its two terms */
            const double w = (spin && l1 == 1 && m1 == 0 && l2 == 1 && m2 == 0) ? 2.0 : 1.0;
            for (size_t j = 0; j < N; ++j)
              src[im][j] += w*G.combine(K, c, j);
          }
        }
      }
    }
  }

  cout << "Saving results" << endl;
  H5F src_h5(spin ? "src_H_S1.h5" : "src_H.h5", H5F_ACC_TRUNC);
  for (int im = 0; im < num_modes; ++im) {
    boost::multi_array<double,2> data(boost::extents[N][2]);
    stringstream srcds;
    srcds << "src i=" << modes[im].i << " l=" << modes[im].l << " m=" << modes[im].m;
    for (size_t j = 0; j < N; ++j) {
      data[j][0] = real(src[im][j]);
      data[j][1] = imag(src[im][j]);
    }
    src_h5.write_dataset(srcds.str(), data);
  }
  boost::multi_array<double,1> rdata(boost::extents[N]);
  for (size_t j = 0; j < N; ++j) rdata[j] = r[j];
  src_h5.write_dataset("r", rdata);
  cout << "Done" << endl;
  return 0;
}
