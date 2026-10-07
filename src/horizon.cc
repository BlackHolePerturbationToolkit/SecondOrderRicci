#include "horizon.h"
#include "field.h"
#include "modeblock.h"
#include "taylor.h"
#include "R2V.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <regex>
#include <sys/stat.h>
#include <hdf5.h>
#include "utils.h"
#ifdef _OPENMP
#  include <omp.h>
#endif

using namespace std;
typedef complex<double> cplx;

static const double M = 1.0;

double rstar(double r) { return r + 2.0*M*log(r/(2.0*M) - 1.0); }

/* kernels: values on grid points, and Taylor series about the horizon */
typedef void (*R2V_grid)(const double, const double, const vector<double>&, const vector<double>&,
                         const vector<double>&, const ModeBlock&, const ModeBlock&, const ModeBlock&,
                         const ModeBlock&, const ModeBlock&, const ModeBlock&,
                         int, int, int, int, double *const*, double *const*, size_t, size_t);
typedef void (*R2V_series)(const double, const double, const Taylor3*, const Taylor3*,
                           const ModeBlockT<Taylor3>&, const ModeBlockT<Taylor3>&, const ModeBlockT<Taylor3>&,
                           const ModeBlockT<Taylor3>&, const ModeBlockT<Taylor3>&, const ModeBlockT<Taylor3>&,
                           int, int, int, int, Taylor3 *const*, Taylor3 *const*);
typedef void (*R2C_function)(int, int, int, int, int, int, double*);
static const R2V_grid KG[11] = {nullptr, R2V_1, R2V_2, R2V_3, R2V_4, R2V_5, R2V_6, R2V_7, R2V_8, R2V_9, R2V_10};
static const R2V_series KS[11] = {nullptr, R2V_1, R2V_2, R2V_3, R2V_4, R2V_5, R2V_6, R2V_7, R2V_8, R2V_9, R2V_10};
static const R2C_function KC[11] = {nullptr, R2C_1, R2C_2, R2C_3, R2C_4, R2C_5, R2C_6, R2C_7, R2C_8, R2C_9, R2C_10};

/* --- polynomials in x = r - 2M ------------------------------------------------------- */
typedef vector<cplx> Poly;
static Poly pmul(const Poly &a, const Poly &b) {
  if (a.empty() || b.empty()) return Poly();
  Poly c(a.size() + b.size() - 1, 0.0);
  for (size_t i = 0; i < a.size(); ++i)
    for (size_t j = 0; j < b.size(); ++j) c[i+j] += a[i]*b[j];
  return c;
}
static Poly psub(const Poly &a, const Poly &b) {
  Poly c(max(a.size(), b.size()), 0.0);
  for (size_t i = 0; i < a.size(); ++i) c[i] += a[i];
  for (size_t i = 0; i < b.size(); ++i) c[i] -= b[i];
  return c;
}
static Poly pscale(Poly a, cplx s) { for (auto &v : a) v *= s; return a; }
/* divide by x^n; the dropped coefficients vanish for a mode regular on the future horizon (to the
 * accuracy of the input: ~1e-12 relative for h1Lorenz) */
static Poly pshift(const Poly &a, size_t n, const string &what) {
  double amax = 0.0, low = 0.0;
  for (size_t i = 0; i < a.size(); ++i) amax = max(amax, abs(a[i]));
  for (size_t i = 0; i < min(n, a.size()); ++i) low = max(low, abs(a[i]));
  if (low > 1e-10*amax && low > 1e-280)
    cerr << "Warning: " << what << " is not regular at the horizon (" << low << " vs " << amax << ")" << endl;
  return a.size() > n ? Poly(a.begin() + n, a.end()) : Poly();
}
/* value and first two derivatives at x */
static void peval(const Poly &p, double x, cplx &v, cplx &d, cplx &dd) {
  v = d = dd = 0.0;
  for (size_t k = p.size(); k-- > 0;) { dd = dd*x + 2.0*d; d = d*x + v; v = v*x + p[k]; }
}
/* Taylor coefficients 0..2 of p, p', p'' */
static void ptaylor(const Poly &p, cplx v[3], cplx d[3], cplx dd[3]) {
  auto c = [&](size_t k) { return k < p.size() ? p[k] : cplx(0.0); };
  for (int k = 0; k < 3; ++k) { v[k] = c(k); d[k] = double(k + 1)*c(k + 1); dd[k] = double((k + 1)*(k + 2))*c(k + 2); }
}

/* file numbering of the fields in an h1 mode file (as h1.cc) */
static vector<int> file_fields(int l, int m)
{
  if ((l + m) % 2 == 0) {
    if (l == 0) return {1, 3, 6, 2};
    if (l == 1) return {1, 3, 5, 6, 2, 4};
    return {1, 3, 5, 6, 7, 2, 4};
  }
  if (l == 1) return {9, 8};
  return {9, 10, 8};
}

void HorizonData::load(const string &source, int lmax)
{
  struct stat st;
  const bool table = stat(source.c_str(), &st) == 0 && S_ISREG(st.st_mode);
  if (!table) {
    /* directories of h1 mode files with a "horizon" dataset (K x 2 nf: row n, Re/Im of each
     * field in file order); a later directory overrides an earlier one */
    lmax_ = -1; ncoef_ = 0;
    stringstream dirs(source);
    string dir;
    int nfiles = 0;
    while (getline(dirs, dir, ':')) {
      if (dir.empty()) continue;
      for (const string &file : list_files(dir, regex("h1-l([0-9]+)m([0-9]+).h5"))) {
        const lm_mode lm = filenameToMode(file);
        if (lm.l > lmax) continue;
        hid_t fid = H5Fopen(file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        if (fid < 0) continue;
        if (H5Lexists(fid, "horizon", H5P_DEFAULT) > 0) {
          hid_t did = H5Dopen2(fid, "horizon", H5P_DEFAULT);
          hid_t sid = H5Dget_space(did);
          hsize_t dims[2];
          H5Sget_simple_extent_dims(sid, dims, nullptr);
          vector<double> data(dims[0]*dims[1]);
          H5Dread(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
          H5Sclose(sid); H5Dclose(did);
          const vector<int> fields = file_fields(lm.l, lm.m);
          if (dims[1] != 2*fields.size())
            throw runtime_error("unexpected shape of the horizon dataset in " + file);
          for (size_t k = 0; k < fields.size(); ++k) {
            Poly c(dims[0]);
            for (size_t n = 0; n < dims[0]; ++n) c[n] = cplx(data[n*dims[1] + 2*k], data[n*dims[1] + 2*k + 1]);
            const int ii = fields[k];
            b_[make_tuple(lm.l, lm.m, (ii == 3) ? 6 : ((ii == 6) ? 3 : ii))] = c;
          }
          lmax_ = max(lmax_, lm.l);
          ncoef_ = max(ncoef_, (int) dims[0]);
          ++nfiles;
        }
        H5Fclose(fid);
      }
    }
    /* every (l, m) up to lmax_ must be present */
    for (int l = 0; l <= lmax_; ++l)
      for (int m = 0; m <= l; ++m)
        if (!b_.count(make_tuple(l, m, (l + m) % 2 == 0 ? 1 : 9)))
          throw runtime_error("no horizon coefficients for (l, m) = (" + to_string(l) + ", " + to_string(m) + ") in " + source);
    cout << "Horizon coefficients: l <= " << lmax_ << ", " << ncoef_ << " terms per mode, " << nfiles << " files (" << source << ")" << endl;
    return;
  }
  ifstream in(source);
  if (!in) throw runtime_error("cannot read horizon coefficients " + source);
  string line;
  lmax_ = -1; ncoef_ = 0;
  while (getline(in, line)) {
    istringstream ss(line);
    double l, m, i;
    if (!(ss >> l >> m >> i)) continue;
    Poly c;
    double re, im;
    while (ss >> re >> im) c.push_back(cplx(re, im));
    const int li = (int) lround(l), mi = (int) lround(m), ii = (int) lround(i);
    if (li > lmax) continue;
    const int ipe = (ii == 3) ? 6 : ((ii == 6) ? 3 : ii);   /* file numbering -> PerturbationEquations */
    b_[make_tuple(li, mi, ipe)] = c;
    lmax_ = max(lmax_, li);
    ncoef_ = max(ncoef_, (int) c.size());
  }
  if (lmax_ < 0) throw runtime_error("no horizon coefficients in " + source);
  cout << "Horizon coefficients: l <= " << lmax_ << ", " << ncoef_ << " terms per mode (" << source << ")" << endl;
}

/* Fill the regular-frame fields of all modes (or only the spin perturbation) into A, dA, ddA
 * with N points per plane; fill(i, l, m, polynomial a_i or nullptr, h1, h3, spin coefficient,
 * value/derivative planes) does the evaluation (on the grid or as Taylor coefficients). */
template<class Fill>
void HorizonData::build(Field &A, Field &dA, Field &ddA, size_t N, bool spin_only, Fill fill) const
{
  if (spin_only) {
    A.init(1, N); dA.init(1, N); ddA.init(1, N);
    const double c = 16.0*sqrt(M_PI/3.0);   /* a8 = -c/r, a9 = +c/r */
    for (int i : {8, 9})
      fill(i, 1, 0, nullptr, nullptr, nullptr, i == 8 ? -c : c, A, dA, ddA);
    return;
  }
  A.init(lmax_, N); dA.init(lmax_, N); ddA.init(lmax_, N);
  const Poly rpoly = {2.0*M, 1.0};
  for (int l = 0; l <= lmax_; ++l)
    for (int m = 0; m <= l; ++m) {
      Poly h[11];
      for (int i = 1; i <= 10; ++i) {
        auto it = b_.find(make_tuple(l, m, i));
        if (it != b_.end()) h[i] = it->second;
      }
      stringstream lm; lm << " (l, m) = (" << l << ", " << m << ")";
      const Poly d12 = psub(h[1], h[2]);
      Poly a[11];
      a[2] = psub(pscale(h[3], -1.0), pmul(rpoly, pshift(d12, 1, "h1 - h2" + lm.str())));
      a[3] = pscale(pmul(pmul(rpoly, rpoly), pshift(d12, 2, "h1 - h2" + lm.str())), 2.0);
      a[4] = h[4];
      a[5] = pmul(rpoly, pshift(psub(h[5], h[4]), 1, "h5 - h4" + lm.str()));
      a[6] = h[6]; a[7] = h[7]; a[8] = h[8];
      a[9] = pmul(rpoly, pshift(psub(h[9], h[8]), 1, "h9 - h8" + lm.str()));
      a[10] = h[10];
      for (int i = 1; i <= 10; ++i)
        if (Field::present(i, l, m))
          fill(i, l, m, i == 1 ? nullptr : &a[i], &h[1], &h[3], 0.0, A, dA, ddA);
    }
}

template<bool Series>
void HorizonData::sum(Piece piece, double r0, const vector<double> &r, const vector<ilm_mode> &modes,
                      vector<vector<cplx>> &out) const
{
  const size_t N = Series ? 3 : r.size();        /* points, or Taylor orders */
  const bool spin = (piece == Piece::S1);
  vector<double> f(r.size()), fp(r.size());
  for (size_t j = 0; j < r.size(); ++j) { f[j] = 1.0 - 2.0*M/r[j]; fp[j] = 2.0*M/(r[j]*r[j]); }

  auto set = [&](Field &F, int i, int l, int m, size_t j, cplx v) {
    F.slot_re(i, l, m)[j] = v.real(); F.slot_im(i, l, m)[j] = v.imag();
  };
  /* evaluation of a_i: polynomial; a1 = h1 + f h3; spin field c/r */
  auto fill = [&](int i, int l, int m, const Poly *a, const Poly *h1, const Poly *h3, double c,
                  Field &A, Field &dA, Field &ddA) {
    if (Series) {
      cplx v[3], d[3], dd[3];
      if (a) {
        ptaylor(*a, v, d, dd);
      } else {
        /* Taylor series of the non-polynomial fields to order 4 (enough for a'' to order 2) */
        Poly s(5, 0.0);
        if (h1) {
          Poly fser(5, 0.0);                       /* f = x/(2M + x) */
          for (int n = 1; n < 5; ++n) fser[n] = ((n % 2) ? 1.0 : -1.0)/pow(2.0*M, n);
          Poly fh3 = pmul(fser, *h3);
          for (int n = 0; n < 5; ++n) s[n] = (n < (int) h1->size() ? (*h1)[n] : 0.0) + (n < (int) fh3.size() ? fh3[n] : 0.0);
        } else {
          for (int n = 0; n < 5; ++n) s[n] = c*((n % 2) ? -1.0 : 1.0)/pow(2.0*M, n + 1);   /* c/(2M + x) */
        }
        ptaylor(s, v, d, dd);
      }
      for (int k = 0; k < 3; ++k) { set(A, i, l, m, k, v[k]); set(dA, i, l, m, k, d[k]); set(ddA, i, l, m, k, dd[k]); }
    } else {
      for (size_t j = 0; j < N; ++j) {
        const double x = r[j] - 2.0*M;
        cplx v, d, dd;
        if (a) {
          peval(*a, x, v, d, dd);
        } else if (h1) {
          cplx v1, d1, dd1, v3, d3, dd3;
          peval(*h1, x, v1, d1, dd1); peval(*h3, x, v3, d3, dd3);
          const double ff = f[j], f1 = fp[j], f2 = -4.0*M/(r[j]*r[j]*r[j]);
          v = v1 + ff*v3; d = d1 + f1*v3 + ff*d3; dd = dd1 + f2*v3 + 2.0*f1*d3 + ff*dd3;
        } else {
          v = c/r[j]; d = -c/(r[j]*r[j]); dd = 2.0*c/(r[j]*r[j]*r[j]);
        }
        set(A, i, l, m, j, v); set(dA, i, l, m, j, d); set(ddA, i, l, m, j, dd);
      }
    }
  };
  Field A, dA, ddA, S, dS, ddS;
  build(A, dA, ddA, N, false, fill);
  if (spin) build(S, dS, ddS, N, true, fill);

  out.assign(modes.size(), vector<cplx>(N, 0.0));
  const int num_modes = modes.size();
#pragma omp parallel
  {
    ModeBuffer bX[3], bY[3];
    for (int k = 0; k < 3; ++k) { bX[k].reserve(N); bY[k].reserve(N); }
    ModeBlock X[3], Y[3];
    GBuffer G(R2_NC_MAX, N);
    double c[R2_NC_MAX];
    /* series: the resolved coefficient planes as Taylor3 values (one point) */
    Taylor3 tv[2][3][2][11], tG[2][R2_NC_MAX];
    ModeBlockT<Taylor3> TX[3], TY[3];
    Taylor3 *pGr[R2_NC_MAX], *pGi[R2_NC_MAX];
    for (int k = 0; k < R2_NC_MAX; ++k) { pGr[k] = &tG[0][k]; pGi[k] = &tG[1][k]; }
    Taylor3 rT(2.0*M);                        /* r = 2M + x */
    rT.c[1] = 1.0;
    const Taylor3 fT = 1.0 - 2.0*M*(1.0/rT);
    auto to_taylor = [&](const ModeBlock *B, int w, ModeBlockT<Taylor3> *T) {
      for (int s = 0; s < 3; ++s) {
        T[s].re[0] = T[s].im[0] = nullptr;
        for (int i = 1; i <= 10; ++i) {
          for (int k = 0; k < 3; ++k) { tv[w][s][0][i].c[k] = B[s].re[i][k]; tv[w][s][1][i].c[k] = B[s].im[i][k]; }
          T[s].re[i] = &tv[w][s][0][i]; T[s].im[i] = &tv[w][s][1][i];
        }
      }
    };
    /* add c_k G_k(X, Y) for output (i3, l3, m3) with weight w into acc */
    auto term = [&](const Field *FX, const Field *dFX, const Field *ddFX, int lx, int mx,
                    const Field *FY, const Field *dFY, const Field *ddFY, int ly, int my,
                    int i3, int l3, int m3, cplx *acc) {
      const int K = R2_NC[i3];
      bX[0].resolve(*FX, lx, mx, 0, N, X[0]); bX[1].resolve(*dFX, lx, mx, 0, N, X[1]); bX[2].resolve(*ddFX, lx, mx, 0, N, X[2]);
      bY[0].resolve(*FY, ly, my, 0, N, Y[0]); bY[1].resolve(*dFY, ly, my, 0, N, Y[1]); bY[2].resolve(*ddFY, ly, my, 0, N, Y[2]);
      KC[i3](l3, m3, lx, mx, ly, my, c);
      if (Series) {
        to_taylor(X, 0, TX); to_taylor(Y, 1, TY);
        for (int k = 0; k < K; ++k) tG[0][k] = tG[1][k] = Taylor3(0.0);
        KS[i3](M, r0, &rT, &fT, TX[0], TX[1], TX[2], TY[0], TY[1], TY[2], lx, mx, ly, my, pGr, pGi);
        for (int n = 0; n < 3; ++n) {
          double vr = 0.0, vi = 0.0;
          for (int k = 0; k < K; ++k) { vr += c[k]*tG[0][k].c[n]; vi += c[k]*tG[1][k].c[n]; }
          acc[n] += cplx(vr, vi);
        }
      } else {
        G.zero(K, 0, N);
        KG[i3](M, r0, r, f, fp, X[0], X[1], X[2], Y[0], Y[1], Y[2], lx, mx, ly, my, G.re(), G.im(), 0, N);
        for (size_t j = 0; j < N; ++j) acc[j] += G.combine(K, c, j);
      }
    };
#pragma omp for schedule(dynamic)
    for (int im = 0; im < num_modes; ++im) {
      const int i3 = modes[im].i, l3 = modes[im].l, m3 = modes[im].m;
      cplx *acc = out[im].data();
      if (!spin) {
        for (int l1 = 0; l1 <= lmax_; ++l1)
          for (int m1 = -l1; m1 <= l1; ++m1)
            for (int l2 = abs(l3 - l1); l2 <= min(lmax_, l3 + l1); ++l2) {
              const int m2 = m3 - m1;
              if (abs(m2) > l2) continue;
              term(&A, &dA, &ddA, l1, m1, &A, &dA, &ddA, l2, m2, i3, l3, m3, acc);
            }
      } else {
        /* as Ricci_S: the spin dipole (1, 0) with every h1 mode l2 = l3 - 1 ... l3 + 1, m2 = m3 */
        for (int l2 = max(abs(l3 - 1), 0); l2 <= min(lmax_, l3 + 1); ++l2) {
          if (abs(m3) > l2) continue;
          term(&S, &dS, &ddS, 1, 0, &A, &dA, &ddA, l2, m3, i3, l3, m3, acc);
          term(&A, &dA, &ddA, l2, m3, &S, &dS, &ddS, 1, 0, i3, l3, m3, acc);
        }
      }
    }
  }
}

vector<vector<cplx>> HorizonData::values(Piece piece, double r0, const vector<double> &r,
                                         const vector<ilm_mode> &modes) const
{
  vector<vector<cplx>> out;
  if (!r.empty()) sum<false>(piece, r0, r, modes, out);
  else out.assign(modes.size(), vector<cplx>());
  return out;
}

vector<array<cplx, 3>> HorizonData::coefficients(Piece piece, double r0, const vector<ilm_mode> &modes) const
{
  vector<vector<cplx>> out;
  sum<true>(piece, r0, vector<double>(), modes, out);
  vector<array<cplx, 3>> c(modes.size());
  for (size_t im = 0; im < modes.size(); ++im)
    for (int n = 0; n < 3; ++n) c[im][n] = out[im][n];
  return c;
}
