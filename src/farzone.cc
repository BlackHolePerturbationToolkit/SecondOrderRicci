#include "farzone.h"
#include "field.h"
#include "modeblock.h"
#include "yseries.h"
#include "taylor.h"
#include "R2U.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <regex>
#include <sys/stat.h>
#include <hdf5.h>
#ifdef _OPENMP
#  include <omp.h>
#endif

using namespace std;
typedef complex<double> cplx;

static const double M = 1.0;

typedef void (*R2U_grid)(const double, const double, const vector<double>&, const vector<double>&,
                         const vector<double>&, const ModeBlock&, const ModeBlock&, const ModeBlock&,
                         const ModeBlock&, const ModeBlock&, const ModeBlock&,
                         int, int, int, int, double *const*, double *const*, size_t, size_t);
typedef void (*R2U_series)(const double, const double, const InfR*, const YLog*,
                           const ModeBlockT<YLog>&, const ModeBlockT<YLog>&, const ModeBlockT<YLog>&,
                           const ModeBlockT<YLog>&, const ModeBlockT<YLog>&, const ModeBlockT<YLog>&,
                           int, int, int, int, YLog *const*, YLog *const*);
typedef void (*R2C_function)(int, int, int, int, int, int, double*);
static const R2U_grid KG[11] = {nullptr, R2U_1, R2U_2, R2U_3, R2U_4, R2U_5, R2U_6, R2U_7, R2U_8, R2U_9, R2U_10};
static const R2U_series KS[11] = {nullptr, R2U_1, R2U_2, R2U_3, R2U_4, R2U_5, R2U_6, R2U_7, R2U_8, R2U_9, R2U_10};
static const R2C_function KC[11] = {nullptr, R2C_1, R2C_2, R2C_3, R2C_4, R2C_5, R2C_6, R2C_7, R2C_8, R2C_9, R2C_10};

/* number of doubles per point in the series planes (YLog coefficients) */
static const size_t NY = YLOG_N*(YLOG_P + 1);

/* --- value and first two r-derivatives ----------------------------------------------- */
template<class S> struct Jet { S v, d, dd; };
template<class S> static Jet<S> jmul(const Jet<S> &a, const Jet<S> &b) {
  return {a.v*b.v, a.d*b.v + a.v*b.d, a.dd*b.v + 2.0*(a.d*b.d) + a.v*b.dd};
}
template<class S> static Jet<S> jadd(const Jet<S> &a, const Jet<S> &b) { return {a.v + b.v, a.d + b.d, a.dd + b.dd}; }
template<class S> static Jet<S> jscale(const Jet<S> &a, double s) { return {s*a.v, s*a.d, s*a.dd}; }
/* complex jets as real and imaginary parts */
template<class S> struct CJet { Jet<S> re, im; };
template<class S> static CJet<S> cmul(const Jet<S> &g, const CJet<S> &h) { return {jmul(g, h.re), jmul(g, h.im)}; }
template<class S> static CJet<S> cadd(const CJet<S> &a, const CJet<S> &b) { return {jadd(a.re, b.re), jadd(a.im, b.im)}; }
template<class S> static CJet<S> cscale(const CJet<S> &a, double s) { return {jscale(a.re, s), jscale(a.im, s)}; }

/* h = sum_n (b_n + c_n log r)/r^n at r (values) */
static CJet<double> field_at(const vector<cplx> &b, const vector<cplx> &c, double r)
{
  const double L = log(r);
  cplx v = 0.0, d = 0.0, dd = 0.0;
  for (size_t n = 0; n < max(b.size(), c.size()); ++n) {
    const cplx bn = n < b.size() ? b[n] : 0.0, cn = n < c.size() ? c[n] : 0.0;
    const cplx q = bn + cn*L;
    const double p0 = pow(r, -double(n)), p1 = p0/r, p2 = p1/r;
    v += q*p0;
    d += -double(n)*q*p1 + cn*p1;
    dd += double(n*(n + 1))*q*p2 - double(2*n + 1)*cn*p2;
  }
  return {{v.real(), d.real(), dd.real()}, {v.imag(), d.imag(), dd.imag()}};
}
/* the same as series in y = 1/r and log r */
static CJet<YLog> field_series(const vector<cplx> &b, const vector<cplx> &c)
{
  CJet<YLog> h;
  for (size_t n = 0; n < (size_t) YLOG_N; ++n) {
    const cplx bn = n < b.size() ? b[n] : 0.0, cn = n < c.size() ? c[n] : 0.0;
    for (int part = 0; part < 2; ++part) {
      Jet<YLog> &J = part ? h.im : h.re;
      const double B = part ? bn.imag() : bn.real(), C = part ? cn.imag() : cn.real();
      J.v.c[n][0] += B; J.v.c[n][1] += C;
      if (n + 1 < (size_t) YLOG_N) { J.d.c[n+1][0] += -double(n)*B + C; J.d.c[n+1][1] += -double(n)*C; }
      if (n + 2 < (size_t) YLOG_N) { J.dd.c[n+2][0] += double(n*(n + 1))*B - double(2*n + 1)*C; J.dd.c[n+2][1] += double(n*(n + 1))*C; }
    }
  }
  return h;
}
/* 1/r, f = 1 - 2M/r, 1/f, 1/f^2 as jets */
static void radial_jets(double r, Jet<double> &rinv, Jet<double> &f, Jet<double> &g, Jet<double> &g2)
{
  rinv = {1.0/r, -1.0/(r*r), 2.0/(r*r*r)};
  f = {1.0 - 2.0*M/r, 2.0*M/(r*r), -4.0*M/(r*r*r)};
  g = {1.0/f.v, -f.d/(f.v*f.v), (2.0*f.d*f.d - f.v*f.dd)/(f.v*f.v*f.v)};
  g2 = jmul(g, g);
}
static void radial_jets(Jet<YLog> &rinv, Jet<YLog> &f, Jet<YLog> &g, Jet<YLog> &g2)
{
  rinv = Jet<YLog>{}; f = Jet<YLog>{}; g = Jet<YLog>{};
  /* d/dr y^n = -n y^(n+1) */
  rinv.v.c[1][0] = 1.0; rinv.d.c[2][0] = -1.0; rinv.dd.c[3][0] = 2.0;
  f.v.c[0][0] = 1.0; f.v.c[1][0] = -2.0*M; f.d.c[2][0] = 2.0*M; f.dd.c[3][0] = -4.0*M;
  for (int n = 0; n < YLOG_N; ++n) {                       /* 1/f = sum (2My)^n */
    const double cn = pow(2.0*M, n);
    g.v.c[n][0] = cn;
    if (n + 1 < YLOG_N) g.d.c[n+1][0] = -n*cn;
    if (n + 2 < YLOG_N) g.dd.c[n+2][0] = n*(n + 1.0)*cn;
  }
  g2 = jmul(g, g);
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

/* read a K x 2nf dataset into per-field coefficient lists; false if absent */
static bool read_coeffs(hid_t fid, const char *name, vector<vector<cplx>> &out, size_t nf, const string &file)
{
  if (H5Lexists(fid, name, H5P_DEFAULT) <= 0) return false;
  hid_t did = H5Dopen2(fid, name, H5P_DEFAULT);
  hid_t sid = H5Dget_space(did);
  hsize_t dims[2];
  H5Sget_simple_extent_dims(sid, dims, nullptr);
  vector<double> data(dims[0]*dims[1]);
  H5Dread(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
  H5Sclose(sid); H5Dclose(did);
  if (dims[1] != 2*nf) throw runtime_error(string("unexpected shape of ") + name + " in " + file);
  out.assign(nf, vector<cplx>(dims[0]));
  for (size_t k = 0; k < nf; ++k)
    for (size_t n = 0; n < dims[0]; ++n) out[k][n] = cplx(data[n*dims[1] + 2*k], data[n*dims[1] + 2*k + 1]);
  return true;
}

void InfinityData::load(const string &source, int lmax)
{
  struct stat st;
  const string first = source.substr(0, source.find(':'));
  if (stat(first.c_str(), &st) == 0 && S_ISDIR(st.st_mode)) {
    /* directories of h1 mode files with datasets "infinity" and "infinity_log" (K x 2 nf: row n,
     * Re/Im of each field in file order); a later directory overrides an earlier one */
    lmax_ = -1;
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
        const vector<int> fields = file_fields(lm.l, lm.m);
        vector<vector<cplx>> bb, cc;
        if (read_coeffs(fid, "infinity", bb, fields.size(), file)) {
          const bool logs = read_coeffs(fid, "infinity_log", cc, fields.size(), file);
          for (size_t k = 0; k < fields.size(); ++k) {
            const int ii = fields[k], ipe = (ii == 3) ? 6 : ((ii == 6) ? 3 : ii);
            b_[make_tuple(lm.l, lm.m, ipe)] = bb[k];
            if (logs) c_[make_tuple(lm.l, lm.m, ipe)] = cc[k];
            else c_.erase(make_tuple(lm.l, lm.m, ipe));
          }
          lmax_ = max(lmax_, lm.l);
          ++nfiles;
        }
        H5Fclose(fid);
      }
    }
    for (int l = 0; l <= lmax_; ++l)
      for (int m = 0; m <= l; ++m)
        if (!b_.count(make_tuple(l, m, (l + m) % 2 == 0 ? 1 : 9)))
          throw runtime_error("no coefficients at infinity for (l, m) = (" + to_string(l) + ", " + to_string(m) + ") in " + source);
    cout << "Coefficients at infinity: l <= " << lmax_ << ", " << nfiles << " files (" << source << ")" << endl;
    return;
  }
  const size_t colon = source.find(':');
  const string files[2] = {source.substr(0, colon), colon == string::npos ? string() : source.substr(colon + 1)};
  lmax_ = -1;
  for (int k = 0; k < 2; ++k) {
    if (files[k].empty()) continue;
    ifstream in(files[k]);
    if (!in) throw runtime_error("cannot read the expansion coefficients at infinity " + files[k]);
    string line;
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
      (k == 0 ? b_ : c_)[make_tuple(li, mi, ipe)] = c;
      if (k == 0) lmax_ = max(lmax_, li);
    }
  }
  if (lmax_ < 0) throw runtime_error("no expansion coefficients at infinity in " + source);
  cout << "Coefficients at infinity: l <= " << lmax_ << " (" << source << ")" << endl;
}

template<bool Series>
void InfinityData::sum(Piece piece, double r0, const vector<double> &r, const vector<ilm_mode> &modes,
                       vector<vector<cplx>> &out) const
{
  typedef typename conditional<Series, YLog, double>::type S;
  const size_t N = Series ? NY : r.size();
  const size_t npts = Series ? 1 : r.size();
  const bool spin = (piece == Piece::S1);
  vector<double> f(r.size()), fp(r.size());
  for (size_t j = 0; j < r.size(); ++j) { f[j] = 1.0 - 2.0*M/r[j]; fp[j] = 2.0*M/(r[j]*r[j]); }

  /* store a jet component at point j (series: the YLog coefficients as NY consecutive values) */
  auto put = [&](double *pv, size_t j, const S &v) {
    if (Series) { const YLog &y = reinterpret_cast<const YLog&>(v); for (size_t q = 0; q < NY; ++q) pv[q] = (&y.c[0][0])[q]; }
    else pv[j] = reinterpret_cast<const double&>(v);
  };
  auto store = [&](Field &A, Field &dA, Field &ddA, int i, int l, int m, size_t j, const CJet<S> &a) {
    put(A.slot_re(i, l, m), j, a.re.v); put(A.slot_im(i, l, m), j, a.im.v);
    put(dA.slot_re(i, l, m), j, a.re.d); put(dA.slot_im(i, l, m), j, a.im.d);
    put(ddA.slot_re(i, l, m), j, a.re.dd); put(ddA.slot_im(i, l, m), j, a.im.dd);
  };
  /* (u, r)-frame fields from the BLS fields h[1..10] */
  auto frame = [&](const CJet<S> *h, const Jet<S> &fj, const Jet<S> &g, const Jet<S> &g2, CJet<S> *a) {
    const CJet<S> h12 = cadd(h[1], h[2]);
    a[1] = cadd(h[1], cmul(fj, h[3]));
    a[2] = cadd(h[3], cmul(g, h12));
    a[3] = cscale(cmul(g2, h12), 2.0);
    a[4] = h[4]; a[5] = cmul(g, cadd(h[5], h[4]));
    a[6] = h[6]; a[7] = h[7]; a[8] = h[8];
    a[9] = cmul(g, cadd(h[9], h[8]));
    a[10] = h[10];
  };
  auto jets_at = [&](size_t j, Jet<S> &ri, Jet<S> &fj, Jet<S> &g, Jet<S> &g2) {
    if (Series) radial_jets(reinterpret_cast<Jet<YLog>&>(ri), reinterpret_cast<Jet<YLog>&>(fj),
                            reinterpret_cast<Jet<YLog>&>(g), reinterpret_cast<Jet<YLog>&>(g2));
    else radial_jets(r[j], reinterpret_cast<Jet<double>&>(ri), reinterpret_cast<Jet<double>&>(fj),
                     reinterpret_cast<Jet<double>&>(g), reinterpret_cast<Jet<double>&>(g2));
  };
  auto hfield = [&](const Poly &b, const Poly &c, size_t j) -> CJet<S> {
    if (Series) { CJet<YLog> h = field_series(b, c); return reinterpret_cast<CJet<S>&>(h); }
    CJet<double> h = field_at(b, c, r[j]); return reinterpret_cast<CJet<S>&>(h);
  };

  Field A, dA, ddA, Sf, dSf, ddSf;
  A.init(lmax_, N); dA.init(lmax_, N); ddA.init(lmax_, N);
  const Poly none;
  for (int l = 0; l <= lmax_; ++l)
    for (int m = 0; m <= l; ++m)
      for (size_t j = 0; j < npts; ++j) {
        Jet<S> ri, fj, g, g2;
        jets_at(j, ri, fj, g, g2);
        CJet<S> h[11], a[11];
        for (int i = 1; i <= 10; ++i) {
          auto ib = b_.find(make_tuple(l, m, i)), ic = c_.find(make_tuple(l, m, i));
          h[i] = hfield(ib != b_.end() ? ib->second : none, ic != c_.end() ? ic->second : none, j);
        }
        frame(h, fj, g, g2, a);
        for (int i = 1; i <= 10; ++i)
          if (Field::present(i, l, m)) store(A, dA, ddA, i, l, m, j, a[i]);
      }
  if (spin) {
    /* spin perturbation (l, m) = (1, 0): h8 = -c/r, h9 = -2Mc/r^2 */
    Sf.init(1, N); dSf.init(1, N); ddSf.init(1, N);
    const double c = 16.0*sqrt(M_PI/3.0);
    for (size_t j = 0; j < npts; ++j) {
      Jet<S> ri, fj, g, g2;
      jets_at(j, ri, fj, g, g2);
      CJet<S> h[11], a[11];
      for (int i = 1; i <= 10; ++i) h[i] = CJet<S>{};
      h[8].re = jscale(ri, -c);
      h[9].re = jscale(jmul(ri, ri), -2.0*M*c);
      frame(h, fj, g, g2, a);
      store(Sf, dSf, ddSf, 8, 1, 0, j, a[8]);
      store(Sf, dSf, ddSf, 9, 1, 0, j, a[9]);
    }
  }

  out.assign(modes.size(), vector<cplx>(N, 0.0));
  const int num_modes = modes.size();
#pragma omp parallel
  {
    ModeBuffer bX[3], bY[3];
    for (int k = 0; k < 3; ++k) { bX[k].reserve(N); bY[k].reserve(N); }
    ModeBlock X[3], Y[3];
    GBuffer G(R2_NC_MAX, N);
    double c[R2_NC_MAX];
    YLog sv[2][3][2][11], sG[2][R2_NC_MAX];
    ModeBlockT<YLog> SX[3], SY[3];
    YLog *pGr[R2_NC_MAX], *pGi[R2_NC_MAX];
    for (int k = 0; k < R2_NC_MAX; ++k) { pGr[k] = &sG[0][k]; pGi[k] = &sG[1][k]; }
    const InfR rS;
    YLog fS(1.0); fS.c[1][0] = -2.0*M;
    auto to_series = [&](const ModeBlock *B, int w, ModeBlockT<YLog> *T) {
      for (int s = 0; s < 3; ++s) {
        T[s].re[0] = T[s].im[0] = nullptr;
        for (int i = 1; i <= 10; ++i) {
          for (size_t q = 0; q < NY; ++q) { (&sv[w][s][0][i].c[0][0])[q] = B[s].re[i][q]; (&sv[w][s][1][i].c[0][0])[q] = B[s].im[i][q]; }
          T[s].re[i] = &sv[w][s][0][i]; T[s].im[i] = &sv[w][s][1][i];
        }
      }
    };
    auto term = [&](const Field *FX, const Field *dFX, const Field *ddFX, int lx, int mx,
                    const Field *FY, const Field *dFY, const Field *ddFY, int ly, int my,
                    int i3, int l3, int m3, cplx *acc) {
      const int K = R2_NC[i3];
      bX[0].resolve(*FX, lx, mx, 0, N, X[0]); bX[1].resolve(*dFX, lx, mx, 0, N, X[1]); bX[2].resolve(*ddFX, lx, mx, 0, N, X[2]);
      bY[0].resolve(*FY, ly, my, 0, N, Y[0]); bY[1].resolve(*dFY, ly, my, 0, N, Y[1]); bY[2].resolve(*ddFY, ly, my, 0, N, Y[2]);
      KC[i3](l3, m3, lx, mx, ly, my, c);
      if (Series) {
        to_series(X, 0, SX); to_series(Y, 1, SY);
        for (int k = 0; k < K; ++k) sG[0][k] = sG[1][k] = YLog(0.0);
        KS[i3](M, r0, &rS, &fS, SX[0], SX[1], SX[2], SY[0], SY[1], SY[2], lx, mx, ly, my, pGr, pGi);
        for (size_t q = 0; q < NY; ++q) {
          double vr = 0.0, vi = 0.0;
          for (int k = 0; k < K; ++k) { vr += c[k]*(&sG[0][k].c[0][0])[q]; vi += c[k]*(&sG[1][k].c[0][0])[q]; }
          acc[q] += cplx(vr, vi);
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
        for (int l2 = max(abs(l3 - 1), 0); l2 <= min(lmax_, l3 + 1); ++l2) {
          if (abs(m3) > l2) continue;
          term(&Sf, &dSf, &ddSf, 1, 0, &A, &dA, &ddA, l2, m3, i3, l3, m3, acc);
          term(&A, &dA, &ddA, l2, m3, &Sf, &dSf, &ddSf, 1, 0, i3, l3, m3, acc);
        }
      }
    }
  }
}

vector<vector<cplx>> InfinityData::values(Piece piece, double r0, const vector<double> &r,
                                          const vector<ilm_mode> &modes) const
{
  vector<vector<cplx>> out;
  if (!r.empty()) sum<false>(piece, r0, r, modes, out);
  else out.assign(modes.size(), vector<cplx>());
  return out;
}

vector<array<cplx, 4>> InfinityData::coefficients(Piece piece, double r0, const vector<ilm_mode> &modes, double *residual) const
{
  vector<vector<cplx>> out;
  sum<true>(piece, r0, vector<double>(), modes, out);
  vector<array<cplx, 4>> A(modes.size());
  double res = 0.0, big = 0.0;
  auto at = [&](size_t im, int n, int p) { return out[im][n*(YLOG_P + 1) + p]; };
  for (size_t im = 0; im < modes.size(); ++im) {
    A[im] = {at(im, 2, 0), at(im, 3, 0), at(im, 4, 0), at(im, 4, 1)};
    for (int n = 0; n < YLOG_N; ++n)
      for (int p = 0; p <= YLOG_P; ++p) {
        const bool kept = (p == 0 && n >= 2) || (n == 4 && p == 1);
        (kept ? big : res) = max(kept ? big : res, abs(at(im, n, p)));
      }
  }
  if (residual) *residual = big > 0 ? res/big : res;
  return A;
}

double InfinityData::attach_radius(double tol) const
{
  /* smallest r beyond which, for every mode, the last coefficient kept, |b_K-1| + |c_K-1| log r
   * over r^(K-1), is below tol times the largest leading coefficient of any radiative mode */
  double big = 0.0;
  for (const auto &kv : b_)
    if (get<1>(kv.first) > 0)
      for (size_t n = 0; n < min<size_t>(3, kv.second.size()); ++n) big = max(big, abs(kv.second[n]));
  double rmin = 0.0;
  for (const auto &kv : b_) {
    const Poly &b = kv.second;
    if (b.size() < 2) continue;
    auto it = c_.find(kv.first);
    const double bl = abs(b.back()), cl = (it != c_.end() && !it->second.empty()) ? abs(it->second.back()) : 0.0;
    const double p = b.size() - 1.0;
    double lo = 10.0, hi = 1e8;
    auto bad = [&](double r) { return (bl + cl*log(r))/pow(r, p) > tol*big; };
    if (!bad(lo)) continue;
    for (int it2 = 0; it2 < 100; ++it2) { const double mid = sqrt(lo*hi); (bad(mid) ? lo : hi) = mid; }
    rmin = max(rmin, hi);
  }
  return rmin;
}

