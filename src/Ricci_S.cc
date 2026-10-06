/*
 * Copyright 2014 Barry Wardell
 *
 * This file is distributed under the University of Illinois/NCSA
 * Open Source License. See LICENSE file for details.
 *
 */

#include <vector>
#include <complex>
#include <iostream>
#include <sstream>
#include <limits>
#include <boost/multi_array.hpp>
#ifdef _OPENMP
#  include <omp.h>
#endif
#include "field.h"
#include "modeblock.h"
#include "h1.h"
#include "h1_S.h"
#include "R2G.h"
#include "utils.h"
#include "h5wrapper.h"

using namespace std;
using boost::multi_array;
typedef boost::multi_array_types::extent_range range;

const double M = 1.0;

typedef void (*R2G_function)(const double, const double, const vector<double>&, const vector<double>&,
                             const vector<double>&, const ModeBlock&, const ModeBlock&, const ModeBlock&,
                             const ModeBlock&, const ModeBlock&, const ModeBlock&,
                             int, int, int, int, complex<double> *const*, size_t, size_t);
typedef void (*R2C_function)(int, int, int, int, int, int, double*);
static const R2G_function R2G[11] = {nullptr, R2G_1, R2G_2, R2G_3, R2G_4, R2G_5, R2G_6, R2G_7, R2G_8, R2G_9, R2G_10};
static const R2C_function R2C[11] = {nullptr, R2C_1, R2C_2, R2C_3, R2C_4, R2C_5, R2C_6, R2C_7, R2C_8, R2C_9, R2C_10};

int main(int argc, char* argv[])
{
  cout << "Version: " << __GIT_VERSION << endl;

  string dir;
  int lmax;
  const int dlmax = 1;
  switch(argc) {
    case 2:
      dir = argv[1];
      lmax = std::numeric_limits<int>::max();
      break;
    case 3:
      dir = argv[1];
      lmax = atoi(argv[2]) + dlmax;
      break;
    default:
      cout << "Usage: " << argv[0] << " <dir> [<lmax>]" << endl;
      exit(1);
  }

  int num_threads = 1;
#ifdef _OPENMP
  num_threads = omp_get_max_threads();
  cout << "Running with " << num_threads << " OpenMP threads" << endl;
#endif

  /* Read in first-order fields (B) and the angular momentum perturbation (A) */
  vector<double> r, f, fp;
  double r0;
  int r0i;
  Field hA, dhA, ddhA, hB, dhB, ddhB;
  int h1lmax = read_h1(dir, r0, r0i, r, f, fp, hB, dhB, ddhB, lmax);
  h1_S(r0, r, f, fp, hA, dhA, ddhA);
  lmax = min(lmax, h1lmax);
  const size_t N = r.size();

  /* Compute the source */
  multi_array<complex<double>,4> src(boost::extents[range(1,11)][lmax+1][lmax+1][N]);
  fill(src.data(), src.data() + src.num_elements(), 0.0);

  /* Loop over i3, l3, m3 */
  vector<ilm_mode> modes;
  for(int l3=0; l3<=lmax; ++l3) {
    for(int m3=0; m3<=l3; ++m3) {
      for(int i3=1; i3<=10; ++i3) {
        if(((i3<=7) && isOdd(l3+m3)) || ((i3>7) && isEven(l3+m3)))
          continue;
        if(((i3==4)||(i3==5)||(i3==8)||(i3==9)) && l3<1)
          continue;
        if(((i3==7)||(i3==10)) && l3<2)
          continue;
        modes.push_back(ilm_mode({i3, l3, m3}));
      }
    }
  }

  cout << "Computing source: " << modes.size() << " (i,l,m) modes, l_max = " << lmax << endl;

  int status = 0;
  const int status_frequency = 4*num_threads;
  const int num_modes = modes.size();
#pragma omp parallel
  {
    ModeBuffer bA[3], bB[3];
    for (int k = 0; k < 3; ++k) { bA[k].reserve(N); bB[k].reserve(N); }
    ModeBlock A[3], B[3];
    vector<complex<double>> tmp1(N), tmp2(N);
    vector<complex<double>> Gbuf(static_cast<size_t>(R2_NC_MAX)*N);
    complex<double> *G[R2_NC_MAX];
    for (int k = 0; k < R2_NC_MAX; ++k) G[k] = Gbuf.data() + static_cast<size_t>(k)*N;
    double c[R2_NC_MAX];
    /* d2R_i3 for one ordering of the two modes, written into out */
    auto eval = [&](const ModeBlock *X, int lx, int mx, const ModeBlock *Y, int ly, int my, int i3, int l3, int m3,
                    complex<double> *out) {
      const int K = R2_NC[i3];
      for (int k = 0; k < K; ++k) fill(G[k], G[k] + N, complex<double>(0.0, 0.0));
      R2G[i3](M, r0, r, f, fp, X[0], X[1], X[2], Y[0], Y[1], Y[2], lx, mx, ly, my, G, 0, N);
      R2C[i3](l3, m3, lx, mx, ly, my, c);
      for (size_t j = 0; j < N; ++j) {
        complex<double> v = 0.0;
        for (int k = 0; k < K; ++k) v += c[k]*G[k][j];
        out[j] = v;
      }
    };
#pragma omp for schedule(dynamic)
    for(int im = 0; im < num_modes; ++im) {
      const int i3 = modes[im].i, l3 = modes[im].l, m3 = modes[im].m;
      complex<double> *out = &src[i3][l3][m3][0];

      /* Sum over l2, m1, m2 = m3-m1 with l1 = 1 fixed (the angular momentum perturbation), both
       * orderings d2R[hA, hB] + d2R[hB, hA]. */
      int l1 = 1;
      for(int l2=max(abs(l3-l1),l3-dlmax); l2<=min(min(l3+l1, h1lmax),l3+dlmax); ++l2) {
        for(int m1=0; m1<=0; ++m1) {
          const int m2 = m3 - m1;
          if(abs(m2) > l2)
            continue;
          bA[0].resolve(hA, l1, m1, 0, N, A[0]); bA[1].resolve(dhA, l1, m1, 0, N, A[1]); bA[2].resolve(ddhA, l1, m1, 0, N, A[2]);
          bB[0].resolve(hB, l2, m2, 0, N, B[0]); bB[1].resolve(dhB, l2, m2, 0, N, B[1]); bB[2].resolve(ddhB, l2, m2, 0, N, B[2]);
          eval(A, l1, m1, B, l2, m2, i3, l3, m3, tmp1.data());
          eval(B, l2, m2, A, l1, m1, i3, l3, m3, tmp2.data());
          for(size_t j = 0; j < N; ++j)
            out[j] += tmp1[j] + tmp2[j];
        }
      }

#pragma omp critical
      {
        if(++status % status_frequency == 0)
          cout << "Finished " << status << " of " << num_modes << " modes" << endl;
      }
    }
  }
  if(num_modes % status_frequency != 0)
    cout << "Finished " << status << " of " << num_modes << " modes" << endl;

  /* Output result */
  cout << "Saving results" << endl;
  H5F src_h5("src_S.h5", H5F_ACC_TRUNC);
  for(const ilm_mode &mode : modes) {
    multi_array<double,2> data(boost::extents[N][2]);
    stringstream srcds;
    srcds << "src i=" << mode.i << " l=" << mode.l << " m=" << mode.m;
    for(size_t j=0; j<N; ++j) {
      data[j][0] = real(src[mode.i][mode.l][mode.m][j]);
      data[j][1] = imag(src[mode.i][mode.l][mode.m][j]);
    }
    src_h5.write_dataset(srcds.str(), data);
  }
  multi_array<double,1> rdata(boost::extents[N]);
  for(size_t j=0; j<N; ++j)
    rdata[j] = r[j];
  src_h5.write_dataset("r", rdata);
  cout << "Done" << endl;
  return 0;
}
