/*
 * Copyright 2014 Barry Wardell (original dense-array version)
 *
 * This file is distributed under the University of Illinois/NCSA
 * Open Source License. See LICENSE file for details.
 *
 * Packed-storage reader: same arithmetic as the original read_h1, but no dense
 * temporaries and no explicitly stored negative-m modes.
 */
#include "h1.h"
#include "utils.h"
#include "h5wrapper.h"
#include <complex>
#include <iostream>
#include <cassert>
#include <cmath>
#include <regex>

using namespace std;

static const double M = 1.0;

int read_h1(const string dir, double &r0, int &r0i, vector<double> &r, vector<double> &f,
            vector<double> &fp, Field &h, Field &dh, Field &ddh, int l_max)
{
  cout << "Reading first order fields from directory: " << dir << endl;

  regex filepattern("h1-l([0-9]+)m([0-9]+).h5");
  vector<string> all_files = list_files(dir, filepattern);
  if (all_files.size() == 0)
    exit(1);

  vector<string> files;
  int lmax_found = 0;
  for (auto file: all_files) {
    lm_mode mode = filenameToMode(file);
    if (mode.l > l_max)
      continue;
    files.push_back(file);
    lmax_found = max(lmax_found, mode.l);
  }

  /* Grid: assume the same for all fields. The worldline point appears twice. */
  size_t N;
  {
    H5F h5_file(files[0], H5F_ACC_RDONLY);
    H5D dataset(h5_file, "grid");
    H5S dataspace(dataset);
    assert(dataspace.getSimpleExtentNDims() == 1);
    vector<hsize_t> gridSize = dataspace.getSimpleExtentDims();
    N = gridSize[0] + 1;
    r.resize(N);
    f.resize(N);
    fp.resize(N);
    H5S memspace(gridSize);
    H5Dread(dataset.getId(), H5T_NATIVE_DOUBLE, memspace.getId(), dataspace.getId(), H5P_DEFAULT, r.data());

    H5D dataset_left(h5_file, "inhom_left");
    H5S dataspace_left(dataset_left);
    r0i = dataspace_left.getSimpleExtentDims()[0] - 1;
    r0  = r[r0i];
    for (int i = N-1; i > r0i; --i)
      r[i] = r[i-1];
    for (size_t i = 0; i < N; ++i) {
      f[i] = 1.0 - 2.0*M/r[i];
      fp[i] = 2.0*M/(r[i]*r[i]);
    }
  }

  cout << "Grid size: [" << r.front() << ", " << r.back() << "] (" << N << " points)" << endl;
  cout << "Worldline: r_0 = " << r0 << " (index " << r0i << ")" << endl;
  cout << "Modes: (" << files.size() << " total (l,m) modes, l_max = " << lmax_found << ")" << endl;

  h.init(lmax_found, N);
  dh.init(lmax_found, N);
  ddh.init(lmax_found, N);

  for (auto file: files) {
    lm_mode lm = filenameToMode(file);
    const int l = lm.l;
    const int m = lm.m;

    H5F h5_file(file, H5F_ACC_RDONLY);
    H5D dataset_left(h5_file, "inhom_left");
    H5D dataset_right(h5_file, "inhom_right");
    H5S dataspace_left(dataset_left);
    H5S dataspace_right(dataset_right);
    assert(dataspace_left.getSimpleExtentNDims() == 2);
    assert(dataspace_right.getSimpleExtentNDims() == 2);
    vector<hsize_t> size_left = dataspace_left.getSimpleExtentDims();
    vector<hsize_t> size_right = dataspace_right.getSimpleExtentDims();
    const size_t N_left = size_left[0];
    const size_t N_right = size_right[0];
    assert(N_left + N_right == N);
    assert(N_left == (size_t)r0i + 1);
    assert(size_left[1] == size_right[1]);
    const size_t num_fields = size_left[1];

    vector<int> fields;
    if (isEven(l+m)) {
      if (l == 1) { assert(num_fields == 6*6); fields = {1, 3, 5, 6, 2, 4}; }
      else if (l == 0) { assert(num_fields == 6*4); fields = {1, 3, 6, 2}; }
      else { assert(num_fields == 6*7); fields = {1, 3, 5, 6, 7, 2, 4}; }
    } else {
      if (l == 1) { assert(num_fields == 6*2); fields = {9, 8}; }
      else { assert(num_fields == 6*3); fields = {9, 10, 8}; }
    }

    vector<double> data_left(N_left*num_fields), data_right(N_right*num_fields);
    H5S memspace_left(size_left);
    H5S memspace_right(size_right);
    H5Dread(dataset_left.getId(), H5T_NATIVE_DOUBLE, memspace_left.getId(), dataspace_left.getId(), H5P_DEFAULT, data_left.data());
    H5Dread(dataset_right.getId(), H5T_NATIVE_DOUBLE, memspace_right.getId(), dataspace_right.getId(), H5P_DEFAULT, data_right.data());

    /* Field values and their first two r-derivatives as stored in the files (the kernels
     * generated from PerturbationEquations work with these directly), with components 3
     * and 6 swapped, as in the dense-array reader. */
    for (size_t it = 0; it != fields.size(); ++it) {
      const int i = fields[it];
      const int ih = (i == 3) ? 6 : ((i == 6) ? 3 : i);
      double *hr = h.slot_re(ih, l, m), *hi = h.slot_im(ih, l, m);
      double *dhr = dh.slot_re(ih, l, m), *dhi = dh.slot_im(ih, l, m);
      double *ddhr = ddh.slot_re(ih, l, m), *ddhi = ddh.slot_im(ih, l, m);
      for (size_t j = 0; j < N; ++j) {
        const double *row = (j < N_left) ? &data_left[j*num_fields] : &data_right[(j - N_left)*num_fields];
        hr[j] = row[6*it];     hi[j] = row[6*it+1];
        dhr[j] = row[6*it+2];  dhi[j] = row[6*it+3];
        ddhr[j] = row[6*it+4]; ddhi[j] = row[6*it+5];
      }
    }
  }
  cout << "Stored " << (h.bytes() + dh.bytes() + ddh.bytes())/1.0e9 << " GB of first-order data" << endl;
  return lmax_found;
}
