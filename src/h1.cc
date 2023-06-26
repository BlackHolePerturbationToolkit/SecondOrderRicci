/*
 * Copyright 2014 Barry Wardell
 *
 * This file is distributed under the University of Illinois/NCSA
 * Open Source License. See LICENSE file for details.
 *
 */

#include "h1.h"
#include "utils.h"
#include "h5wrapper.h"
#include <complex>
#include <iostream>
#include <assert.h>
#include <regex>

extern "C"
{
#include <gsl/gsl_sf_ellint.h>
}

using namespace std;
using boost::multi_array;
typedef boost::multi_array_types::extent_range range;
typedef boost::multi_array_types::index_range irange;
typedef boost::multi_array<complex<double>,4> field_type;

const double M = 1.0;

/* Coefficients appearing in the Barack-Lousto basis of tensor harmonics */
double a_il(int i, int l) {
  double a = pow(2, -0.5);
  switch(i) {
    case 1:
    case 2:
    case 3:
    case 6:
      break;
    case 4:
    case 5:
    case 8:
    case 9:
      a *= pow(l*(l+1), -0.5);
      break;
    case 7:
    case 10:
      a *= pow((l-1)*l*(l+1)*(l+2), -0.5);
      break;
  }
  return a;
}

/*
 * Given the field and its first derivative, use the field equations to compute
 * the second derivative. This assumes a circular orbit and requires the
 * orbital radius r0, spacetime mass M, field indices (i, l, m), radial grid
 * points r, f = 1-2M/r and fp = df/dr.
 */
complex<double> d2hdr2(int i, int l, int m, double r, double f, double fp, double r0,
   multi_array<complex<double>, 1> hbar, multi_array<complex<double>, 1> dhbar)
{
  double Omega = sqrt(M)*pow(r0, -1.5);
  complex<double> dt = - complex<double>(0.0, m)*Omega;
  complex<double> ddhbar;
  switch(i) {
    case 1:
      ddhbar = -(fp/f)*dhbar[1] + hbar[1]*(-Omega*Omega*m*m)/(f*f)
               + hbar[1]*((2.0*M)/r + l*(l+1))/(f*r*r)
               + 4.0/(f*f)*(0.5*f*fp*dhbar[1] - 0.5*fp*dt*hbar[2]
               + f*f/(2.0*r*r)*(hbar[1] - f*hbar[3] - hbar[5] - f*hbar[6]));
      break;
    case 2:
      ddhbar = - ((dhbar[2]*fp)/f) + (4.*(-(dt*hbar[1]*fp)/2. + (dhbar[2]*f*fp)/2.
                 + (pow(f,2)*(hbar[2] - hbar[4]))/(2.*pow(r,2))))/pow(f,2)
               + (hbar[2]*((2*M)/pow(r,3) + (l*(1 + l))/pow(r,2)))/f
               + hbar[2]*(-Omega*Omega*m*m)/(f*f);
      break;
    case 3:
      ddhbar = -((dhbar[3]*fp)/f) + (hbar[3]*((2*M)/pow(r,3) + (l*(1 + l))/pow(r,2)))/f
               + hbar[3]*(-Omega*Omega*m*m)/(f*f)
               - (2.*(hbar[1] - hbar[5] - (hbar[3] + hbar[6])*(1 - (4*M)/r)))/(f*pow(r,2));
      break;
    case 4:
      ddhbar = -((dhbar[4]*fp)/f) + (hbar[4]*((2*M)/pow(r,3) + (l*(1 + l))/pow(r,2)))/f
               + (4.*(-(dt*hbar[5]*fp)/4. + (dhbar[4]*f*fp)/4. - (f*hbar[2]*(l*(1. + l)))/(2.*pow(r,2))
               - (3.*f*fp*hbar[4])/(4.*r)))/pow(f,2) + hbar[4]*(-Omega*Omega*m*m)/(f*f);
      break;
    case 5:
      ddhbar = -((dhbar[5]*fp)/f) + (hbar[5]*((2*M)/pow(r,3) + (l*(1 + l))/pow(r,2)))/f
               + (4.*(-(dt*hbar[4]*fp)/4. + (dhbar[5]*f*fp)/4. - (pow(f,2)*hbar[7])/(2.*pow(r,2))
               - (f*(hbar[1] - f*hbar[3] - f*hbar[6])*(l*(1. + l)))/(2.*pow(r,2))
               + (f*hbar[5]*(1 - (7*M)/(2.*r)))/pow(r,2)))/pow(f,2) + hbar[5]*(-Omega*Omega*m*m)/(f*f);
      break;
    case 6:
      ddhbar = -((dhbar[6]*fp)/f) + (hbar[6]*((2*M)/pow(r,3) + (l*(1. + l))/pow(r,2)))/f
               + hbar[6]*(-Omega*Omega*m*m)/(f*f)
               - (2.*(hbar[1] - hbar[5] - (hbar[3] + hbar[6])*(1 - (4*M)/r)))/(f*pow(r,2));
      break;
    case 7:
      ddhbar = - ((dhbar[7]*fp)/f) + (hbar[7]*((2*M)/pow(r,3) + (l*(1. + l))/pow(r,2)))/f
               + hbar[7]*(-Omega*Omega*m*m)/(f*f)
               - (2.*(hbar[7] + hbar[5]*(-1. + l)*(2. + l)))/(f*pow(r,2));
      break;
    case 8:
      ddhbar = - ((dhbar[8]*fp)/f) + (hbar[8]*((2*M)/pow(r,3) + (l*(1. + l))/pow(r,2)))/f
               + (4.*(-(dt*hbar[9]*fp)/4. + (dhbar[8]*f*fp)/4. - (3*f*fp*hbar[8])/(4.*r)))/pow(f,2)
               + hbar[8]*(-Omega*Omega*m*m)/(f*f);
      break;
    case 9:
      ddhbar = - ((dhbar[9]*fp)/f) + (hbar[9]*((2*M)/pow(r,3) + (l*(1. + l))/pow(r,2)))/f
               + (4.*(((-dt*hbar[8] + dhbar[9]*f)*fp)/4.
               + (f*(-(f*hbar[10])/2. + hbar[9]*(1 - (7*M)/(2.*r))))/pow(r,2)))/pow(f,2)
               + hbar[9]*(-Omega*Omega*m*m)/(f*f);
      break;
    case 10:
      ddhbar = - ((dhbar[10]*fp)/f) + (hbar[10]*((2*M)/pow(r,3) + (l*(1. + l))/pow(r,2)))/f
               + hbar[10]*(-Omega*Omega*m*m)/(f*f)
               - (2.*(hbar[10] + hbar[9]*(-1. + l)*(2. + l)))/(f*pow(r,2));
      break;
  }
  return ddhbar;
}

/* Read HDF5 files containing first order fields */
void read_h1(const string dir, double &r0, vector<double> &r, vector<double> &f, vector<double> &fp,
             field_type &h, field_type &dh, field_type &ddh, int l_max)
{
  cout << "Reading first order fields from directory: " << dir << endl;

  /* Determine all available modes */
  std::regex filepattern("h1-l([0-9]+)m([0-9]+).h5");
  vector<string> files = list_files(dir, filepattern);
  if(files.size() == 0)
    exit(1);

  vector<lm_mode> modes;
  for (auto file: files) {
    lm_mode mode = filenameToMode(file);
    if(mode.l > l_max)
      continue; 
    modes.push_back(mode);
  }

  l_max = 0;
  for (auto mode: modes) {
    if (mode.l > l_max)
      l_max = mode.l;
  }

  /* Read the grid structure - assume the same for all fields */
  size_t N;  /* Number of grid points */
  int r0i;   /* Index of worldline point */
  {
    H5F h5_file(files[0], H5F_ACC_RDONLY);
    H5D dataset(h5_file, "grid");
    H5S dataspace(dataset);

    const int rank = dataspace.getSimpleExtentNDims();
    assert(rank == 1);

    vector<hsize_t> gridSize = dataspace.getSimpleExtentDims();
    N = gridSize[0] + 1;
    r.resize(N);
    f.resize(N);
    fp.resize(N);

    H5S memspace(gridSize);

    H5Dread(dataset.getId(), H5T_NATIVE_DOUBLE, memspace.getId(), dataspace.getId(), H5P_DEFAULT, r.data());

    /* Determine the radius of the worldline point */
    H5D dataset_left(h5_file, "inhom_left");
    H5S dataspace_left(dataset_left);
    r0i = dataspace_left.getSimpleExtentDims()[0] - 1;
    r0  = r[r0i];

    /* The actual grid has r0 twice, once for the left solutions and once for the right */
    for(int i=N-1; i>r0i; --i) {
      r[i] = r[i-1];
    }

    /* Compute f and f' on the grid */
    for(size_t i=0; i<N; ++i) {
      f[i] = 1.0 - 2.0*M/r[i];
      fp[i] = 2.0*M/(r[i]*r[i]);
    }
  }

  cout << "Grid size: [" << r.front() << ", " << r.back() << "] (" << N << " points)" << endl;
  cout << "Worldline: r_0 = " << r0 << " (index " << r0i << ")" << endl;
  cout << "Modes: (" << modes.size() << " total (l,m) modes, l_max = " << l_max << ")" << endl;

  /* Read the first order data */
  field_type hbar(boost::extents[range(1,11)][l_max+1][range(-l_max,l_max+1)][N]);
  field_type dhbar(boost::extents[range(1,11)][l_max+1][range(-l_max,l_max+1)][N]);
  field_type ddhbar(boost::extents[range(1,11)][l_max+1][range(-l_max,l_max+1)][N]);
  h.resize(boost::extents[range(1,11)][l_max+1][range(-l_max,l_max+1)][N]);
  dh.resize(boost::extents[range(1,11)][l_max+1][range(-l_max,l_max+1)][N]);
  ddh.resize(boost::extents[range(1,11)][l_max+1][range(-l_max,l_max+1)][N]);

  fill(hbar.data(), hbar.data() + hbar.num_elements(), 0.0);
  fill(dhbar.data(), dhbar.data() + dhbar.num_elements(), 0.0);
  fill(ddhbar.data(), ddhbar.data() + ddhbar.num_elements(), 0.0);
  fill(h.data(), h.data() + h.num_elements(), 0.0);
  fill(dh.data(), dh.data() + dh.num_elements(), 0.0);
  fill(ddh.data(), ddh.data() + ddh.num_elements(), 0.0);

  for(auto file: files)
  {
    lm_mode lm = filenameToMode(file);
    const int l = lm.l;
    const int m = lm.m;

    if(l > l_max)
      continue; 

    cout << "[" << l << ", " << m << "] ";

    H5F h5_file(file, H5F_ACC_RDONLY);
    H5D dataset_left(h5_file, "inhom_left");
    H5D dataset_right(h5_file, "inhom_right");
    H5S dataspace_left(dataset_left);
    H5S dataspace_right(dataset_right);

    /* Check the dataset is a 2D array or the right size */
    const int rank_left = dataspace_left.getSimpleExtentNDims();
    const int rank_right = dataspace_right.getSimpleExtentNDims();
    assert(rank_left == 2);
    assert(rank_right == 2);

    vector<hsize_t> size_left = dataspace_left.getSimpleExtentDims();
    vector<hsize_t> size_right = dataspace_right.getSimpleExtentDims();
    size_t N_left = size_left[0];
    size_t N_right = size_right[0];
    hsize_t num_fields_left = size_left[1];
    hsize_t num_fields_right = size_right[1];

    assert(N_left + N_right == N);
    assert(N_left == (size_t)r0i + 1);
    assert(num_fields_left == num_fields_right);
    hsize_t num_fields = num_fields_left;

    vector<int> fields;
    if (isEven(l+m)) {
      if (l==1) {
        assert(num_fields == 6*6);
        fields = {1, 3, 5, 6, 2, 4};
      } else if (l==0) {
        assert(num_fields == 6*4);
        fields = {1, 3, 6, 2};
      } else {
        assert(num_fields == 6*7);
        fields = {1, 3, 5, 6, 7, 2, 4};
      }
    } else {
      if (l==1) {
        assert(num_fields == 6*2);
        fields = {9, 8};
      } else {
        assert(num_fields == 6*3);
        fields = {9, 10, 8};
      }
    }

    /* Read data for a single dataset */
    H5S memspace_left(size_left);
    H5S memspace_right(size_right);
    multi_array<double, 2> data_left(boost::extents[N_left][num_fields_left]);
    multi_array<double, 2> data_right(boost::extents[N_right][num_fields_right]);
    H5Dread(dataset_left.getId(), H5T_NATIVE_DOUBLE, memspace_left.getId(), dataspace_left.getId(), H5P_DEFAULT, data_left.data());
    H5Dread(dataset_right.getId(), H5T_NATIVE_DOUBLE, memspace_right.getId(), dataspace_right.getId(), H5P_DEFAULT, data_right.data());

    /* Transfer the data to the appropriate locations */

    /* The trace-reversed field and its first derivative.
     * We don't include the a or 1/r factors here.
     */
    for(vector<int>::size_type it=0; it!=fields.size(); ++it) {
      int i = fields[it];
      for(size_t j=0; j<N_left; ++j) {
        hbar[i][l][m][j]  = complex<double>(data_left[j][6*it], data_left[j][6*it+1]);
        dhbar[i][l][m][j] = complex<double>(data_left[j][6*it+2], data_left[j][6*it+3]);
        ddhbar[i][l][m][j] = complex<double>(data_left[j][6*it+4], data_left[j][6*it+5]);
      }
      for(size_t j=0; j<N_right; ++j) {
        hbar[i][l][m][N_left+j]  = complex<double>(data_right[j][6*it], data_right[j][6*it+1]);
        dhbar[i][l][m][N_left+j] = complex<double>(data_right[j][6*it+2], data_right[j][6*it+3]);
        ddhbar[i][l][m][N_left+j] = complex<double>(data_right[j][6*it+4], data_right[j][6*it+5]);
      }
    }

    /* Second derivative of the trace-reversed field */
    // for(vector<int>::size_type it=0; it!=fields.size(); ++it) {
    //   int i = fields[it];
    //   for(size_t j=0; j<N; ++j) {
    //     const double rj = r[j], fj = f[j], fpj = fp[j];
    //     field_type::array_view<1>::type hbarj = hbar[boost::indices[irange()][l][m][j]];
    //     field_type::array_view<1>::type dhbarj = dhbar[boost::indices[irange()][l][m][j]];
    //     hbarj.reindex(1);
    //     dhbarj.reindex(1);
    //     ddhbar[i][l][m][j] = d2hdr2(i, l, m, rj, fj, fpj, r0, hbarj, dhbarj);
    //   }
    // }

    /* The non-trace-reversed field and its first and second derivatives.
     * Here we do include the factor of a and 1/r.
     */
    for(vector<int>::size_type it=0; it!=fields.size(); ++it) {
      int i = fields[it];
      /* Trace-reversal corresponds to swapping 3 and 6 */
      int ih = (i==3) ? 6 : ((i==6) ? 3: i);

      double a = a_il(ih, l);
      for(size_t j=0; j<N; ++j) {
        const complex<double> hbarj = hbar[i][l][m][j];
        const complex<double> dhbarj = dhbar[i][l][m][j];
        const complex<double> ddhbarj = ddhbar[i][l][m][j];
        const double rj = r[j];
        h[ih][l][m][j] = a*hbarj/rj;
        dh[ih][l][m][j] = a*(dhbarj - hbarj/rj)/rj;
        ddh[ih][l][m][j] = a*(ddhbarj - 2.0*(dhbarj - hbarj/rj)/rj)/rj;
      }
    }

    /* Negative m modes */
    double sign = isOdd(m) ? -1.0 : 1.0;
    for(vector<int>::size_type it=0; it!=fields.size(); ++it) {
      int i = fields[it];
      for(size_t j=0; j<N; ++j) {
        h[i][l][-m][j] = sign*conj(h[i][l][m][j]);
        dh[i][l][-m][j] = sign*conj(dh[i][l][m][j]);
        ddh[i][l][-m][j] = sign*conj(ddh[i][l][m][j]);
      }
    }
  }
  cout << endl;
}
