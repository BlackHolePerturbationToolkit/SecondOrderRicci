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
               + hbar[1]*((2.0*M)/r0 + l*(l+1))/(f*r0*r0)
               + 4.0/(f*f)*(0.5*f*fp*dhbar[1] - 0.5*fp*dt*hbar[2]
               + f*f/(2.0*r0*r0)*(hbar[1] - f*hbar[3] - hbar[5] - f*hbar[6]));
      break;
    case 2:
      ddhbar = - ((dhbar[2]*fp)/f) + (4.*(-(dt*hbar[1]*fp)/2. + (dhbar[2]*f*fp)/2.
                 + (pow(f,2)*(hbar[2] - hbar[4]))/(2.*pow(r0,2))))/pow(f,2)
               + (hbar[2]*((2*M)/pow(r0,3) + (l*(1 + l))/pow(r0,2)))/f
               - (hbar[2]*pow(m,2)*M)/(pow(f,2)*pow(r0,3));
      break;
    case 3:
      ddhbar = -((dhbar[3]*fp)/f) + (hbar[3]*((2*M)/pow(r0,3) + (l*(1 + l))/pow(r0,2)))/f
               -(hbar[3]*pow(m,2)*M)/(pow(f,2)*pow(r0,3))
               - (2.*(hbar[1] - hbar[5] - (hbar[3] + hbar[6])*(1 - (4*M)/r0)))/(f*pow(r0,2));
      break;
    case 4:
      ddhbar = -((dhbar[4]*fp)/f) + (hbar[4]*((2*M)/pow(r0,3) + (l*(1 + l))/pow(r0,2)))/f
               + (4.*(-(dt*hbar[5]*fp)/4. + (dhbar[4]*f*fp)/4. - (f*hbar[2]*(l*(1. + l)))/(2.*pow(r0,2))
               - (3.*f*fp*hbar[4])/(4.*r0)))/pow(f,2) - (hbar[4]*pow(m,2)*M)/(pow(f,2)*pow(r0,3));
      break;
    case 5:
      ddhbar = -((dhbar[5]*fp)/f) + (hbar[5]*((2*M)/pow(r0,3) + (l*(1 + l))/pow(r0,2)))/f
               + (4.*(-(dt*hbar[4]*fp)/4. + (dhbar[5]*f*fp)/4. - (pow(f,2)*hbar[7])/(2.*pow(r0,2))
               - (f*(hbar[1] - f*hbar[3] - f*hbar[6])*(l*(1. + l)))/(2.*pow(r0,2))
               + (f*hbar[5]*(1 - (7*M)/(2.*r0)))/pow(r0,2)))/pow(f,2) - (hbar[5]*pow(m,2)*M)/(pow(f,2)*pow(r0,3));
      break;
    case 6:
      ddhbar = -((dhbar[6]*fp)/f) + (hbar[6]*((2*M)/pow(r0,3) + (l*(1. + l))/pow(r0,2)))/f
               - (hbar[6]*pow(m,2)*M)/(pow(f,2)*pow(r0,3))
               - (2.*(hbar[1] - hbar[5] - (hbar[3] + hbar[6])*(1 - (4*M)/r0)))/(f*pow(r0,2));
      break;
    case 7:
      ddhbar = - ((dhbar[7]*fp)/f) + (hbar[7]*((2*M)/pow(r0,3) + (l*(1. + l))/pow(r0,2)))/f
               - (hbar[7]*pow(m,2)*M)/(pow(f,2)*pow(r0,3))
               - (2.*(hbar[7] + hbar[5]*(-1. + l)*(2. + l)))/(f*pow(r0,2));
      break;
    case 8:
      ddhbar = - ((dhbar[8]*fp)/f) + (hbar[8]*((2*M)/pow(r0,3) + (l*(1. + l))/pow(r0,2)))/f
               + (4.*(-(dt*hbar[9]*fp)/4. + (dhbar[8]*f*fp)/4. - (3*f*fp*hbar[8])/(4.*r0)))/pow(f,2)
               - (hbar[8]*pow(m,2)*M)/(pow(f,2)*pow(r0,3));
      break;
    case 9:
      ddhbar = - ((dhbar[9]*fp)/f) + (hbar[9]*((2*M)/pow(r0,3) + (l*(1. + l))/pow(r0,2)))/f
               + (4.*(((-dt*hbar[8] + dhbar[9]*f)*fp)/4.
               + (f*(-(f*hbar[10])/2. + hbar[9]*(1 - (7*M)/(2.*r0))))/pow(r0,2)))/pow(f,2)
               - (hbar[9]*pow(m,2)*M)/(pow(f,2)*pow(r0,3));
      break;
    case 10:
      ddhbar = - ((dhbar[10]*fp)/f) + (hbar[10]*((2*M)/pow(r0,3) + (l*(1. + l))/pow(r0,2)))/f
               - (hbar[10]*pow(m,2)*M)/(pow(f,2)*pow(r0,3))
               - (2.*(hbar[10] + hbar[9]*(-1. + l)*(2. + l)))/(f*pow(r0,2));
      break;
  }
  return ddhbar;
}

/* Read HDF5 files containing first order fields */
void read_h1(vector<double> &r, vector<double> &f, vector<double> &fp,
             field_type &h, field_type &dh, field_type &ddh)
{
  cout << "Reading first order fields" << endl;

  /* Determine all available modes */
  vector<string> files = list_files("../h1_fields");

  vector<lm_mode> modes;
  for (auto file: files) {
    modes.push_back(filenameToMode(file));
  }

  int l_max = 0;
  for (auto mode: modes) {
    if (mode.l > l_max)
      l_max = mode.l;
  }

  /* Read the grid structure - assume the same for all fields */
  int N;     /* Number of grid points */
  int r0i;   /* Index of worldline point */
  double r0; /* Radius of the worldline point */
  {
    H5F h5_file(files[0], H5F_ACC_RDONLY);
    H5D dataset(h5_file, "grid");
    H5S dataspace(dataset);

    const int rank = dataspace.getSimpleExtentNDims();
    assert(rank == 1);

    vector<hsize_t> gridSize = dataspace.getSimpleExtentDims();
    N = gridSize[0];
    r.resize(N);
    f.resize(N);
    fp.resize(N);

    H5S memspace(gridSize);

    H5Dread(dataset.getId(), H5T_NATIVE_DOUBLE, memspace.getId(), dataspace.getId(), H5P_DEFAULT, r.data());
    for(int i=0; i<N; ++i) {
      f[i] = 1.0 - 2.0*M/r[i];
      fp[i] = 2.0*M/(r[i]*r[i]);
    }

    /* Determine the radius of the worldline point */
    H5A r0index(dataset.getId(), "r0index");

    H5S r0_dataspace(r0index);
    assert(r0_dataspace.getSimpleExtentNDims() == 1);

    vector<hsize_t> r0_dims = r0_dataspace.getSimpleExtentDims();
    assert(r0_dims[0] == 1);

    H5Aread(r0index.getId(), H5T_NATIVE_INT, &r0i);
    r0 = r[r0i];
  }

  cout << "Grid size: [" << r.front() << ", " << r.back() << "] (" << N << " points)" << endl;
  cout << "Worldline: r_0 = " << r0 << " (index " << r0i << ")" << endl;
  cout << "Modes: (" << modes.size() << " total, l_max = " << l_max << ")" << endl;

  /* Read the first order data */
  h.resize(boost::extents[range(1,11)][l_max+1][range(-l_max,l_max+1)][N]);
  dh.resize(boost::extents[range(1,11)][l_max+1][range(-l_max,l_max+1)][N]);
  ddh.resize(boost::extents[range(1,11)][l_max+1][range(-l_max,l_max+1)][N]);

  fill(h.data(), h.data() + h.num_elements(), 0.0);
  fill(dh.data(), dh.data() + dh.num_elements(), 0.0);
  fill(ddh.data(), ddh.data() + ddh.num_elements(), 0.0);

  for(auto file: files)
  {
    lm_mode lm = filenameToMode(file);
    const int l = lm.l;
    const int m = lm.m;
    cout << "[" << l << ", " << m << "] ";

    H5F h5_file(file, H5F_ACC_RDONLY);
    H5D dataset(h5_file, "inhom");
    H5S dataspace(dataset);

    /* Check the dataset is a 2D array or the right size */
    const int rank = dataspace.getSimpleExtentNDims();
    assert(rank == 2);

    vector<hsize_t> size = dataspace.getSimpleExtentDims();
    assert(size[0] == N);
    hsize_t num_fields = size[1];

    vector<int> fields;
    if (isEven(l+m)) {
      if (l==1) {
        assert(num_fields == 4*6);
        fields = {1, 3, 5, 6, 2, 4};
      } else {
        assert(num_fields == 4*7);
        fields = {1, 3, 5, 6, 7, 2, 4};
      }
    } else {
      if (l==1) {
        assert(num_fields == 4*2);
        fields = {9, 8};
      } else {
        assert(num_fields == 4*3);
        fields = {9, 10, 8};
      }
    }

    /* Read data for a single dataset */
    H5S memspace(size);
    multi_array<double, 2> data(boost::extents[N][num_fields]);
    H5Dread(dataset.getId(), H5T_NATIVE_DOUBLE, memspace.getId(), dataspace.getId(), H5P_DEFAULT, data.data());

    /* Transfer the data to the appropriate locations */
    /* The field and its first derivative */
    for(vector<int>::size_type it=0; it!=fields.size(); ++it) {
      int i = fields[it];
      double a = a_il(i, l);
      for(int j=0; j<N; ++j) {
        const complex<double> hbar(data[j][2*it], data[j][2*it+1]);
        const complex<double> dhbar(data[j][2*it+2], data[j][2*it+3]);
        const double rj = r[j];
        h[i][l][m][j] = a*hbar/rj;
        dh[i][l][m][j] = a*(dhbar - hbar/rj)/rj;
      }
    }

    /* Second derivative of the field */
    for(vector<int>::size_type it=0; it!=fields.size(); ++it) {
      int i = fields[it];
      double a = a_il(i, l);
      for(int j=0; j<N; ++j) {
        field_type::array_view<1>::type hbarj = h[boost::indices[irange()][l][m][j]];
        field_type::array_view<1>::type dhbarj = dh[boost::indices[irange()][l][m][j]];
        hbarj.reindex(1);
        dhbarj.reindex(1);
        const double rj = r[j], fj = f[j], fpj = fp[j];
        const complex<double> ddhbar = d2hdr2(i, l, m, rj, fj, fpj, r0, hbarj, dhbarj);
        const complex<double> hbar  = hbarj[i];
        const complex<double> dhbar = dhbarj[i];
        ddh[i][l][m][j] = a*(ddhbar - 2.0*(dhbar - hbar/rj)/rj)/rj;
      }
    }

    /* Negative m modes */
    double sign = isOdd(m) ? -1.0 : 1.0;
    for(vector<int>::size_type it=0; it!=fields.size(); ++it) {
      int i = fields[it];
      for(int j=0; j<N; ++j) {
        h[i][l][-m][j] = sign*conj(h[i][l][m][j]);
        dh[i][l][-m][j] = sign*conj(dh[i][l][m][j]);
        ddh[i][l][-m][j] = sign*conj(ddh[i][l][m][j]);
      }
    }
  }
  cout << endl;
}