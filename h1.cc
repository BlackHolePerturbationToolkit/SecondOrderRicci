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
#include "WignerDMatrix.h"

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
             field_type &h, field_type &dh, field_type &ddh)
{
  cout << "Reading first order fields from directory: " << dir << endl;

  /* Determine all available modes */
  vector<string> files = list_files(dir);
  if(files.size() == 0)
    exit(1);

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
  size_t N;  /* Number of grid points */
  int r0i;   /* Index of worldline point */
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
    for(size_t i=0; i<N; ++i) {
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

    /* The trace-reversed field and its first derivative.
     * We don't include the a or 1/r factors here.
     */
    for(vector<int>::size_type it=0; it!=fields.size(); ++it) {
      int i = fields[it];
      for(size_t j=0; j<N; ++j) {
        hbar[i][l][m][j]  = complex<double>(data[j][4*it], data[j][4*it+1]);
        dhbar[i][l][m][j] = complex<double>(data[j][4*it+2], data[j][4*it+3]);
      }
    }

    /* Second derivative of the trace-reversed field */
    for(vector<int>::size_type it=0; it!=fields.size(); ++it) {
      int i = fields[it];
      for(size_t j=0; j<N; ++j) {
        const double rj = r[j], fj = f[j], fpj = fp[j];
        field_type::array_view<1>::type hbarj = hbar[boost::indices[irange()][l][m][j]];
        field_type::array_view<1>::type dhbarj = dhbar[boost::indices[irange()][l][m][j]];
        hbarj.reindex(1);
        dhbarj.reindex(1);
        ddhbar[i][l][m][j] = d2hdr2(i, l, m, rj, fj, fpj, r0, hbarj, dhbarj);
      }
    }

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

  /* Include analytic expressions for the l=0 mode */
  for(size_t j=0; j<N; ++j) {
    if(r[j]<=r0) {
      h[1][0][0][j] = (-1.6710855164206668*(4. + (3. - 1.*r0)*log(1. - 2./r0))*pow(2. - 1.*r[j],3)*(2. + r[j]))/(sqrt((-3. + r0)*r0)*pow(r[j],4));
      h[2][0][0][j] = 0.0;
      h[3][0][0][j] = (-1.6710855164206668*(64. - 2.*pow(r[j],3) + (3. - 1.*r0)*log(1. - 2./r0)*(16. + pow(r[j],3))))/(sqrt((-3. + r0)*r0)*pow(r[j],3));
      h[6][0][0][j] = (-3.3421710328413337*(32. - 1.*pow(r[j],3) + (3. - 1.*r0)*log(1. - 2./r0)*(8. - 1.*pow(r[j],3))))/(sqrt((-3. + r0)*r0)*pow(r[j],3));
      dh[1][0][0][j] = (-6.684342065682667*(-4. + (-3. + r0)*log((-2. + r0)/r0))*pow(-2. + r[j],2)*(4. + r[j]))/(sqrt((-3. + r0)*r0)*pow(r[j],5));
      dh[2][0][0][j] = 0.0;
      dh[3][0][0][j] = (-80.21210478819201*(-4. + (-3. + r0)*log((-2. + r0)/r0)))/(sqrt((-3. + r0)*r0)*pow(r[j],4));
      dh[6][0][0][j] = (-80.21210478819201*(-4. + (-3. + r0)*log((-2. + r0)/r0)))/(sqrt((-3. + r0)*r0)*pow(r[j],4));
      ddh[1][0][0][j] = (13.368684131365335*(-2. + r[j])*(-20. + 2.*r[j] + pow(r[j],2))*(-4. + (-3. + r0)*log((-2. + r0)/r0)))/(pow(r[j],6)*sqrt((-3. + r0)*r0));
      ddh[2][0][0][j] = 0.0;
      ddh[3][0][0][j] = (320.84841915276803*(-4. + (-3. + r0)*log((-2. + r0)/r0)))/(pow(r[j],5)*sqrt((-3. + r0)*r0));
      ddh[6][0][0][j] = (320.84841915276803*(-4. + (-3. + r0)*log((-2. + r0)/r0)))/(pow(r[j],5)*sqrt((-3. + r0)*r0));
    } else {
      h[1][0][0][j] = (-1.6710855164206668*((3. - 1.*r0)*(16.*log(r[j]/r0)*(1. - 1.*r[j]) + log(1. - 2./r[j])*pow(2. - 1.*r[j],3)*(2. + r[j])) + 2.*(32. - 2.*r0*pow(r[j],3) - 4.*(3.*r0 + 5.*r[j]) - 1.*(-1.*r0 + r[j])*(r0 + 9.*r[j]) + r[j]*(-1.*pow(r0,2) + 3.*r0*r[j] + 6.*pow(r[j],2)))))/(sqrt((-3. + r0)*r0)*pow(r[j],4));
      h[2][0][0][j] = 0.0;
      h[3][0][0][j] = (-1.6710855164206668*((3. - 1.*r0)*(16.*log((-2. + r[j])/r0) + log(1. - 2./r[j])*pow(r[j],3)) + 2.*(32. + pow(r0,2) - 4.*r0*r[j] + 3.*pow(r[j],2) - 1.*r0*pow(r[j],2) + 12.*(-1.*r0 + r[j]))))/(sqrt((-3. + r0)*r0)*pow(r[j],3));
      h[6][0][0][j] = (-3.3421710328413337*(32. + pow(r0,2) - 4.*r0*r[j] + 3.*pow(r[j],2) - 1.*r0*pow(r[j],2) + 12.*(-1.*r0 + r[j]) + (3. - 1.*r0)*(8.*log((-2. + r[j])/r0) - 1.*log(1. - 2./r[j])*pow(r[j],3))))/(sqrt((-3. + r0)*r0)*pow(r[j],3));
      dh[1][0][0][j] = (-3.3421710328413337*(-4.*(-8. + r0)*(-4. + r0) + 32.*(-3. + r0)*log(1. - 2./r[j]) + 32.*(-3. + r0)*log(r[j]/r0) + r[j]*(48. + r0*(-20. + 3.*r0) - 24.*(-3. + r0)*log(1. - 2./r[j]) - 24.*(-3. + r0)*log(r[j]/r0) + (-3. + r0)*r[j]*(-8. + (3. + 2.*log(1. - 2./r[j]))*r[j]))))/(sqrt((-3. + r0)*r0)*pow(r[j],5));
      dh[2][0][0][j] = 0.0;
      dh[3][0][0][j] = (10.026513098524001*(-2.*(-8. + r0)*(-4. + r0) + 16.*(-3. + r0)*log((-2. + r[j])/r0) + r[j]*(8. + (-4. + r0)*r0 - 8.*(-3. + r0)*log((-2. + r[j])/r0) - 2.*(-3. + r0)*r[j])))/(sqrt((-3. + r0)*r0)*(-2. + r[j])*pow(r[j],4));
      dh[6][0][0][j] = (10.026513098524001*((-8. + r0)*(-4. + r0) - 8.*(-3. + r0)*log((-2. + r[j])/r0) - 1.*(-3. + r0)*r[j]*(4. + r[j])))/(sqrt((-3. + r0)*r0)*pow(r[j],4));
      ddh[1][0][0][j] = (6.684342065682667*(-320. + 72.*r[j] + 42.*pow(r[j],2) - 9.*pow(r[j],3) + 120.*r0 - 32.*r[j]*r0 - 14.*pow(r[j],2)*r0 + 3.*pow(r[j],3)*r0 - 10.*pow(r0,2) + 6.*r[j]*pow(r0,2) + 2.*(40. - 24.*r[j] + pow(r[j],3))*(-3. + r0)*log((-2. + r[j])/r[j]) - 16.*(-5. + 3.*r[j])*(-3. + r0)*log(r[j]/r0)))/(pow(r[j],6)*sqrt((-3. + r0)*r0));
      ddh[2][0][0][j] = 0.0;
      ddh[3][0][0][j] = (20.053026197048002*(3.*pow(r[j],3)*(-3. + r0) - 2.*pow(r[j],2)*(-4. + pow(r0,2)) - 8.*(32. - 12.*r0 + pow(r0,2)) + 8.*r[j]*(20. - 8.*r0 + pow(r0,2)) + 16.*pow(-2. + r[j],2)*(-3. + r0)*log((-2. + r[j])/r0)))/(pow(-2. + r[j],2)*pow(r[j],5)*sqrt((-3. + r0)*r0));
      ddh[6][0][0][j] = (20.053026197048002*(4.*pow(r[j],2)*(-3. + r0) + pow(r[j],3)*(-3. + r0) + 4.*(32. - 12.*r0 + pow(r0,2)) - 2.*r[j]*(8. - 4.*r0 + pow(r0,2)) + 16.*(-2. + r[j])*(-3. + r0)*log((-2. + r[j])/r0)))/((-2. + r[j])*pow(r[j],5)*sqrt((-3. + r0)*r0));
    }
  }
}

/* signum function - returns -1, 0, 1 depending on the sign of val*/
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

/* Compute the first-order puncture fields */
void compute_h1P(const double r0, const vector<double> &r, const int l_max, field_type &hP, field_type &dhP, field_type &ddhP)
{
  cout << "Computing first order punctures " << endl;

  if (l_max > 100) {
    cout << "Punctures can only be computed up to l_max=100, but l_max is " << l_max << endl;
  }

  const size_t N = r.size();  /* Number of grid points */
  
  hP.resize(boost::extents[range(1,11)][l_max+1][range(-l_max,l_max+1)][N]);
  dhP.resize(boost::extents[range(1,11)][l_max+1][range(-l_max,l_max+1)][N]);
  ddhP.resize(boost::extents[range(1,11)][l_max+1][range(-l_max,l_max+1)][N]);

  fill(hP.data(), hP.data() + hP.num_elements(), 0.0);
  fill(dhP.data(), dhP.data() + dhP.num_elements(), 0.0);
  fill(ddhP.data(), ddhP.data() + ddhP.num_elements(), 0.0);

  /* Complete elliptic integrals */
  const double ellE = gsl_sf_ellint_Ecomp(sqrt(M/(r0-2.0*M)), GSL_PREC_DOUBLE);
  const double ellK = gsl_sf_ellint_Kcomp(sqrt(M/(r0-2.0*M)), GSL_PREC_DOUBLE);

  /* Wigner-D matrix */
  WignerDMatrix WignerD;

  /* Window function size */
  const double sigma = 1.0*M;

  const complex<double> I(0.0, 1.0);

  /* Loop over all l, m modes */
#pragma omp parallel for
  for(int l=0; l<=l_max; ++l) {
    const double hP1bar000 = (2.5464790894703255*ellK*pow(-2.*M + r0,1.5))/(pow(r0,2)*sqrt(-3.*M + r0)) + (1.2732395447351628*pow(-2.*M + r0,1.5)*(ellE*(17.*M - 7.*r0)*(-2.*M + r0) + ellK*(-3.*M + r0)*(-32.*M + 7.*r0)))/((-1. + 2.*l)*(3. + 2.*l)*pow(r0,3)*pow(-3.*M + r0,1.5));
    const double hP1bar001 = (1.2732395447351628*sqrt(-2.*M + r0)*(-2.*ellK*(-4.*M + r0) + ellE*(-2.*M + r0)))/(pow(r0,3)*sqrt(-3.*M + r0)) + (0.15915494309189535*sqrt(-2.*M + r0)*(ellK*(-9246.*pow(M,4) + 9973.*pow(M,3)*r0 - 3767.*pow(M,2)*pow(r0,2) + 559.*M*pow(r0,3) - 23.*pow(r0,4)) + ellE*(3952.*pow(M,4) - 5402.*pow(M,3)*r0 + 2689.*pow(M,2)*pow(r0,2) - 510.*M*pow(r0,3) + 23.*pow(r0,4))))/((-1. + 2.*l)*(3. + 2.*l)*pow(r0,4)*pow(-3.*M + r0,2.5));
    const double hP1bar002 = (0.3183098861837907*ellE*pow(1. + 2.*l,2)*sqrt(-2.*M + r0))/(pow(r0,3)*sqrt(-3.*M + r0)) + (0.15915494309189535*sqrt(-2.*M + r0)*(5.*ellK*(-7.*M + r0)*(-3.*M + r0) + ellE*(14.*pow(M,2) - 5.*M*r0 + pow(r0,2))))/(pow(r0,4)*pow(-3.*M + r0,1.5));
    const double hP1bar003 = 0.;
    const double hP1bar004 = (0.0022104853207207684*(-1. + 2.*l)*pow(1. + 2.*l,2)*(3. + 2.*l)*(-1.*ellK*(-3.*M + r0) + 2.*ellE*(-5.*M + 2.*r0)))/(pow(r0,4)*sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5));
    const double hP1bar010 = (-2.*(1. + 2.*l)*(-2.*M + r0))/(pow(r0,2.5)*sqrt(-3.*M + r0));
    const double hP1bar011 = (2.*(1. + 2.*l)*sqrt(-3.*M + r0))/pow(r0,3.5);
    const double hP1bar012 = (0.16666666666666666*l*(1. + l)*(1. + 2.*l)*(5.*M - 2.*r0))/(pow(r0,3.5)*sqrt(-3.*M + r0)*(-2.*M + r0));
    const double hP1bar013 = 0.;

    const double hP1bar200 = 0.15915494309189535*((-64.*pow(-2.*M + r0,1.5)*(-2.*ellE*(-2.*M + r0) + ellK*(-5.*M + 2.*r0)))/((-1. + 2.*l)*(3. + 2.*l)*M*pow(r0,2)*sqrt(-3.*M + r0)) + (32.*pow((-2.*M + r0)/(-3.*M + r0),1.5)*(ellE*(1022.*pow(M,3) - 913.*pow(M,2)*r0 + 217.*M*pow(r0,2) - 8.*pow(r0,3)) + ellK*(-1284.*pow(M,3) + 1019.*pow(M,2)*r0 - 221.*M*pow(r0,2) + 8.*pow(r0,3))))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*M*pow(r0,3)));
    const double hP1bar201 = 0.15915494309189535*((-32.*sqrt(-2.*M + r0)*(-2.*ellK*(26.*pow(M,2) - 18.*M*r0 + 3.*pow(r0,2)) + ellE*(42.*pow(M,2) - 33.*M*r0 + 6.*pow(r0,2))))/((-1. + 2.*l)*(3. + 2.*l)*M*pow(r0,3)*sqrt(-3.*M + r0)) - (4.*sqrt(-2.*M + r0)*(ellK*(-3.*M + r0)*(43678.*pow(M,4) - 45091.*pow(M,3)*r0 + 14858.*pow(M,2)*pow(r0,2) - 1501.*M*pow(r0,3) - 20.*pow(r0,4)) + ellE*(104536.*pow(M,5) - 153254.*pow(M,4)*r0 + 81693.*pow(M,3)*pow(r0,2) - 18606.*pow(M,2)*pow(r0,3) + 1451.*M*pow(r0,4) + 20.*pow(r0,5))))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*M*pow(r0,4)*pow(-3.*M + r0,2.5)));
    const double hP1bar202 = 0.15915494309189535*((2.6666666666666665*sqrt(-2.*M + r0)*(-2.*ellK*(-3.*M + r0) + ellE*(-5.*M + 2.*r0)))/(M*pow(r0,3)*sqrt(-3.*M + r0)) - (4.*sqrt(-2.*M + r0)*(ellE*(590.*pow(M,3) - 731.*pow(M,2)*r0 + 287.*M*pow(r0,2) - 36.*pow(r0,3)) + ellK*(-741.*pow(M,3) + 838.*pow(M,2)*r0 - 305.*M*pow(r0,2) + 36.*pow(r0,3))))/((-1. + 2.*l)*(3. + 2.*l)*M*pow(r0,4)*pow(-3.*M + r0,1.5)));
    const double hP1bar203 = 0.;
    const double hP1bar204 = (0.02122065907891938*(-1. + l + pow(l,2))*(ellK*(-15.*pow(M,2) + 11.*M*r0 - 2.*pow(r0,2)) + 2.*ellE*(7.*pow(M,2) - 5.*M*r0 + pow(r0,2))))/(M*pow(r0,4)*sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5));
    const double hP1bar212 = (-0.08333333333333333*(1. + 2.*l)*M)/(pow(r0,3.5)*sqrt(-3.*M + r0)*(-2.*M + r0));
    const double hP1bar213 = 0.;

    const double hP2bar100 = (-2.5464790894703255*(2.*M - 1.*r0)*(-2.*ellK*(3.*M - 1.*r0)*(71.*pow(M,3) - 48.*(-11. + 3.*l + 3.*pow(l,2))*pow(M,2)*r0 + 2.*(-187. + 48.*l + 48.*pow(l,2))*M*pow(r0,2) + (63. - 16.*l - 16.*pow(l,2))*pow(r0,3)) + ellE*(2.*M - 1.*r0)*(143.*pow(M,3) - 2.*(-529. + 144.*l + 144.*pow(l,2))*pow(M,2)*r0 + (-751. + 192.*l + 192.*pow(l,2))*M*pow(r0,2) - 2.*(-63. + 16.*l + 16.*pow(l,2))*pow(r0,3))))/((45. - 72.*l - 56.*pow(l,2) + 32.*pow(l,3) + 16.*pow(l,4))*pow(r0*(-3.*M + r0),2.5)*sqrt(M*(-2.*M + r0)));
    const double hP2bar101 = (10.185916357881302*(-3.*ellK*(12.*pow(M,2) - 7.*M*r0 + pow(r0,2)) + ellE*(27.*pow(M,2) - 19.*M*r0 + 3.*pow(r0,2))))/((-1. + 2.*l)*(3. + 2.*l)*(3.*M - 1.*r0)*pow(r0,2.5)*sqrt((M*(-3.*M + r0))/(-2.*M + r0)));
    const double hP2bar102 = (0.3183098861837907*(2.*ellE*(222.*pow(M,3) - 3.*(67. + 48.*l + 48.*pow(l,2))*pow(M,2)*r0 + 3.*(19. + 32.*l + 32.*pow(l,2))*M*pow(r0,2) - 4.*pow(1. + 2.*l,2)*pow(r0,3)) - 1.*ellK*(3.*M - 1.*r0)*(189.*pow(M,2) - 1.*(109. + 96.*l + 96.*pow(l,2))*M*r0 + 8.*pow(r0 + 2.*l*r0,2))))/((-1. + 2.*l)*(3. + 2.*l)*(3.*M - 1.*r0)*pow(r0,3.5)*sqrt(M*(-3.*M + r0)*(-2.*M + r0)));

    const double hP3bar000 = (2.5464790894703255*ellK*pow(-2.*M + r0,1.5))/(pow(r0,2)*sqrt(-3.*M + r0)) + (1.2732395447351628*pow(-2.*M + r0,1.5)*(ellK*(16.*M - 9.*r0)*(-3.*M + r0) + ellE*(14.*pow(M,2) - 17.*M*r0 + 9.*pow(r0,2))))/((-1. + 2.*l)*(3. + 2.*l)*pow(r0,3)*pow(-3.*M + r0,1.5));
    const double hP3bar001 = (1.2732395447351628*sqrt(-2.*M + r0)*(-2.*ellK*(-4.*M + r0) + ellE*(-2.*M + r0)))/(pow(r0,3)*sqrt(-3.*M + r0)) - (0.15915494309189535*sqrt(-2.*M + r0)*(ellK*(-7602.*pow(M,4) + 8603.*pow(M,3)*r0 - 3433.*pow(M,2)*pow(r0,2) + 545.*M*pow(r0,3) - 25.*pow(r0,4)) + ellE*(2960.*pow(M,4) - 4294.*pow(M,3)*r0 + 2367.*pow(M,2)*pow(r0,2) - 498.*M*pow(r0,3) + 25.*pow(r0,4))))/((-1. + 2.*l)*(3. + 2.*l)*pow(r0,4)*pow(-3.*M + r0,2.5));
    const double hP3bar002 = (0.3183098861837907*ellE*pow(1. + 2.*l,2)*sqrt(-2.*M + r0))/(pow(r0,3)*sqrt(-3.*M + r0)) + (0.15915494309189535*(ellE*(2.*M - 1.*r0)*(130.*pow(M,2) - 91.*M*r0 + 15.*pow(r0,2)) + ellK*(-3.*M + r0)*(214.*pow(M,2) - 157.*M*r0 + 21.*pow(r0,2))))/(pow(r0,4)*pow(-3.*M + r0,1.5)*sqrt(-2.*M + r0));
    const double hP3bar003 = 0.;
    const double hP3bar004 = (0.0022104853207207684*(-1. + 2.*l)*pow(1. + 2.*l,2)*(3. + 2.*l)*(-1.*ellK*(-3.*M + r0) + 2.*ellE*(-5.*M + 2.*r0)))/(pow(r0,4)*sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5));
    const double hP3bar010 = (-2.*(1. + 2.*l)*(-2.*M + r0))/(pow(r0,2.5)*sqrt(-3.*M + r0));
    const double hP3bar011 = (2.*(1. + 2.*l)*sqrt(-3.*M + r0))/pow(r0,3.5);
    const double hP3bar012 = (0.16666666666666666*l*(1. + l)*(1. + 2.*l)*(5.*M - 2.*r0))/(pow(r0,3.5)*sqrt(-3.*M + r0)*(-2.*M + r0));
    const double hP3bar013 = 0.;

    const double hP3bar200 = (-10.185916357881302*pow(-2.*M + r0,1.5)*(-2.*ellE*(-2.*M + r0) + ellK*(-5.*M + 2.*r0)))/((-1. + 2.*l)*(3. + 2.*l)*M*pow(r0,2)*sqrt(-3.*M + r0)) - (5.092958178940651*pow((-2.*M + r0)/(-3.*M + r0),1.5)*(ellE*(466.*pow(M,3) - 767.*pow(M,2)*r0 + 375.*M*pow(r0,2) - 56.*pow(r0,3)) + ellK*(-588.*pow(M,3) + 901.*pow(M,2)*r0 - 403.*M*pow(r0,2) + 56.*pow(r0,3))))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*M*pow(r0,3));
    const double hP3bar201 = (-5.092958178940651*sqrt(-2.*M + r0)*(3.*ellE*(-2.*M + r0)*(-7.*M + 2.*r0) - 2.*ellK*(26.*pow(M,2) - 18.*M*r0 + 3.*pow(r0,2))))/((-1. + 2.*l)*(3. + 2.*l)*M*pow(r0,3)*sqrt(-3.*M + r0)) + (0.6366197723675814*sqrt(-2.*M + r0)*(ellK*(-3.*M + r0)*(33650.*pow(M,4) - 44381.*pow(M,3)*r0 + 21110.*pow(M,2)*pow(r0,2) - 4371.*M*pow(r0,3) + 340.*pow(r0,4)) + ellE*(80360.*pow(M,5) - 141274.*pow(M,4)*r0 + 96259.*pow(M,3)*pow(r0,2) - 31970.*pow(M,2)*pow(r0,3) + 5221.*M*pow(r0,4) - 340.*pow(r0,5))))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*M*pow(r0,4)*pow(-3.*M + r0,2.5));
    const double hP3bar202 = (0.4244131815783876*sqrt(-2.*M + r0)*(-2.*ellK*(-3.*M + r0) + ellE*(-5.*M + 2.*r0)))/(M*pow(r0,3)*sqrt(-3.*M + r0)) - (0.6366197723675814*(-1.*ellK*(-3.*M + r0)*(1790.*pow(M,3) - 2161.*pow(M,2)*r0 + 829.*M*pow(r0,2) - 100.*pow(r0,3)) + ellE*(-2.*M + r0)*(2174.*pow(M,3) - 2459.*pow(M,2)*r0 + 879.*M*pow(r0,2) - 100.*pow(r0,3))))/((-1. + 2.*l)*(3. + 2.*l)*M*pow(r0,4)*pow(-3.*M + r0,1.5)*sqrt(-2.*M + r0));
    const double hP3bar203 = 0.;
    const double hP3bar204 = (0.02122065907891938*(-1. + l + pow(l,2))*(ellK*(-15.*pow(M,2) + 11.*M*r0 - 2.*pow(r0,2)) + 2.*ellE*(7.*pow(M,2) - 5.*M*r0 + pow(r0,2))))/(M*pow(r0,4)*sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5));
    const double hP3bar212 = (-0.08333333333333333*(1. + 2.*l)*M)/(pow(r0,3.5)*sqrt(-3.*M + r0)*(-2.*M + r0));
    const double hP3bar213 = 0.;

    const double hP4bar100 = (-5.092958178940651*(ellE - 1.*ellK)*pow(-2.*M + r0,1.5))/(sqrt(M)*sqrt(r0)*sqrt(-3.*M + r0)) - (12.732395447351628*sqrt(M)*sqrt(-2.*M + r0)*(ellE*(82.*pow(M,2) - 39.*M*r0 - 1.*pow(r0,2)) - 2.*ellK*(51.*pow(M,2) - 20.*M*r0 + pow(r0,2))))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*pow(r0,1.5)*pow(-3.*M + r0,1.5)) - (2.5464790894703255*sqrt(-2.*M + r0)*(2.*ellE*(12.*pow(M,3) - 8.*pow(M,2)*r0 + 3.*M*pow(r0,2) - 1.*pow(r0,3)) + ellK*(-42.*pow(M,3) + 35.*pow(M,2)*r0 - 13.*M*pow(r0,2) + 2.*pow(r0,3))))/((-1. + 2.*l)*(3. + 2.*l)*sqrt(M)*pow(r0,1.5)*pow(-3.*M + r0,1.5));
    const double hP4bar101 = (2.5464790894703255*sqrt(-2.*M + r0)*(-1.*ellK*(-5.*M + r0) + ellE*(-4.*M + r0)))/(sqrt(M)*pow(r0,1.5)*sqrt(-3.*M + r0)) + (2.5464790894703255*sqrt(M)*(2.*ellK*M + ellE*(-2.*M + r0)))/((-1. + 2.*l)*(3. + 2.*l)*pow(r0,1.5)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) + (38.197186342054884*(ellE*(10.*pow(M,2) - 7.*M*r0 + pow(r0,2)) - 1.*ellK*(12.*pow(M,2) - 7.*M*r0 + pow(r0,2))))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*sqrt(M)*pow(r0,1.5)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0));
    const double hP4bar102 = (-0.2122065907891938*pow(1. + 2.*l,2)*(ellE*(-4.*M + r0) - 1.*ellK*(-3.*M + r0)))/(sqrt(M)*pow(r0,1.5)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) + (0.954929658551372*(2.*ellE*(60.*pow(M,3) - 80.*pow(M,2)*r0 + 31.*M*pow(r0,2) - 3.*pow(r0,3)) + ellK*(-162.*pow(M,3) + 191.*pow(M,2)*r0 - 65.*M*pow(r0,2) + 6.*pow(r0,3))))/((-1. + 2.*l)*(3. + 2.*l)*(2.*M - 1.*r0)*pow(r0,2.5)*sqrt(M*(-3.*M + r0)*(-2.*M + r0))) + (0.3183098861837907*(2.*ellK*(141.*pow(M,3) - 146.*pow(M,2)*r0 + 45.*M*pow(r0,2) - 4.*pow(r0,3)) + ellE*(-230.*pow(M,3) + 251.*pow(M,2)*r0 - 83.*M*pow(r0,2) + 8.*pow(r0,3))))/((3.*M - 1.*r0)*pow(r0,2.5)*sqrt(M*(-3.*M + r0)*(-2.*M + r0)));
    const double hP4bar103 = 0.;
    const double hP4bar104 = (-0.001768388256576615*pow(1. + 2.*l,2)*(-9. + 2.*l + 2.*pow(l,2))*(ellE*(M - 1.*r0) + ellK*(-3.*M + r0)))/(pow(r0,2.5)*sqrt((M*(-2.*M + r0))/(-3.*M + r0))*(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) + (0.0008841941282883075*(45. - 15.*pow(1. + 2.*l,2) + 4.*(-1. + 2.*l)*pow(1. + 2.*l,2)*(3. + 2.*l))*sqrt(M*r0*(-3.*M + r0))*(ellE*(4.*M - 1.*r0) + ellK*(-3.*M + r0)))/(M*pow(r0,3)*pow(-2.*M + r0,2.5));
    const double hP4bar110 = (-2.*(1. + 2.*l)*M)/(r0*sqrt(M*(-3.*M + r0)));
    const double hP4bar111 = 0.;
    const double hP4bar112 = (0.08333333333333333*(1. + 2.*l)*(2. - 3.*l - 3.*pow(l,2))*sqrt(M)*sqrt(-3.*M + r0))/(pow(r0,2)*pow(-2.*M + r0,2)) + (0.08333333333333333*(-1. - 2.*l)*sqrt(M))/(pow(r0,2)*sqrt(-3.*M + r0)*(-2.*M + r0));
    const double hP4bar113 = 0.;

    const double hP4bar300 = (-6.790610905254201*pow(-2.*M + r0,1.5)*(ellE*(17.*M - 8.*r0) + ellK*(-21.*M + 8.*r0)))/((-1. + 2.*l)*(3. + 2.*l)*pow(M,1.5)*sqrt(r0)*sqrt(-3.*M + r0)) - (1283.425461093044*sqrt(-2.*M + r0)*(ellE*(3802.*pow(M,3) - 4271.*pow(M,2)*r0 + 1541.*M*pow(r0,2) - 178.*pow(r0,3)) + 2.*ellK*(-2355.*pow(M,3) + 2429.*pow(M,2)*r0 - 815.*M*pow(r0,2) + 89.*pow(r0,3))))/((-7. + 2.*l)*(-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(9. + 2.*l)*sqrt(M)*pow(r0,1.5)*pow(-3.*M + r0,1.5)) + (10.185916357881302*sqrt(-2.*M + r0)*(ellE*(1886.*pow(M,4) + 4907.*pow(M,3)*r0 - 6097.*pow(M,2)*pow(r0,2) + 1906.*M*pow(r0,3) - 160.*pow(r0,4)) - 2.*ellK*(1155.*pow(M,4) + 3158.*pow(M,3)*r0 - 3440.*pow(M,2)*pow(r0,2) + 993.*M*pow(r0,3) - 80.*pow(r0,4))))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*pow(M,1.5)*pow(r0,1.5)*pow(-3.*M + r0,1.5)) + (10.185916357881302*sqrt(-2.*M + r0)*(ellK*(1182.*pow(M,4) + 1475.*pow(M,3)*r0 - 2891.*pow(M,2)*pow(r0,2) + 1284.*M*pow(r0,3) - 176.*pow(r0,4)) - 4.*ellE*(238.*pow(M,4) + 277.*pow(M,3)*r0 - 620.*pow(M,2)*pow(r0,2) + 299.*M*pow(r0,3) - 44.*pow(r0,4))))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*pow(M,1.5)*pow(r0,1.5)*pow(-3.*M + r0,1.5));
    const double hP4bar301 = (3.3953054526271007*((-210.*(ellE*(398.*pow(M,3) - 477.*pow(M,2)*r0 + 187.*M*pow(r0,2) - 24.*pow(r0,3)) + ellK*(-492.*pow(M,3) + 545.*pow(M,2)*r0 - 199.*M*pow(r0,2) + 24.*pow(r0,3))))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)) + (ellE*(344.*pow(M,3) - 438.*pow(M,2)*r0 + 181.*M*pow(r0,2) - 24.*pow(r0,3)) + ellK*(-426.*pow(M,3) + 503.*pow(M,2)*r0 - 193.*M*pow(r0,2) + 24.*pow(r0,3)))/((-1. + 2.*l)*(3. + 2.*l)) + (3.*(ellE*(1666.*pow(M,3) - 2151.*pow(M,2)*r0 + 899.*M*pow(r0,2) - 120.*pow(r0,3)) + ellK*(-2064.*pow(M,3) + 2473.*pow(M,2)*r0 - 959.*M*pow(r0,2) + 120.*pow(r0,3))))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l))))/(pow(M,1.5)*pow(r0,1.5)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0));
    const double hP4bar302 = (0.16976527263135505*(ellK*(-45.*pow(M,2) + 39.*M*r0 - 8.*pow(r0,2)) + ellE*(36.*pow(M,2) - 35.*M*r0 + 8.*pow(r0,2))))/(pow(M,1.5)*pow(r0,1.5)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) - (89.12676813146139*(ellK*(-9222.*pow(M,4) + 15665.*pow(M,3)*r0 - 9633.*pow(M,2)*pow(r0,2) + 2536.*M*pow(r0,3) - 240.*pow(r0,4)) + 4.*ellE*(1862.*pow(M,4) - 3335.*pow(M,3)*r0 + 2170.*pow(M,2)*pow(r0,2) - 604.*M*pow(r0,3) + 60.*pow(r0,4))))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*pow(M,1.5)*pow(r0,2.5)*sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) + (0.4244131815783876*(-2.*ellK*(8289.*pow(M,4) - 13299.*pow(M,3)*r0 + 7847.*pow(M,2)*pow(r0,2) - 2021.*M*pow(r0,3) + 192.*pow(r0,4)) + ellE*(13398.*pow(M,4) - 22723.*pow(M,3)*r0 + 14177.*pow(M,2)*pow(r0,2) - 3850.*M*pow(r0,3) + 384.*pow(r0,4))))/((-1. + 2.*l)*(3. + 2.*l)*pow(M,1.5)*pow(r0,2.5)*pow(-3.*M + r0,1.5)*sqrt(-2.*M + r0)) - (1.2732395447351628*(2.*ellE*(69216.*pow(M,5) - 139190.*pow(M,4)*r0 + 110825.*pow(M,3)*pow(r0,2) - 43615.*pow(M,2)*pow(r0,3) + 8470.*M*pow(r0,4) - 648.*pow(r0,5)) + ellK*(-171270.*pow(M,5) + 328731.*pow(M,4)*r0 - 250216.*pow(M,3)*pow(r0,2) + 94323.*pow(M,2)*pow(r0,3) - 17588.*M*pow(r0,4) + 1296.*pow(r0,5))))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*pow(M,1.5)*pow(r0,2.5)*pow(-3.*M + r0,1.5)*pow(-2.*M + r0,1.5));
    const double hP4bar303 = 0.;
    const double hP4bar304 = (0.0015157613627799558*(-13. + 4.*l + 4.*pow(l,2))*sqrt(-3.*M + r0)*(-1.*ellK*(3.*M - 1.*r0)*(15.*(57. + 8.*l + 8.*pow(l,2))*pow(M,2) - 12.*(60. + 11.*l + 11.*pow(l,2))*M*r0 + 16.*(9. + 2.*l + 2.*pow(l,2))*pow(r0,2)) + 2.*ellE*(3.*(345. + 44.*l + 44.*pow(l,2))*pow(M,3) - 6.*(219. + 37.*l + 37.*pow(l,2))*pow(M,2)*r0 + 2.*(270. + 53.*l + 53.*pow(l,2))*M*pow(r0,2) - 8.*(9. + 2.*l + 2.*pow(l,2))*pow(r0,3))))/((-1. + 2.*l)*(3. + 2.*l)*pow(M,1.5)*(3.*M - 1.*r0)*pow(r0*(-2.*M + r0),2.5));
    const double hP4bar312 = (-0.08333333333333333*(1. + 2.*l)*pow(M,1.5))/(pow(r0,2)*sqrt(-3.*M + r0)*pow(-2.*M + r0,2));
    const double hP4bar313 = 0.;

    const double hP5bar000 = (20.371832715762604*(ellE*(2.*M - 1.*r0) + ellK*(-3.*M + r0)))/((-1. + 2.*l)*(3. + 2.*l)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0));
    const double hP5bar001 = (10.185916357881302*(ellK*(18.*pow(M,3) - 24.*pow(M,2)*r0 + 9.*M*pow(r0,2) - 1.*pow(r0,3)) + ellE*(-18.*pow(M,3) + 23.*pow(M,2)*r0 - 9.*M*pow(r0,2) + pow(r0,3))))/((-1. + 2.*l)*(3. + 2.*l)*r0*pow(-3.*M + r0,1.5)*pow(-2.*M + r0,1.5));
    const double hP5bar002 = (-2.5464790894703255*(ellE - 1.*ellK)*sqrt(-3.*M + r0))/(r0*pow(-2.*M + r0,1.5));
    
    const double hP5bar200 = (-122.23099629457562*(ellE*(2.*M - 1.*r0) + ellK*(-3.*M + r0)))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) + (713.0141450516911*(ellK*(-33.*pow(M,2) + 23.*M*r0 - 4.*pow(r0,2)) + ellE*(26.*pow(M,2) - 21.*M*r0 + 4.*pow(r0,2))))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*M*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) + (20.371832715762604*(ellK*(-27.*pow(M,2) + 21.*M*r0 - 4.*pow(r0,2)) + ellE*(22.*pow(M,2) - 19.*M*r0 + 4.*pow(r0,2))))/((-1. + 2.*l)*(3. + 2.*l)*M*sqrt(-3.*M + r0)*sqrt(-2.*M + r0));
    const double hP5bar201 = (50.92958178940651*(ellK*(-594.*pow(M,4) + 942.*pow(M,3)*r0 - 539.*pow(M,2)*pow(r0,2) + 133.*M*pow(r0,3) - 12.*pow(r0,4)) + ellE*(486.*pow(M,4) - 809.*pow(M,3)*r0 + 489.*pow(M,2)*pow(r0,2) - 127.*M*pow(r0,3) + 12.*pow(r0,4))))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*M*r0*pow(-3.*M + r0,1.5)*pow(-2.*M + r0,1.5)) - (10.185916357881302*(-1.*ellK*(558.*pow(M,4) - 894.*pow(M,3)*r0 + 521.*pow(M,2)*pow(r0,2) - 131.*M*pow(r0,3) + 12.*pow(r0,4)) + ellE*(450.*pow(M,4) - 763.*pow(M,3)*r0 + 471.*pow(M,2)*pow(r0,2) - 125.*M*pow(r0,3) + 12.*pow(r0,4))))/((-1. + 2.*l)*(3. + 2.*l)*M*r0*pow(-3.*M + r0,1.5)*pow(-2.*M + r0,1.5));
    const double hP5bar202 = (5.092958178940651*(ellE - 1.*ellK)*sqrt(-3.*M + r0))/((-1. + 2.*l)*(3. + 2.*l)*r0*pow(-2.*M + r0,1.5)) + (38.197186342054884*sqrt(-3.*M + r0)*(ellE*(9.*M - 4.*r0) + ellK*(-11.*M + 4.*r0)))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*M*r0*pow(-2.*M + r0,1.5)) - (0.8488263631567752*sqrt(-3.*M + r0)*(ellE*(7.*M - 4.*r0) + ellK*(-9.*M + 4.*r0)))/(M*r0*pow(-2.*M + r0,1.5));

    const double hP6bar000 = (2.5464790894703255*ellK*M*r0)/(sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) - (1.2732395447351628*(ellK*(-3.*M + r0)*(32.*pow(M,2) - 31.*M*r0 + 8.*pow(r0,2)) + ellE*(2.*M - 1.*r0)*(17.*pow(M,2) - 15.*M*r0 + 8.*pow(r0,2))))/((-1. + 2.*l)*(3. + 2.*l)*pow(-3.*M + r0,1.5)*sqrt(-2.*M + r0));
    const double hP6bar001 = (1.2732395447351628*(ellE + 2.*ellK)*M)/(sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) + (0.15915494309189535*(-1.*ellK*(-3.*M + r0)*(197.*pow(M,3) - 164.*pow(M,2)*r0 + 11.*M*pow(r0,2) + 8.*pow(r0,3)) + ellE*(-596.*pow(M,4) + 917.*pow(M,3)*r0 - 412.*pow(M,2)*pow(r0,2) + 35.*M*pow(r0,3) + 8.*pow(r0,4))))/((-1. + 2.*l)*(3. + 2.*l)*r0*pow(-3.*M + r0,2.5)*sqrt(-2.*M + r0));
    const double hP6bar002 = (0.3183098861837907*ellE*pow(1. + 2.*l,2)*M)/(sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) + (0.15915494309189535*(ellK*(-3.*M + r0)*(61.*pow(M,2) - 51.*M*r0 + 8.*pow(r0,2)) + ellE*(158.*pow(M,3) - 173.*pow(M,2)*r0 + 65.*M*pow(r0,2) - 8.*pow(r0,3))))/(r0*pow(-3.*M + r0,1.5)*pow(-2.*M + r0,1.5));
    const double hP6bar003 = 0.;
    const double hP6bar004 = (0.08333333333333333*(1. + 2.*l)*(2. - 3.*l - 3.*pow(l,2))*sqrt(M)*sqrt(-3.*M + r0))/(pow(r0,2)*pow(-2.*M + r0,2)) + (0.08333333333333333*(-1. - 2.*l)*sqrt(M))/(pow(r0,2)*sqrt(-3.*M + r0)*(-2.*M + r0));
    const double hP6bar010 = (-2.*(1. + 2.*l)*M*sqrt(r0))/(sqrt(-3.*M + r0)*(-2.*M + r0));
    const double hP6bar011 = (-2.*(1. + 2.*l)*M*sqrt(-3.*M + r0))/(sqrt(r0)*pow(-2.*M + r0,2));
    const double hP6bar012 = (-0.16666666666666666*l*(1. + l)*(1. + 2.*l)*M*(-5.*M + 2.*r0))/(sqrt(r0)*sqrt(-3.*M + r0)*pow(-2.*M + r0,3));
    const double hP6bar013 = 0.;

    const double hP6bar200 = (10.185916357881302*r0*(2.*ellE*(-2.*M + r0) - 1.*ellK*(-5.*M + 2.*r0)))/((-1. + 2.*l)*(3. + 2.*l)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) + (5.092958178940651*(ellK*(-1284.*pow(M,4) + 2003.*pow(M,3)*r0 - 1197.*pow(M,2)*pow(r0,2) + 320.*M*pow(r0,3) - 32.*pow(r0,4)) + ellE*(1022.*pow(M,4) - 1697.*pow(M,3)*r0 + 1073.*pow(M,2)*pow(r0,2) - 304.*M*pow(r0,3) + 32.*pow(r0,4))))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*M*pow(-3.*M + r0,1.5)*sqrt(-2.*M + r0));
    const double hP6bar201 = (5.092958178940651*(1/((-1. + 2.*l)*(3. + 2.*l)) + 3./((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)))*(-2.*ellK*(-2.*M + r0) + ellE*(-3.*M + 2.*r0)))/(sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) + (0.6366197723675814*(ellE*(10064.*pow(M,5) - 18017.*pow(M,4)*r0 + 12332.*pow(M,3)*pow(r0,2) - 3983.*pow(M,2)*pow(r0,3) + 596.*M*pow(r0,4) - 32.*pow(r0,5)) + ellK*(-12453.*pow(M,5) + 21071.*pow(M,4)*r0 - 13683.*pow(M,3)*pow(r0,2) + 4229.*pow(M,2)*pow(r0,3) - 612.*M*pow(r0,4) + 32.*pow(r0,5))))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*M*r0*pow(-3.*M + r0,2.5)*sqrt(-2.*M + r0));
    const double hP6bar202 = (-0.4244131815783876*(ellE*(5.*M - 2.*r0) + 2.*ellK*(-3.*M + r0)))/(sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) + (0.6366197723675814*(ellK*(-1563.*pow(M,4) + 2354.*pow(M,3)*r0 - 1343.*pow(M,2)*pow(r0,2) + 340.*M*pow(r0,3) - 32.*pow(r0,4)) + ellE*(1282.*pow(M,4) - 2029.*pow(M,3)*r0 + 1217.*pow(M,2)*pow(r0,2) - 324.*M*pow(r0,3) + 32.*pow(r0,4))))/((-1. + 2.*l)*(3. + 2.*l)*M*r0*pow(-3.*M + r0,1.5)*pow(-2.*M + r0,1.5));
    const double hP6bar203 = 0.;
    const double hP6bar204 = (0.02122065907891938*(-1. + l + pow(l,2))*(ellK*(-15.*pow(M,2) + 11.*M*r0 - 2.*pow(r0,2)) + 2.*ellE*(7.*pow(M,2) - 5.*M*r0 + pow(r0,2))))/(r0*sqrt(-3.*M + r0)*pow(-2.*M + r0,3.5));
    const double hP6bar212 = (-0.08333333333333333*(1. + 2.*l)*pow(M,2))/(sqrt(r0)*sqrt(-3.*M + r0)*pow(-2.*M + r0,3));
    const double hP6bar213 = 0.;

    const double hP7bar000 = (10.185916357881302*(1/((-1. + 2.*l)*(3. + 2.*l)) + 3./((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)) + 30./((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)))*r0*(2.*ellE*(-2.*M + r0) - 1.*ellK*(-5.*M + 2.*r0)))/(sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) + (5.092958178940651*(ellK*(3.*M - 1.*r0)*(44.*pow(M,2) - 43.*M*r0 + 12.*pow(r0,2)) - 1.*ellE*(2.*M - 1.*r0)*(31.*pow(M,2) - 37.*M*r0 + 12.*pow(r0,2))))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (25.464790894703256*(ellK*(3.*M - 1.*r0)*(44.*pow(M,2) - 43.*M*r0 + 12.*pow(r0,2)) - 1.*ellE*(2.*M - 1.*r0)*(31.*pow(M,2) - 37.*M*r0 + 12.*pow(r0,2))))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (1069.5212175775366*(ellK*(3.*M - 1.*r0)*(44.*pow(M,2) - 43.*M*r0 + 12.*pow(r0,2)) - 1.*ellE*(2.*M - 1.*r0)*(31.*pow(M,2) - 37.*M*r0 + 12.*pow(r0,2))))/((-7. + 2.*l)*(-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(9. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2)));
    const double hP7bar001 = (5.092958178940651*(1/((-1. + 2.*l)*(3. + 2.*l)) + 3./((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)) + 30./((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)))*(-2.*ellK*(-2.*M + r0) + ellE*(-3.*M + 2.*r0)))/(sqrt(-3.*M + r0)*sqrt(-2.*M + r0));
    const double hP7bar002 = (0.4244131815783876*(1. + 3./((-1. + 2.*l)*(3. + 2.*l)) - 18./((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)))*(-2.*ellK*(-3.*M + r0) + ellE*(-5.*M + 2.*r0)))/(sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) + (0.6366197723675814*(ellE*(130.*pow(M,3) - 139.*pow(M,2)*r0 + 55.*M*pow(r0,2) - 8.*pow(r0,3)) + ellK*(-123.*pow(M,3) + 134.*pow(M,2)*r0 - 55.*M*pow(r0,2) + 8.*pow(r0,3))))/((-1. + 2.*l)*(3. + 2.*l)*r0*pow(6.*pow(M,2) - 5.*M*r0 + pow(r0,2),1.5)) + (1.909859317102744*(ellE*(130.*pow(M,3) - 139.*pow(M,2)*r0 + 55.*M*pow(r0,2) - 8.*pow(r0,3)) + ellK*(-123.*pow(M,3) + 134.*pow(M,2)*r0 - 55.*M*pow(r0,2) + 8.*pow(r0,3))))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*r0*pow(6.*pow(M,2) - 5.*M*r0 + pow(r0,2),1.5)) + (19.098593171027442*(ellE*(130.*pow(M,3) - 139.*pow(M,2)*r0 + 55.*M*pow(r0,2) - 8.*pow(r0,3)) + ellK*(-123.*pow(M,3) + 134.*pow(M,2)*r0 - 55.*M*pow(r0,2) + 8.*pow(r0,3))))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*r0*pow(6.*pow(M,2) - 5.*M*r0 + pow(r0,2),1.5));
    const double hP7bar003 = 0.;
    const double hP7bar004 = (0.005305164769729845*(-5. + pow(1. + 2.*l,2) + 30./((-1. + 2.*l)*(3. + 2.*l)))*(-1.*ellK*(-3.*M + r0)*(-5.*M + 2.*r0) + 2.*ellE*(7.*pow(M,2) - 5.*M*r0 + pow(r0,2))))/(r0*sqrt(-3.*M + r0)*pow(-2.*M + r0,3.5));
    const double hP7bar012 = (-0.08333333333333333*(1. + 2.*l)*pow(M,2))/(sqrt(r0)*sqrt(-3.*M + r0)*pow(-2.*M + r0,3));
    const double hP7bar013 = 0.;

    const double hP7bar200 = (10.185916357881302*ellK*M*r0)/((-1. + 2.*l)*(3. + 2.*l)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) + (53.47606087887684*r0*(ellK*(-247.*pow(M,2) + 200.*M*r0 - 40.*pow(r0,2)) + 20.*ellE*(10.*pow(M,2) - 9.*M*r0 + 2.*pow(r0,2))))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*M*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) - (0.8488263631567752*r0*(4.*ellE*(-2.*M + r0)*(-5.*M + 2.*r0) - 1.*ellK*(51.*pow(M,2) - 40.*M*r0 + 8.*pow(r0,2))))/(M*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) + (1604.281826366305*r0*(4.*ellE*(-2.*M + r0)*(-5.*M + 2.*r0) - 1.*ellK*(51.*pow(M,2) - 40.*M*r0 + 8.*pow(r0,2))))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*M*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) + (15.278874536821952*(ellE*(34.*pow(M,3) - 79.*pow(M,2)*r0 + 47.*M*pow(r0,2) - 8.*pow(r0,3)) + ellK*(-96.*pow(M,3) + 125.*pow(M,2)*r0 - 55.*M*pow(r0,2) + 8.*pow(r0,3))))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (1.2732395447351628*(3.*ellE*(154.*pow(M,3) - 155.*pow(M,2)*r0 + 47.*M*pow(r0,2) - 4.*pow(r0,3)) + ellK*(-552.*pow(M,3) + 511.*pow(M,2)*r0 - 145.*M*pow(r0,2) + 12.*pow(r0,3))))/((-1. + 2.*l)*(3. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (5294.130027008807*(ellE*(734.*pow(M,3) - 1097.*pow(M,2)*r0 + 517.*M*pow(r0,2) - 76.*pow(r0,3)) + ellK*(-1320.*pow(M,3) + 1511.*pow(M,2)*r0 - 585.*M*pow(r0,2) + 76.*pow(r0,3))))/((-7. + 2.*l)*(-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(9. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (34.3774677078494*(ellE*(3778.*pow(M,3) - 4519.*pow(M,2)*r0 + 1739.*M*pow(r0,2) - 212.*pow(r0,3)) + ellK*(-5400.*pow(M,3) + 5577.*pow(M,2)*r0 - 1895.*M*pow(r0,2) + 212.*pow(r0,3))))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2)));
    const double hP7bar201 = (5.092958178940651*(ellE + 2.*ellK)*M)/((-1. + 2.*l)*(3. + 2.*l)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) + (0.4244131815783876*(-2.*ellK*(39.*pow(M,2) - 26.*M*r0 + 4.*pow(r0,2)) + ellE*(67.*pow(M,2) - 48.*M*r0 + 8.*pow(r0,2))))/(M*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) - (802.1409131831525*(-2.*ellK*(39.*pow(M,2) - 26.*M*r0 + 4.*pow(r0,2)) + ellE*(67.*pow(M,2) - 48.*M*r0 + 8.*pow(r0,2))))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*M*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) - (26.73803043943842*(-2.*ellK*(203.*pow(M,2) - 130.*M*r0 + 20.*pow(r0,2)) + ellE*(327.*pow(M,2) - 240.*M*r0 + 40.*pow(r0,2))))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*M*sqrt(-3.*M + r0)*sqrt(-2.*M + r0));
    const double hP7bar202 = (-1.2732395447351628*ellE*M)/(sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) - (0.02122065907891938*pow(1. + 2.*l,2)*(-4.*ellK*(15.*pow(M,2) - 11.*M*r0 + 2.*pow(r0,2)) + ellE*(41.*pow(M,2) - 40.*M*r0 + 8.*pow(r0,2))))/(M*sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) - (66.84507609859604*(-4.*ellK*(-3.*M + r0)*(-5.*M + 2.*r0) + ellE*(49.*pow(M,2) - 40.*M*r0 + 8.*pow(r0,2))))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*M*sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) + (1.5915494309189535*(-36.*ellK*(15.*pow(M,2) - 11.*M*r0 + 2.*pow(r0,2)) + ellE*(433.*pow(M,2) - 360.*M*r0 + 72.*pow(r0,2))))/((-1. + 2.*l)*(3. + 2.*l)*M*sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) + (0.6366197723675814*(ellE*(158.*pow(M,3) - 173.*pow(M,2)*r0 + 65.*M*pow(r0,2) - 8.*pow(r0,3)) + ellK*(-183.*pow(M,3) + 166.*pow(M,2)*r0 - 59.*M*pow(r0,2) + 8.*pow(r0,3))))/((-1. + 2.*l)*(3. + 2.*l)*r0*pow(6.*pow(M,2) - 5.*M*r0 + pow(r0,2),1.5)) + (0.05305164769729845*(ellE*(-2094.*pow(M,4) + 3757.*pow(M,3)*r0 - 2301.*pow(M,2)*pow(r0,2) + 572.*M*pow(r0,3) - 48.*pow(r0,4)) + ellK*(2619.*pow(M,4) - 4458.*pow(M,3)*r0 + 2551.*pow(M,2)*pow(r0,2) - 596.*M*pow(r0,3) + 48.*pow(r0,4))))/(M*(2.*M - 1.*r0)*(3.*M - 1.*r0)*r0*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) + (100.26761414789407*(ellK*(-2619.*pow(M,4) + 4458.*pow(M,3)*r0 - 2551.*pow(M,2)*pow(r0,2) + 596.*M*pow(r0,3) - 48.*pow(r0,4)) + ellE*(2094.*pow(M,4) - 3757.*pow(M,3)*r0 + 2301.*pow(M,2)*pow(r0,2) - 572.*M*pow(r0,3) + 48.*pow(r0,4))))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*M*r0*pow(6.*pow(M,2) - 5.*M*r0 + pow(r0,2),1.5)) + (3.3422538049298023*(ellK*(-14559.*pow(M,4) + 23618.*pow(M,3)*r0 - 13227.*pow(M,2)*pow(r0,2) + 3044.*M*pow(r0,3) - 240.*pow(r0,4)) + ellE*(11734.*pow(M,4) - 20169.*pow(M,3)*r0 + 12025.*pow(M,2)*pow(r0,2) - 2924.*M*pow(r0,3) + 240.*pow(r0,4))))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*M*(2.*M - 1.*r0)*(3.*M - 1.*r0)*r0*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2)));
    const double hP7bar203 = 0.;
    const double hP7bar204 = (-0.0010105075751866372*(513. - 452.*l - 388.*pow(l,2) + 128.*pow(l,3) + 64.*pow(l,4))*sqrt(-3.*M + r0)*(ellE*(5.*M - 2.*r0) + 2.*ellK*(-3.*M + r0)))/(M*r0*pow(-2.*M + r0,2.5)) - (0.00006315672344916482*(1575. - 496.*pow(1. + 2.*l,2) + 19.*pow(1. + 2.*l,4) + 9450./((-1. + 2.*l)*(3. + 2.*l)))*(ellK*(3.*M - 1.*r0)*(47.*pow(M,2) - 40.*M*r0 + 8.*pow(r0,2)) + 2.*ellE*(-55.*pow(M,3) + 72.*pow(M,2)*r0 - 30.*M*pow(r0,2) + 4.*pow(r0,3))))/(M*r0*sqrt(-3.*M + r0)*pow(-2.*M + r0,3.5));
    const double hP7bar210 = (-1.*(1. + 2.*l)*M*sqrt(r0))/(sqrt(-3.*M + r0)*(-2.*M + r0));
    const double hP7bar211 = (-1.*(1. + 2.*l)*M*sqrt(-3.*M + r0))/(sqrt(r0)*pow(-2.*M + r0,2));
    const double hP7bar212 = (-0.08333333333333333*(1. + 2.*l)*(-4. + l + pow(l,2))*M*(-5.*M + 2.*r0))/(sqrt(r0)*sqrt(-3.*M + r0)*pow(-2.*M + r0,3));
    const double hP7bar213 = 0.;

    const double hP7bar400 = (0.6790610905254202*r0*(ellK*(5.*M - 2.*r0)*(-21.*M + 8.*r0)*(-19.*M + 8.*r0) + 2.*ellE*(-2.*M + r0)*(403.*pow(M,2) - 320.*M*r0 + 64.*pow(r0,2))))/((-1. + 2.*l)*(3. + 2.*l)*pow(M,2)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) - (550589.5228089159*r0*(ellK*(5.*M - 2.*r0)*(-21.*M + 8.*r0)*(-19.*M + 8.*r0) + 2.*ellE*(-2.*M + r0)*(403.*pow(M,2) - 320.*M*r0 + 64.*pow(r0,2))))/((-9. + 2.*l)*(-7. + 2.*l)*(-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(9. + 2.*l)*(11. + 2.*l)*pow(M,2)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) + (91.67324722093171*r0*(ellK*(14165.*pow(M,3) - 16866.*pow(M,2)*r0 + 6720.*M*pow(r0,2) - 896.*pow(r0,3)) - 2.*ellE*(5722.*pow(M,3) - 7341.*pow(M,2)*r0 + 3136.*M*pow(r0,2) - 448.*pow(r0,3))))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*pow(M,2)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) + (3.3953054526271007*r0*(ellK*(14085.*pow(M,3) - 16834.*pow(M,2)*r0 + 6720.*M*pow(r0,2) - 896.*pow(r0,3)) - 2.*ellE*(5690.*pow(M,3) - 7325.*pow(M,2)*r0 + 3136.*M*pow(r0,2) - 448.*pow(r0,3))))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*pow(M,2)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) + (7058.840036011742*r0*(2.*ellE*(7094.*pow(M,3) - 9307.*pow(M,2)*r0 + 4032.*M*pow(r0,2) - 576.*pow(r0,3)) + ellK*(-17555.*pow(M,3) + 21422.*pow(M,2)*r0 - 8640.*M*pow(r0,2) + 1152.*pow(r0,3))))/((-7. + 2.*l)*(-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(9. + 2.*l)*pow(M,2)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) - (1.6976527263135504*(-1.*ellK*(35964.*pow(M,5) + 18033.*pow(M,4)*r0 - 85379.*pow(M,3)*pow(r0,2) + 61604.*pow(M,2)*pow(r0,3) - 17536.*M*pow(r0,4) + 1792.*pow(r0,5)) + ellE*(29082.*pow(M,5) + 11909.*pow(M,4)*r0 - 71377.*pow(M,3)*pow(r0,2) + 55188.*pow(M,2)*pow(r0,3) - 16640.*M*pow(r0,4) + 1792.*pow(r0,5))))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*pow(M,2)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (1.7697520376000868e6*(ellK*(-74940.*pow(M,5) + 40031.*pow(M,4)*r0 + 51347.*pow(M,3)*pow(r0,2) - 52196.*pow(M,2)*pow(r0,3) + 16512.*M*pow(r0,4) - 1792.*pow(r0,5)) + ellE*(60794.*pow(M,5) - 37899.*pow(M,4)*r0 - 40705.*pow(M,3)*pow(r0,2) + 46292.*pow(M,2)*pow(r0,3) - 15616.*M*pow(r0,4) + 1792.*pow(r0,5))))/((-11. + 2.*l)*(-9. + 2.*l)*(-7. + 2.*l)*(-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(9. + 2.*l)*(11. + 2.*l)*(13. + 2.*l)*pow(M,2)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (15.278874536821952*(-1.*ellK*(127380.*pow(M,5) + 25067.*pow(M,4)*r0 - 239121.*pow(M,3)*pow(r0,2) + 180108.*pow(M,2)*pow(r0,3) - 52096.*M*pow(r0,4) + 5376.*pow(r0,5)) + ellE*(103102.*pow(M,5) + 10823.*pow(M,4)*r0 - 198795.*pow(M,3)*pow(r0,2) + 161116.*pow(M,2)*pow(r0,3) - 49408.*M*pow(r0,4) + 5376.*pow(r0,5))))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*pow(M,2)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (1787.6283208081686*(-1.*ellK*(1.10886e6*pow(M,5) + 711217.*pow(M,4)*r0 - 2.885571e6*pow(M,3)*pow(r0,2) + 2.051748e6*pow(M,2)*pow(r0,3) - 580736.*M*pow(r0,4) + 59136.*pow(r0,5)) + ellE*(896282.*pow(M,5) + 492613.*pow(M,4)*r0 - 2.416785e6*pow(M,3)*pow(r0,2) + 1.838996e6*pow(M,2)*pow(r0,3) - 551168.*M*pow(r0,4) + 59136.*pow(r0,5))))/((-9. + 2.*l)*(-7. + 2.*l)*(-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(9. + 2.*l)*(11. + 2.*l)*pow(M,2)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (15.278874536821952*(ellK*(-1.79094e6*pow(M,5) + 304903.*pow(M,4)*r0 + 2.290011e6*pow(M,3)*pow(r0,2) - 1.887108e6*pow(M,2)*pow(r0,3) + 562816.*M*pow(r0,4) - 59136.*pow(r0,5)) + ellE*(1.451242e6*pow(M,5) - 379027.*pow(M,4)*r0 - 1.880025e6*pow(M,3)*pow(r0,2) + 1.683316e6*pow(M,2)*pow(r0,3) - 533248.*M*pow(r0,4) + 59136.*pow(r0,5))))/((-7. + 2.*l)*(-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(9. + 2.*l)*pow(M,2)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2)));
    const double hP7bar401 = (0.3395305452627101*(-2.*ellK*(19.*M - 8.*r0)*(21.*M - 8.*r0)*(10.*M - 3.*r0) + ellE*(6451.*pow(M,3) - 7698.*pow(M,2)*r0 + 3008.*M*pow(r0,2) - 384.*pow(r0,3))))/((-1. + 2.*l)*(3. + 2.*l)*pow(M,2)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) - (275294.76140445797*(-2.*ellK*(19.*M - 8.*r0)*(21.*M - 8.*r0)*(10.*M - 3.*r0) + ellE*(6451.*pow(M,3) - 7698.*pow(M,2)*r0 + 3008.*M*pow(r0,2) - 384.*pow(r0,3))))/((-9. + 2.*l)*(-7. + 2.*l)*(-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(9. + 2.*l)*(11. + 2.*l)*pow(M,2)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) + (45.83662361046586*(ellE*(45037.*pow(M,3) - 53806.*pow(M,2)*r0 + 21056.*M*pow(r0,2) - 2688.*pow(r0,3)) + 2.*ellK*(-27850.*pow(M,3) + 30739.*pow(M,2)*r0 - 11200.*M*pow(r0,2) + 1344.*pow(r0,3))))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*pow(M,2)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) + (1.6976527263135504*(ellE*(45085.*pow(M,3) - 53838.*pow(M,2)*r0 + 21056.*M*pow(r0,2) - 2688.*pow(r0,3)) + 2.*ellK*(-27882.*pow(M,3) + 30755.*pow(M,2)*r0 - 11200.*M*pow(r0,2) + 1344.*pow(r0,3))))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*pow(M,2)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) - (3529.420018005871*(ellE*(58299.*pow(M,3) - 69442.*pow(M,2)*r0 + 27072.*M*pow(r0,2) - 3456.*pow(r0,3)) + 2.*ellK*(-36070.*pow(M,3) + 39653.*pow(M,2)*r0 - 14400.*M*pow(r0,2) + 1728.*pow(r0,3))))/((-7. + 2.*l)*(-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(9. + 2.*l)*pow(M,2)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0));
    const double hP7bar402 = (-26470.650135044034*(2.*ellK*(-3.*M + r0)*(-21.*M + 8.*r0)*(-19.*M + 8.*r0) + ellE*(5.*M - 2.*r0)*(387.*pow(M,2) - 320.*M*r0 + 64.*pow(r0,2))))/((-7. + 2.*l)*(-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(9. + 2.*l)*pow(M,2)*sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) - (0.4244131815783876*(ellE*(1975.*pow(M,3) - 2390.*pow(M,2)*r0 + 960.*M*pow(r0,2) - 128.*pow(r0,3)) + 2.*ellK*(-1221.*pow(M,3) + 1367.*pow(M,2)*r0 - 512.*M*pow(r0,2) + 64.*pow(r0,3))))/((-1. + 2.*l)*(3. + 2.*l)*pow(M,2)*sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) - (0.012126090902239647*(ellE*(2015.*pow(M,3) - 2406.*pow(M,2)*r0 + 960.*M*pow(r0,2) - 128.*pow(r0,3)) + 2.*ellK*(-1245.*pow(M,3) + 1375.*pow(M,2)*r0 - 512.*M*pow(r0,2) + 64.*pow(r0,3))))/(pow(M,2)*sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) + (515.662015617741*(ellE*(13465.*pow(M,3) - 16586.*pow(M,2)*r0 + 6720.*M*pow(r0,2) - 896.*pow(r0,3)) + 2.*ellK*(-8331.*pow(M,3) + 9497.*pow(M,2)*r0 - 3584.*M*pow(r0,2) + 448.*pow(r0,3))))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*pow(M,2)*sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) - (0.4244131815783876*(ellE*(68805.*pow(M,3) - 83522.*pow(M,2)*r0 + 33600.*M*pow(r0,2) - 4480.*pow(r0,3)) + 2.*ellK*(-42543.*pow(M,3) + 47781.*pow(M,2)*r0 - 17920.*M*pow(r0,2) + 2240.*pow(r0,3))))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*pow(M,2)*sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) + (34411.845175557246*(ellE*(390642.*pow(M,5) - 862915.*pow(M,4)*r0 + 740903.*pow(M,3)*pow(r0,2) - 310320.*pow(M,2)*pow(r0,3) + 63616.*M*pow(r0,4) - 5120.*pow(r0,5)) + ellK*(-483435.*pow(M,5) + 1.023438e6*pow(M,4)*r0 - 840151.*pow(M,3)*pow(r0,2) + 336688.*pow(M,2)*pow(r0,3) - 66176.*M*pow(r0,4) + 5120.*pow(r0,5))))/((-9. + 2.*l)*(-7. + 2.*l)*(-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(9. + 2.*l)*(11. + 2.*l)*pow(M,2)*r0*pow(6.*pow(M,2) - 5.*M*r0 + pow(r0,2),1.5)) - (0.04244131815783876*(ellE*(390642.*pow(M,5) - 862915.*pow(M,4)*r0 + 740903.*pow(M,3)*pow(r0,2) - 310320.*pow(M,2)*pow(r0,3) + 63616.*M*pow(r0,4) - 5120.*pow(r0,5)) + ellK*(-483435.*pow(M,5) + 1.023438e6*pow(M,4)*r0 - 840151.*pow(M,3)*pow(r0,2) + 336688.*pow(M,2)*pow(r0,3) - 66176.*M*pow(r0,4) + 5120.*pow(r0,5))))/((-1. + 2.*l)*(3. + 2.*l)*pow(M,2)*(2.*M - 1.*r0)*(3.*M - 1.*r0)*r0*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (0.2122065907891938*(ellE*(2.676078e6*pow(M,5) - 5.954125e6*pow(M,4)*r0 + 5.135849e6*pow(M,3)*pow(r0,2) - 2.158224e6*pow(M,2)*pow(r0,3) + 443776.*M*pow(r0,4) - 35840.*pow(r0,5)) + ellK*(-3.311973e6*pow(M,5) + 7.064082e6*pow(M,4)*r0 - 5.825209e6*pow(M,3)*pow(r0,2) + 2.342032e6*pow(M,2)*pow(r0,3) - 461696.*M*pow(r0,4) + 35840.*pow(r0,5))))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*pow(M,2)*r0*pow(6.*pow(M,2) - 5.*M*r0 + pow(r0,2),1.5)) - (5.729577951308232*(ellE*(2.637134e6*pow(M,5) - 5.896605e6*pow(M,4)*r0 + 5.102201e6*pow(M,3)*pow(r0,2) - 2.14888e6*pow(M,2)*pow(r0,3) + 442752.*M*pow(r0,4) - 35840.*pow(r0,5)) + ellK*(-3.263925e6*pow(M,5) + 6.997426e6*pow(M,4)*r0 - 5.787977e6*pow(M,3)*pow(r0,2) + 2.332176e6*pow(M,2)*pow(r0,3) - 460672.*M*pow(r0,4) + 35840.*pow(r0,5))))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*pow(M,2)*r0*pow(6.*pow(M,2) - 5.*M*r0 + pow(r0,2),1.5)) + (441.1775022507339*(ellE*(3.710498e6*pow(M,5) - 8.053835e6*pow(M,4)*r0 + 6.836367e6*pow(M,3)*pow(r0,2) - 2.8396e6*pow(M,2)*pow(r0,3) + 577664.*M*pow(r0,4) - 46080.*pow(r0,5)) + ellK*(-4.591155e6*pow(M,5) + 9.544222e6*pow(M,4)*r0 - 7.747519e6*pow(M,3)*pow(r0,2) + 3.079472e6*pow(M,2)*pow(r0,3) - 600704.*M*pow(r0,4) + 46080.*pow(r0,5))))/((-7. + 2.*l)*(-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(9. + 2.*l)*pow(M,2)*r0*pow(6.*pow(M,2) - 5.*M*r0 + pow(r0,2),1.5));
    const double hP7bar403 = 0.;
    const double hP7bar404 = (0.005389373734328732*(-81. + 4.*pow(1. + 2.*l,2) - 945./((-1. + 2.*l)*(3. + 2.*l)) - 39690./((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)))*sqrt(-3.*M + r0)*(-8.*ellK*(5.*M - 2.*r0)*(3.*M - 1.*r0) + ellE*(97.*pow(M,2) - 80.*M*r0 + 16.*pow(r0,2))))/(pow(M,2)*r0*pow(-2.*M + r0,2.5)) - (0.00008420896459888643*(-675. + 31.*pow(1. + 2.*l,2) - 8505./((-1. + 2.*l)*(3. + 2.*l)) - 19845./((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)) - 5.35815e6/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)))*(-1.*ellK*(5.*M - 2.*r0)*(3.*M - 1.*r0)*(383.*pow(M,2) - 320.*M*r0 + 64.*pow(r0,2)) + 2.*ellE*(2321.*pow(M,4) - 3855.*pow(M,3)*r0 + 2371.*pow(M,2)*pow(r0,2) - 640.*M*pow(r0,3) + 64.*pow(r0,4))))/(pow(M,2)*r0*sqrt(-3.*M + r0)*pow(-2.*M + r0,3.5));
    const double hP7bar412 = (-0.041666666666666664*(1. + 2.*l)*pow(M,2))/(sqrt(r0)*sqrt(-3.*M + r0)*pow(-2.*M + r0,3));
    const double hP7bar413 = 0;

    const double hP8bar100 = (-5.092958178940651*sqrt(-2.*M + r0)*(ellK*(3.*M - 1.*r0) + ellE*(-2.*M + r0)))/(sqrt(M)*sqrt(r0)*sqrt(-3.*M + r0)) - (12.732395447351628*pow(-2.*M + r0,1.5)*(ellK*(-39.*pow(M,2) + 49.*M*r0 - 12.*pow(r0,2)) + 2.*ellE*(23.*pow(M,2) - 23.*M*r0 + 6.*pow(r0,2))))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*sqrt(M)*pow(r0,1.5)*pow(-3.*M + r0,1.5)) - (2.5464790894703255*sqrt(-2.*M + r0)*(ellE*(34.*pow(M,3) - 43.*pow(M,2)*r0 + 17.*M*pow(r0,2) - 2.*pow(r0,3)) + 2.*ellK*(-9.*pow(M,3) + 21.*pow(M,2)*r0 - 9.*M*pow(r0,2) + pow(r0,3))))/((-1. + 2.*l)*(3. + 2.*l)*sqrt(M)*(3.*M - 1.*r0)*pow(r0,1.5)*sqrt(-3.*M + r0));
    const double hP8bar101 = (38.197186342054884*sqrt(-2.*M + r0)*(-1.*ellK*(-5.*M + r0) + ellE*(-4.*M + r0)))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*sqrt(M)*pow(r0,1.5)*sqrt(-3.*M + r0)) - (2.5464790894703255*sqrt(M)*(2.*ellK*M + ellE*(-2.*M + r0)))/((-1. + 2.*l)*(3. + 2.*l)*pow(r0,1.5)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) + (2.5464790894703255*(ellE*(10.*pow(M,2) - 7.*M*r0 + pow(r0,2)) - 1.*ellK*(12.*pow(M,2) - 7.*M*r0 + pow(r0,2))))/(sqrt(M)*pow(r0,1.5)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0));
    const double hP8bar102 = (-0.2122065907891938*pow(1. + 2.*l,2)*(-1.*ellK*(-3.*M + r0) + ellE*(-1.*M + r0)))/(sqrt(M)*pow(r0,1.5)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) + (0.3183098861837907*(-1.*ellK*(162.*pow(M,3) - 179.*pow(M,2)*r0 + 55.*M*pow(r0,2) - 4.*pow(r0,3)) + 2.*ellE*(60.*pow(M,3) - 76.*pow(M,2)*r0 + 27.*M*pow(r0,2) - 2.*pow(r0,3))))/(sqrt(M)*pow(r0,2.5)*sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) - (0.954929658551372*(ellE*(230.*pow(M,3) - 269.*pow(M,2)*r0 + 95.*M*pow(r0,2) - 10.*pow(r0,3)) + 2.*ellK*(-141.*pow(M,3) + 155.*pow(M,2)*r0 - 51.*M*pow(r0,2) + 5.*pow(r0,3))))/((-1. + 2.*l)*(3. + 2.*l)*sqrt(M)*pow(r0,2.5)*pow(-3.*M + r0,1.5)*sqrt(-2.*M + r0));
    const double hP8bar103 = 0.;
    const double hP8bar104 = (-0.0008841941282883075*(-4.*ellE*((-63. - 158.*l - 142.*pow(l,2) + 32.*pow(l,3) + 16.*pow(l,4))*pow(M,2) + 5.*(9. + 10.*l + 14.*pow(l,2) + 8.*pow(l,3) + 4.*pow(l,4))*M*r0 - 3.*(3. - 2.*l + 2.*pow(l,2) + 8.*pow(l,3) + 4.*pow(l,4))*pow(r0,2)) - 2.*ellK*(15.*(9. + 2.*l + 10.*pow(l,2) + 16.*pow(l,3) + 8.*pow(l,4))*pow(M,2) - 1.*(99. - 26.*l + 86.*pow(l,2) + 224.*pow(l,3) + 112.*pow(l,4))*M*r0 + 6.*(3. - 2.*l + 2.*pow(l,2) + 8.*pow(l,3) + 4.*pow(l,4))*pow(r0,2))))/(sqrt(M)*pow(r0,2.5)*sqrt(-3.*M + r0)*pow(-2.*M + r0,2.5));
    const double hP8bar110 = (2.*(1. + 2.*l)*sqrt(M))/(r0*sqrt(-3.*M + r0));
    const double hP8bar111 = 0.;
    const double hP8bar112 = (0.08333333333333333*(1. + 2.*l)*sqrt(M)*(M - 1.*r0 + 3.*l*(1. + l)*(-2.*M + r0)))/(pow(r0,2)*sqrt(-3.*M + r0)*pow(-2.*M + r0,2));
    const double hP8bar113 = 0.;

    const double hP8bar300 = (6.790610905254201*sqrt(-2.*M + r0)*(ellK*(-57.*pow(M,2) + 43.*M*r0 - 8.*pow(r0,2)) + ellE*(46.*pow(M,2) - 39.*M*r0 + 8.*pow(r0,2))))/((-1. + 2.*l)*(3. + 2.*l)*pow(M,1.5)*sqrt(r0)*sqrt(-3.*M + r0)) - (1283.425461093044*pow(-2.*M + r0,1.5)*(2.*ellE*(237.*pow(M,3) - 660.*pow(M,2)*r0 + 431.*M*pow(r0,2) - 80.*pow(r0,3)) + ellK*(-585.*pow(M,3) + 1581.*pow(M,2)*r0 - 942.*M*pow(r0,2) + 160.*pow(r0,3))))/((-7. + 2.*l)*(-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(9. + 2.*l)*pow(M*r0,1.5)*pow(-3.*M + r0,1.5)) - (30.557749073643905*sqrt(-2.*M + r0)*(ellE*(634.*pow(M,4) + 45.*pow(M,3)*r0 - 661.*pow(M,2)*pow(r0,2) + 336.*M*pow(r0,3) - 48.*pow(r0,4)) - 2.*ellK*(393.*pow(M,4) + 64.*pow(M,3)*r0 - 389.*pow(M,2)*pow(r0,2) + 180.*M*pow(r0,3) - 24.*pow(r0,4))))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*pow(M*r0*(-3.*M + r0),1.5)) - (30.557749073643905*sqrt(-2.*M + r0)*(ellK*(-9030.*pow(M,4) + 14347.*pow(M,3)*r0 - 9245.*pow(M,2)*pow(r0,2) + 2782.*M*pow(r0,3) - 320.*pow(r0,4)) + 2.*ellE*(3644.*pow(M,4) - 6132.*pow(M,3)*r0 + 4137.*pow(M,2)*pow(r0,2) - 1311.*M*pow(r0,3) + 160.*pow(r0,4))))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*pow(M*r0*(-3.*M + r0),1.5));
    const double hP8bar301 = (3.3953054526271007*((ellE*(398.*pow(M,3) - 477.*pow(M,2)*r0 + 187.*M*pow(r0,2) - 24.*pow(r0,3)) + ellK*(-492.*pow(M,3) + 545.*pow(M,2)*r0 - 199.*M*pow(r0,2) + 24.*pow(r0,3)))/((-1. + 2.*l)*(3. + 2.*l)) + (210.*(ellK*(426.*pow(M,3) - 503.*pow(M,2)*r0 + 193.*M*pow(r0,2) - 24.*pow(r0,3)) + ellE*(-344.*pow(M,3) + 438.*pow(M,2)*r0 - 181.*M*pow(r0,2) + 24.*pow(r0,3))))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)) + (3.*(ellE*(2044.*pow(M,3) - 2424.*pow(M,2)*r0 + 941.*M*pow(r0,2) - 120.*pow(r0,3)) + ellK*(-2526.*pow(M,3) + 2767.*pow(M,2)*r0 - 1001.*M*pow(r0,2) + 120.*pow(r0,3))))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l))))/(pow(M,1.5)*pow(r0,1.5)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0));
    const double hP8bar302 = (0.16976527263135505*(ellK*(-75.*pow(M,2) + 49.*M*r0 - 8.*pow(r0,2)) + ellE*(61.*pow(M,2) - 45.*M*r0 + 8.*pow(r0,2))))/(pow(M,1.5)*pow(r0,1.5)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) - (0.4244131815783876*(ellK*(-9222.*pow(M,4) + 16493.*pow(M,3)*r0 - 10611.*pow(M,2)*pow(r0,2) + 2914.*M*pow(r0,3) - 288.*pow(r0,4)) + 2.*ellE*(3724.*pow(M,4) - 7006.*pow(M,3)*r0 + 4766.*pow(M,2)*pow(r0,2) - 1385.*M*pow(r0,3) + 144.*pow(r0,4))))/((-1. + 2.*l)*(3. + 2.*l)*pow(M,1.5)*(2.*M - 1.*r0)*pow(r0,2.5)*sqrt((-3.*M + r0)*(-2.*M + r0))) - (89.12676813146139*(-2.*ellK*(8289.*pow(M,4) - 12840.*pow(M,3)*r0 + 7325.*pow(M,2)*pow(r0,2) - 1826.*M*pow(r0,3) + 168.*pow(r0,4)) + ellE*(13398.*pow(M,4) - 21985.*pow(M,3)*r0 + 13265.*pow(M,2)*pow(r0,2) - 3484.*M*pow(r0,3) + 336.*pow(r0,4))))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*pow(M,1.5)*pow(r0,2.5)*pow(-3.*M + r0,1.5)*sqrt(-2.*M + r0)) - (1.2732395447351628*(ellE*(107268.*pow(M,5) - 210880.*pow(M,4)*r0 + 159925.*pow(M,3)*pow(r0,2) - 57995.*pow(M,2)*pow(r0,3) + 9920.*M*pow(r0,4) - 624.*pow(r0,5)) - 4.*ellK*(33210.*pow(M,5) - 62211.*pow(M,4)*r0 + 44956.*pow(M,3)*pow(r0,2) - 15573.*pow(M,2)*pow(r0,3) + 2558.*M*pow(r0,4) - 156.*pow(r0,5))))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*pow(M,1.5)*(2.*M - 1.*r0)*(3.*M - 1.*r0)*pow(r0,2.5)*sqrt((-3.*M + r0)*(-2.*M + r0)));
    const double hP8bar303 = 0.;
    const double hP8bar304 = (0.0015157613627799558*(-13. + 4.*l + 4.*pow(l,2))*(9. + 16.*l + 16.*pow(l,2))*(ellK*(-75.*pow(M,2) + 49.*M*r0 - 8.*pow(r0,2)) + ellE*(61.*pow(M,2) - 45.*M*r0 + 8.*pow(r0,2))))/((-1. + 2.*l)*(3. + 2.*l)*pow(M,1.5)*pow(r0,2.5)*sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) - (0.004547284088339867*(-13. + 4.*l + 4.*pow(l,2))*M*sqrt(-3.*M + r0)*(ellK*(-45.*pow(M,2) + 39.*M*r0 - 8.*pow(r0,2)) + ellE*(36.*pow(M,2) - 35.*M*r0 + 8.*pow(r0,2))))/pow(M*r0*(-2.*M + r0),2.5);
    const double hP8bar312 = (0.08333333333333333*(1. + 2.*l)*pow(M,1.5))/(pow(r0,2)*sqrt(-3.*M + r0)*pow(-2.*M + r0,2));
    const double hP8bar313 = 0;

    const double hP9bar200 = (-122.23099629457562*(-1.*ellK*(-3.*M + r0) + ellE*(-2.*M + r0)))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) + (1426.0282901033822*sqrt(-2.*M + r0)*(-2.*ellK*(-3.*M + r0) + ellE*(-5.*M + 2.*r0)))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*M*sqrt(-3.*M + r0)) + (40.74366543152521*sqrt(-3.*M + r0)*(2.*ellE*(-2.*M + r0) - 1.*ellK*(-5.*M + 2.*r0)))/((-1. + 2.*l)*(3. + 2.*l)*M*sqrt(-2.*M + r0));
    const double hP9bar201 = (-101.85916357881302*(ellE*(108.*pow(M,3) - 131.*pow(M,2)*r0 + 50.*M*pow(r0,2) - 6.*pow(r0,3)) + ellK*(-135.*pow(M,3) + 150.*pow(M,2)*r0 - 53.*M*pow(r0,2) + 6.*pow(r0,3))))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*M*r0*pow(-3.*M + r0,1.5)*sqrt(-2.*M + r0)) + (20.371832715762604*(3.*ellE*(26.*pow(M,3) - 35.*pow(M,2)*r0 + 15.*M*pow(r0,2) - 2.*pow(r0,3)) + ellK*(-96.*pow(M,3) + 121.*pow(M,2)*r0 - 48.*M*pow(r0,2) + 6.*pow(r0,3))))/((-1. + 2.*l)*(3. + 2.*l)*M*r0*sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5));
    const double hP9bar202 = (-5.092958178940651*(ellE - 1.*ellK)*sqrt(-3.*M + r0))/((-1. + 2.*l)*(3. + 2.*l)*r0*pow(-2.*M + r0,1.5)) - (1.6976527263135504*sqrt(-3.*M + r0)*(ellE*(5.*M - 2.*r0) + 2.*ellK*(-3.*M + r0)))/(M*r0*pow(-2.*M + r0,1.5)) + (76.39437268410977*sqrt(-3.*M + r0)*(ellE*(4.*M - 2.*r0) + ellK*(-5.*M + 2.*r0)))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*M*r0*pow(-2.*M + r0,1.5));

    const double hP10bar200 = (-10.185916357881302*ellK*M*r0)/((-1. + 2.*l)*(3. + 2.*l)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) + (10695.212175775367*ellE*M*r0)/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) - (13262.063097961456*ellK*M*r0)/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) - (9625.69095819783*ellE*pow(r0,2))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) + (10695.212175775367*ellK*pow(r0,2))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) + (2139.0424351550732*ellE*pow(r0,3))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*M*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) - (2139.0424351550732*ellK*pow(r0,3))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*M*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) + (16.976527263135505*ellE*r0*sqrt(-2.*M + r0))/sqrt(-3.*M + r0) - (20.371832715762604*ellK*r0*sqrt(-2.*M + r0))/sqrt(-3.*M + r0) - (32085.6365273261*ellE*r0*sqrt(-2.*M + r0))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*sqrt(-3.*M + r0)) + (38502.76383279132*ellK*r0*sqrt(-2.*M + r0))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*sqrt(-3.*M + r0)) - (6.790610905254201*ellE*pow(r0,2)*sqrt(-2.*M + r0))/(M*sqrt(-3.*M + r0)) + (6.790610905254201*ellK*pow(r0,2)*sqrt(-2.*M + r0))/(M*sqrt(-3.*M + r0)) + (12834.25461093044*ellE*pow(r0,2)*sqrt(-2.*M + r0))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*M*sqrt(-3.*M + r0)) - (12834.25461093044*ellK*pow(r0,2)*sqrt(-2.*M + r0))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*M*sqrt(-3.*M + r0)) - (631.5268141886407*ellE*pow(M,3))/((-1. + 2.*l)*(3. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) + (825.0592249883855*ellK*pow(M,3))/((-1. + 2.*l)*(3. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (519.4817342519464*ellE*pow(M,3))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) + (1466.7719555349074*ellK*pow(M,3))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (100657.22544858303*ellE*pow(M,3))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) + (103132.40312354818*ellK*pow(M,3))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (1.1858851260499726e6*ellE*pow(M,3))/((-7. + 2.*l)*(-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(9. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (635295.6032410568*ellK*pow(M,3))/((-7. + 2.*l)*(-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(9. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) + (692.6423123359285*ellE*pow(M,2)*r0)/((-1. + 2.*l)*(3. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (809.7803504515635*ellK*pow(M,2)*r0)/((-1. + 2.*l)*(3. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) + (1207.0310884089342*ellE*pow(M,2)*r0)/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (1909.8593171027442*ellK*pow(M,2)*r0)/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) + (87456.27784876886*ellE*pow(M,2)*r0)/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (84293.55081964671*ellK*pow(M,2)*r0)/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (465883.44237677497*ellE*pow(M,2)*r0)/((-7. + 2.*l)*(-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(9. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) + (1.9270633298312058e6*ellK*pow(M,2)*r0)/((-7. + 2.*l)*(-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(9. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (239.3690344102106*ellE*M*pow(r0,2))/((-1. + 2.*l)*(3. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) + (254.64790894703256*ellK*M*pow(r0,2))/((-1. + 2.*l)*(3. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (718.1071032306318*ellE*M*pow(r0,2))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) + (840.3380995252074*ellK*M*pow(r0,2))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (19388.891787227058*ellE*M*pow(r0,2))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) + (17876.283208081684*ellK*M*pow(r0,2))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) + (995296.4450776557*ellE*M*pow(r0,2))/((-7. + 2.*l)*(-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(9. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (1.2705912064821136e6*ellK*M*pow(r0,2))/((-7. + 2.*l)*(-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(9. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) + (25.464790894703256*ellE*pow(r0,3))/((-1. + 2.*l)*(3. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (25.464790894703256*ellK*pow(r0,3))/((-1. + 2.*l)*(3. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) + (122.23099629457562*ellE*pow(r0,3))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (122.23099629457562*ellK*pow(r0,3))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) + (412.52961249419275*ellE*pow(r0,3))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (412.52961249419275*ellK*pow(r0,3))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (232941.72118838748*ellE*pow(r0,3))/((-7. + 2.*l)*(-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(9. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) + (232941.72118838748*ellK*pow(r0,3))/((-7. + 2.*l)*(-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(9. + 2.*l)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2)));
    const double hP10bar201 = (27.162443621016806*ellE*M)/(sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) - (35.65070725258456*ellK*M)/(sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) - (5.092958178940651*ellE*M)/((-1. + 2.*l)*(3. + 2.*l)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) - (10.185916357881302*ellK*M)/((-1. + 2.*l)*(3. + 2.*l)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) - (8770.0739841358*ellE*M)/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) + (10802.164297533121*ellK*M)/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) - (51337.01844372176*ellE*M)/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) + (67379.83670738482*ellK*M)/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) - (20.371832715762604*ellE*r0)/(sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) + (22.069485442076154*ellK*r0)/(sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) + (6417.12730546522*ellE*r0)/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) - (6951.887914253989*ellK*r0)/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) + (38502.76383279132*ellE*r0)/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) - (41711.32748552393*ellK*r0)/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) + (3.3953054526271007*ellE*pow(r0,2))/(M*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) - (3.3953054526271007*ellK*pow(r0,2))/(M*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) - (1069.5212175775366*ellE*pow(r0,2))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*M*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) + (1069.5212175775366*ellK*pow(r0,2))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*M*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) - (6417.12730546522*ellE*pow(r0,2))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*M*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) + (6417.12730546522*ellK*pow(r0,2))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*M*sqrt(-3.*M + r0)*sqrt(-2.*M + r0));
    const double hP10bar202 = (0.08488263631567752*ellE*M)/(sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) + (1.2732395447351628*ellK*M)/(sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) - (4.753427633677941*ellE*l*M)/(sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) + (5.092958178940651*ellK*l*M)/(sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) - (4.753427633677941*ellE*pow(l,2)*M)/(sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) + (5.092958178940651*ellK*pow(l,2)*M)/(sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) + (700.2817496043396*ellE*M)/((-1. + 2.*l)*(3. + 2.*l)*sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) - (859.4366926962349*ellK*M)/((-1. + 2.*l)*(3. + 2.*l)*sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) + (0.8488263631567752*ellE*r0)/(sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) - (0.9337089994724527*ellK*r0)/(sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) + (3.3953054526271007*ellE*l*r0)/(sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) - (3.7348359978898107*ellK*l*r0)/(sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) + (3.3953054526271007*ellE*pow(l,2)*r0)/(sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) - (3.7348359978898107*ellK*pow(l,2)*r0)/(sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) - (572.9577951308232*ellE*r0)/((-1. + 2.*l)*(3. + 2.*l)*sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) + (630.2535746439056*ellK*r0)/((-1. + 2.*l)*(3. + 2.*l)*sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) - (0.16976527263135505*ellE*pow(r0,2))/(M*sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) + (0.16976527263135505*ellK*pow(r0,2))/(M*sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) - (0.6790610905254202*ellE*l*pow(r0,2))/(M*sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) + (0.6790610905254202*ellK*l*pow(r0,2))/(M*sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) - (0.6790610905254202*ellE*pow(l,2)*pow(r0,2))/(M*sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) + (0.6790610905254202*ellK*pow(l,2)*pow(r0,2))/(M*sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) + (114.59155902616465*ellE*pow(r0,2))/((-1. + 2.*l)*(3. + 2.*l)*M*sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) - (114.59155902616465*ellK*pow(r0,2))/((-1. + 2.*l)*(3. + 2.*l)*M*sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) + (1069.5212175775366*ellE*sqrt(-3.*M + r0))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*pow(-2.*M + r0,1.5)) - (1336.9015219719208*ellK*sqrt(-3.*M + r0))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*pow(-2.*M + r0,1.5)) - (534.7606087887683*ellE*r0*sqrt(-3.*M + r0))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*M*pow(-2.*M + r0,1.5)) + (534.7606087887683*ellK*r0*sqrt(-3.*M + r0))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*M*pow(-2.*M + r0,1.5)) + (110.13522061959158*ellE*pow(M,2))/((-1. + 2.*l)*(3. + 2.*l)*pow(6.*pow(M,2) - 5.*M*r0 + pow(r0,2),1.5)) - (105.67888221301851*ellK*pow(M,2))/((-1. + 2.*l)*(3. + 2.*l)*pow(6.*pow(M,2) - 5.*M*r0 + pow(r0,2),1.5)) - (66831.70708337633*ellE*pow(M,2))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*pow(6.*pow(M,2) - 5.*M*r0 + pow(r0,2),1.5)) + (78382.53623321372*ellK*pow(M,2))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*pow(6.*pow(M,2) - 5.*M*r0 + pow(r0,2),1.5)) - (100.58592403407786*ellE*pow(M,3))/((-1. + 2.*l)*(3. + 2.*l)*r0*pow(6.*pow(M,2) - 5.*M*r0 + pow(r0,2),1.5)) + (116.50141834326739*ellK*pow(M,3))/((-1. + 2.*l)*(3. + 2.*l)*r0*pow(6.*pow(M,2) - 5.*M*r0 + pow(r0,2),1.5)) + (38689.93004586739*ellE*pow(M,3))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*r0*pow(6.*pow(M,2) - 5.*M*r0 + pow(r0,2),1.5)) - (48048.24069967084*ellK*pow(M,3))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*r0*pow(6.*pow(M,2) - 5.*M*r0 + pow(r0,2),1.5)) - (41.38028520389279*ellE*M*r0)/((-1. + 2.*l)*(3. + 2.*l)*pow(6.*pow(M,2) - 5.*M*r0 + pow(r0,2),1.5)) + (37.5605665696873*ellK*M*r0)/((-1. + 2.*l)*(3. + 2.*l)*pow(6.*pow(M,2) - 5.*M*r0 + pow(r0,2),1.5)) + (39973.35550696043*ellE*M*r0)/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*pow(6.*pow(M,2) - 5.*M*r0 + pow(r0,2),1.5)) - (44010.79810331564*ellK*M*r0)/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*pow(6.*pow(M,2) - 5.*M*r0 + pow(r0,2),1.5)) + (5.092958178940651*ellE*pow(r0,2))/((-1. + 2.*l)*(3. + 2.*l)*pow(6.*pow(M,2) - 5.*M*r0 + pow(r0,2),1.5)) - (5.092958178940651*ellK*pow(r0,2))/((-1. + 2.*l)*(3. + 2.*l)*pow(6.*pow(M,2) - 5.*M*r0 + pow(r0,2),1.5)) - (9746.012095175303*ellE*pow(r0,2))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*pow(6.*pow(M,2) - 5.*M*r0 + pow(r0,2),1.5)) + (10147.08255176688*ellK*pow(r0,2))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*pow(6.*pow(M,2) - 5.*M*r0 + pow(r0,2),1.5)) + (802.1409131831525*ellE*pow(r0,3))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*M*pow(6.*pow(M,2) - 5.*M*r0 + pow(r0,2),1.5)) - (802.1409131831525*ellK*pow(r0,3))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*M*pow(6.*pow(M,2) - 5.*M*r0 + pow(r0,2),1.5)) + (60.47887837492023*ellE*M)/((2.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (68.96714200648798*ellK*M)/((2.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (114305.08012859923*ellE*M)/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(2.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) + (130347.89839226229*ellK*M)/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(2.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (45.41221042888747*ellE*pow(M,2))/((2.*M - 1.*r0)*r0*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) + (56.02253996834716*ellK*pow(M,2))/((2.*M - 1.*r0)*r0*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) + (85829.07771059732*ellE*pow(M,2))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(2.*M - 1.*r0)*r0*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (105882.60054017614*ellK*pow(M,2))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(2.*M - 1.*r0)*r0*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (23.979344759178897*ellE*r0)/((2.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) + (25.25258430391406*ellK*r0)/((2.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) + (45320.961594848115*ellE*r0)/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(2.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (47727.384334397575*ellK*r0)/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(2.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) + (2.5464790894703255*ellE*pow(r0,2))/(M*(2.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (2.5464790894703255*ellK*pow(r0,2))/(M*(2.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (4812.845479098915*ellE*pow(r0,2))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*M*(2.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) + (4812.845479098915*ellK*pow(r0,2))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*M*(2.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2)));
    const double hP10bar203 = 0.;
    const double hP10bar204 = (0.0002526268937966593*(ellK*(15.*(29241. - 8304.*l - 9280.*pow(l,2) - 1440.*pow(l,3) + 560.*pow(l,4) + 1536.*pow(l,5) + 512.*pow(l,6))*pow(M,3) - 5.*(103113. - 30768.*l - 33184.*pow(l,2) - 3168.*pow(l,3) + 2576.*pow(l,4) + 4992.*pow(l,5) + 1664.*pow(l,6))*pow(M,2)*r0 + 384.*(513. - 156.*l - 166.*pow(l,2) - 12.*pow(l,3) + 14.*pow(l,4) + 24.*pow(l,5) + 8.*pow(l,6))*M*pow(r0,2) - 48.*(513. - 156.*l - 166.*pow(l,2) - 12.*pow(l,3) + 14.*pow(l,4) + 24.*pow(l,5) + 8.*pow(l,6))*pow(r0,3)) - 2.*ellE*(5.*(35397. - 7824.*l - 10544.*pow(l,2) - 4608.*pow(l,3) - 224.*pow(l,4) + 2496.*pow(l,5) + 832.*pow(l,6))*pow(M,3) - 2.*(112347. - 31224.*l - 35444.*pow(l,2) - 6408.*pow(l,3) + 1876.*pow(l,4) + 6096.*pow(l,5) + 2032.*pow(l,6))*pow(M,2)*r0 + 180.*(513. - 156.*l - 166.*pow(l,2) - 12.*pow(l,3) + 14.*pow(l,4) + 24.*pow(l,5) + 8.*pow(l,6))*M*pow(r0,2) - 24.*(513. - 156.*l - 166.*pow(l,2) - 12.*pow(l,3) + 14.*pow(l,4) + 24.*pow(l,5) + 8.*pow(l,6))*pow(r0,3))))/((-1. + 2.*l)*(3. + 2.*l)*M*pow(2.*M - 1.*r0,3)*r0*sqrt((-3.*M + r0)*(-2.*M + r0)));
    const double hP10bar210 = (M*sqrt(r0))/(sqrt(-3.*M + r0)*(-2.*M + r0)) + (2.*l*M*sqrt(r0))/(sqrt(-3.*M + r0)*(-2.*M + r0));
    const double hP10bar211 = (M*sqrt(-3.*M + r0))/(sqrt(r0)*pow(-2.*M + r0,2)) + (2.*l*M*sqrt(-3.*M + r0))/(sqrt(r0)*pow(-2.*M + r0,2));
    const double hP10bar212 = (0.6666666666666666*M*sqrt(r0*(-3.*M + r0)))/(pow(2.*M - 1.*r0,3)*(-3.*M + r0)) + (1.1666666666666667*l*M*sqrt(r0*(-3.*M + r0)))/(pow(2.*M - 1.*r0,3)*(-3.*M + r0)) - (0.5*pow(l,2)*M*sqrt(r0*(-3.*M + r0)))/(pow(2.*M - 1.*r0,3)*(-3.*M + r0)) - (0.3333333333333333*pow(l,3)*M*sqrt(r0*(-3.*M + r0)))/(pow(2.*M - 1.*r0,3)*(-3.*M + r0)) - (1.6666666666666667*pow(M,2)*sqrt(r0*(-3.*M + r0)))/(pow(2.*M - 1.*r0,3)*r0*(-3.*M + r0)) - (2.9166666666666665*l*pow(M,2)*sqrt(r0*(-3.*M + r0)))/(pow(2.*M - 1.*r0,3)*r0*(-3.*M + r0)) + (1.25*pow(l,2)*pow(M,2)*sqrt(r0*(-3.*M + r0)))/(pow(2.*M - 1.*r0,3)*r0*(-3.*M + r0)) + (0.8333333333333334*pow(l,3)*pow(M,2)*sqrt(r0*(-3.*M + r0)))/(pow(2.*M - 1.*r0,3)*r0*(-3.*M + r0));
    const double hP10bar213 = 0.;

    const double hP10bar400 = (5.4324887242033615*r0*sqrt(-2.*M + r0)*(-8.*ellK*(-3.*M + r0)*(-5.*M + 2.*r0) + ellE*(97.*pow(M,2) - 80.*M*r0 + 16.*pow(r0,2))))/((-1. + 2.*l)*(3. + 2.*l)*pow(M,2)*sqrt(-3.*M + r0)) - (4.4047161824713275e6*r0*sqrt(-2.*M + r0)*(-8.*ellK*(-3.*M + r0)*(-5.*M + 2.*r0) + ellE*(97.*pow(M,2) - 80.*M*r0 + 16.*pow(r0,2))))/((-9. + 2.*l)*(-7. + 2.*l)*(-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(9. + 2.*l)*(11. + 2.*l)*pow(M,2)*sqrt(-3.*M + r0)) + (56470.72028809394*r0*(ellE*(1786.*pow(M,3) - 2333.*pow(M,2)*r0 + 1008.*M*pow(r0,2) - 144.*pow(r0,3)) + 2.*ellK*(-1105.*pow(M,3) + 1342.*pow(M,2)*r0 - 540.*M*pow(r0,2) + 72.*pow(r0,3))))/((-7. + 2.*l)*(-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(9. + 2.*l)*pow(M,2)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) + (733.3859777674537*r0*(ellK*(1655.*pow(M,3) - 2062.*pow(M,2)*r0 + 840.*M*pow(r0,2) - 112.*pow(r0,3)) + ellE*(-1338.*pow(M,3) + 1789.*pow(M,2)*r0 - 784.*M*pow(r0,2) + 112.*pow(r0,3))))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*pow(M,2)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) + (27.162443621016806*r0*(ellK*(1665.*pow(M,3) - 2066.*pow(M,2)*r0 + 840.*M*pow(r0,2) - 112.*pow(r0,3)) + ellE*(-1346.*pow(M,3) + 1793.*pow(M,2)*r0 - 784.*M*pow(r0,2) + 112.*pow(r0,3))))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*pow(M,2)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) + (1.4158016300800694e7*(2.*M - 1.*r0)*(ellK*(-645.*pow(M,4) + 5114.*pow(M,3)*r0 - 5305.*pow(M,2)*pow(r0,2) + 1896.*M*pow(r0,3) - 224.*pow(r0,4)) + ellE*(536.*pow(M,4) - 4173.*pow(M,3)*r0 + 4651.*pow(M,2)*pow(r0,2) - 1784.*M*pow(r0,3) + 224.*pow(r0,4))))/((-11. + 2.*l)*(-9. + 2.*l)*(-7. + 2.*l)*(-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(9. + 2.*l)*(11. + 2.*l)*(13. + 2.*l)*pow(M,2)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (13.581221810508403*(-1.*ellK*(3582.*pow(M,5) + 3615.*pow(M,4)*r0 - 11470.*pow(M,3)*pow(r0,2) + 7921.*pow(M,2)*pow(r0,3) - 2216.*M*pow(r0,4) + 224.*pow(r0,5)) + ellE*(2892.*pow(M,5) + 2656.*pow(M,4)*r0 - 9641.*pow(M,3)*pow(r0,2) + 7107.*pow(M,2)*pow(r0,3) - 2104.*M*pow(r0,4) + 224.*pow(r0,5))))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*pow(M,2)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (122.23099629457562*(-1.*ellK*(8310.*pow(M,5) + 14474.*pow(M,4)*r0 - 36537.*pow(M,3)*pow(r0,2) + 24351.*pow(M,2)*pow(r0,3) - 6712.*M*pow(r0,4) + 672.*pow(r0,5)) + ellE*(6694.*pow(M,5) + 11081.*pow(M,4)*r0 - 30840.*pow(M,3)*pow(r0,2) + 21877.*pow(M,2)*pow(r0,3) - 6376.*M*pow(r0,4) + 672.*pow(r0,5))))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*pow(M,2)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (122.23099629457562*(-1.*ellK*(42690.*pow(M,5) + 231794.*pow(M,4)*r0 - 444447.*pow(M,3)*pow(r0,2) + 279621.*pow(M,2)*pow(r0,3) - 75112.*M*pow(r0,4) + 7392.*pow(r0,5)) + ellE*(33994.*pow(M,5) + 184151.*pow(M,4)*r0 - 377580.*pow(M,3)*pow(r0,2) + 251767.*pow(M,2)*pow(r0,3) - 71416.*M*pow(r0,4) + 7392.*pow(r0,5))))/((-7. + 2.*l)*(-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(9. + 2.*l)*pow(M,2)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (14301.026566465349*(-1.*ellK*(127950.*pow(M,5) + 104779.*pow(M,4)*r0 - 370002.*pow(M,3)*pow(r0,2) + 259041.*pow(M,2)*pow(r0,3) - 72872.*M*pow(r0,4) + 7392.*pow(r0,5)) + ellE*(103364.*pow(M,5) + 75196.*pow(M,4)*r0 - 310485.*pow(M,3)*pow(r0,2) + 232307.*pow(M,2)*pow(r0,3) - 69176.*M*pow(r0,4) + 7392.*pow(r0,5))))/((-9. + 2.*l)*(-7. + 2.*l)*(-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(9. + 2.*l)*(11. + 2.*l)*pow(M,2)*(3.*M - 1.*r0)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2)));
    const double hP10bar401 = (2.7162443621016807*((-1.*(2.*ellE*(-2.*M + r0)*(203.*pow(M,2) - 140.*M*r0 + 24.*pow(r0,2)) - 1.*ellK*(-3.*M + r0)*(335.*pow(M,2) - 256.*M*r0 + 48.*pow(r0,2))))/((-1. + 2.*l)*(3. + 2.*l)) + (810810.*(2.*ellE*(-2.*M + r0)*(203.*pow(M,2) - 140.*M*r0 + 24.*pow(r0,2)) - 1.*ellK*(-3.*M + r0)*(335.*pow(M,2) - 256.*M*r0 + 48.*pow(r0,2))))/((-9. + 2.*l)*(-7. + 2.*l)*(-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(9. + 2.*l)*(11. + 2.*l)) + (5.*(ellE*(5693.*pow(M,3) - 6768.*pow(M,2)*r0 + 2632.*M*pow(r0,2) - 336.*pow(r0,3)) + ellK*(-7047.*pow(M,3) + 7727.*pow(M,2)*r0 - 2800.*M*pow(r0,2) + 336.*pow(r0,3))))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)) + (135.*(ellE*(5699.*pow(M,3) - 6772.*pow(M,2)*r0 + 2632.*M*pow(r0,2) - 336.*pow(r0,3)) + ellK*(-7055.*pow(M,3) + 7731.*pow(M,2)*r0 - 2800.*M*pow(r0,2) + 336.*pow(r0,3))))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)) - (10395.*(ellE*(7278.*pow(M,3) - 8674.*pow(M,2)*r0 + 3384.*M*pow(r0,2) - 432.*pow(r0,3)) + ellK*(-9005.*pow(M,3) + 9907.*pow(M,2)*r0 - 3600.*M*pow(r0,2) + 432.*pow(r0,3))))/((-7. + 2.*l)*(-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(9. + 2.*l))))/(pow(M,2)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0));
    const double hP10bar402 = (211765.20108035227*sqrt(-3.*M + r0)*(8.*ellE*(-2.*M + r0)*(-5.*M + 2.*r0) - 1.*ellK*(-11.*M + 4.*r0)*(-9.*M + 4.*r0)))/((-7. + 2.*l)*(-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(9. + 2.*l)*pow(M,2)*pow(-2.*M + r0,1.5)) - (0.09700872721791717*(2.*ellE*(115.*pow(M,3) - 146.*pow(M,2)*r0 + 60.*M*pow(r0,2) - 8.*pow(r0,3)) + ellK*(-285.*pow(M,3) + 335.*pow(M,2)*r0 - 128.*M*pow(r0,2) + 16.*pow(r0,3))))/(pow(M,2)*sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) - (3.3953054526271007*(ellE*(235.*pow(M,3) - 294.*pow(M,2)*r0 + 120.*M*pow(r0,2) - 16.*pow(r0,3)) + ellK*(-291.*pow(M,3) + 337.*pow(M,2)*r0 - 128.*M*pow(r0,2) + 16.*pow(r0,3))))/((-1. + 2.*l)*(3. + 2.*l)*pow(M,2)*sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) + (4125.296124941928*(2.*ellE*(845.*pow(M,3) - 1038.*pow(M,2)*r0 + 420.*M*pow(r0,2) - 56.*pow(r0,3)) + ellK*(-2091.*pow(M,3) + 2377.*pow(M,2)*r0 - 896.*M*pow(r0,2) + 112.*pow(r0,3))))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*pow(M,2)*sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) - (3.3953054526271007*(ellE*(8265.*pow(M,3) - 10306.*pow(M,2)*r0 + 4200.*M*pow(r0,2) - 560.*pow(r0,3)) + ellK*(-10233.*pow(M,3) + 11811.*pow(M,2)*r0 - 4480.*M*pow(r0,2) + 560.*pow(r0,3))))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*pow(M,2)*sqrt(-3.*M + r0)*pow(-2.*M + r0,1.5)) - (0.3395305452627101*(ellK*(-22020.*pow(M,4) + 37907.*pow(M,3)*r0 - 23825.*pow(M,2)*pow(r0,2) + 6472.*M*pow(r0,3) - 640.*pow(r0,4)) + ellE*(17798.*pow(M,4) - 32269.*pow(M,3)*r0 + 21429.*pow(M,2)*pow(r0,2) - 6152.*M*pow(r0,3) + 640.*pow(r0,4))))/((-1. + 2.*l)*(3. + 2.*l)*pow(M,2)*(2.*M - 1.*r0)*r0*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) + (275294.76140445797*(ellK*(-22020.*pow(M,4) + 37907.*pow(M,3)*r0 - 23825.*pow(M,2)*pow(r0,2) + 6472.*M*pow(r0,3) - 640.*pow(r0,4)) + ellE*(17798.*pow(M,4) - 32269.*pow(M,3)*r0 + 21429.*pow(M,2)*pow(r0,2) - 6152.*M*pow(r0,3) + 640.*pow(r0,4))))/((-9. + 2.*l)*(-7. + 2.*l)*(-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(9. + 2.*l)*(11. + 2.*l)*pow(M,2)*(2.*M - 1.*r0)*r0*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - (45.83662361046586*(ellE*(385928.*pow(M,5) - 820210.*pow(M,4)*r0 + 686407.*pow(M,3)*pow(r0,2) - 282115.*pow(M,2)*pow(r0,3) + 56824.*M*pow(r0,4) - 4480.*pow(r0,5)) + ellK*(-477435.*pow(M,5) + 971017.*pow(M,4)*r0 - 777309.*pow(M,3)*pow(r0,2) + 305767.*pow(M,2)*pow(r0,3) - 59064.*M*pow(r0,4) + 4480.*pow(r0,5))))/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*pow(M,2)*r0*pow(6.*pow(M,2) - 5.*M*r0 + pow(r0,2),1.5)) - (1.6976527263135504*(ellE*(381060.*pow(M,5) - 813020.*pow(M,4)*r0 + 682201.*pow(M,3)*pow(r0,2) - 280947.*pow(M,2)*pow(r0,3) + 56696.*M*pow(r0,4) - 4480.*pow(r0,5)) + ellK*(-471429.*pow(M,5) + 962685.*pow(M,4)*r0 - 772655.*pow(M,3)*pow(r0,2) + 304535.*pow(M,2)*pow(r0,3) - 58936.*M*pow(r0,4) + 4480.*pow(r0,5))))/((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*pow(M,2)*r0*pow(6.*pow(M,2) - 5.*M*r0 + pow(r0,2),1.5)) + (3529.420018005871*(ellE*(456206.*pow(M,5) - 995495.*pow(M,4)*r0 + 847974.*pow(M,3)*pow(r0,2) - 353125.*pow(M,2)*pow(r0,3) + 72008.*M*pow(r0,4) - 5760.*pow(r0,5)) + ellK*(-564510.*pow(M,5) + 1.180009e6*pow(M,4)*r0 - 961168.*pow(M,3)*pow(r0,2) + 383009.*pow(M,2)*pow(r0,3) - 74888.*M*pow(r0,4) + 5760.*pow(r0,5))))/((-7. + 2.*l)*(-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)*(9. + 2.*l)*pow(M,2)*r0*pow(6.*pow(M,2) - 5.*M*r0 + pow(r0,2),1.5));
    const double hP10bar403 = 0.;
    const double hP10bar404 = (0.0006736717167910915*(675. - 31.*pow(1. + 2.*l,2) + 8505./((-1. + 2.*l)*(3. + 2.*l)) + 19845./((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)) + 5.35815e6/((-5. + 2.*l)*(-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)*(7. + 2.*l)))*sqrt(-3.*M + r0)*(-8.*ellK*(5.*M - 2.*r0)*(3.*M - 1.*r0) + ellE*(97.*pow(M,2) - 80.*M*r0 + 16.*pow(r0,2))))/(pow(M,2)*r0*pow(-2.*M + r0,2.5)) + (0.0006736717167910915*(-81. + 4.*pow(1. + 2.*l,2) - 945./((-1. + 2.*l)*(3. + 2.*l)) - 39690./((-3. + 2.*l)*(-1. + 2.*l)*(3. + 2.*l)*(5. + 2.*l)))*(-1.*ellK*(5.*M - 2.*r0)*(3.*M - 1.*r0)*(383.*pow(M,2) - 320.*M*r0 + 64.*pow(r0,2)) + 2.*ellE*(2321.*pow(M,4) - 3855.*pow(M,3)*r0 + 2371.*pow(M,2)*pow(r0,2) - 640.*M*pow(r0,3) + 64.*pow(r0,4))))/(pow(M,2)*r0*sqrt(-3.*M + r0)*pow(-2.*M + r0,3.5));
    const double hP10bar412 = (0.041666666666666664*(1. + 2.*l)*pow(M,2))/(sqrt(r0)*sqrt(-3.*M + r0)*pow(-2.*M + r0,3));
    const double hP10bar413 = 0;

    for(int m=0; m<=l; ++m) {

      const double ld = l;

      /* The Wigner-D matrices */
      const complex<double> w0 = WignerD(l, 0, m);
      const complex<double> w1p = l>=1 ? WignerD(l, 1, m) : 0.0;
      const complex<double> w1m = l>=1 ? WignerD(l, -1, m) : 0.0;
      const complex<double> w2p = l>=2 ? WignerD(l, 2, m) : 0.0;
      const complex<double> w2m = l>=2 ? WignerD(l, -2, m) : 0.0;
      const complex<double> w3p = l>=3 ? WignerD(l, 3, m) : 0.0;
      const complex<double> w3m = l>=3 ? WignerD(l, -3, m) : 0.0;
      const complex<double> w4p = l>=4 ? WignerD(l, 4, m) : 0.0;
      const complex<double> w4m = l>=4 ? WignerD(l, -4, m) : 0.0;
      
      /* The trace-reversed puncture and its first derivative.
       * We don't include the a or 1/r factors here.
       */

      for(size_t j=0; j<N; ++j) {
        const double dr = r[j]-r0;
        const double adr = std::abs(dr);
        const double sdr = sgn(dr);

        const complex<double> hP1bar  = 2.*r[j]*sqrt(M_PI/(2.*ld+1.))*
          (w0*(hP1bar000 + adr*hP1bar010
              + dr*(hP1bar001 + dr*hP1bar002 + adr*hP1bar011
                    + dr*dr*hP1bar003 + dr*dr*dr*hP1bar004
                    + adr*dr*hP1bar012 + adr*dr*dr*hP1bar013))
           +(w2p+w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(hP1bar200
              + dr*(hP1bar201 + dr*hP1bar202
                    + dr*dr*hP1bar203 + dr*dr*dr*hP1bar204
                    + adr*dr*hP1bar212 + adr*dr*dr*hP1bar213)));

        const complex<double> hP2bar  = 2.*sqrt(M_PI/(2.*ld+1.))*
          (w1p-w1m)*sqrt(ld*(ld+1.))*(
            (r[j]-2.*M)*(hP2bar100 + dr*(hP2bar101 + dr*hP2bar102)));

        const complex<double> hP3bar  = 2.*r[j]*r[j]/(r[j]-2.*M)*sqrt(M_PI/(2.*ld+1.))*
          (w0*(hP3bar000 + adr*hP3bar010
              + dr*(hP3bar001 + dr*hP3bar002 + adr*hP3bar011
                    + dr*dr*hP3bar003 + dr*dr*dr*hP3bar004
                    + adr*dr*hP3bar012 + adr*dr*dr*hP3bar013))
           +(w2p+w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(hP3bar200
              + dr*(hP3bar201 + dr*hP3bar202
                    + dr*dr*hP3bar203 + dr*dr*dr*hP3bar204
                    + adr*dr*hP3bar212 + adr*dr*dr*hP3bar213)));

        const complex<double> hP4bar  = 2.*sqrt(M_PI/(2.*ld+1.))*
          ((w1p-w1m)*sqrt(ld*(ld+1.))*(hP4bar100 + adr*hP4bar110
              + dr*(hP4bar101 + dr*hP4bar102 + adr*hP4bar111
                    + dr*dr*hP4bar103 + dr*dr*dr*hP4bar104
                    + adr*dr*hP4bar112 + adr*dr*dr*hP4bar113))
           +(w3p-w3m)*sqrt((ld-2.)*(ld-1.)*ld*(ld+1.)*(ld+2.)*(ld+3.))*(hP4bar300
              + dr*(hP4bar301 + dr*hP4bar302
                    + dr*dr*hP4bar303 + dr*dr*dr*hP4bar304
                    + adr*dr*hP4bar312 + adr*dr*dr*hP4bar313)));

        const complex<double> hP5bar  = 2.*sqrt(M_PI/(2.*ld+1.))*(1.-2.*M/r[j])*
          (w0*(ld*(ld+1.))*(hP5bar000 + dr*(hP5bar001 + dr*hP5bar002))
           +(w2p+w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(hP5bar200
              + dr*(hP5bar201 + dr*hP5bar202)));

        const complex<double> hP6bar  = 2./r[j]*sqrt(M_PI/(2.*ld+1.))*
          (w0*(hP6bar000 + adr*hP6bar010
              + dr*(hP6bar001 + dr*hP6bar002 + adr*hP6bar011
                    + dr*dr*hP6bar003 + dr*dr*dr*hP6bar004
                    + adr*dr*hP6bar012 + adr*dr*dr*hP6bar013))
           +(w2p+w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(hP6bar200
              + dr*(hP6bar201 + dr*hP6bar202
                    + dr*dr*hP6bar203 + dr*dr*dr*hP6bar204
                    + adr*dr*hP6bar212 + adr*dr*dr*hP6bar213)));

        const complex<double> hP7bar  = 2./r[j]*sqrt(M_PI/(2.*ld+1.))*
          (((ld-1.)*ld*(ld+1.)*(ld+2.))*w0*(hP7bar000
              + dr*(hP7bar001 + dr*hP7bar002
                    + dr*dr*hP7bar003 + dr*dr*dr*hP7bar004
                    + adr*dr*hP7bar012 + adr*dr*dr*hP7bar013))
           +(w2p+w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(hP7bar200
              + adr*hP7bar210 + dr*(hP7bar201 + dr*hP7bar202 + adr*hP7bar211
                    + dr*dr*hP7bar203 + dr*dr*dr*hP7bar204
                    + adr*dr*hP7bar212 + adr*dr*dr*hP7bar213))
           +(w4p+w4m)*
            sqrt((ld-3.)*(ld-2.)*(ld-1.)*ld*(ld+1.)*(ld+2.)*(ld+3.)*(ld+4.))*
            (hP7bar400 + dr*(hP7bar401 + dr*hP7bar402
                    + dr*dr*hP7bar403 + dr*dr*dr*hP7bar404
                    + adr*dr*hP7bar412 + adr*dr*dr*hP7bar413)));

        const complex<double> hP8bar  = 2.*I*sqrt(M_PI/(2.*ld+1.))*
          ((w1p+w1m)*sqrt(ld*(ld+1.))*(hP8bar100 + adr*hP8bar110
              + dr*(hP8bar101 + dr*hP8bar102 + adr*hP8bar111
                    + dr*dr*hP8bar103 + dr*dr*dr*hP8bar104
                    + adr*dr*hP8bar112 + adr*dr*dr*hP8bar113))
           +(w3p+w3m)*sqrt((ld-2.)*(ld-1.)*ld*(ld+1.)*(ld+2.)*(ld+3.))*(hP8bar300
              + dr*(hP8bar301 + dr*hP8bar302
                    + dr*dr*hP8bar303 + dr*dr*dr*hP8bar304
                    + adr*dr*hP8bar312 + adr*dr*dr*hP8bar313)));

        const complex<double> hP9bar  = 2.*I*sqrt(M_PI/(2.*ld+1.))*(1.-2.*M/r[j])*
          ((w2p-w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(hP9bar200
              + dr*(hP9bar201 + dr*hP9bar202)));

        const complex<double> hP10bar  = 2.*I/r[j]*sqrt(M_PI/(2.*ld+1.))*
          ((w2p-w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(hP10bar200
              + adr*hP10bar210 + dr*(hP10bar201 + dr*hP10bar202 + adr*hP10bar211
                    + dr*dr*hP10bar203 + dr*dr*dr*hP10bar204
                    + adr*dr*hP10bar212 + adr*dr*dr*hP10bar213))
           +(w4p-w4m)*
            sqrt((ld-3.)*(ld-2.)*(ld-1.)*ld*(ld+1.)*(ld+2.)*(ld+3.)*(ld+4.))*
            (hP10bar400 + dr*(hP10bar401 + dr*hP10bar402
                    + dr*dr*hP10bar403 + dr*dr*dr*hP10bar404
                    + adr*dr*hP10bar412 + adr*dr*dr*hP10bar413)));

        const complex<double> dhP1bar  = 2.*sqrt(M_PI/(2.*ld+1.))*(
          (w0*(hP1bar000 + adr*hP1bar010
              + dr*(hP1bar001 + dr*hP1bar002 + adr*hP1bar011
                    + dr*dr*hP1bar003 + dr*dr*dr*hP1bar004
                    + adr*dr*hP1bar012 + adr*dr*dr*hP1bar013))
            +(w2p+w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(hP1bar200
              + dr*(hP1bar201 + dr*hP1bar202
                    + dr*dr*hP1bar203 + dr*dr*dr*hP1bar204
                    + adr*dr*hP1bar212 + adr*dr*dr*hP1bar213)))
          +r[j]*(w0*(sdr*hP1bar010 + hP1bar001
              + 2.0*dr*hP1bar002 + 2.0*adr*hP1bar011
                    + 3.0*dr*dr*hP1bar003 + 4.0*dr*dr*dr*hP1bar004
                    + 3.0*adr*dr*hP1bar012 + 4.0*adr*dr*dr*hP1bar013)
            +(w2p+w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(hP1bar201
              + 2.0*dr*hP1bar202
                    + 3.0*dr*dr*hP1bar203 + 4.0*dr*dr*dr*hP1bar204
                    + 3.0*adr*dr*hP1bar212 + 4.0*adr*dr*dr*hP1bar213)));

        const complex<double> dhP2bar  = 2.*sqrt(M_PI/(2.*ld+1.))*
          (w1p-w1m)*sqrt(ld*(ld+1.))*(
            (r[j]-2.*M)*(hP2bar101 + 2.*dr*hP2bar102)
            + (hP2bar100 + dr*(hP2bar101 + dr*hP2bar102)));

        const complex<double> dhP3bar  = 2.*sqrt(M_PI/(2.*ld+1.))*(
          r[j]*(r[j]-4.*M)/pow(r[j]-2.*M,2)*(w0*(hP3bar000 + adr*hP3bar010
              + dr*(hP3bar001 + dr*hP3bar002 + adr*hP3bar011
                    + dr*dr*hP3bar003 + dr*dr*dr*hP3bar004
                    + adr*dr*hP3bar012 + adr*dr*dr*hP3bar013))
            +(w2p+w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(hP3bar200
              + dr*(hP3bar201 + dr*hP3bar202
                    + dr*dr*hP3bar203 + dr*dr*dr*hP3bar204
                    + adr*dr*hP3bar212 + adr*dr*dr*hP3bar213)))
          +r[j]*r[j]/(r[j]-2.*M)*(w0*(sdr*hP3bar010 + hP3bar001
              + 2.0*dr*hP3bar002 + 2.0*adr*hP3bar011
                    + 3.0*dr*dr*hP3bar003 + 4.0*dr*dr*dr*hP3bar004
                    + 3.0*adr*dr*hP3bar012 + 4.0*adr*dr*dr*hP3bar013)
            +(w2p+w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(hP3bar201
              + 2.0*dr*hP3bar202
                    + 3.0*dr*dr*hP3bar203 + 4.0*dr*dr*dr*hP3bar204
                    + 3.0*adr*dr*hP3bar212 + 4.0*adr*dr*dr*hP3bar213)));

        const complex<double> dhP4bar  = 2.*sqrt(M_PI/(2.*ld+1.))*(
          +(w1p-w1m)*sqrt(ld*(ld+1.))*(sdr*hP4bar110 + hP4bar101
              + 2.0*dr*hP4bar102 + 2.0*adr*hP4bar111
                    + 3.0*dr*dr*hP4bar103 + 4.0*dr*dr*dr*hP4bar104
                    + 3.0*adr*dr*hP4bar112 + 4.0*adr*dr*dr*hP4bar113)
          +(w3p-w3m)*sqrt((ld-2.)*(ld-1.)*ld*(ld+1.)*(ld+2.)*(ld+3.))*(
              hP4bar301 + 2.0*dr*hP4bar302
                    + 3.0*dr*dr*hP4bar303 + 4.0*dr*dr*dr*hP4bar304
                    + 3.0*adr*dr*hP4bar312 + 4.0*adr*dr*dr*hP4bar313));

        const complex<double> dhP5bar  = 2.*sqrt(M_PI/(2.*ld+1.))*(
          2.*M/pow(r[j],2)*(w0*(ld*(ld+1.))*(
              hP5bar000 + dr*(hP5bar001 + dr*hP5bar002))
           +(w2p+w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(hP5bar200
              + dr*(hP5bar201 + dr*hP5bar202)))
          +(1.-2.*M/r[j])*(w0*(ld*(ld+1.))*(hP5bar001 + 2.*dr*hP5bar002)
           +(w2p+w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(
              hP5bar201 + 2.*dr*hP5bar202)));

        const complex<double> dhP6bar  = 2.*sqrt(M_PI/(2.*ld+1.))*(
          -1./pow(r[j],2)*(w0*(hP6bar000 + adr*hP6bar010
              + dr*(hP6bar001 + dr*hP6bar002 + adr*hP6bar011
                    + dr*dr*hP6bar003 + dr*dr*dr*hP6bar004
                    + adr*dr*hP6bar012 + adr*dr*dr*hP6bar013))
            +(w2p+w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(hP6bar200
              + dr*(hP6bar201 + dr*hP6bar202
                    + dr*dr*hP6bar203 + dr*dr*dr*hP6bar204
                    + adr*dr*hP6bar212 + adr*dr*dr*hP6bar213)))
          +1./r[j]*(w0*(sdr*hP6bar010 + hP6bar001
              + 2.0*dr*hP6bar002 + 2.0*adr*hP6bar011
                    + 3.0*dr*dr*hP6bar003 + 4.0*dr*dr*dr*hP6bar004
                    + 3.0*adr*dr*hP6bar012 + 4.0*adr*dr*dr*hP6bar013)
            +(w2p+w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(hP6bar201
              + 2.0*dr*hP6bar202
                    + 3.0*dr*dr*hP6bar203 + 4.0*dr*dr*dr*hP6bar204
                    + 3.0*adr*dr*hP6bar212 + 4.0*adr*dr*dr*hP6bar213)));

        const complex<double> dhP7bar  = 2.*sqrt(M_PI/(2.*ld+1.))*(
          -1./pow(r[j],2)*(((ld-1.)*ld*(ld+1.)*(ld+2.))*w0*(hP7bar000
              + dr*(hP7bar001 + dr*hP7bar002
                    + dr*dr*hP7bar003 + dr*dr*dr*hP7bar004
                    + adr*dr*hP7bar012 + adr*dr*dr*hP7bar013))
           +(w2p+w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(hP7bar200
              + adr*hP7bar210 + dr*(hP7bar201 + dr*hP7bar202 + adr*hP7bar211
                    + dr*dr*hP7bar203 + dr*dr*dr*hP7bar204
                    + adr*dr*hP7bar212 + adr*dr*dr*hP7bar213))
           +(w4p+w4m)*
            sqrt((ld-3.)*(ld-2.)*(ld-1.)*ld*(ld+1.)*(ld+2.)*(ld+3.)*(ld+4.))*
            (hP7bar400 + dr*(hP7bar401 + dr*hP7bar402
                    + dr*dr*hP7bar403 + dr*dr*dr*hP7bar404
                    + adr*dr*hP7bar412 + adr*dr*dr*hP7bar413)))
          +1./r[j]*(((ld-1.)*ld*(ld+1.)*(ld+2.))*w0*(
              hP7bar001 + 2.0*dr*hP7bar002
                    + 3.0*dr*dr*hP7bar003 + 4.0*dr*dr*dr*hP7bar004
                    + 3.0*adr*dr*hP7bar012 + 4.0*adr*dr*dr*hP7bar013)
           +(w2p+w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(sdr*hP7bar210
              + hP7bar201 + 2.0*dr*hP7bar202 + 2.0*adr*hP7bar211
                    + 3.0*dr*dr*hP7bar203 + 4.0*dr*dr*dr*hP7bar204
                    + 3.0*adr*dr*hP7bar212 + 4.0*adr*dr*dr*hP7bar213)
           +(w4p+w4m)*
            sqrt((ld-3.)*(ld-2.)*(ld-1.)*ld*(ld+1.)*(ld+2.)*(ld+3.)*(ld+4.))*
              (hP7bar401 + 2.*dr*hP7bar402
                    + 3.0*dr*dr*hP7bar403 + 4.0*dr*dr*dr*hP7bar404
                    + 3.0*adr*dr*hP7bar412 + 4.0*adr*dr*dr*hP7bar413)));

        const complex<double> dhP8bar  = 2.*I*sqrt(M_PI/(2.*ld+1.))*(
          +(w1p+w1m)*sqrt(ld*(ld+1.))*(sdr*hP8bar110 + hP8bar101
              + 2.0*dr*hP8bar102 + 2.0*adr*hP8bar111
                    + 3.0*dr*dr*hP8bar103 + 4.0*dr*dr*dr*hP8bar104
                    + 3.0*adr*dr*hP8bar112 + 4.0*adr*dr*dr*hP8bar113)
          +(w3p+w3m)*sqrt((ld-2.)*(ld-1.)*ld*(ld+1.)*(ld+2.)*(ld+3.))*(
              hP8bar301 + 2.0*dr*hP8bar302
                    + 3.0*dr*dr*hP8bar303 + 4.0*dr*dr*dr*hP8bar304
                    + 3.0*adr*dr*hP8bar312 + 4.0*adr*dr*dr*hP8bar313));

        const complex<double> dhP9bar  = 2.*I*sqrt(M_PI/(2.*ld+1.))*(
          2.*M/pow(r[j],2)*((w2p-w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(
            hP9bar200 + dr*(hP9bar201 + dr*hP9bar202)))
          +(1.-2.*M/r[j])*((w2p-w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(
              hP9bar201 + 2.*dr*hP9bar202)));

        const complex<double> dhP10bar  = 2.*I*sqrt(M_PI/(2.*ld+1.))*(
          -1./pow(r[j],2)*(
           (w2p-w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(hP10bar200
              + adr*hP10bar210 + dr*(hP10bar201 + dr*hP10bar202 + adr*hP10bar211
                    + dr*dr*hP10bar203 + dr*dr*dr*hP10bar204
                    + adr*dr*hP10bar212 + adr*dr*dr*hP10bar213))
           +(w4p-w4m)*
            sqrt((ld-3.)*(ld-2.)*(ld-1.)*ld*(ld+1.)*(ld+2.)*(ld+3.)*(ld+4.))*
            (hP10bar400 + dr*(hP10bar401 + dr*hP10bar402
                    + dr*dr*hP10bar403 + dr*dr*dr*hP10bar404
                    + adr*dr*hP10bar412 + adr*dr*dr*hP10bar413)))
          +1./r[j]*(
           (w2p-w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(sdr*hP10bar210
              + hP10bar201 + 2.0*dr*hP10bar202 + 2.0*adr*hP10bar211
                    + 3.0*dr*dr*hP10bar203 + 4.0*dr*dr*dr*hP10bar204
                    + 3.0*adr*dr*hP10bar212 + 4.0*adr*dr*dr*hP10bar213)
           +(w4p-w4m)*
            sqrt((ld-3.)*(ld-2.)*(ld-1.)*ld*(ld+1.)*(ld+2.)*(ld+3.)*(ld+4.))*
              (hP10bar401 + 2.*dr*hP10bar402
                    + 3.0*dr*dr*hP10bar403 + 4.0*dr*dr*dr*hP10bar404
                    + 3.0*adr*dr*hP10bar412 + 4.0*adr*dr*dr*hP10bar413)));

        const complex<double> ddhP1bar  = 2.*sqrt(M_PI/(2.*ld+1.))*(
          2.0*(w0*(sdr*hP1bar010 + hP1bar001
              + 2.0*dr*hP1bar002 + 2.0*adr*hP1bar011
                    + 3.0*dr*dr*hP1bar003 + 4.0*dr*dr*dr*hP1bar004
                    + 3.0*adr*dr*hP1bar012 + 4.0*adr*dr*dr*hP1bar013)
            +(w2p+w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(
              + hP1bar201 + 2.0*dr*hP1bar202
                    + 3.0*dr*dr*hP1bar203 + 4.0*dr*dr*dr*hP1bar204
                    + 3.0*adr*dr*hP1bar212 + 4.0*adr*dr*dr*hP1bar213))
          +r[j]*(w0*(2.0*hP1bar002 + 2.0*sdr*hP1bar011
                    + 6.0*dr*hP1bar003 + 12.0*dr*dr*hP1bar004
                    + 6.0*adr*hP1bar012 + 12.0*adr*dr*hP1bar013)
            +(w2p+w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(2.0*hP1bar202
                    + 6.0*dr*hP1bar203 + 12.0*dr*dr*hP1bar204
                    + 6.0*adr*hP1bar212 + 12.0*adr*dr*hP1bar213)));

        const complex<double> ddhP2bar  = 2.*sqrt(M_PI/(2.*ld+1.))*
          (w1p-w1m)*sqrt(ld*(ld+1.))*(
            2.*(hP2bar101 + 2.*dr*hP2bar102)
            +(r[j]-2.*M)*(2.*hP2bar102));

        const complex<double> ddhP3bar  = 2.*sqrt(M_PI/(2.*ld+1.))*(
          2.0*r[j]*(r[j]-4.*M)/pow(r[j]-2.*M,2)*(w0*(sdr*hP3bar010 + hP3bar001
              + 2.0*dr*hP3bar002 + 2.0*adr*hP3bar011
                    + 3.0*dr*dr*hP3bar003 + 4.0*dr*dr*dr*hP3bar004
                    + 3.0*adr*dr*hP3bar012 + 4.0*adr*dr*dr*hP3bar013)
            +(w2p+w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(
              + hP3bar201 + 2.0*dr*hP3bar202
                    + 3.0*dr*dr*hP3bar203 + 4.0*dr*dr*dr*hP3bar204
                    + 3.0*adr*dr*hP3bar212 + 4.0*adr*dr*dr*hP3bar213))
          +r[j]*r[j]/(r[j]-2.*M)*(w0*(2.0*hP3bar002 + 2.0*sdr*hP3bar011
                    + 6.0*dr*hP3bar003 + 12.0*dr*dr*hP3bar004
                    + 6.0*adr*hP3bar012 + 12.0*adr*dr*hP3bar013)
            +(w2p+w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(2.0*hP3bar202
                    + 6.0*dr*hP3bar203 + 12.0*dr*dr*hP3bar204
                    + 6.0*adr*hP3bar212 + 12.0*adr*dr*hP3bar213))
          +8.*M*M/pow(r[j]-2.*M,3)*(w0*(hP3bar000 + adr*hP3bar010
              + dr*(hP3bar001 + dr*hP3bar002 + adr*hP3bar011
                    + dr*dr*hP3bar003 + dr*dr*dr*hP3bar004
                    + adr*dr*hP3bar012 + adr*dr*dr*hP3bar013))
           +(w2p+w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(hP3bar200
              + dr*(hP3bar201 + dr*hP3bar202
                    + dr*dr*hP3bar203 + dr*dr*dr*hP3bar204
                    + adr*dr*hP3bar212 + adr*dr*dr*hP3bar213))));

        const complex<double> ddhP4bar  = 2.*sqrt(M_PI/(2.*ld+1.))*(
          +(w1p-w1m)*sqrt(ld*(ld+1.))*(2.0*hP4bar102 + 2.0*sdr*hP4bar111
                    + 6.0*dr*hP4bar103 + 12.0*dr*dr*hP4bar104
                    + 6.0*adr*hP4bar112 + 12.0*adr*dr*hP4bar113)
          +(w3p-w3m)*sqrt((ld-2.)*(ld-1.)*ld*(ld+1.)*(ld+2.)*(ld+3.))*(
            2.0*hP4bar302
                    + 6.0*dr*hP4bar303 + 12.0*dr*dr*hP4bar304
                    + 6.0*adr*hP4bar312 + 12.0*adr*dr*hP4bar313));

        const complex<double> ddhP5bar  = 2.*sqrt(M_PI/(2.*ld+1.))*(
          -4.*M/pow(r[j],3)*(w0*(ld*(ld+1.))*(
              hP5bar000 + dr*(hP5bar001 + dr*hP5bar002))
           +(w2p+w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(hP5bar200
              + dr*(hP5bar201 + dr*hP5bar202)))
          +4.*M/pow(r[j],2)*(w0*(ld*(ld+1.))*(
              hP5bar001 + 2.*dr*hP5bar002)
           +(w2p+w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(
              hP5bar201 + 2.*dr*hP5bar202))
          +(1.-2.*M/r[j])*(w0*(ld*(ld+1.))*(2.*hP5bar002)
           +(w2p+w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(
              2.*hP5bar202)));

        const complex<double> ddhP6bar  = 2.*sqrt(M_PI/(2.*ld+1.))*(
          -2./pow(r[j],2)*(w0*(sdr*hP6bar010 + hP6bar001
              + 2.0*dr*hP6bar002 + 2.0*adr*hP6bar011
                    + 3.0*dr*dr*hP6bar003 + 4.0*dr*dr*dr*hP6bar004
                    + 3.0*adr*dr*hP6bar012 + 4.0*adr*dr*dr*hP6bar013)
            +(w2p+w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(
              + hP6bar201 + 2.0*dr*hP6bar202
                    + 3.0*dr*dr*hP6bar203 + 4.0*dr*dr*dr*hP6bar204
                    + 3.0*adr*dr*hP6bar212 + 4.0*adr*dr*dr*hP6bar213))
          +1./r[j]*(w0*(2.0*hP6bar002 + 2.0*sdr*hP6bar011
                    + 6.0*dr*hP6bar003 + 12.0*dr*dr*hP6bar004
                    + 6.0*adr*hP6bar012 + 12.0*adr*dr*hP6bar013)
            +(w2p+w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(2.0*hP6bar202
                    + 6.0*dr*hP6bar203 + 12.0*dr*dr*hP6bar204
                    + 6.0*adr*hP6bar212 + 12.0*adr*dr*hP6bar213))
          +2./pow(r[j],3)*(w0*(hP6bar000 + adr*hP6bar010
              + dr*(hP6bar001 + dr*hP6bar002 + adr*hP6bar011
                    + dr*dr*hP6bar003 + dr*dr*dr*hP6bar004
                    + adr*dr*hP6bar012 + adr*dr*dr*hP6bar013))
           +(w2p+w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(hP6bar200
              + dr*(hP6bar201 + dr*hP6bar202
                    + dr*dr*hP6bar203 + dr*dr*dr*hP6bar204
                    + adr*dr*hP6bar212 + adr*dr*dr*hP6bar213))));

        const complex<double> ddhP7bar  = 2.*sqrt(M_PI/(2.*ld+1.))*(
          -2./pow(r[j],2)*(((ld-1.)*ld*(ld+1.)*(ld+2.))*w0*(
              hP7bar001 + 2.0*dr*hP7bar002
                    + 3.0*dr*dr*hP7bar003 + 4.0*dr*dr*dr*hP7bar004
                    + 3.0*adr*dr*hP7bar012 + 4.0*adr*dr*dr*hP7bar013)
           +(w2p+w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(sdr*hP7bar210
              + hP7bar201 + 2.0*dr*hP7bar202 + 2.0*adr*hP7bar211
                    + 3.0*dr*dr*hP7bar203 + 4.0*dr*dr*dr*hP7bar204
                    + 3.0*adr*dr*hP7bar212 + 4.0*adr*dr*dr*hP7bar213)
           +(w4p+w4m)*
            sqrt((ld-3.)*(ld-2.)*(ld-1.)*ld*(ld+1.)*(ld+2.)*(ld+3.)*(ld+4.))*
            (hP7bar401 + 2.*dr*hP7bar402
                    + 3.0*dr*dr*hP7bar403 + 4.0*dr*dr*dr*hP7bar404
                    + 3.0*adr*dr*hP7bar412 + 4.0*adr*dr*dr*hP7bar413))
          +1./r[j]*(((ld-1.)*ld*(ld+1.)*(ld+2.))*w0*(2.*hP7bar002
                    + 6.0*dr*hP7bar003 + 12.0*dr*dr*hP7bar004
                    + 6.0*adr*hP7bar012 + 12.0*adr*dr*hP7bar013)
           +(w2p+w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(
              2.*hP7bar202 + 2.*sdr*hP7bar211
                    + 6.0*dr*hP7bar203 + 12.0*dr*dr*hP7bar204
                    + 6.0*adr*hP7bar212 + 12.0*adr*dr*hP7bar213)
           +(w4p+w4m)*
            sqrt((ld-3.)*(ld-2.)*(ld-1.)*ld*(ld+1.)*(ld+2.)*(ld+3.)*(ld+4.))*
            (2.*hP7bar402
                    + 6.0*dr*hP7bar403 + 12.0*dr*dr*hP7bar404
                    + 6.0*adr*hP7bar412 + 12.0*adr*dr*hP7bar413))
          +2./pow(r[j],3)*(((ld-1.)*ld*(ld+1.)*(ld+2.))*w0*(hP7bar000
              + dr*(hP7bar001 + dr*hP7bar002
                    + dr*dr*hP7bar003 + dr*dr*dr*hP7bar004
                    + adr*dr*hP7bar012 + adr*dr*dr*hP7bar013))
           +(w2p+w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(hP7bar200
              + adr*hP7bar210 + dr*(hP7bar201 + dr*hP7bar202 + adr*hP7bar211
                    + dr*dr*hP7bar203 + dr*dr*dr*hP7bar204
                    + adr*dr*hP7bar212 + adr*dr*dr*hP7bar213))
           +(w4p+w4m)*
            sqrt((ld-3.)*(ld-2.)*(ld-1.)*ld*(ld+1.)*(ld+2.)*(ld+3.)*(ld+4.))*
            (hP7bar400 + dr*(hP7bar401 + dr*hP7bar402
                    + dr*dr*hP7bar403 + dr*dr*dr*hP7bar404
                    + adr*dr*hP7bar412 + adr*dr*dr*hP7bar413))));

        const complex<double> ddhP8bar  = 2.*I*sqrt(M_PI/(2.*ld+1.))*(
          +(w1p+w1m)*sqrt(ld*(ld+1.))*(2.0*hP8bar102 + 2.0*sdr*hP8bar111
                    + 6.0*dr*hP8bar103 + 12.0*dr*dr*hP8bar104
                    + 6.0*adr*hP8bar112 + 12.0*adr*dr*hP8bar113)
          +(w3p+w3m)*sqrt((ld-2.)*(ld-1.)*ld*(ld+1.)*(ld+2.)*(ld+3.))*(
            2.0*hP8bar302
                    + 6.0*dr*hP8bar303 + 12.0*dr*dr*hP8bar304
                    + 6.0*adr*hP8bar312 + 12.0*adr*dr*hP8bar313));

        const complex<double> ddhP9bar  = 2.*I*sqrt(M_PI/(2.*ld+1.))*(
          -4.*M/pow(r[j],3)*((w2p-w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(
            hP9bar200 + dr*(hP9bar201 + dr*hP9bar202)))
          +4.*M/pow(r[j],2)*((w2p-w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(
            hP9bar201 + 2.*dr*hP9bar202))
          +(1.-2.*M/r[j])*((w2p-w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(
            2.*hP9bar202)));

        const complex<double> ddhP10bar  = 2.*I*sqrt(M_PI/(2.*ld+1.))*(
          -2./pow(r[j],2)*(
           (w2p-w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(sdr*hP10bar210
              + hP10bar201 + 2.0*dr*hP10bar202 + 2.0*adr*hP10bar211
                    + 3.0*dr*dr*hP10bar203 + 4.0*dr*dr*dr*hP10bar204
                    + 3.0*adr*dr*hP10bar212 + 4.0*adr*dr*dr*hP10bar213)
           +(w4p-w4m)*
            sqrt((ld-3.)*(ld-2.)*(ld-1.)*ld*(ld+1.)*(ld+2.)*(ld+3.)*(ld+4.))*
            (hP10bar401 + 2.*dr*hP10bar402
                    + 3.0*dr*dr*hP10bar403 + 4.0*dr*dr*dr*hP10bar404
                    + 3.0*adr*dr*hP10bar412 + 4.0*adr*dr*dr*hP10bar413))
          +1./r[j]*(
           (w2p-w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(
              2.*hP10bar202 + 2.*sdr*hP10bar211
                    + 6.0*dr*hP10bar203 + 12.0*dr*dr*hP10bar204
                    + 6.0*adr*hP10bar212 + 12.0*adr*dr*hP10bar213)
           +(w4p-w4m)*
            sqrt((ld-3.)*(ld-2.)*(ld-1.)*ld*(ld+1.)*(ld+2.)*(ld+3.)*(ld+4.))*
            (2.*hP10bar402
                    + 6.0*dr*hP10bar403 + 12.0*dr*dr*hP10bar404
                    + 6.0*adr*hP10bar412 + 12.0*adr*dr*hP10bar413))
          +2./pow(r[j],3)*(
           (w2p-w2m)*sqrt((ld-1.)*ld*(ld+1.)*(ld+2.))*(hP10bar200
              + adr*hP10bar210 + dr*(hP10bar201 + dr*hP10bar202 + adr*hP10bar211
                    + dr*dr*hP10bar203 + dr*dr*dr*hP10bar204
                    + adr*dr*hP10bar212 + adr*dr*dr*hP10bar213))
           +(w4p-w4m)*
            sqrt((ld-3.)*(ld-2.)*(ld-1.)*ld*(ld+1.)*(ld+2.)*(ld+3.)*(ld+4.))*
            (hP10bar400 + dr*(hP10bar401 + dr*hP10bar402
                    + dr*dr*hP10bar403 + dr*dr*dr*hP10bar404
                    + adr*dr*hP10bar412 + adr*dr*dr*hP10bar413))));
        
        /* The non-trace-reversed field and its first and second derivatives.
         * Here we do include the factor of a and 1/r, and the window function.
         * Trace-reversal corresponds to swapping 3 and 6. */

        /* Step-function window */
        const double ra = r0 - sigma;
        const double rb = r0 + sigma;
        const double W = (r[j]<ra || r[j]>rb) ? 0.0 : 1.0;
        const double dW = (r[j]<ra || r[j]>rb) ? 0.0 : 0.0;
        const double ddW = (r[j]<ra || r[j]>rb) ? 0.0 : 0.0;

        /* Gaussian window, order n */
        // const int n = 4;
        // const double W = exp(-0.5*pow((r[j]-r0)/sigma,n));
        // const double dW = - 0.5 * W * n * pow((r[j]-r0)/sigma,n-1) / sigma;
        // const double ddW = 0.25 * W * n * (2.0 + n * (-2.0 + pow((r[j]-r0)/sigma,n))) * pow((r[j]-r0),n-2) * pow(sigma,-n);

        hP[1][l][m][j] = a_il(1,l)*W*hP1bar/r[j];
        hP[2][l][m][j] = a_il(2,l)*W*hP2bar/r[j];
        hP[3][l][m][j] = a_il(3,l)*W*hP6bar/r[j];
        hP[4][l][m][j] = a_il(4,l)*W*hP4bar/r[j];
        hP[5][l][m][j] = a_il(5,l)*W*hP5bar/r[j];
        hP[6][l][m][j] = a_il(6,l)*W*hP3bar/r[j];
        hP[7][l][m][j] = a_il(7,l)*W*hP7bar/r[j];
        hP[8][l][m][j] = a_il(8,l)*W*hP8bar/r[j];
        hP[9][l][m][j] = a_il(9,l)*W*hP9bar/r[j];
        hP[10][l][m][j] = a_il(10,l)*W*hP10bar/r[j];

        dhP[1][l][m][j] = a_il(1,l)*(W*(dhP1bar - hP1bar/r[j])/r[j] + dW*hP1bar/r[j]);
        dhP[2][l][m][j] = a_il(2,l)*(W*(dhP2bar - hP2bar/r[j])/r[j] + dW*hP2bar/r[j]);
        dhP[3][l][m][j] = a_il(3,l)*(W*(dhP6bar - hP6bar/r[j])/r[j] + dW*hP6bar/r[j]);
        dhP[4][l][m][j] = a_il(4,l)*(W*(dhP4bar - hP4bar/r[j])/r[j] + dW*hP4bar/r[j]);
        dhP[5][l][m][j] = a_il(5,l)*(W*(dhP5bar - hP5bar/r[j])/r[j] + dW*hP5bar/r[j]);
        dhP[6][l][m][j] = a_il(6,l)*(W*(dhP3bar - hP3bar/r[j])/r[j] + dW*hP3bar/r[j]);
        dhP[7][l][m][j] = a_il(7,l)*(W*(dhP7bar - hP7bar/r[j])/r[j] + dW*hP7bar/r[j]);
        dhP[8][l][m][j] = a_il(8,l)*(W*(dhP8bar - hP8bar/r[j])/r[j] + dW*hP8bar/r[j]);
        dhP[9][l][m][j] = a_il(9,l)*(W*(dhP9bar - hP9bar/r[j])/r[j] + dW*hP9bar/r[j]);
        dhP[10][l][m][j] = a_il(10,l)*(W*(dhP10bar - hP10bar/r[j])/r[j] + dW*hP10bar/r[j]);

        ddhP[1][l][m][j] = a_il(1,l)*(W*(ddhP1bar - 2.0*(dhP1bar - hP1bar/r[j])/r[j])/r[j] + 2.0*dW*(dhP1bar - hP1bar/r[j])/r[j] + ddW*hP1bar/r[j]);
        ddhP[2][l][m][j] = a_il(2,l)*(W*(ddhP2bar - 2.0*(dhP2bar - hP2bar/r[j])/r[j])/r[j] + 2.0*dW*(dhP2bar - hP2bar/r[j])/r[j] + ddW*hP2bar/r[j]);
        ddhP[3][l][m][j] = a_il(3,l)*(W*(ddhP6bar - 2.0*(dhP6bar - hP6bar/r[j])/r[j])/r[j] + 2.0*dW*(dhP6bar - hP6bar/r[j])/r[j] + ddW*hP6bar/r[j]);
        ddhP[4][l][m][j] = a_il(4,l)*(W*(ddhP4bar - 2.0*(dhP4bar - hP4bar/r[j])/r[j])/r[j] + 2.0*dW*(dhP4bar - hP4bar/r[j])/r[j] + ddW*hP4bar/r[j]);
        ddhP[5][l][m][j] = a_il(5,l)*(W*(ddhP5bar - 2.0*(dhP5bar - hP5bar/r[j])/r[j])/r[j] + 2.0*dW*(dhP5bar - hP5bar/r[j])/r[j] + ddW*hP5bar/r[j]);
        ddhP[6][l][m][j] = a_il(6,l)*(W*(ddhP3bar - 2.0*(dhP3bar - hP3bar/r[j])/r[j])/r[j] + 2.0*dW*(dhP3bar - hP3bar/r[j])/r[j] + ddW*hP3bar/r[j]);
        ddhP[7][l][m][j] = a_il(7,l)*(W*(ddhP7bar - 2.0*(dhP7bar - hP7bar/r[j])/r[j])/r[j] + 2.0*dW*(dhP7bar - hP7bar/r[j])/r[j] + ddW*hP7bar/r[j]);
        ddhP[8][l][m][j] = a_il(8,l)*(W*(ddhP8bar - 2.0*(dhP8bar - hP8bar/r[j])/r[j])/r[j] + 2.0*dW*(dhP8bar - hP8bar/r[j])/r[j] + ddW*hP8bar/r[j]);
        ddhP[9][l][m][j] = a_il(9,l)*(W*(ddhP9bar - 2.0*(dhP9bar - hP9bar/r[j])/r[j])/r[j] + 2.0*dW*(dhP9bar - hP9bar/r[j])/r[j] + ddW*hP9bar/r[j]);
        ddhP[10][l][m][j] = a_il(10,l)*(W*(ddhP10bar - 2.0*(dhP10bar - hP10bar/r[j])/r[j])/r[j] + 2.0*dW*(dhP10bar - hP10bar/r[j])/r[j] + ddW*hP10bar/r[j]);

        /* Negative m modes */
        const double sign = isOdd(m) ? -1.0 : 1.0;
        for(int i=1; i<=10; ++i) {
          hP[i][l][-m][j] = sign*conj(hP[i][l][m][j]);
          dhP[i][l][-m][j] = sign*conj(dhP[i][l][m][j]);
          ddhP[i][l][-m][j] = sign*conj(ddhP[i][l][m][j]);
        }
      }
    }
  }
}
