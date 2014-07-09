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
#include "WignerDMatrices.hpp"

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

/* Compute the first-order puncture fields */
void compute_h1P(const double r0, const vector<double> &r, const int l_max, field_type &hP, field_type &dhP, field_type &ddhP)
{
  cout << "Computing first order punctures " << endl;
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
  SphericalFunctions::WignerDMatrix WignerD(Quaternions::Quaternion(M_PI, M_PI_2, M_PI_2));

  /* Worldtube size */
  const double ra = r0 - 2.0*M;
  const double rb = r0 + 2.0*M;

  const complex<double> I(0.0, 1.0);

  /* Loop over all l, m modes */
  for(int l=0; l<=l_max; ++l) {
    for(int m=0; m<=l; ++m) {

      const double ld = l;

      /* The Wigner-D matrices - note that our convention is different from that
       * provided by the WignerDMatrix class. The difference corresponds to a
       * complex conjugation and a change in the sign of m and m'. */
      const complex<double> w0 = conj(WignerD(l, -m, -0));
      const complex<double> w1p = l>=1 ? conj(WignerD(l, -m, -1)) : 0.0;
      const complex<double> w1m = l>=1 ? conj(WignerD(l, -m, -(-1))) : 0.0;
      const complex<double> w2p = l>=2 ? conj(WignerD(l, -m, -2)) : 0.0;
      const complex<double> w2m = l>=2 ? conj(WignerD(l, -m, -(-2))) : 0.0;
      
      /* The trace-reversed puncture and its first derivative.
       * We don't include the a or 1/r factors here.
       */

      for(size_t j=0; j<N; ++j) {
        const double pm = r[j] <= r0 ? -1.0 : 1.0;
        const complex<double> hP1bar = (4.*w0*r[j]*(2.*ellE*pow(-2.*M + r0,1.5)*(-1.*r0 + r[j]) - 3.141592653589793*pm*sqrt(r0)*(-2.*(1. + 2.*ld)*M + r0 + 2.*ld*r0)*(-1.*r0 + r[j]) - 4.*ellK*sqrt(-2.*M + r0)*(r0*(-2.*r0 + r[j]) + 2.*M*(r0 - 2.*(-1.*r0 + r[j])))))/(sqrt(3.141592653589793 + 6.283185307179586*ld)*pow(r0,3)*sqrt(-3.*M + r0));
        const complex<double> hP3bar = (4.*w0*pow(r[j],2)*(2.*ellE*pow(-2.*M + r0,1.5)*(-1.*r0 + r[j]) - 3.141592653589793*pm*sqrt(r0)*(-2.*(1. + 2.*ld)*M + r0 + 2.*ld*r0)*(-1.*r0 + r[j]) - 4.*ellK*sqrt(-2.*M + r0)*(r0*(-2.*r0 + r[j]) + 2.*M*(r0 - 2.*(-1.*r0 + r[j])))))/(sqrt(3.141592653589793 + 6.283185307179586*ld)*pow(r0,3)*sqrt(-3.*M + r0)*(-2.*M + r[j]));
        const complex<double> hP4bar = (2.256758334191025*sqrt(ld*(1. + ld))*(-1.*w1m + w1p)*((4.*sqrt((M*(2.*M - 1.*r0)*r0)/(3.*M - 1.*r0))*((ellK*M)/(ld*(1. + ld)) + 2.*(-1.*ellE + ellK)*(-2.*M + r0)))/(sqrt(1. + 2.*ld)*M) - 3.141592653589793*pm*sqrt((-1.*(M + 2.*ld*M))/(3.*M - 1.*r0))*(-1.*r0 + r[j]) + (2.*(2.*(-2.*M + r0)*(ellK*(5.*M - 1.*r0) + ellE*(-4.*M + r0)) + (M*(2.*ellK*M + ellE*(-2.*M + r0)))/(ld*(1. + ld)))*(-1.*r0 + r[j]))/(sqrt(1. + 2.*ld)*sqrt(M*r0*(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))))))/r0;
        const complex<double> hP6bar = (4.*w0*((4.*ellK*M*r0)/(sqrt(-3.*M + r0)*sqrt(-2.*M + r0)) - (3.141592653589793*(1. + 2.*ld)*M*pm*r0*(-1.*r0 + r[j]))/(sqrt(r0*(-3.*M + r0))*(-2.*M + r0)) + (2.*(ellE + 2.*ellK)*M*(-1.*r0 + r[j]))/(sqrt(-3.*M + r0)*sqrt(-2.*M + r0))))/(sqrt(3.141592653589793 + 6.283185307179586*ld)*r[j]);
        const complex<double> hP7bar = (0.3761263890318375*sqrt((-1. + ld)*ld*(1. + ld)*(2. + ld))*(w2m + w2p)*((4.*r0*((2.4*ellK*pow(M + 2.*ld*M,2))/((-1. + ld)*ld*(1. + ld)*(2. + ld)) - 8.*ellE*(-2.*M + r0)*(-5.*M + 2.*r0) + 2.*ellK*(51.*pow(M,2) - 40.*M*r0 + 8.*pow(r0,2))))/(sqrt(1. + 2.*ld)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) + (9.42477796076938*pow(M,2)*pm*sqrt((-1.*(r0 + 2.*ld*r0))/(3.*M - 1.*r0))*(-1.*r0 + r[j]))/(2.*M - 1.*r0) + (2.*((2.4*(ellE + 2.*ellK)*pow(M + 2.*ld*M,2))/((-1. + ld)*ld*(1. + ld)*(2. + ld)) + 2.*(-2.*ellK*(39.*pow(M,2) - 26.*M*r0 + 4.*pow(r0,2)) + ellE*(67.*pow(M,2) - 48.*M*r0 + 8.*pow(r0,2))))*(-1.*r0 + r[j]))/(sqrt(1. + 2.*ld)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0))))/(M*r[j]);
        const complex<double> hP8bar = 2.256758334191025*I*sqrt(ld*(1. + ld))*r0*(w1m + w1p)*((4.*sqrt((M*r0*(-2.*M + r0))/(-3.*M + r0))*(ellE*(4.*M - 2.*r0) - (1.*ellK*((1. + 6.*ld + 6.*pow(ld,2))*M - 2.*ld*(1. + ld)*r0))/(ld*(1. + ld))))/(sqrt(1. + 2.*ld)*M*pow(r0,2)) + (3.141592653589793*pm*sqrt(((1. + 2.*ld)*M)/(-3.*M + r0))*(-1.*r0 + r[j]))/pow(r0,2) + (2.*((-1.*M*(2.*ellK*M + ellE*(-2.*M + r0)))/(ld*(1. + ld)) + 2.*(-1.*ellK*(-4.*M + r0)*(-3.*M + r0) + ellE*(-5.*M + r0)*(-2.*M + r0)))*(-1.*r0 + r[j]))/sqrt((1. + 2.*ld)*M*pow(r0,5)*(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))));
        const complex<double> hP10bar =(-0.3761263890318375*I*sqrt((-1. + ld)*ld*(1. + ld)*(2. + ld))*(-1.*w2m + w2p)*((-4.*r0*((-2.4*ellK*pow(M + 2.*ld*M,2))/((-1. + ld)*ld*(1. + ld)*(2. + ld)) + 16.*ellK*(-3.*M + r0)*(-2.*M + r0) - 8.*ellE*(-2.*M + r0)*(-5.*M + 2.*r0)))/(sqrt(1. + 2.*ld)*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) + (9.42477796076938*pow(M,2)*pm*sqrt((-1.*(r0 + 2.*ld*r0))/(3.*M - 1.*r0))*(-1.*r0 + r[j]))/(2.*M - 1.*r0) - (2.*((-2.4*(ellE + 2.*ellK)*pow(M + 2.*ld*M,2))/((-1. + ld)*ld*(1. + ld)*(2. + ld)) + 8.*(2.*ellE*(-4.*M + r0)*(-2.*M + r0) - 1.*ellK*(-3.*M + r0)*(-7.*M + 2.*r0)))*(-1.*r0 + r[j]))/(sqrt(1. + 2.*ld)*sqrt(-3.*M + r0)*sqrt(-2.*M + r0))))/(M*r[j]);
        
        const complex<double> dhP1bar = (4.*w0*(-3.141592653589793*(1. + 2.*ld)*pm*sqrt(r0)*(-2.*M + r0)*(r0 + 2.*(-1.*r0 + r[j])) + 2.*ellE*pow(-2.*M + r0,1.5)*(r0 + 2.*(-1.*r0 + r[j])) + 8.*ellK*sqrt(-2.*M + r0)*(-1.*r0*(-1.*r0 + r[j]) + M*(r0 + 4.*(-1.*r0 + r[j])))))/(sqrt(3.141592653589793 + 6.283185307179586*ld)*pow(r0,3)*sqrt(-3.*M + r0));
        const complex<double> dhP3bar = (-4.*w0*r[j]*(8.*ellK*sqrt(-2.*M + r0)*(-1.*r0 + r[j])*(12.*pow(M,2) + r0*r[j] - 4.*M*(r0 + r[j])) + 3.141592653589793*(1. + 2.*ld)*pm*sqrt(r0)*(4.*pow(M,2)*(r0 + 3.*(-1.*r0 + r[j])) - 4.*M*(pow(r0,2) + 3.*r0*(-1.*r0 + r[j]) + pow(-1.*r0 + r[j],2)) + r0*(pow(r0,2) + 3.*r0*(-1.*r0 + r[j]) + 2.*pow(-1.*r0 + r[j],2))) - 2.*ellE*sqrt(-2.*M + r0)*(4.*pow(M,2)*(r0 + 3.*(-1.*r0 + r[j])) - 4.*M*(pow(r0,2) + 3.*r0*(-1.*r0 + r[j]) + pow(-1.*r0 + r[j],2)) + r0*(pow(r0,2) + 3.*r0*(-1.*r0 + r[j]) + 2.*pow(-1.*r0 + r[j],2)))))/(sqrt(3.141592653589793 + 6.283185307179586*ld)*pow(r0,3)*sqrt(-3.*M + r0)*pow(-2.*M + r[j],2));
        const complex<double> dhP4bar = (4.*sqrt(ld*(1. + ld))*(-1.7724538509055159*pm*sqrt(((1. + 2.*ld)*M)/(-3.*M + r0)) + (2.*(2.*(-2.*M + r0)*(ellK*(5.*M - 1.*r0) + ellE*(-4.*M + r0)) + (M*(2.*ellK*M + ellE*(-2.*M + r0)))/(ld*(1. + ld))))/(sqrt(3.141592653589793 + 6.283185307179586*ld)*sqrt(M*r0*(6.*pow(M,2) - 5.*M*r0 + pow(r0,2)))))*(-1.*w1m + w1p))/r0;
        const complex<double> dhP6bar = (-4.*M*sqrt(r0*(-3.*M + r0))*(2.*ellE*(2.*M - 1.*r0)*sqrt(r0*(-3.*M + r0)) + 3.141592653589793*(1. + 2.*ld)*pm*r0*sqrt(-3.*M + r0)*sqrt(-2.*M + r0))*w0)/(sqrt(3.141592653589793 + 6.283185307179586*ld)*pow(-3.*M + r0,1.5)*pow(-2.*M + r0,1.5)*pow(r[j],2));
        const complex<double> dhP7bar = (1.4142135623730951*sqrt(ld*(-2. - 1.*ld + 2.*pow(ld,2) + pow(ld,3)))*(w2m + w2p)*((-0.752252778063675*r0*((2.4*ellK*pow(M + 2.*ld*M,2))/((-1. + ld)*ld*(1. + ld)*(2. + ld)) - 8.*ellE*(-2.*M + r0)*(-5.*M + 2.*r0) + 2.*ellK*(51.*pow(M,2) - 40.*M*r0 + 8.*pow(r0,2))))/(sqrt(0.5 + ld)*M*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) + ((2.5066282746310002*M*pm*sqrt(((1. + 2.*ld)*r0)/(-3.*M + r0)))/(2.*M - 1.*r0) + (0.3761263890318375*((2.4*(ellE + 2.*ellK)*pow(M + 2.*ld*M,2))/((-1. + ld)*ld*(1. + ld)*(2. + ld)) + 2.*(-2.*ellK*(39.*pow(M,2) - 26.*M*r0 + 4.*pow(r0,2)) + ellE*(67.*pow(M,2) - 48.*M*r0 + 8.*pow(r0,2)))))/(sqrt(0.5 + ld)*M*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)))*r[j] - (2.5066282746310002*M*pm*sqrt(((1. + 2.*ld)*r0)/(-3.*M + r0))*(-1.*r0 + r[j]))/(2.*M - 1.*r0) - (0.3761263890318375*((2.4*(ellE + 2.*ellK)*pow(M + 2.*ld*M,2))/((-1. + ld)*ld*(1. + ld)*(2. + ld)) + 2.*(-2.*ellK*(39.*pow(M,2) - 26.*M*r0 + 4.*pow(r0,2)) + ellE*(67.*pow(M,2) - 48.*M*r0 + 8.*pow(r0,2))))*(-1.*r0 + r[j]))/(sqrt(0.5 + ld)*M*sqrt(-3.*M + r0)*sqrt(-2.*M + r0))))/pow(r[j],2);
        const complex<double> dhP8bar = 1.4142135623730951*sqrt(ld*(1. + ld))*r0*((5.0132565492620005*I*pm*sqrt(((1. + 2.*ld)*M)/(-3.*M + r0)))/pow(r0,2) + (3.1915382432114616*I*((-1.*M*(2.*ellK*M + ellE*(-2.*M + r0)))/(ld*(1. + ld)) + 2.*(-1.*ellK*(-4.*M + r0)*(-3.*M + r0) + ellE*(-5.*M + r0)*(-2.*M + r0))))/sqrt((1. + 2.*ld)*M*pow(r0,5)*(-3.*M + r0)*(-2.*M + r0)))*(w1m + w1p);
        const complex<double> dhP10bar =(1.4142135623730951*sqrt(ld*(-2. - 1.*ld + 2.*pow(ld,2) + pow(ld,3)))*(-1.*w2m + w2p)*((-0.752252778063675*I*r0*((-2.4*ellK*pow(M + 2.*ld*M,2))/((-1. + ld)*ld*(1. + ld)*(2. + ld)) + 16.*ellK*(-3.*M + r0)*(-2.*M + r0) - 8.*ellE*(-2.*M + r0)*(-5.*M + 2.*r0)))/(sqrt(0.5 + ld)*M*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) + ((-2.5066282746310002*I*M*pm*sqrt(((1. + 2.*ld)*r0)/(-3.*M + r0)))/(2.*M - 1.*r0) + (0.3761263890318375*I*((-2.4*(ellE + 2.*ellK)*pow(M + 2.*ld*M,2))/((-1. + ld)*ld*(1. + ld)*(2. + ld)) + 8.*(2.*ellE*(-4.*M + r0)*(-2.*M + r0) - 1.*ellK*(-3.*M + r0)*(-7.*M + 2.*r0))))/(sqrt(0.5 + ld)*M*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)))*r[j] + (2.5066282746310002*I*M*pm*sqrt(((1. + 2.*ld)*r0)/(-3.*M + r0))*(-1.*r0 + r[j]))/(2.*M - 1.*r0) - (0.3761263890318375*I*((-2.4*(ellE + 2.*ellK)*pow(M + 2.*ld*M,2))/((-1. + ld)*ld*(1. + ld)*(2. + ld)) + 8.*(2.*ellE*(-4.*M + r0)*(-2.*M + r0) - 1.*ellK*(-3.*M + r0)*(-7.*M + 2.*r0)))*(-1.*r0 + r[j]))/(sqrt(0.5 + ld)*M*sqrt(-3.*M + r0)*sqrt(-2.*M + r0))))/pow(r[j],2);
        
        const complex<double> ddhP1bar = (8.*(4.*ellK*(4.*M - 1.*r0)*sqrt(-2.*M + r0) + 2.*ellE*pow(-2.*M + r0,1.5) - 3.141592653589793*pm*sqrt(r0)*(-2.*M - 4.*ld*M + r0 + 2.*ld*r0))*w0)/(sqrt(3.141592653589793 + 6.283185307179586*ld)*pow(r0,3)*sqrt(-3.*M + r0));
        const complex<double> ddhP3bar = (8.*w0*(2.*M*(4.*ellK*(4.*M - 1.*r0)*sqrt(-2.*M + r0) + 2.*ellE*pow(-2.*M + r0,1.5) - 3.141592653589793*pm*sqrt(r0)*(-2.*M - 4.*ld*M + r0 + 2.*ld*r0))*(2.*M - 1.*r[j])*r[j] - 1.*(4.*ellK*(-4.*M + r0)*sqrt(-2.*M + r0) - 2.*ellE*pow(-2.*M + r0,1.5) + 3.141592653589793*pm*sqrt(r0)*(-2.*M - 4.*ld*M + r0 + 2.*ld*r0))*r[j]*pow(-2.*M + r[j],2) + 4.*pow(M,2)*(2.*ellE*pow(-2.*M + r0,1.5)*(-1.*r0 + r[j]) - 3.141592653589793*pm*sqrt(r0)*(-2.*(1. + 2.*ld)*M + r0 + 2.*ld*r0)*(-1.*r0 + r[j]) - 4.*ellK*sqrt(-2.*M + r0)*(r0*(-2.*r0 + r[j]) + 2.*M*(r0 - 2.*(-1.*r0 + r[j]))))))/(sqrt(3.141592653589793 + 6.283185307179586*ld)*pow(r0,3)*sqrt(-3.*M + r0)*pow(-2.*M + r[j],3));
        const complex<double> ddhP4bar = 0.;
        const complex<double> ddhP6bar = (8.*M*sqrt(r0*(-3.*M + r0))*(2.*ellE*(2.*M - 1.*r0)*sqrt(r0*(-3.*M + r0)) + 3.141592653589793*(1. + 2.*ld)*pm*r0*sqrt(-3.*M + r0)*sqrt(-2.*M + r0))*w0)/(sqrt(3.141592653589793 + 6.283185307179586*ld)*pow(-3.*M + r0,1.5)*pow(-2.*M + r0,1.5)*pow(r[j],3));
        const complex<double> ddhP7bar = (2.8284271247461903*sqrt(ld*(-2. - 1.*ld + 2.*pow(ld,2) + pow(ld,3)))*(w2m + w2p)*((0.752252778063675*r0*((2.4*ellK*pow(M + 2.*ld*M,2))/((-1. + ld)*ld*(1. + ld)*(2. + ld)) - 8.*ellE*(-2.*M + r0)*(-5.*M + 2.*r0) + 2.*ellK*(51.*pow(M,2) - 40.*M*r0 + 8.*pow(r0,2))))/(sqrt(0.5 + ld)*M*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - 1.*((2.5066282746310002*M*pm*sqrt(((1. + 2.*ld)*r0)/(-3.*M + r0)))/(2.*M - 1.*r0) + (0.3761263890318375*((2.4*(ellE + 2.*ellK)*pow(M + 2.*ld*M,2))/((-1. + ld)*ld*(1. + ld)*(2. + ld)) + 2.*(-2.*ellK*(39.*pow(M,2) - 26.*M*r0 + 4.*pow(r0,2)) + ellE*(67.*pow(M,2) - 48.*M*r0 + 8.*pow(r0,2)))))/(sqrt(0.5 + ld)*M*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)))*r[j] + (2.5066282746310002*M*pm*sqrt(((1. + 2.*ld)*r0)/(-3.*M + r0))*(-1.*r0 + r[j]))/(2.*M - 1.*r0) + (0.3761263890318375*((2.4*(ellE + 2.*ellK)*pow(M + 2.*ld*M,2))/((-1. + ld)*ld*(1. + ld)*(2. + ld)) + 2.*(-2.*ellK*(39.*pow(M,2) - 26.*M*r0 + 4.*pow(r0,2)) + ellE*(67.*pow(M,2) - 48.*M*r0 + 8.*pow(r0,2))))*(-1.*r0 + r[j]))/(sqrt(0.5 + ld)*M*sqrt(-3.*M + r0)*sqrt(-2.*M + r0))))/pow(r[j],3);
        const complex<double> ddhP8bar = 0.;
        const complex<double> ddhP10bar =(2.8284271247461903*sqrt(ld*(-2. - 1.*ld + 2.*pow(ld,2) + pow(ld,3)))*(-1.*w2m + w2p)*((0.752252778063675*I*r0*((-2.4*ellK*pow(M + 2.*ld*M,2))/((-1. + ld)*ld*(1. + ld)*(2. + ld)) + 16.*ellK*(-3.*M + r0)*(-2.*M + r0) - 8.*ellE*(-2.*M + r0)*(-5.*M + 2.*r0)))/(sqrt(0.5 + ld)*M*sqrt(6.*pow(M,2) - 5.*M*r0 + pow(r0,2))) - 1.*((-2.5066282746310002*I*M*pm*sqrt(((1. + 2.*ld)*r0)/(-3.*M + r0)))/(2.*M - 1.*r0) + (0.3761263890318375*I*((-2.4*(ellE + 2.*ellK)*pow(M + 2.*ld*M,2))/((-1. + ld)*ld*(1. + ld)*(2. + ld)) + 8.*(2.*ellE*(-4.*M + r0)*(-2.*M + r0) - 1.*ellK*(-3.*M + r0)*(-7.*M + 2.*r0))))/(sqrt(0.5 + ld)*M*sqrt(-3.*M + r0)*sqrt(-2.*M + r0)))*r[j] - (2.5066282746310002*I*M*pm*sqrt(((1. + 2.*ld)*r0)/(-3.*M + r0))*(-1.*r0 + r[j]))/(2.*M - 1.*r0) + (0.3761263890318375*I*((-2.4*(ellE + 2.*ellK)*pow(M + 2.*ld*M,2))/((-1. + ld)*ld*(1. + ld)*(2. + ld)) + 8.*(2.*ellE*(-4.*M + r0)*(-2.*M + r0) - 1.*ellK*(-3.*M + r0)*(-7.*M + 2.*r0)))*(-1.*r0 + r[j]))/(sqrt(0.5 + ld)*M*sqrt(-3.*M + r0)*sqrt(-2.*M + r0))))/pow(r[j],3);
        
        /* The non-trace-reversed field and its first and second derivatives.
         * Here we do include the factor of a and 1/r, and the window function.
         * Trace-reversal corresponds to swapping 3 and 6. */
        const double W = (r[j]<ra || r[j]>rb) ? 0.0 : 1.0;
        const double dW = (r[j]<ra || r[j]>rb) ? 0.0 : 0.0;
        const double ddW = (r[j]<ra || r[j]>rb) ? 0.0 : 0.0;
        hP[1][l][m][j] = a_il(1,l)*W*hP1bar/r[j];
        hP[3][l][m][j] = a_il(3,l)*W*hP6bar/r[j];
        hP[4][l][m][j] = a_il(4,l)*W*hP4bar/r[j];
        hP[6][l][m][j] = a_il(6,l)*W*hP3bar/r[j];
        hP[7][l][m][j] = a_il(7,l)*W*hP7bar/r[j];
        hP[8][l][m][j] = a_il(8,l)*W*hP8bar/r[j];
        hP[10][l][m][j] = a_il(10,l)*W*hP10bar/r[j];

        dhP[1][l][m][j] = a_il(1,l)*(W*(dhP1bar - hP1bar/r[j])/r[j] + dW*hP1bar/r[j]);
        dhP[3][l][m][j] = a_il(3,l)*(W*(dhP6bar - hP6bar/r[j])/r[j] + dW*hP6bar/r[j]);
        dhP[4][l][m][j] = a_il(4,l)*(W*(dhP4bar - hP4bar/r[j])/r[j] + dW*hP4bar/r[j]);
        dhP[6][l][m][j] = a_il(6,l)*(W*(dhP3bar - hP3bar/r[j])/r[j] + dW*hP3bar/r[j]);
        dhP[7][l][m][j] = a_il(7,l)*(W*(dhP7bar - hP7bar/r[j])/r[j] + dW*hP7bar/r[j]);
        dhP[8][l][m][j] = a_il(8,l)*(W*(dhP8bar - hP8bar/r[j])/r[j] + dW*hP8bar/r[j]);
        dhP[10][l][m][j] = a_il(10,l)*(W*(dhP10bar - hP10bar/r[j])/r[j] + dW*hP10bar/r[j]);

        ddhP[1][l][m][j] = a_il(1,l)*(W*(ddhP1bar - 2.0*(dhP1bar - hP1bar/r[j])/r[j])/r[j] + 2.0*dW*(dhP1bar - hP1bar/r[j])/r[j] + ddW*hP1bar/r[j]);
        ddhP[3][l][m][j] = a_il(3,l)*(W*(ddhP6bar - 2.0*(dhP6bar - hP6bar/r[j])/r[j])/r[j] + 2.0*dW*(dhP6bar - hP6bar/r[j])/r[j] + ddW*hP6bar/r[j]);
        ddhP[4][l][m][j] = a_il(4,l)*(W*(ddhP4bar - 2.0*(dhP4bar - hP4bar/r[j])/r[j])/r[j] + 2.0*dW*(dhP4bar - hP4bar/r[j])/r[j] + ddW*hP4bar/r[j]);
        ddhP[6][l][m][j] = a_il(6,l)*(W*(ddhP3bar - 2.0*(dhP3bar - hP3bar/r[j])/r[j])/r[j] + 2.0*dW*(dhP3bar - hP3bar/r[j])/r[j] + ddW*hP3bar/r[j]);
        ddhP[7][l][m][j] = a_il(7,l)*(W*(ddhP7bar - 2.0*(dhP7bar - hP7bar/r[j])/r[j])/r[j] + 2.0*dW*(dhP7bar - hP7bar/r[j])/r[j] + ddW*hP7bar/r[j]);
        ddhP[8][l][m][j] = a_il(8,l)*(W*(ddhP8bar - 2.0*(dhP8bar - hP8bar/r[j])/r[j])/r[j] + 2.0*dW*(dhP8bar - hP8bar/r[j])/r[j] + ddW*hP8bar/r[j]);
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
