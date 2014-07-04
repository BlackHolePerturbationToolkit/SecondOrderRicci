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
#include <regex>
#include <dirent.h>
#include <boost/multi_array.hpp>
#include "h5wrapper.h"

using namespace std;
using boost::multi_array;

/* Return a list of all files in a directory */
vector<string> list_files(string dirname)
{
  DIR *dir;
  struct dirent *ent;
  vector<string> files;

  dir = opendir(dirname.c_str());
  if (dir != NULL) {
    while ((ent = readdir(dir)) != NULL) {
      if(ent->d_type == DT_REG) {
        files.push_back(dirname + '/' + ent->d_name);
      }
    }
    closedir(dir);
  } else {
    cerr << "Cannot open directory: " << dirname << endl;
  }

  sort(files.begin(), files.end());

  return files;
}

struct lm_mode
{
  int l;
  int m;
};

struct lm_mode filenameToMode(string filename)
{
  smatch match;
  regex  pattern("h1-l([\\d]+)m([\\d]+).h5");

  regex_search(filename, match, pattern);
  assert(match.size() == 3);
  int l = stoi(match[1]);
  int m = stoi(match[2]);
  struct lm_mode lm = {l, m};
  return lm;
}

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

complex<double> d2hdr2(int i, int l, int m, double r, double f, double fp, vector<complex<double>> hbar, vector<complex<double>> dhbar, double M, double r0)
{
  double Omega = sqrt(M)*pow(r0, -1.5);
  complex<double> dt = - complex<double>(0.0, m)*Omega;
  complex<double> ddhbar;
  switch(i) {
    case 1:
      ddhbar = -(fp/f)*dhbar[1-1] + hbar[1-1]*(-Omega*Omega*m*m)/(f*f)
               + hbar[1-1]*((2.0*M)/r0 + l*(l+1))/(f*r0*r0)
               + 4.0/(f*f)*(0.5*f*fp*dhbar[1-1] - 0.5*fp*dt*hbar[2-1]
               + f*f/(2.0*r0*r0)*(hbar[1-1] - f*hbar[3-1] - hbar[5-1] - f*hbar[6-1]));
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

int main(int argc, char* argv[])
{
  const double M = 1.0;

  /* Read input HDF5 files containing first order fields */
  vector<string> files = list_files("../h1_fields");

  /* FIXME: don't hardcode these */
  const int l_max = 23;   /* Number of l modes */
  const int N     = 210;  /* Number of grid points */
  const int r0i   = 23;   /* Index of worldline point */

  /* Read the grid coordinates - assume the same for all fields */
  vector<double> r(N), f(N), fp(N);
  double Omega;
  {
    H5F h5_file(files[0], H5F_ACC_RDONLY);
    H5D dataset(h5_file, "grid");
    H5S dataspace(dataset);
    vector<hsize_t> gridSize = {N};
    H5S memspace(gridSize);

    H5Dread(dataset.getId(), H5T_NATIVE_DOUBLE, memspace.getId(), dataspace.getId(), H5P_DEFAULT, r.data());
    Omega = sqrt(M) * pow(r[r0i], -1.5);
    for(int i=0; i<N; ++i) {
      f[i] = 1.0 - 2.0*M/r[i];
      fp[i] = 2.0*M/(r[i]*r[i]);
    }
  }
  const double r0 = r[r0i];

  /* Read the first order data */
  /* FIXME: This is wasteful of memory by up to a factor of 4, but it probably doesn't matter */
  multi_array<complex<double>, 4> h(boost::extents[10][l_max+1][2*l_max+1][N]);
  multi_array<complex<double>, 4> dh(boost::extents[10][l_max+1][2*l_max+1][N]);
  multi_array<complex<double>, 4> ddh(boost::extents[10][l_max+1][2*l_max+1][N]);
  fill(h.origin(), h.origin() + h.size(), 0.0);
  fill(dh.origin(), dh.origin() + dh.size(), 0.0);
  fill(ddh.origin(), ddh.origin() + ddh.size(), 0.0);

  for(auto file: files)
  {
    struct lm_mode lm = filenameToMode(file);
    int l = lm.l;
    int m = lm.m;

    H5F h5_file(file, H5F_ACC_RDONLY);
    H5D dataset(h5_file, "inhom");
    H5S dataspace(dataset);

    /* Check the dataset is a 2D array */
    const int rank = dataspace.getSimpleExtentNDims();
    assert(rank == 2);

    vector<hsize_t> size = dataspace.getSimpleExtentDims();
    assert(size[0] == N);

    H5S memspace(size);
    multi_array<double, 2> data(boost::extents[size[0]][size[1]]);
    multi_array<complex<double>, 2> hbar(boost::extents[size[1]][N]);
    multi_array<complex<double>, 2> dhbar(boost::extents[size[1]][N]);
    multi_array<complex<double>, 2> ddhbar(boost::extents[size[1]][N]);
    H5Dread(dataset.getId(), H5T_NATIVE_DOUBLE, memspace.getId(), dataspace.getId(), H5P_DEFAULT, data.data());
    vector<int> fields;
    if( (l+m) % 2 == 0) {
      if (l==1) {
        assert(size[1] == 4*6);
        fields = {1, 3, 5, 6, 2, 4};
      } else {
        assert(size[1] == 4*7);
        fields = {1, 3, 5, 6, 7, 2, 4};
      }
    } else {
      if (l==1) {
        assert(size[1] == 4*2);
        fields = {9, 8};
      } else {
        assert(size[1] == 4*3);
        fields = {9, 10, 8};
      }
    }
    for(vector<int>::size_type i=0; i!=fields.size(); ++i) {
      double a = a_il(fields[i], l);
      for(int j=0; j<size[0]; ++j) {
        hbar[i][j] = complex<double>(data[j][2*i], data[j][2*i+1]);
        dhbar[i][j] = complex<double>(data[j][2*i+2], data[j][2*i+3]);
        h[fields[i]-1][l][m][j] = a*hbar[i][j]/r[j];
        dh[fields[i]-1][l][m][j] = a*(dhbar[i][j] - hbar[i][j]/r[j])/r[j];
      }
    }

    for(vector<int>::size_type i=0; i!=fields.size(); ++i) {
      double a = a_il(fields[i], l);
      for(int j=0; j<size[0]; ++j) {
        vector<complex<double>> hbarj(10), dhbarj(10);
        for(int k=0; k<10; ++k) {
          hbarj[k] = h[k][l][m][j];
          dhbarj[k] = dh[k][l][m][j];
        }
        ddhbar[i][j] = d2hdr2(i, l, m, r[j], f[j], fp[j], hbarj, dhbarj, M, r0);
        ddh[fields[i]-1][l][m][j] = a*(ddhbar[i][j] - 2.0*(dhbar[i][j] - hbar[i][j]/r[j])/r[j])/r[j];
      }
    }
  }

  return 0;
}
