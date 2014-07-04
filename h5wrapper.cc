/*
 * Copyright 2014 Barry Wardell
 *
 * This file is distributed under the University of Illinois/NCSA
 * Open Source License. See LICENSE file for details.
 *
 */

#include "h5wrapper.h"

using namespace std;

/* Exception handling */
H5Exception::H5Exception(const string& message)
  : message(message) {}

const string& H5Exception::getMessage()
{
   return message;
}

/* Base wrapper class */
hid_t H5Base::getId() const
{
  return id;
}

/* H5A wrapper */
H5A::H5A(hid_t location_id, const string& attr_name, const H5T& type, const H5S& dataspace)
{
  id = H5Acreate(location_id, attr_name.c_str(), type.getId(), dataspace.getId(), H5P_DEFAULT, H5P_DEFAULT);
  if( id<0 )
    throw H5Exception("H5Acreate failed");
}

H5A::H5A(hid_t location_id, const string& attr_name)
{
  id = H5Aopen(location_id, attr_name.c_str(), H5P_DEFAULT);
  if( id<0 )
    throw H5Exception("H5Aopen failed");
}

H5A::~H5A()
{
  H5Aclose(id);
}

/* H5D wrapper */
H5D::H5D(const H5F& file, const string& dataset, const H5T& type, const H5S& dataspace, const H5P& prop)
{
  id = H5Dcreate(file.getId(), dataset.c_str(), type.getId(), dataspace.getId(), H5P_DEFAULT, prop.getId(), H5P_DEFAULT);
  if (id<0)
    throw H5Exception("H5Dcreate failed");
}

H5D::H5D(const H5F& file, const string& dataset)
{
    id = H5Dopen(file.getId(), dataset.c_str(), H5P_DEFAULT);
    if (id<0)
      throw H5Exception("H5Dopen failed");
}

H5D::~H5D()
{
  H5Dclose(id);
}

void H5D::resize(vector<hsize_t> newSize) const
{
  if( H5Dset_extent(id, newSize.data())<0 )
    throw H5Exception("H5Dset_extent failed");
}

int H5D::getNumAttrs() const
{
  H5O_info_t object_info;
  H5Oget_info(id, &object_info);
  return object_info.num_attrs;
}

void H5D::write(H5S& memspace, H5S& filespace, vector<double> data) const
{
  if( H5Dwrite(id, H5T_NATIVE_DOUBLE, memspace.getId(), filespace.getId(), H5P_DEFAULT, data.data())<0 )
    throw H5Exception("H5Dwrite failed");
}

void H5D::write(H5S& memspace, H5S& filespace, const double* data) const
{
  if( H5Dwrite(id, H5T_NATIVE_DOUBLE, memspace.getId(), filespace.getId(), H5P_DEFAULT, data)<0 )
    throw H5Exception("H5Dwrite failed");
}

/* H5F wrapper */
H5F::H5F(const string& filename, const unsigned int flags)
{
  if (flags & (H5F_ACC_CREAT|H5F_ACC_EXCL|H5F_ACC_TRUNC|H5F_ACC_DEBUG))
  {
    id = H5Fcreate(filename.c_str(), flags, H5P_DEFAULT, H5P_DEFAULT);
    if (id<0)
    {
      throw H5Exception("H5Fcreate failed");
    }
  } else {
    id = H5Fopen(filename.c_str(), flags, H5P_DEFAULT);
    if (id<0)
    {
      throw H5Exception("H5Fopen failed");
    }
  }
}

H5F::~H5F()
{
  H5Fclose(id);
}

template<class T, size_t N>
void H5F::write_dataset(const string& dsName, const boost::multi_array<T, N>& d)
{
  H5T type(H5T_NATIVE_DOUBLE);
  vector<hsize_t> dims(d.shape(), d.shape() + d.num_dimensions());
  H5S dataspace(dims);
  H5P prop(H5P_DATASET_CREATE);
  H5D dataset(*this, dsName, type, dataspace, prop);
  dataset.write(dataspace, dataspace, d.data());
}

template void H5F::write_dataset(const string& dsName, const boost::multi_array<double, 1>& d);
template void H5F::write_dataset(const string& dsName, const boost::multi_array<double, 2>& d);
template void H5F::write_dataset(const string& dsName, const boost::multi_array<double, 3>& d);

/* H5P wrapper for property lists */
H5P::H5P(hid_t cls_id)
{
  id = H5Pcreate(cls_id);
  if (id<0)
    throw H5Exception("H5Pcreate failed");
}

void H5P::set_chunk(vector<hsize_t> chunk_dims)
{
  if(H5Pset_chunk(id, chunk_dims.size(), chunk_dims.data())<0)
    throw H5Exception("H5Pset_chunk failed");
}

void H5P::set_deflate(const unsigned int level)
{
  if(H5Pset_deflate(id, level)<0)
    throw H5Exception("H5Pset_deflate failed");
}

void H5P::set_shuffle()
{
  if(H5Pset_shuffle(id)<0)
    throw H5Exception("H5Pset_shuffle failed");
}

H5P::~H5P()
{
  H5Pclose(id);
}

/* H5S dataspace wrapper  */
H5S::H5S(const H5S_class_t type)
{
  id = H5Screate(type);
  if( id<0 )
    throw H5Exception("H5Screate failed");
}

H5S::H5S(const vector<hsize_t> dims, const vector<hsize_t> maxdims)
{
  id = H5Screate_simple(dims.size(), dims.data(), maxdims.data());
  if( id<0 )
    throw H5Exception("H5Screate failed");
}

H5S::H5S(const vector<hsize_t> dims)
{
  id = H5Screate_simple(dims.size(), dims.data(), NULL);
  if( id<0 )
    throw H5Exception("H5Screate failed");
}

H5S::H5S(const H5A& attr)
{
  id = H5Aget_space(attr.getId());
  if( id<0 )
    throw H5Exception("H5Aget_space failed");
}

H5S::H5S(const H5D& ds)
{
  id = H5Dget_space(ds.getId());
  if( id<0 )
    throw H5Exception("H5Dget_space failed");
}

H5S::~H5S()
{
  H5Sclose(id);
}

vector<hsize_t> H5S::getSimpleExtentDims() const
{
  int ndims = getSimpleExtentNDims();

  vector<hsize_t> dims(ndims);
  if( H5Sget_simple_extent_dims(id, dims.data(), NULL) < 0 )
    throw H5Exception("H5Sget_simple_extent_dims failed");

  return dims;
}

int H5S::getSimpleExtentNDims() const
{
  int ndims = H5Sget_simple_extent_ndims(id);
  if( ndims < 0 )
    throw H5Exception("H5Sget_simple_extent_ndims failed");
  return ndims;
}

void H5S::select_hyperslab(H5S_seloper_t op, const hsize_t *start, const hsize_t *stride, const hsize_t *count, const hsize_t *block)
{
  H5Sselect_hyperslab(id, op, start, stride, count, block);
}

/* H5T wrapper */
H5T::H5T(const H5T_class_t cls, const size_t size)
 : isBuiltin(false)
{
  id = H5Tcreate(cls, size);
  if( id<0 )
    throw H5Exception("H5Tcreate failed");
}

H5T::H5T(const hid_t type)
 : isBuiltin(true)
{
  id = type;
}

H5T::H5T(const H5D& ds)
 : isBuiltin(false)
{
  id = H5Dget_type(ds.getId());
  if( id<0 )
    throw H5Exception("H5Dget_type failed");
}

H5T::H5T(const H5A& attr)
 : isBuiltin(false)
{
  id = H5Aget_type(attr.getId());
  if( id<0 )
    throw H5Exception("H5Aget_type failed");
}

H5T::~H5T()
{
  if( !isBuiltin )
    H5Tclose(id);
}

size_t H5T::getSize() const
{
  return H5Tget_size(id);
}

/* H5PT - class for storing time series data in a dataset */
H5TS::H5TS(const H5F& file, const string& name, const hsize_t size, const hsize_t chunk_size)
{
  /* Data type */
  H5T type(H5T_NATIVE_DOUBLE);

  /* Dataspace - we have 'size' columns and an extensible number (initially 0) of rows */
  vector<hsize_t> dims(2);
  dims[0] = 0;
  dims[1] = size;

  vector<hsize_t> maxdims(2);
  maxdims[0] = H5S_UNLIMITED;
  maxdims[1] = size;
  H5S dataspace(dims, maxdims);

  /* Creation property list */
  H5P prop(H5P_DATASET_CREATE);
  vector<hsize_t> chunk_dims(2);
  chunk_dims[0] = chunk_size;
  chunk_dims[1] = size;
  prop.set_chunk(chunk_dims);

  dataset = new H5D(file, name, type, dataspace, prop);
}

H5TS::H5TS(const H5F& file, const string& name)
{
  dataset = new H5D(file, name);
}

H5TS::~H5TS()
{
  delete dataset;
}

void H5TS::AppendData(vector<double> data) const
{
  H5S oldDataspace(*dataset);
  vector<hsize_t> oldSize = oldDataspace.getSimpleExtentDims();

  vector<hsize_t> extSize(2);
  extSize[0] = 1;
  extSize[1] = oldSize[1];

  vector<hsize_t> newSize(2);
  newSize[0] = oldSize[0] + extSize[0];
  newSize[1] = oldSize[1];

  dataset->resize(newSize);

  /* Select a hyperslab in the extended portion of the dataset  */
  H5S newDataspace(*dataset);
  vector<hsize_t> offset(2);
  offset[0] = oldSize[0];
  offset[1] = 0;
  newDataspace.select_hyperslab(H5S_SELECT_SET, offset.data(), NULL, extSize.data(), NULL);

  /* Define memory dataspace */
  H5S memspace(extSize);

  /* Write the data to the extended portion of dataset  */
  dataset->write(memspace, newDataspace, data);
}
