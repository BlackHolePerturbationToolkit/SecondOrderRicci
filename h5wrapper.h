/*
 * Copyright 2014 Barry Wardell
 *
 * This file is distributed under the University of Illinois/NCSA
 * Open Source License. See LICENSE file for details.
 *
 */

#ifndef H5WRAPPER_H
#define H5WRAPPER_H

#include <hdf5.h>
#include <string>
#include <vector>
#include <boost/multi_array.hpp>

class H5A;
class H5D;
class H5F;
class H5S;
class H5P;
class H5T;

class H5Exception
{
public:
  H5Exception(const std::string& message);
  const std::string& getMessage();

private:
  std::string message;
};

class H5Base
{
public:
  hid_t getId() const;
  virtual ~H5Base() { };

protected:
  hid_t id;
};

class H5A : public H5Base
{
public:
  H5A(hid_t location_id, const std::string& attr_name, const H5T& type, const H5S& dataspace);
  H5A(hid_t location_id, const std::string& attr_name);
  ~H5A();
};

class H5D : public H5Base
{
public:
  H5D(const H5F& file, const std::string& dataset, const H5T& type, const H5S& dataspace, const H5P& prop);
  H5D(const H5F& file, const std::string& dataset);
  ~H5D();

  void resize(std::vector<hsize_t> newSize) const;
  int getNumAttrs() const;
  void write(H5S& memspace, H5S& filespace, std::vector<double> data) const;
  void write(H5S& memspace, H5S& filespace, const double* data) const;
};

class H5F : public H5Base
{
public:
  H5F(const std::string& filename, const unsigned int flags);
  ~H5F();
  template<class T, size_t N>
  void write_dataset(const std::string& dsName, const boost::multi_array<T, N>& d);
};

class H5P : public H5Base
{
public:
  H5P(hid_t cls_id);
  void set_chunk(std::vector<hsize_t> chunk_dims);
  void set_deflate(const unsigned int level);
  void set_shuffle();
  ~H5P();
};

class H5S : public H5Base
{
public:
  H5S(const H5S_class_t type);
  H5S(const std::vector<hsize_t> dims, const std::vector<hsize_t> maxdims);
  H5S(const std::vector<hsize_t> dims);
  H5S(const H5A& attribute);
  H5S(const H5D& dataset);
  ~H5S();

  std::vector<hsize_t> getSimpleExtentDims() const;
  int getSimpleExtentNDims() const;
  void select_hyperslab(H5S_seloper_t op, const hsize_t *start, const hsize_t *stride, const hsize_t *count, const hsize_t *block);
};

class H5T : public H5Base
{
public:
  H5T(const H5T_class_t cls, const size_t size);
  H5T(const hid_t type);
  H5T(const H5A& dataset);
  H5T(const H5D& dataset);
  ~H5T();

  size_t getSize() const;

private:
  const bool isBuiltin;
};

class H5TS
{
public:
  H5TS(const H5F& file, const std::string& name, const hsize_t size, const hsize_t chunk_size = 200);
  H5TS(const H5F& file, const std::string& name);
  ~H5TS();

  void AppendData(std::vector<double> data) const;

private:
  H5D *dataset;
};

#endif // H5WRAPPER_H
