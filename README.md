# Second Order Ricci

Compute the source for the second order Lorenz gauge equation


### Requirements

The code requires the following libraries:

* GNU Scientific Library: This provides functions for computing Wigner 3j symbols.
* HDF5 1.8: Data is read in and output in HDF5 format.
* C++11 conformant C++ compiler with support for OpenMP.
* Boost: Multi-dimensional arrays are provided by Boost's multi_array type.
* SCons: The build system requires SCons

It has been tested with the following versions:
GSL 1.16
HDF5 1.8.13
GCC 4.9.0
Boost 1.55.0
SCons 2.3.1


### Building

To compile the software, run 'scons' from inside the source directory. If you have
libraries in non-standard paths, edit the src/SConscript file and change the settings
in the "Build options" section.


### Installation

No installation is required.


### Usage

Make sure data for the first order field is available in a directory by itself
(i.e. no other files present), then run:

./Ricci <directory>

where <directory> is the path to the directory of first order data.


### License

This code is distributed under the University of Illinois/NCSA
Open Source License. Details can be found in the LICENSE file.

### Authors

Barry Wardell
