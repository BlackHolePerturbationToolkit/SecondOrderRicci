/*
 * Copyright 2014 Barry Wardell (original dense-array version)
 *
 * This file is distributed under the University of Illinois/NCSA
 * Open Source License. See LICENSE file for details.
 */
#pragma once

#include <string>
#include <vector>
#include "field.h"

/* Read the first-order field modes h1-l*m*.h5 in `dir` (modes with l <= l_max) into packed
 * storage: h, dh, ddh are the field values and their first two r-derivatives as stored in
 * the files (components 3 and 6 swapped), as expected by the generated kernels.
 * Returns the largest l read. */
int read_h1(const std::string dir, double &r0, int &r0i,
            std::vector<double> &r, std::vector<double> &f, std::vector<double> &fp,
            Field &h, Field &dh, Field &ddh, int l_max);
