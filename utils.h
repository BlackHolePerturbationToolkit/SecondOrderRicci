/*
 * Copyright 2014 Barry Wardell
 *
 * This file is distributed under the University of Illinois/NCSA
 * Open Source License. See LICENSE file for details.
 *
 */

struct lm_mode
{
  int l;
  int m;
};

lm_mode filenameToMode(std::string filename);
std::vector<std::string> list_files(std::string dirname);
