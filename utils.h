/*
 * Copyright 2014 Barry Wardell
 *
 * This file is distributed under the University of Illinois/NCSA
 * Open Source License. See LICENSE file for details.
 *
 */

#include <string>
#include <vector>
#include <regex>

struct ilm_mode
{
  int i;
  int l;
  int m;
};

struct lm_mode
{
  int l;
  int m;
};

lm_mode filenameToMode(std::string filename);
std::vector<std::string> list_files(std::string dirname, std::regex pattern);

bool isOdd(int n);
bool isEven(int n);
