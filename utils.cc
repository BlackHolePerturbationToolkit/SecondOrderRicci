/*
 * Copyright 2014 Barry Wardell
 *
 * This file is distributed under the University of Illinois/NCSA
 * Open Source License. See LICENSE file for details.
 *
 */

#include <assert.h>
#include <iostream>
#include <string>
#include <vector>
#include <regex>
#include <dirent.h>
#include "utils.h"
#include <stdio.h>
#include <string.h>

using namespace std;

/* Return a list of all files in a directory matching a pattern */
vector<string> list_files(string dirname, regex pattern)
{
  DIR *dir;
  struct dirent *ent;
  vector<string> files;

  dir = opendir(dirname.c_str());
  if (dir != NULL) {
    while ((ent = readdir(dir)) != NULL) {
      if(!regex_match(ent->d_name, pattern))
        continue;
      files.push_back(dirname + '/' + ent->d_name);
    }
    closedir(dir);
  } else {
    cerr << "Cannot open directory: " << dirname << endl;
  }

  sort(files.begin(), files.end());

  return files;
}

/* Convert a filename to an lm_mode object */
lm_mode filenameToMode(string filename)
{
  smatch match;
  regex  pattern("h1-l([0-9]+)m([0-9]+).h5");

  regex_search(filename, match, pattern);
  assert(match.size() == 3);
  int l = stoi(match[1]);
  int m = stoi(match[2]);
  lm_mode lm = {l, m};
  return lm;
}

bool isOdd(int n)
{
  return n % 2;
}

bool isEven(int n)
{
  return !isOdd(n);
}
