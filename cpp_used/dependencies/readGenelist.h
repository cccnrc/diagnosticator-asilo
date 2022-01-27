#ifndef READGENELIST_H
#define READGENELIST_H

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <set>
#include <sstream>
#include <algorithm>

#include "basicFunctions.h"
#include "atavClasses.h"

std::vector<std::string> get_gene_lists(
  std::string &nomeFile,
  std::vector<Gene> &vGeneClass
);

#endif
