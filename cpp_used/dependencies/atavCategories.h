#ifndef ATAVCATEGORIES_H
#define ATAVCATEGORIES_H

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <set>
#include <sstream>
#include <map>
#include <algorithm>

#include "atavClasses.h"

bool recsingleable( Variant &var );
bool rechomable( Variant &var );
bool reccompable( Variant &var );
bool domable( Variant &var );
bool dom_pass( Variant &var );
bool patho_var( Variant &var );
bool patho_acmg_var( Variant &var );

std::map<std::string, int> create_headers_map( std::map<std::vector<std::string>, std::string> &mapVsS, std::vector<std::string> &vString);

void getVariants(
  std::string &nomeFile,
  std::map<std::string, int> &mapAtavHeader,
  std::vector<std::string> &variant_abbreviations,
  std::vector<std::string> &patient_abbreviations,
  std::vector<std::string> &gene_abbreviations,
  std::vector<std::string>  &patient_var_abbreviations,
  std::vector<std::string> &vVarClassNames,
  std::vector<Variant> &vVarClass,
  std::vector<std::string> &vPatientClassNames,
  std::vector<Patient> &vPatientClass,
  std::vector<std::string> &vGeneClassNames,
  std::vector<Gene> &vGeneClass
);

#endif
