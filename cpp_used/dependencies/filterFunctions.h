#ifndef FILTERFUNCTIONS_H
#define FILTERFUNCTIONS_H

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <set>
#include <map>
#include <sstream>
#include <algorithm>
#include <iomanip> // setprecision

#include "atavClasses.h"
#include "basicFunctions.h"
#include "atavCategories.h"

/* ----------
  funzione iterativa per fargli passare tutte le varianti
  ----------------- modello DOMINANTE -------------------
--------- */
void filter_variants (
  std::vector<Variant> &vVarClass,
  std::vector<std::string> &vPopulationSelected,
  std::vector<Gene> &vGeneClass,
  std::vector<std::string> &vGeneClassNames,
  float rare_threshold,
  int casesMaxAC
);

void filter_variants_recessive (
  std::vector<Patient> &vPatientClass,
  std::vector<Gene> &vGeneClass,
  std::vector<std::string> &vGeneClassNames,
  std::vector<Variant> &vVarClass,
  std::vector<std::string> &vVarClassNames,
  int comphet_analysis
);

void filter_variants_acmg (
  std::vector<Variant> &vVarClass,
  std::vector<Gene> &vGeneClass,
  std::vector<std::string> &vGeneClassNames,
  int &parents_aff
);

void assign_priority(
  std::vector<Variant> &vVarClass
);

void check_and_add_patients( Variant &var );

void acmg_priority_screening( std::vector<Variant> &vVarClass, std::vector<std::string> &vVarClassNames );

void acmg_pm3_screening( std::vector<Variant> &vVarClass, std::vector<std::string> &vVarClassNames );

#endif
