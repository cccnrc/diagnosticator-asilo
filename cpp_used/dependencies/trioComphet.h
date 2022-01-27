#ifndef TRIOCOMPHET_H
#define TRIOCOMPHET_H

#include "basicFunctions.h"
#include "atavCategories.h"
#include "atavClasses.h"

void get_comphet_trio_file
(
    std::map<std::vector<std::string>, std::string> &atav_categories_abbreviation_map,
    std::vector< std::string > &variant_abbreviations,
    std::vector< std::string > &patient_abbreviations,
    std::vector< std::string > &patient_var_abbreviations,
    std::string nomeFileTrioComphet,
    std::vector<Patient> &vPatientClass,
    std::vector<std::string> &vPatientClassNames,
    std::vector<Gene> &vGeneClass,
    std::vector<std::string> &vGeneClassNames,
    std::vector<Variant> &vVarClass,
    std::vector<std::string> &vVarClassNames,
    std::map < std::string, std::map < std::string, std::vector<std::string >>> &mPatientsComphetVars
);

#endif
