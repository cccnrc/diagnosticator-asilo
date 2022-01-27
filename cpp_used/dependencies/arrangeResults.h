#ifndef ARRANGERESULTS_H
#define ARRANGERESULTS_H

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <set>
#include <sstream>
#include <map>
#include <algorithm>

#include "basicFunctions.h"
#include "atavCategories.h"
#include "filterFunctions.h"

void arrange_and_write_results
(
    std::vector<Variant> vVarClass,
    std::vector<std::string> vVarClassNames,
    int selected_only_patho,
    std::string nomeFileResVar,
    std::vector<Patient> vPatientClass,
    std::vector<std::string> vPatientClassNames,
    std::string nomeFileResPz,
    std::string nomeFileResVarPz,
    std::vector<Gene> vGeneClass,
    std::vector<std::string> vGeneClassNames,
    std::string nomeFileResGene
);


#endif
