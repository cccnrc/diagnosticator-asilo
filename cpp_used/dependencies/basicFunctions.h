#ifndef BASICFUNCTIONS_H
#define BASICFUNCTIONS_H

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <map>
#include <set>
#include <sstream>
#include <algorithm>
#include <math.h>
#include <functional>
#include <numeric>
#include <sys/stat.h>
#include <sys/types.h>
#include <limits>

void remove_duplicates_from_V (std::vector<std::string> &v);
std::string stringize_float( float f );
void write_map_VV
(
    std::map<std::string, std::vector<std::vector<std::string>>> mVV,
    std::string nomeFile,
    std::string separator,
    std::string suffix
);
std::vector<std::string> convert_map_keys_to_V ( std::map<std::string, std::vector<std::string>> m );
std::string append_str_to_str ( std::string a, std::string b );
bool create_dir( std::string nomeDirResults );
std::string add_results_dir_path( std::string nameDir, std::string nameFile );
std::string stringize_vector (const std::vector<std::string> & vec, std::string delimiter);
bool floatable( std::string &f_string );
bool cmpf(float A, float B, float epsilon = 10E-10f);
bool cmpf_greater(float A, float B, float epsilon = 10E-10f);
void vectorize_column( std::string sCol, std::vector<std::string> &v, char separator );
void read_VEP_cin( std::string &nomeFile, std::vector<std::vector<std::string>> &vvVar );
std::string create_var_name( std::string &chr, std::string &pos, std::string &ref, std::string &alt );
void write_cin_file_vv( std::string &nomeFile, std::vector< std::vector<std::string>> &vvGnomadfounds );
std::vector<std::string> getHeader( std::string &nomeFile );
int add_absent_to_vector( std::vector<std::string> &v, std::string &value );
std::vector<std::string> getVectorFile( std::string &nomeFile );
std::vector<std::vector<std::vector<std::string>>> get_vvFile( std::string &nomeFile, char sep1, char sep2 );
std::vector<std::string> take_vvv_column( std::vector<std::vector<std::vector<std::string>>> &vvvHGVSPgenes, int numCol );
std::map<std::string, std::vector<std::string>> create_vvv_map( std::vector<std::vector<std::vector<std::string>>> &vvv, int col_key, int col_values );
void replaceAll(std::string& str, const std::string& from, const std::string& to);
bool replace_substring (std::string& str, const std::string& from, const std::string& to);
std::string stringize_map_str_str( std::map<std::string, std::string> m, std::string sep1, std::string sep2 );
std::string stringize_map_str_V( std::map<std::string, std::vector<std::string>> m, std::string sep1, std::string sep2, std::string sep3 );
std::vector<std::string> vectorize_column_v2( std::string sCol, char separator );
std::vector<std::string> vectorize_column_comma_quote( std::string sCol, char separator );
std::vector<std::vector<std::string>> getVVectorFile( std::string nomeFile, char sep );
bool inRange(unsigned low, unsigned high, unsigned x);
std::string stringize_map_map_str_V
(
    std::map< std::string, std::map<std::string, std::vector<std::string>>> mm,
    std::string sep1,
    std::string sep2,
    std::string sep3,
    std::string sep4,
    std::string sep5
);

#endif
