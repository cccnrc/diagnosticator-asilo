#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <set>
#include <sstream>
#include <algorithm>
#include <math.h>
#include <functional>
#include <numeric>

/*
  questo molto importante PRE-asilo.1.x poichè prende le varie N liste geniche passate (gene_name \t inh_modes)
  e le merga in un file unico merged_gene_lists.txt organizzato come geni unici e header:
  nome_gene - inh_modes - genelists
*/

/*
time gcc -O3 -lstdc++ --std=c++11 -g \
    mergeGenelists.0.cpp\
    -o mergeGenelists.0 \
    && ./mergeGenelists.0 \
    /home/enrico/columbia/asilo2/sito/app/users/enrico5/genelist/acmg.59.txt \
    /home/enrico/columbia/asilo2/sito/app/users/enrico5/genelist/new_kidney_genelist.due_col.txt \
    sticazzi.txt

*/

/* DEBUG:
G_SLICE=always-malloc \
  G_DEBUG=gc-friendly \
  valgrind -v --tool=memcheck --leak-check=full \
  --num-callers=40 \
  --log-file=valgrind.log \
  $(which ./mergeGenelists.0) \
    /home/enrico/columbia/asilo2/sito/app/users/enrico5/genelist/acmg.59.txt \
    /home/enrico/columbia/asilo2/sito/app/users/enrico5/genelist/new_kidney_genelist.due_col.txt \
    sticazzi.txt

*/


// add a std::string to a std::vector<std::string> only if absent e mi ritorna 1/0 se fatto/non
int add_if_absent( std::vector<std::string> &v, std::string &s ) {
  int result = 0;
  if (std::find (v.begin(), v.end(), s) == v.end()) { v.push_back(s); result = 1; }
  return result;
};


class Gene {
  std::string geneName;
  std::vector<std::string> vInhMode;
  std::vector<std::string> vGeneLists;
  public :
    void defineName( std::string gene_name ) { geneName = gene_name; }
    std::string getGeneName() { return geneName; }
    void addInhMode( std::string inh_mode ) { int added = add_if_absent( vInhMode, inh_mode ); }
    std::vector<std::string> getInhModes() { return vInhMode; }
    void addGeneList( std::string gene_list ) { int added = add_if_absent( vGeneLists, gene_list ); }
    std::vector<std::string> getGeneLists() { return vGeneLists; }
};

//vettorializzo la colonna dalla stringa
void vectorize_column( std::string &sCol, std::vector<std::string> &v, char &separator ) {
  std::istringstream ss(sCol); // string stream delle posizioni
  while (getline ( ss, sCol, separator) ) v.push_back(sCol); //riempo array con la pos convertita a int
}


//definisco funzione per CREARE VAR NAME da sottostringhe
std::string create_var_name( std::string &chr, std::string &pos, std::string &ref, std::string &alt ) {
  std::stringstream streamVarName;
  streamVarName << chr << "-" << pos << "-" << ref << "-" << alt;
  std::string varName = streamVarName.str();
  return varName;
}

//funzione per WRITE FILE di vec of vec
void write_cin_file_vv( std::string &nomeFile, std::vector< std::vector<std::string>> &vvGnomadfounds ) {

    std::ofstream out(nomeFile.c_str());
    std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
    std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!

    for (auto i : vvGnomadfounds) {
      for (auto vi : i) std::cout << vi << "\t";
      std::cout << '\n';
    }

    std::cout.rdbuf(coutbuf); //reset to standard output again
}


/*--------------- READ HEADER -------------------*/
void extract_genes(
  std::vector<Gene> &vGeneClass,
  std::vector<std::string> &vGeneClassNames,
  std::string &nomeList
) {

    std::string line;
    int lineCount = 0;
    // apro loop sulle righe
    while(std::getline(std::cin, line)) {

      if (line[0] == '#') continue;

      Gene tmp;

      char fSep = '\t'; std::vector<std::string> vRiga;
      vectorize_column(line, vRiga, fSep);

      std::string gene_name = vRiga[0];

      std::string inh_mode;
      if ( vRiga.size() > 1 ) {
          inh_mode = vRiga[1];
      }
      else {
          inh_mode = "U";
      }

      char iSep = ','; std::vector<std::string> vInhMode;
      vectorize_column(inh_mode, vInhMode, iSep);

      int geneOk = 1;

      if ( vInhMode.size() > 0 ) {
          for ( std::string inh : vInhMode )
          {
              if ( inh != "AD" && inh != "AR" && inh != "U" && inh != "XL" && inh != "XLD" && inh != "XLR" && inh != "MT"  ) {
                  std::cout << " !!! ERROR !!! genelist: " << nomeList << " invalid inh mode: " << inh << " in gene " << gene_name << ": gene REMOVED! !!!" << '\n';
                  geneOk = 0;
              }
          }
      }

      if ( geneOk == 0 ) {
        continue;
      }

      int added = add_if_absent( vGeneClassNames, gene_name );

      // se immesso (quindi nuovo)
      if ( added == 1 ) {
        tmp.defineName(gene_name);
        for ( std::string single_inh_mode : vInhMode )
        {
            tmp.addInhMode(single_inh_mode);
        }
        tmp.addGeneList(nomeList);
        vGeneClass.push_back(tmp);
      // se era gia presente riprendo il numero e ci metto i valori
      } else {
        std::ptrdiff_t posGene = find( vGeneClassNames.begin(), vGeneClassNames.end(), gene_name) - vGeneClassNames.begin();
        vGeneClass[ posGene ].addInhMode(inh_mode);
        vGeneClass[ posGene ].addGeneList(nomeList);
      }

      lineCount++;

    //chiudo LOOP sulle linee
    }

    std::cout << "---> genelist: " << nomeList << ", n genes: " << lineCount << '\n';
}

//questa cambia STD::CIN con il risultato della lettura del file
//poi chiama funzione su quella
// e mi restituisce un map di map di vettori std::string organizzata come: gene_name<vInhMode, vGeneLists>
void get_genes(
  std::string &nomeFile,
  std::string &nomeList,
  std::vector<Gene> &vGeneClass,
  std::vector<std::string> &vGeneClassNames
) {

    std::ifstream in(nomeFile.c_str());
    std::streambuf *cinbuf = std::cin.rdbuf(); //save old buf
    std::cin.rdbuf(in.rdbuf()); //redirect std::cin to in.txt!

    extract_genes( vGeneClass, vGeneClassNames, nomeList ); //call function

    std::cin.rdbuf(cinbuf);   //reset to standard input again
}
/*-----------------------------------------------------*/


/*
 * Erase all Occurrences of given substring from main string.
 */
void eraseAllSubStr(std::string & mainStr, const std::string & toErase)
{
	size_t pos = std::string::npos;

	// Search for the substring in string in a loop untill nothing is found
	while ((pos  = mainStr.find(toErase) )!= std::string::npos)
	{
		// If found then erase it from string
		mainStr.erase(pos, toErase.length());
	}
}

/* Erase all Occurrences of all given substrings from main string using C++11 stuff
*/
std::string eraseSubStrings(std::string & mainStr, const std::vector<std::string> & strList) {
 // Iterate over the given list of substrings. For each substring call eraseAllSubStr() to
 // remove its all occurrences from main string.
 std::string newName = mainStr;
 std::vector<std::string> vS; char vSep = '/'; vectorize_column(newName, vS, vSep);
 newName = vS.back();
 std::for_each(strList.begin(), strList.end(), std::bind(eraseAllSubStr, std::ref(newName), std::placeholders::_1));
 return newName;
}

// trasforma vettore in stringa con delimiter specificato
std::string stringize_vector (const std::vector<std::string> & vec, std::string delimiter) {
    return std::accumulate(std::next(vec.begin()), vec.end(),
        vec[0],
        [&delimiter](std::string& a, std::string b) {
        return a + delimiter + b;
    });
}


int main(int argc, char const *argv[]) {

  std::vector<std::string> vGeneListsPassed;

  //std::string nomeFileOutput = "merged_gene_lists.txt";
  std::string nomeFileOutput = argv[argc-1];

  if ( argc > 2 ) {
    std::cout << "num arg correct" << '\n';
    std::string counter = "--";
    for (int i = 1; i < argc-1; ++i)
    {
        counter = counter + "--";
        std::cout << counter << "> list " << i << ": " << argv[i] << '\n';
        vGeneListsPassed.push_back(argv[i]);
    }
    std::cout << "(I/O) --> writing results in: " << nomeFileOutput << '\n';
  }

  // inizializzo vettore classe e nomi della stessa
  std::vector<Gene> vGeneClass;
  std::vector<std::string> vGeneClassNames;

  // vettore per tenere le substring da eliminare (estensioni dei file)
  std::vector<std::string> vExtensionToRemove = {".txt", ".list"};

  // assegno alla classe Gene i vari geni con relativi vettori di inh_mode e genelists di appartenenza
  for ( std::string gene_list : vGeneListsPassed ) {
    std::string nameList = eraseSubStrings(gene_list, vExtensionToRemove);
    get_genes( gene_list, nameList, vGeneClass, vGeneClassNames );
  }

  // metto in VV per poi scrivere su file
  std::vector<std::vector<std::string>> vvGenes;
  for ( Gene gene : vGeneClass ) {
    std::vector<std::string> vGene;
    vGene.push_back(gene.getGeneName());
    vGene.push_back(stringize_vector(gene.getInhModes(), ","));
    vGene.push_back(stringize_vector(gene.getGeneLists(), ","));
    vvGenes.push_back(vGene);
  }

  write_cin_file_vv( nomeFileOutput, vvGenes );

  return 0;
}
