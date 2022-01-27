#include "readGenelist.h"

/*
  questa serve a leggere il file delle merged_gene_lists generato da ./mergeGenelists.0
  e trasformarlo in una classe Gene e relativo vettore nomi
*/



/*--------------- READ HEADER -------------------*/
void extractGenelistFile(
  std::vector<Gene> &vGeneClass,
  std::vector<std::string> &vGeneClassNames
) {

    std::string line;
    int lineCount = 0;
    // apro loop sulle righe
    while(std::getline(std::cin, line)) {

      Gene tmp;

      // vettorializzo la riga
      char fSep = '\t'; std::vector<std::string> vRiga;
      vectorize_column(line, vRiga, fSep);

      // vettorializzo gli inh_mode
      char iSep = ','; std::vector<std::string> vInhMode;
      vectorize_column(vRiga[1], vInhMode, iSep);
      // vettorializzo le genelists
      std::vector<std::string> vGenelist;
      vectorize_column(vRiga[2], vGenelist, iSep);

      vGeneClassNames.push_back( vRiga[0] );

      tmp.addAttr("gene_name_correct", vRiga[0]);
      remove_duplicates_from_V(vInhMode);
      remove_duplicates_from_V(vGenelist);
      tmp.addInhModes( vInhMode );
      tmp.addGenelists( vGenelist );

      vGeneClass.push_back(tmp);

      lineCount++;

    //chiudo LOOP sulle linee
    }
    std::cout << "--> read n: " << lineCount << " genes from lists" << '\n';
}

//questa cambia STD::CIN con il risultato della lettura del file
//poi chiama funzione su quella
std::vector<std::string> getGenelistFile( std::string &nomeFile, std::vector<Gene> &vGeneClass ) {

    std::vector<std::string> vGeneClassNames;

    std::ifstream in(nomeFile.c_str());
    std::streambuf *cinbuf = std::cin.rdbuf(); //save old buf
    std::cin.rdbuf(in.rdbuf()); //redirect std::cin to in.txt!

    extractGenelistFile( vGeneClass, vGeneClassNames ); //call function

    std::cin.rdbuf(cinbuf);   //reset to standard input again

    return vGeneClassNames;
}
/*-----------------------------------------------------*/


std::vector<std::string> get_gene_lists(
  std::string &nomeFile,
  std::vector<Gene> &vGeneClass
){
  std::vector<std::string> vGeneClassNames = getGenelistFile( nomeFile, vGeneClass );
  return vGeneClassNames;
}
