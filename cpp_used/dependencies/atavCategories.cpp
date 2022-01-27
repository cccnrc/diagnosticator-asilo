#include "atavCategories.h"
#include "basicFunctions.h"

// funzione per capire se VAR in alemno qualche modello DOM
bool patho_var( Variant &var )
{
    // check for DOM model
    if ( domable(var) ) {
        if ( var.getValue("dom") != "NA" ) {
            return true;
        }
    }
    // check for REC_HOM model
    if ( rechomable(var) ) {
        if ( var.getValue("rec_hom") != "NA" ) {
            return true;
        }
    }
    // check for REC COMPHET model
    if ( reccompable(var) ) {
        if ( var.getValue("rec") != "NA" ) {
            return true;
        }
    }
    return false;
}

// funzione per capire se VAR in alemno qualche modello oppure se ACMG P/LP
bool patho_acmg_var( Variant &var )
{
    if ( patho_var(var) || var.getValue("ACMG") == "P" || var.getValue("ACMG") == "LP" ) {
        return true;
    }
    return false;
}

// questa mi ritorna un TRUE se la var rimasta in SINGLE-HET in almeo un paziente
bool recsingleable( Variant &var ) {
  if ( var.getPzSingleHet().size() > 0 ) { return true; }
  else { return false; }
}

// questa mi ritorna un TRUE se la var finisce in almeno un pz COMPHET
bool reccompable( Variant &var ) {
  try {
    var.getValue("rec");
    return true;
  } catch(...) { return false; }
}

// questa mi ritorna un TRUE se la var finisce in almeno un pz HOM
bool rechomable( Variant &var ) {
  try {
    var.getValue("rec_hom");
    return true;
  } catch(...) { return false; }
}

// questa mi ritorna un TRUE se la var finisce in filtraggio per modelli
bool domable( Variant &var ) {
  try {
    var.getValue("dom");
    return true;
  } catch(...) { return false; }
}

// questa mi dice se la var finisce in almeno un modello
bool dom_pass( Variant &var ) {
  if ( domable(var) ) {
    if (var.getValue("dom") != "NA" ) { return true; }
    else { return false; }
  } else { return false; }
}



/*
  questa crea il map di abbreviazione come FIRST e numero colonna corrispondente come SECOND
  dati una map di vector_string : string dove ci sono i possibili header_atav e corrispondente abbreviazione
  e un vettore string che corrisponde a header_atav del file in analisi
*/

std::map<std::string, int> create_headers_map(std::map<std::vector<std::string>, std::string> &atav_categories_abbreviation_map, std::vector<std::string> &atav_header) {
  std::map<std::string, int> atav_categories_map;
  // apro loop su vector_string : string
  for (auto const& abbr_pair : atav_categories_abbreviation_map) {
    std::ptrdiff_t posHeader;
    // apro loop su vector_string
    for (auto abbr : abbr_pair.first) {
      // se il vettore i viene trovato in header_atav
      if (std::find (atav_header.begin(), atav_header.end(), abbr) != atav_header.end()) {
        // assegna alla posizione il valore corrispondente
        posHeader = find( atav_header.begin(), atav_header.end(), abbr) - atav_header.begin();
        break;
      // se ASSENTE in HEADER mette valore -1
      } else posHeader = -1;
    }
  // inserisco in nuova map
  atav_categories_map.insert(std::make_pair(abbr_pair.second,posHeader));
  }
  return atav_categories_map;
}


/*--------------- READ ATAV FILE -------------------*/
// questa funzione guarda se la var sia nuovo o no e aggiunge la var all apposito vettore
void add_var_pz(
  int &var_new,
  std::string &pz_name,
  std::string &var_name,
  Variant &tmpVar,
  std::vector<Variant> &vVarClass,
  std::vector<std::string> &vVarClassNames
) {
  if (var_new == 0) {
    std::ptrdiff_t posVar = find( vVarClassNames.begin(), vVarClassNames.end(), var_name) - vVarClassNames.begin();
    vVarClass[ posVar ].addPz(pz_name);
  } else {
    tmpVar.addPz(pz_name);
  }
}

// questa funzione guarda se il gene sia nuovo o no e aggiunge il gene all apposito vettore
void add_gene_var(
  int &gene_new,
  std::string &gene_name,
  Gene &tmpGene,
  std::vector<Gene> &vGeneClass,
  std::vector<std::string> &vGeneClassNames,
  std::vector<std::string> &vRiga,
  std::map<std::string, int> &mapAtavHeader
) {
  if (gene_new == 0) {
    std::ptrdiff_t posGene = find( vGeneClassNames.begin(), vGeneClassNames.end(), gene_name) - vGeneClassNames.begin();
    vGeneClass[ posGene ].addVar(vRiga[mapAtavHeader["var_name"]]);
  } else {
    tmpGene.addVar(vRiga[mapAtavHeader["var_name"]]);
  }
}

// questa corregge il field di ATAV del nome gene togliendone le virgolette
std::string correct_gene_name ( std::string &gene_name ) {
  char gSep = '\''; std::vector<std::string> vGeneName; vectorize_column( gene_name, vGeneName, gSep );
  std::string gene_name_correct = gene_name;
  if ( gene_name.find(gSep) != std::string::npos ) {
      gene_name_correct = vGeneName[1];
  }
  return gene_name_correct;
}


void check_empty_char_add_pzVar(
  std::string &abbr,
  int &pos,
  std::vector<std::string> &vRiga,
  pzVars &tmpPzVar
) {
  std::string result;
  if ( pos != -1 && !vRiga[pos].empty() ) result = vRiga[pos];
  else result = "NA";
  tmpPzVar.addAttr( abbr, result );
}

// questa funzione aggiunge le variabili vettoriali da mettere nella classe del paziente
// capendo se pz_new di metterla in tmpPz oppure in corretta posizione del vettore se pz noto
void add_pz_vectors(
  int &pz_new,
  std::string &pz_name,
  Patient &tmpPz,
  std::vector<Patient> &vPatientClass,
  std::vector<std::string> &vPatientClassNames,
  std::vector<std::string> &vRiga,
  std::map<std::string, int> &mapAtavHeader,
  std::vector<std::string>  &patient_var_abbreviations
) {
  // prendo le caratteristiche che sono parte della var del paziente e la metto nel tmpPzVar
  // che mi serve come trucchetto per non stare a scrivere ogni sotto categoria che voglio salvre della variante
  pzVars tmpPzVar;
  for ( auto var_char : patient_var_abbreviations ) {
    // se nome gene aggiungo sia la versione sfigata del nome (ATAV) che quella corretta
    check_empty_char_add_pzVar( var_char, mapAtavHeader[var_char], vRiga, tmpPzVar );
  }
  // metto anche nome corretto del gene
  std::string gene_name = vRiga[mapAtavHeader["gene_name"]];
  if ( gene_name == "NA" ) {
    gene_name = vRiga[mapAtavHeader["updated_gene_name"]];
  }
  tmpPzVar.addAttr( "gene_name_correct", correct_gene_name(gene_name) );
  if (pz_new == 0) {
    std::ptrdiff_t posPz = find( vPatientClassNames.begin(), vPatientClassNames.end(), pz_name) - vPatientClassNames.begin();
    // se VAR NON gia presente in paziente
    std::vector<std::string> vPzVars = vPatientClass[ posPz ].getVars();
    if ( std::find ( vPzVars.begin(), vPzVars.end(), vRiga[mapAtavHeader["var_name"]] ) == vPzVars.end() ) {
        vPatientClass[ posPz ].addVar( vRiga[mapAtavHeader["var_name"]] );
        vPatientClass[ posPz ].addPzVarVector( tmpPzVar );
    }
  } else {
    tmpPz.addVar( vRiga[mapAtavHeader["var_name"]] );
    tmpPz.addPzVarVector( tmpPzVar );
  }
}



void check_empty_char_add_var(
  std::string &abbr,
  int &pos,
  std::vector<std::string> &vRiga,
  Variant &tmpVar
) {
  std::string result;
  if ( pos != -1 && !vRiga[pos].empty() ) result = vRiga[pos];
  else result = "NA";
  tmpVar.addAttr( abbr, result );
}

void check_empty_char_add_pz(
  std::string &abbr,
  int &pos,
  std::vector<std::string> &vRiga,
  Patient &tmpPz
) {
  std::string result;
  if ( pos != -1 && !vRiga[pos].empty() ) {
      result = vRiga[pos];
  }
  else result = "NA";
  tmpPz.addAttr( abbr, result );
}

void check_empty_char_add_gene(
  std::string &abbr,
  int &pos,
  std::vector<std::string> &vRiga,
  Gene &tmpGene
) {
  std::string result;
  /*
  try {
    result = vRiga[pos];
  } catch(...) { result = "NA"; }
  */
  if ( pos != -1 && !vRiga[pos].empty() ) result = vRiga[pos];
  else result = "NA";
  tmpGene.addAttr( abbr, result );
}

// questa funzione per ogni riga del file vettorializzata
// assegna il valore corrispondente all header atav
// alla classe di pertinenza (con NA se vuoto o header assente)
void check_empty_and_add_chars(
  std::vector<Variant> &vVarClass,
  std::vector<std::string> &vVarClassNames,
  std::vector<Patient> &vPatientClass,
  std::vector<std::string> &vPatientClassNames,
  std::vector<Gene> &vGeneClass,
  std::vector<std::string> &vGeneClassNames,
  std::map<std::string, int> &mapAtavHeader,
  std::vector<std::string> &variant_abbreviations,
  std::vector<std::string> &patient_abbreviations,
  std::vector<std::string> &gene_abbreviations,
  std::vector<std::string> &patient_var_abbreviations,
  std::vector<std::string> &vRiga
) {

  // prendo i nomi di var, pz e gene attuale
  std::string var_name = vRiga[mapAtavHeader["var_name"]];
  std::string pz_name = vRiga[mapAtavHeader["pz_name"]];
  std::string gene_name = vRiga[mapAtavHeader["gene_name"]];

  if ( gene_name == "NA" ) {
      gene_name = vRiga[mapAtavHeader["updated_gene_name"]];
  }

  // prendo il nome del gene senza virgolette di ATAV
  std::string gene_name_correct = correct_gene_name( gene_name );
  //std::cout << gene_name << " > " << gene_name_correct << '\n';

  // se GENE NON uno di quelli di interesse esce dalla funzione
  if (std::find (vGeneClassNames.begin(), vGeneClassNames.end(), gene_name_correct) == vGeneClassNames.end()) {
      return;
  }
  // if ( mapGeniInteresse.find(gene_name_correct) == mapGeniInteresse.end() ) { return; }

  // se NON ancora nel vettore delle VAR aggiunge il nome della var al suddetto vettore e idem per gene e paziente
  // e mi ricordo se new (1) o no (0)
  int var_new = add_absent_to_vector( vVarClassNames, var_name );
  int pz_new = add_absent_to_vector( vPatientClassNames, pz_name );
  int gene_new = add_absent_to_vector( vGeneClassNames, gene_name_correct );

  // inizializzo oggetti temporanei
  Variant tmpVar;
  Patient tmpPz;
  Gene tmpGene;

  // aggiungo la var e le vettoriali a relativo pz (new or not)
  add_pz_vectors( pz_new, pz_name, tmpPz, vPatientClass, vPatientClassNames, vRiga, mapAtavHeader, patient_var_abbreviations );
  // aggiungo la var a relativo gene (new or not)
  add_gene_var( gene_new, gene_name_correct, tmpGene, vGeneClass, vGeneClassNames, vRiga, mapAtavHeader );
  // aggiungo il paziente al vettore della variante
  add_var_pz( var_new, pz_name, var_name, tmpVar, vVarClass, vVarClassNames );

  // apro loop su MAP di abbr_header : posizione_header
  for (auto const& vi : mapAtavHeader) {

    // se nuova VAR
    if (var_new == 1) {
      // se abbreviation tra quelle della VAR ne aggiunge i valori se esistenti
      if (std::find (variant_abbreviations.begin(), variant_abbreviations.end(), vi.first) != variant_abbreviations.end()) {
        std::string corresponding_abbr = vi.first; int corresponding_pos = vi.second;
        check_empty_char_add_var( corresponding_abbr, corresponding_pos, vRiga, tmpVar );
      }
    }

    // se PZ nuovo
    if (pz_new == 1) {
      // se char del paziente la aggiungo a tmpPz
      if (std::find (patient_abbreviations.begin(), patient_abbreviations.end(), vi.first) != patient_abbreviations.end()) {
        std::string corresponding_abbr = vi.first; int corresponding_pos = vi.second;
        check_empty_char_add_pz( corresponding_abbr, corresponding_pos, vRiga, tmpPz );
      }
    }

    if (gene_new == 1) {
      if (std::find (gene_abbreviations.begin(), gene_abbreviations.end(), vi.first) != gene_abbreviations.end()) {
        std::string corresponding_abbr = vi.first; int corresponding_pos = vi.second;
        check_empty_char_add_gene( corresponding_abbr, corresponding_pos, vRiga, tmpGene );
      }
    }

  // chiudo loop su header atav
  }
  if ( var_new == 1 ) {
      tmpVar.addAttr("gene_name_correct", gene_name_correct);
      vVarClass.push_back(tmpVar);
  }
  if ( pz_new == 1 ) {
      vPatientClass.push_back(tmpPz);
  }
  if ( gene_new == 1 ) {
      tmpGene.addAttr("gene_name_correct", gene_name_correct);
      vGeneClass.push_back(tmpGene);
  }
}

void extract_variants(
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
) {

    std::string line;
    int lineCount = 0;
    // apro loop sulle righe
    while(std::getline(std::cin, line)) {

      // non prendo la prima riga (header)
      if (lineCount > 0) {

        // vectorize taking in consideration commas as separator IF outside quotes
        std::vector<std::string> vRiga = vectorize_column_comma_quote(line, ',');

        /*
        for (auto v : vRiga) {
            if ( v[0] == '"' ) {
                std::cout << v << '\n';
            }
        }
        */

        // check posizioni e aggiungo a classi rispettive (con NA se vuoti)
        check_empty_and_add_chars( vVarClass, vVarClassNames, vPatientClass, vPatientClassNames, vGeneClass, vGeneClassNames, mapAtavHeader, variant_abbreviations, patient_abbreviations, gene_abbreviations, patient_var_abbreviations, vRiga );

      // chiudo IF non header (lineCount > 0)
      }
      lineCount++;
    //chiudo LOOP sulle linee
  }
}

//questa cambia STD::CIN con il risultato della lettura del file
//poi chiama funzione su quella
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
) {
//  std::string &nomeFile, std::map<std::string, int> &mapAtavHeader, std::vector<std::string> &vVarClassNames, std::vector<Variant> &vVarClass) {

    std::ifstream in(nomeFile.c_str());
    std::streambuf *cinbuf = std::cin.rdbuf(); //save old buf
    std::cin.rdbuf(in.rdbuf()); //redirect std::cin to in.txt!

    extract_variants(
      mapAtavHeader,
      variant_abbreviations,
      patient_abbreviations,
      gene_abbreviations,
      patient_var_abbreviations,
      vVarClassNames,
      vVarClass,
      vPatientClassNames,
      vPatientClass,
      vGeneClassNames,
      vGeneClass
    ); //call function

    std::cin.rdbuf(cinbuf);   //reset to standard input again
}
/*-----------------------------------------------------*/
