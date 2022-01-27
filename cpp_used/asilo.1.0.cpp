#include "dependencies/atavCategories.h"
#include "dependencies/filterFunctions.h"
#include "dependencies/atavClasses.h"
#include "dependencies/atavHeaderNames.h"
#include "dependencies/basicFunctions.h"
#include "dependencies/readGenelist.h"
#include "dependencies/arrangeResults.h"
#include "dependencies/trioComphet.h"

/*
  PRE:
    ./mergeGenelists.0 <gene-list-0> <gene-list-1> <...> <gene-list-N>
*/

/*
trio_file_denovo_lab="/Users/cccnrc/Documents/columbia/cardio/cardio0/cardio_analisi0/2019-06-04_14-31-00_cardio_analisi0_denovoandhom.csv"
trio_file_comphet_lab="/Users/cccnrc/Documents/columbia/cardio/cardio0/cardio_analisi0/2019-06-04_14-31-00_cardio_analisi0_comphet.csv"
trio_file_denovo_home="/home/enrico/columbia/cardio/cardio0/cardio_analisi0/2019-06-04_14-31-00_cardio_analisi0_denovoandhom.csv"
trio_file_comphet_home="/home/enrico/columbia/cardio/cardio0/cardio_analisi0/2019-06-04_14-31-00_cardio_analisi0_comphet.csv"
ckd_file_home="/home/enrico/columbia/ckd/ckd_analisi0/2019-03-27_11-13-23_ckd_analisi0_genotypes.csv"
prova_file="ckd.500k.csv"
gcc -O3 -lstdc++ --std=c++11 -g -lm \
    asilo.1.0.cpp\
    dependencies/readGenelist.cpp \
    dependencies/filterFunctions.cpp \
    dependencies/atavCategories.cpp \
    dependencies/basicFunctions.cpp \
    dependencies/arrangeResults.cpp \
    dependencies/trioComphet.cpp \
    -o asilo.1.0 \
    && time ./asilo.1.0 \
      /Users/enrico/columbia2/home_prove_asilo/cpp_used/prove_asilo/ckd_head10k.csv \
      NONE \
      /Users/enrico/columbia2/home_prove_asilo/cpp_used/prove_asilo/new_kidney_genelist.gl \
      analisi_results \
      prova.var \
      prova.pz \
      prova.pzvar \
      prova.gene \
      "10E-5" \
      "0" \
      "0" \
      "gnomad_ex_global_af,gnomad_ex_afr_af,gnomad_gen_global_af,gnomad_gen_afr_af" \
      "5"
*/

/* DEBUG:
G_SLICE=always-malloc \
  G_DEBUG=gc-friendly \
  valgrind -v --tool=memcheck --leak-check=full \
  --num-callers=40 \
  --log-file=valgrind.log \
  $(which ./asilo.1.0) \
      "$trio_file_denovo_home" \
      "$trio_file_comphet_home" \
      merged_gene_lists.AR.txt \
      results3 \
      prova.var \
      prova.pz \
      prova.pzvar \
      prova.gene

      /home/enrico/columbia/asilo2/sito/app/users/enrico7/dm0/upload/atav.csv \
      NONE \
      /home/enrico/columbia/asilo2/sito/app/users/enrico7/dm0/result/merged_gene_lists.txt \
      /home/enrico/columbia/asilo2/sito/app/users/enrico7/dm0/result/analisi_result prova.var \
      prova.pz \
      prova.pzvar \
      prova.gene \
      10E-4 \
      0 \
      0 \
      gnomad_ex_controls_afr_af,gnomad_gen_controls_afr_af,gnomad_ex_controls_amr_af,gnomad_gen_controls_amr_af,gnomad_ex_controls_asj_af,gnomad_gen_controls_asj_af,gnomad_ex_controls_eas_af,gnomad_gen_controls_eas_af,gnomad_ex_controls_fin_af,gnomad_gen_controls_fin_af \
      5
*/

int check_arg_num ( int &argc, char const *argv[], int &numArg );

struct stat info;


/* -------- MAIN FUNCTION -------- */
int main(int argc, char const *argv[]) {

  // metto i geni da carcare (dovro prenderlo da FILE poi)
  //std::map<std::string, std::string> vGeniInteresse = { {"SYMPK", "AD"}, {"COL26A1", "AR"}, {"RHD", "U"}};


  //controllo argomenti
  int numArg = 14; int checkRes = check_arg_num ( argc, argv, numArg ); if ( checkRes != 0 ) return 1;
  //aggiudico a nome file e output
  std::string nomeFile = argv[1];
  std::string nomeFileComphet = argv[2]; // questo OTIONAL (serve per TRIO analisi) quindi devo poi controllare che non sia NONE
  std::string nomeFileGenelist = argv[3];
  std::string nomeDirResults = argv[4];
  //std::string nomeFileResVar = add_results_dir_path( nomeDirResults, argv[4] );
  std::string nomeFileResPz = add_results_dir_path( nomeDirResults, argv[6] );
  //std::string nomeFileResVarPz = add_results_dir_path( nomeDirResults, argv[6] );
  std::string nomeFileResGene = add_results_dir_path( nomeDirResults, argv[8] );

  // da far inserire a User quando lancia analisi
  //float rare_threshold = 10E-4;
  float rare_threshold = std::stoi( argv[9] );
  //int parents_aff = 0; // (da aggiungere tra selezioni utente)
  int parents_aff = std::stoi( argv[10] );
  //int selected_only_patho = 0; // (da aggiungere tra selezioni utente) --> questo decide se stampare SOLO le patho oppure tutte nei geni selezionati
  int selected_only_patho = std::stoi( argv[11] );
  //std::vector<std::string> vPopulationSelected = {"gnomad_ex_global_af", "gnomad_ex_afr_af","gnomad_gen_global_af", "gnomad_gen_afr_af"};
  std::string populationSelectedInput = argv[12];
  // num di AC MAX dei casi da accettare nel modello DOM
  int casesMaxAC = std::stoi( argv[13] );

  std::vector<std::string> vPopulationSelected = vectorize_column_v2( populationSelectedInput, ',' );





  // creo main DIR results
  if ( ! create_dir( nomeDirResults ) ) {
      return 1;
  }
  // creo sub-DIR VAR per tenerci i risultati di Variant divisi per chromosoma
  std::string nomeFileResVar = add_results_dir_path( nomeDirResults, add_results_dir_path( "var_results", argv[5] ));
  if ( !create_dir( add_results_dir_path( nomeDirResults, "var_results" ) ) ) {
      return 1;
  }
  // creo sub-DIR var_pz_results per tenerci i risultati di varPz divisi per chromosoma
  std::string nomeFileResVarPz = add_results_dir_path( nomeDirResults, add_results_dir_path( "var_pz_results", argv[7] ));
  if ( !create_dir( add_results_dir_path( nomeDirResults, "var_pz_results" ) ) ) {
      return 1;
  }


  // prendo header di ATAV file
  std::vector<std::string> atav_header = getHeader( nomeFile );
  // assegno ad ogni abbreviazione il numero colonna corrispondente in file ATAV (tramite atav_categories_abbreviation_map che si trova in atavHeaderNames.h)
  std::map<std::string, int> mapAtavHeader = create_headers_map( atav_categories_abbreviation_map, atav_header );
  //for (auto const& vi : mapAtavHeader) std::cout << vi.first << ": " << vi.second << '\n';

  //controllo che ci siano i valori minimi per fare analisi decentemente
  std::vector<std::string> vHeaderNecessari = {
      "effect",
      "hgvs_p",
      "gene_name",
      "updated_gene_name",
      "pz_name",
      "gt",
      "clinvar_clinrevstars",
      "clinvar_clinsig",
      "clinvar_disease",
      "humdiv_cat",
      "humvar_cat",
      "gnomad_ex_global_af",
      "gnomad_ex_global_an",
      "gnomad_ex_controls_af",
      "gnomad_ex_controls_an",
      "gnomad_gen_global_af",
      "gnomad_gen_global_an",
      "gnomad_gen_controls_af",
      "gnomad_gen_controls_an"
  };
  // se ALMENO uno non trovato manda in stop analisi
  int passedCheck = 1;
  for ( std::string header : vHeaderNecessari )
  {
      if ( mapAtavHeader[ header ] == -1 ) {
          std::cout << "" << '\n';
          std::cout << "==> !!! ERROR !!! <==" << '\n';
          std::cout << "======> header necessary missing: " << header << " <======" << '\n';
          std::cout << "=========> aborting analysis <=========" << '\n';
          std::cout << "" << '\n';
          passedCheck = 0;
      }
  }
  if ( passedCheck == 0 ) {
      return 1;
  }

    std::vector<std::string> vVarClassNames;
    std::vector<Variant> vVarClass;
    std::vector<std::string> vPatientClassNames;
    std::vector<Patient> vPatientClass;
    // inizializzo vettore dei geni e lo riempio con quello delle genelists dal file merged_gene_lists
    std::vector<Gene> vGeneClass;
    std::vector<std::string> vGeneClassNames = get_gene_lists(nomeFileGenelist, vGeneClass);

    // decido se deve fare o meno analisi COMPHET
    // (in funzione del fatto che abbia messo un file delle vere COMPHET (TRIO) oppure no)
    int comphet_analysis = 1;
    std::map < std::string, std::map < std::string, std::vector<std::string >>> mPatientsComphetVars;

    if ( nomeFileComphet != "NONE" ) {
        comphet_analysis = 0;
        get_comphet_trio_file (
            atav_categories_abbreviation_map,
            variant_abbreviations,
            patient_abbreviations,
            patient_var_abbreviations,
            nomeFileComphet,
            vPatientClass,
            vPatientClassNames,
            vGeneClass,
            vGeneClassNames,
            vVarClass,
            vVarClassNames,
            mPatientsComphetVars
        );
    }

  // riempio classi da ATAV file
  // NB: prende solo le VAR dei GENI in vGeneClassNames (riempita da genelist file)
  getVariants( nomeFile, mapAtavHeader, variant_abbreviations, patient_abbreviations, gene_abbreviations, patient_var_abbreviations, vVarClassNames, vVarClass, vPatientClassNames, vPatientClass, vGeneClassNames, vGeneClass );

  std::cout << "----> patients (n): " << vPatientClassNames.size() << '\n';
  std::cout << "----> gene (n): " << vGeneClassNames.size() << '\n';
  std::cout << "----> variants (n): " << vVarClassNames.size() << '\n';

  // questa assegna i filtri e le relative categorie DOM a TUTTE le varianti passate
  // cosi con questi risultati gia dispoibili posso fare la seconda parte comphet
  filter_variants( vVarClass, vPopulationSelected, vGeneClass, vGeneClassNames, rare_threshold, casesMaxAC );
  filter_variants_acmg ( vVarClass, vGeneClass, vGeneClassNames, parents_aff );

  // ultimo campo indica se richiesta o meno l analisi COMPHET
  filter_variants_recessive ( vPatientClass, vGeneClass, vGeneClassNames, vVarClass, vVarClassNames, comphet_analysis );
  assign_priority( vVarClass );


  // questa la funzione per fare lo screening di PM3 (che essendo relativo a essere in trans con altra variante PATHO deve essere fatto DOPO lo screening di tutte le vars e collezionate le COMPHET)
  acmg_pm3_screening( vVarClass, vVarClassNames );

  // questa per togliere le varianti venute P/LP ma con priority troppo bassa
  acmg_priority_screening( vVarClass, vVarClassNames );

  // funzione riassuntiva di arrangeResults.h per mettere a posto poi scrivere i risultati
  arrange_and_write_results( vVarClass, vVarClassNames, selected_only_patho, nomeFileResVar, vPatientClass, vPatientClassNames, nomeFileResPz, nomeFileResVarPz, vGeneClass, vGeneClassNames, nomeFileResGene );

  return 0;

}





/* ----------- CHECK ARGUMENTS NUMBER ----------- */
int check_arg_num ( int &argc, char const *argv[], int &numArg ) {
  //CONTROLLO ARGOMENTI PASSATI
  if ( argc == numArg ) std::cout << "num arg correct: " << '\n';
  else {
    std::cout << '\n';
    std::cout << "PLEASE CALL with: ./asilo <input-file> <genelist-file> <output-file>" << '\n';
    return 1;
  }
  std::cout << "-> INPUT file: " << argv[1] << '\n';
  std::cout << "-> INPUT file COMPHET (optional): " << argv[2] << '\n';
  std::cout << "-> GENELIST file: " << argv[3] << '\n';
  std::cout << "-> RESULTS DIR: " << argv[4] << '\n';
  std::cout << "-> OUTPUT file VAR: " << argv[5] << '\n';
  std::cout << "-> OUTPUT file PZ: " << argv[6] << '\n';
  std::cout << "-> OUTPUT file VAR PZ: " << argv[7] << '\n';
  std::cout << "-> OUTPUT file GENE: " << argv[8] << '\n';
  return 0;
}
