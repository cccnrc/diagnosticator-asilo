#include "arrangeResults.h"


/*
    prima le 3 funzioni per mettere a posto i risultati per Variant - Patient - Gene
    dopodiche quella finale che prende i VV creati da queste e li scrive su file
    con i nomi forniti da main()
*/

// funzione di aggiustamento dei risultati della classe Gene
std::vector<std::vector<std::string>> arrange_gene_results
(
    std::vector<Gene> vGeneClass,
    std::vector<std::string> vGeneClassNames,
    std::vector<Variant> vVarClass,
    std::vector<std::string> vVarClassNames,
    int selected_only_patho
)
{
    std::vector<std::vector<std::string>> vvLineGene;
    // apro loop sui geni
    for ( Gene &gene : vGeneClass )
    {
        std::vector<std::string> vLineGene;
        vLineGene.push_back( gene.getValue("gene_name_correct") );
        gene.removeAttr("gene_name_correct");
        // anche per le VAR da stampare nel file del GEEN capisco la var corrispondente come classe Variant
        // per poter scremare e metterla o no sul fatto che sia patho
        std::vector<std::string> vGeneVarPatho;
        for ( std::string var_name : gene.getVars() )
        {
            if ( selected_only_patho == 1 ) {
                std::ptrdiff_t posVar = find( vVarClassNames.begin(), vVarClassNames.end(), var_name ) - vVarClassNames.begin();
                if ( !patho_var(vVarClass[ posVar ]) ) {
                    continue;
                };
            }
            // aggiunge la VAR solo se nuova
            int varPathoGeneNew = add_absent_to_vector( vGeneVarPatho, var_name );
            //vGeneVarPatho.push_back( var_name );
        }
        // e stampo poi il gene solo se trovata almeno una var Pathologica
        if ( vGeneVarPatho.size() > 0 ) {
            gene.addAttr("gene_var_patho",stringize_vector(vGeneVarPatho,","));
            gene.addAttr("inh_mode_all",stringize_vector(gene.getInhModes(),","));
            gene.addAttr("genelists",stringize_vector(gene.getGenelists(),","));
            vLineGene.push_back( stringize_map_str_str(gene.getAll(), "=", ";" ));
            vvLineGene.push_back(vLineGene);
        }
        vGeneVarPatho.clear();
        // questo sotto se le volessi TUTTE (indipendentemente da status patho)
        // gene.addAttr("gene_var_all",stringize_vector(gene.getVars(),","));
        vLineGene.clear();
    } // chiudo loop sui geni
    return vvLineGene;
}

// funzione di aggiustamento dei risultati della classe Patient
// questa ha bisogno di VV negli arg perchè deve riempire anche il VV per il file delle varianti del paziente
std::vector<std::vector<std::string>> arrange_pz_results
(
    std::vector<Patient> vPatientClass,
    std::vector<std::string> vPatientClassNames,
    std::vector<Variant> vVarClass,
    std::vector<std::string> vVarClassNames,
    std::map<std::string, std::vector<std::vector<std::string>>> &mVV,
    int selected_only_patho
)
{
    std::vector<std::vector<std::string>> vvLinePz;
    // apro loop sui pazienti
    for ( Patient &pz : vPatientClass )
    {
        std::vector<std::string> vLinePz;
        // apro loop su vettore della classe var del paziente
        for ( pzVars &varPZ : pz.getPzVarVector() ) {
            std::vector<std::string> vLinePzVar;
            if ( selected_only_patho == 1 ) {
                // capisco la var corrispondente come classe Variant
                // per poter scremare e metterla o no sul fatto che sia patho
                std::ptrdiff_t posVar = find( vVarClassNames.begin(), vVarClassNames.end(), varPZ.getValue("var_name")) - vVarClassNames.begin();
                if ( !patho_var(vVarClass[ posVar ]) ) {
                    continue;
                };
            }

            // e qui aggiungo le caratteristiche della variante del paziente per stamparla poi in file
            std::vector<std::string> vVarName; vectorize_column( varPZ.getValue("var_name"), vVarName, '-' );
            varPZ.removeAttr("var_name");
            vLinePzVar.push_back( vVarName[0] );
            vLinePzVar.push_back( vVarName[1] );
            vLinePzVar.push_back( "." );
            vLinePzVar.push_back( vVarName[2] );
            vLinePzVar.push_back( vVarName[3] );
            vLinePzVar.push_back( "." );
            vLinePzVar.push_back( "." );
            varPZ.removeAttr("gene_name");
            varPZ.addAttr( "pz_name", pz.getValue("pz_name") );
            vLinePzVar.push_back( stringize_map_str_str(varPZ.getAll(), "=", ";" ));
            // se CHR non ancora tra le key della map inserisce VV come elemento
            if ( mVV.find( vVarName[0] ) == mVV.end() ) {
                std::vector<std::vector<std::string>> vvLinePzVar = { vLinePzVar };
                mVV.insert(std::make_pair(vVarName[0], vvLinePzVar));
            }
            // altrimenti appende la linea al VV rispettivo
            else {
                mVV[ vVarName[0] ].push_back(vLinePzVar);
            }
            vLinePzVar.clear();
        } // chiudo loop su VAR del PAZIENTE

        // col 1 patient name che poi rimuovo (per non ristamparlo)
        vLinePz.push_back( pz.getValue("pz_name") );
        //pz.removeAttr("pz_name");

        // anche per le VAR da stampare nel file del PAZIENTE capisco la var corrispondente come classe Variant
        // per poter scremare e metterla o no sul fatto che sia patho
        std::vector<std::string> vPzVar;
        // e apro vettore per tenere SOLO le pathogeniche
        std::map<std::string, std::vector<std::string>> mPzVarPathoOnly;
        for ( std::string var_name : pz.getVars() )
        {
            std::ptrdiff_t posVar = find( vVarClassNames.begin(), vVarClassNames.end(), var_name ) - vVarClassNames.begin();

            // Aug 16, 2019
            // questa parte mi serve per sistemare bene le caratt della Var per ogni paziente
            // perche se si basa solo su quelle scritte in classe Variant allora se ha trovato
            // almeno 1 COMPHET segnerebbe questa come comphet in ogni paziente in cui trovata la var
            // visto che gli facevo scrivere var.getValue("rec")
            // invece cosi scrive il valore reale in quel paziente anche nel file .pz
            std::string real_comphet_this_patient;
            if ( vVarClass[posVar].check_pzComphet( pz.getValue("pz_name") ) ) {
                real_comphet_this_patient = "comphet";
            }
            else {
                // se NON tra i comphet guarda se tra gli HOM o tra i SINGLE_HET
                if ( vVarClass[posVar].check_pzHom( pz.getValue("pz_name") ) ) {
                    real_comphet_this_patient = "hom";
                }
                else {
                    if ( vVarClass[posVar].check_pzSingleHet( pz.getValue("pz_name") ) ) {
                        real_comphet_this_patient = "single_het";
                    }
                }
            }


            if ( patho_acmg_var(vVarClass[ posVar ]) ) {
                // aggiungo alla KEY = var_name, un V composto da: categoria ACMG, ACMG_bayesian, priority, dom, rec_hom, rec
                std::vector<std::string> vChar = { vVarClass[posVar].getValue("ACMG"), vVarClass[posVar].getValue("ACMG_bayesian"), vVarClass[posVar].getValue("priority"), vVarClass[posVar].getValue("dom"), vVarClass[posVar].getValue("rec_hom"), real_comphet_this_patient };
                mPzVarPathoOnly.insert(std::make_pair(var_name, vChar));
            }
            // se chieste SOLO patho sono le corrispondenti var come KEY in map
            if ( selected_only_patho == 1 ) {
                vPzVar = convert_map_keys_to_V( mPzVarPathoOnly );
            }
            // altrimenti gli mette i nomi di TUTTE
            else {
                vPzVar.push_back( var_name );
            }
        }
        // poi li aggiungo alla Classe Variant cosi poi li stamperà
        if ( vPzVar.size() > 0 ) {
            pz.addAttr("var_all",stringize_vector(vPzVar,","));
            vPzVar.clear();
        }
        else {
            pz.addAttr("var_all","NA");
        }
        if ( mPzVarPathoOnly.size() > 0 ) {
            pz.addAttr("var_only_patho",stringize_map_str_V(mPzVarPathoOnly,"|",":",","));
        }
        else {
            pz.addAttr( "var_only_patho", "NA" );
        }
        // questa INVECE se le volessi stampare TUTTE (indipendentemente da Pathogenicity)
        // pz.addAttr("var_all",stringize_vector(pz.getVars(),","));
        // col 2 tutti gli attributi della classe con un = davanti
        vLinePz.push_back( stringize_map_str_str(pz.getAll(), "=", ";" ));

        vvLinePz.push_back(vLinePz);
        vLinePz.clear();

    } // chiudo loop sui pazienti
    return vvLinePz;
}

// funzione di aggiustamento dei risultati della classe Variant
std::map<std::string, std::vector<std::vector<std::string>>> arrange_var_results
(
    std::vector<Variant> vVarClass,
    std::vector<std::string> vVarClassNames,
    int selected_only_patho
)
{
      std::map<std::string, std::vector<std::vector<std::string>>> mVVvar;
      // converto la classe Var in VV per stamparlo
      for (auto &var : vVarClass)
      {
          // inserisco i pazienti nei relativi campi come attribute della Var
//          check_and_add_patients( var );

          //std::cout << var.getValue("var_name") << '\n';
          // questo fa stampare la var SOLO se finita in almeno una ctegria patho se selected_only_patho == 1
          // altrimenti stampa TUTTE le var
          if ( selected_only_patho == 1 ) {
            if ( patho_var(var) ) {
                ;
            }
            else {
                continue;
            }
          }
          else {
              ;
          }

          std::string nome_var = var.getValue("var_name");

          // qui metto la KEY degli attributi delle VAR che NON voglio siano stampati nel file e poi li rimuovo nel loop sotto
          std::vector<std::string> vAttrToRemove = {
              "ref",
              "alt",
              "gene_name",
              "updated_gene_name"
          };
          for ( std::string x : vAttrToRemove )
          {
              var.removeAttr(x);
          }

          // apro vettore vuoto per tenere i risultati della singola VAR da push_back poi in VV
          std::vector<std::string> vLineVar;

          if ( var.getAcmg().size() > 0 ) {
              var.addAttr("ACMG_categories", stringize_vector(var.getAcmg(), ","));
          }
          else {
              var.addAttr("ACMG_categories", "NA" );
          }
          // deframmento il nome var
          std::vector<std::string> vNome; vectorize_column(nome_var, vNome, '-');
          // creo linea della VAR nel file finale
          vLineVar.push_back( vNome[0] );
          vLineVar.push_back( vNome[1] );
          if ( var.getValue("rs_name") == "NA" ) {
              vLineVar.push_back(".");
          }
          else {
              vLineVar.push_back( var.getValue("rs_name") );
          }
          vLineVar.push_back( vNome[2] );
          vLineVar.push_back( vNome[3] );
          vLineVar.push_back( "." );
          vLineVar.push_back( "." );
          vLineVar.push_back( stringize_map_str_str(var.getAll(), "=", ";" ));
          // se CHR non ancora tra le key della map inserisce VV come elemento
          if ( mVVvar.find( vNome[0] ) == mVVvar.end() ) {
              std::vector<std::vector<std::string>> vvLineVar = { vLineVar };
              mVVvar.insert(std::make_pair(vNome[0], vvLineVar));
          }
          // altrimenti appende la linea al VV rispettivo
          else {
              mVVvar[ vNome[0] ].push_back(vLineVar);
          }
          vLineVar.clear();
      } //chiudo LOOP sulle Variant

      return mVVvar;

}




/*
    questa finale per mettere insieme le 3 sopra e scrivere nei file
*/

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
)
{

    // questa prima riempie riempie la MAP <chr, VV> con le linee di ogni Variant
    std::map<std::string, std::vector<std::vector<std::string>>> mVVvar = arrange_var_results ( vVarClass, vVarClassNames, selected_only_patho );
    // questa seconda riempie la MAP <chr, VV> con le linee di ogni Var del Paziente nel chromosoma (nome paziente nei campi info)
    std::map<std::string, std::vector<std::vector<std::string>>> mVV;
    std::vector<std::vector<std::string>> vvLinePz = arrange_pz_results( vPatientClass, vPatientClassNames, vVarClass, vVarClassNames, mVV, selected_only_patho );
    // questa terza riempie le linee VV per ogni gene
    std::vector<std::vector<std::string>> vvLineGene = arrange_gene_results( vGeneClass, vGeneClassNames, vVarClass, vVarClassNames, selected_only_patho );

    write_map_VV ( mVVvar, nomeFileResVar, ".", ".vcf" );
    write_map_VV ( mVV, nomeFileResVarPz, ".", ".vcf" );

    write_cin_file_vv( nomeFileResPz, vvLinePz );
    write_cin_file_vv( nomeFileResGene, vvLineGene );

}
