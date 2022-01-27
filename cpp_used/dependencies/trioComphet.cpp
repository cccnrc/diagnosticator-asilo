#include "trioComphet.h"

/*
    questa la funzione per prendere le VAR dal file COMPHET del TRIO analysis
    e aggiungerle alla classe Variant e riepire la map< str, map< str, vec>>>
    di PZ: GENE: vars-in-comphet
*/

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
)
{

    std::vector<std::string> vHeaderTrioComphet = getHeader( nomeFileTrioComphet );

    // e la MAP< int, str > dove tenere il numero della colonna relativa e il nome
    // NB metto il NUM come KEY per non dover fare 3 loop diversi quando apro poi loop sulle righe del file (che accedo per POS)
    std::map<int, std::string> mHeaderTrioComphetVar1;
    std::map<int, std::string> mHeaderTrioComphetVar2;
    std::map<int, std::string> mHeaderTrioComphetVar12;
    // e ci metto i relativi header
    int headerCount = 0;
    for ( auto field : vHeaderTrioComphet )
    {
        if ( field.find("(#1)") != std::string::npos ) {
            replaceAll( field, " (#1)", "");
            if ( !mHeaderTrioComphetVar1.insert( std::make_pair( headerCount, field ) ).second ) {
                std::cout << "!!! ATTENTION !!! Field already presenty when trying to insert VAR1 header column name-pos in VAR1 MAP: " << field << '\n';
            }
        }
        else if ( field.find("(#2)") != std::string::npos ) {
            replaceAll( field, " (#2)", "");
            if ( !mHeaderTrioComphetVar2.insert( std::make_pair( headerCount, field ) ).second ) {
                std::cout << "!!! ATTENTION !!! Field already presenty when trying to insert VAR2 header column name-pos in VAR2 MAP: " << field << '\n';
            }
        }
        else {
            // aggiungo anche family, mother and father alle map delle VAR1 e VAR2 che mi servono per la Classe Variant
            if ( field == "Family ID" || field == "Mother" || field == "Father" ) {

                if ( !mHeaderTrioComphetVar1.insert( std::make_pair( headerCount, field ) ).second ) {
                    std::cout << "!!! ATTENTION !!! Field already presenty when trying to insert VAR1 header column name-pos in VAR1 MAP: " << field << '\n';
                }
                if ( !mHeaderTrioComphetVar2.insert( std::make_pair( headerCount, field ) ).second ) {
                    std::cout << "!!! ATTENTION !!! Field already presenty when trying to insert VAR2 header column name-pos in VAR2 MAP: " << field << '\n';
                }
            }
            // e gli altri li metto invece in VAR12
            else {
                if ( !mHeaderTrioComphetVar12.insert( std::make_pair( headerCount, field ) ).second ) {
                    std::cout << "!!! ATTENTION !!! Field already presenty when trying to insert VAR12 header column name-pos in VAR12 MAP: " << field << '\n';
                }
            }

        } // chiudo ELSE (VAR12)
        headerCount++;
    } // chiude loop su HEADER

    // apro la MAP per le ABBREVIAZIONI
    std::map<int, std::string> mHeaderTrioComphetAbbreviationVar1;
    // e un vettore per tenere le KEY (int) delle COLONNE che non mi interessano (da rimuvoere da MAP)
    std::vector<int> vKeysToRemoveVar1;
    for ( auto const &field : mHeaderTrioComphetVar1 )
    {
        int found  = 0;
        // poi apro loop su MAP <std::vector<std::string>, std::string> delle atav_categories_abbreviation_map
        for ( auto const &header : atav_categories_abbreviation_map )
        {
            // se trova il field in uno dei valori della lista possibilita HEADER allora aggiorna il nome del campo con ABBREVIAZIONE
            if (std::find (header.first.begin(), header.first.end(), field.second) != header.first.end()) {
                mHeaderTrioComphetVar1[field.first] = header.second;
                found = 1;
            }
        }
        // se NON trovato tra le ABBR richieste ( quindi NON mi interessa di base) lo scarto
        if (found == 0) {
            vKeysToRemoveVar1.push_back( field.first );
        }
    }
    // quindi RIMUOVO le KEYS che non mi interessano
    for (int i : vKeysToRemoveVar1)
    {
        mHeaderTrioComphetVar1.erase( i );
    }

    std::map<int, std::string> mHeaderTrioComphetAbbreviationVar2;
    std::vector<int> vKeysToRemoveVar2;
    for ( auto const &field : mHeaderTrioComphetVar2 )
    {
        int found = 0;
        // poi apro loop su MAP <std::vector<std::string>, std::string> delle atav_categories_abbreviation_map
        for ( auto const &header : atav_categories_abbreviation_map )
        {
            // se trova il field in uno dei valori della lista possibilita HEADER allora aggiorna il nome del campo con ABBREVIAZIONE
            if (std::find (header.first.begin(), header.first.end(), field.second) != header.first.end()) {
                mHeaderTrioComphetVar2[field.first] = header.second;
                found = 1;
            }
        }
        // se NON trovato tra le ABBR richieste ( quindi NON mi interessa di base) lo scarto
        if (found == 0) {
            vKeysToRemoveVar2.push_back( field.first );
        }
    }
    // quindi RIMUOVO le KEYS che non mi interessano
    for (int i : vKeysToRemoveVar2)
    {
        mHeaderTrioComphetVar2.erase( i );
    }

    // poi fa un CHECK degli HEADER NON trovati VAR1
    std::vector<std::string> vHeaderFoundVar1;
    for ( auto const &field : mHeaderTrioComphetVar1 )
    {
        vHeaderFoundVar1.push_back( field.second );
    }
    std::vector<std::string> vHeaderNotFoundVar1;
    for ( auto const &header : atav_categories_abbreviation_map )
    {
        if (std::find (vHeaderFoundVar1.begin(), vHeaderFoundVar1.end(), header.second) == vHeaderFoundVar1.end()) {
            vHeaderNotFoundVar1.push_back( header.second );
        }
    }
    // poi fa un CHECK degli HEADER NON trovati VAR2
    std::vector<std::string> vHeaderFoundVar2;
    for ( auto const &field : mHeaderTrioComphetVar2 )
    {
        vHeaderFoundVar2.push_back( field.second );
    }
    std::vector<std::string> vHeaderNotFoundVar2;
    for ( auto const &header : atav_categories_abbreviation_map )
    {
        if (std::find (vHeaderFoundVar2.begin(), vHeaderFoundVar2.end(), header.second) == vHeaderFoundVar2.end()) {
            vHeaderNotFoundVar2.push_back( header.second );
        }
    }


    /*
        !!! DA IMPLEMENTARE !!!
        se non ne trova alcuni chiave andare a bloccare analisi
        (idem se il valore della map quando riprende gli header == -1 in asilo main)
    */

    // se qualcuno NON trovato lo scrive in OUTPUT
    if ( vHeaderNotFoundVar1.size() > 0 ) {
        for ( std::string x : vHeaderNotFoundVar1 )
        {
            std::cout << " --!!!---> NOT FOUND COMPHET in VAR1, headers : " << x << '\n';
        }
    }
    if ( vHeaderNotFoundVar2.size() > 0 ) {
        for ( std::string x : vHeaderNotFoundVar2 )
        {
            std::cout << " --!!!---> NOT FOUND COMPHET in VAR2, headers : " << x << '\n';
        }
    }


    // prendo quindi il VV del file
    std::vector<std::vector<std::string>> vvComphetTrio = getVVectorFile( nomeFileTrioComphet, ',' );
    // toglo HEADER
    vvComphetTrio.erase(vvComphetTrio.begin());
    // e metto in un V la MAP < str, str > per tenere nome col e relativo risultato per ogni VAR
    std::vector<std::map<std::string, std::string>> vmVar1;
    std::vector<std::map<std::string, std::string>> vmVar2;
    std::vector<std::map<std::string, std::string>> vmVar12;
    for ( std::vector<std::string> v : vvComphetTrio )
    {
        // apro la MAP per i field VAR1
        std::map<std::string, std::string> mVar1;
        std::map<std::string, std::string> mVar2;
        std::map<std::string, std::string> mVar12;
        // apro contatore e loop sulle colonne del file
        int fCount = 0;
        for ( std::string f : v )
        {

            if ( mHeaderTrioComphetVar1.find(fCount) != mHeaderTrioComphetVar1.end() ) {
                if ( !mHeaderTrioComphetVar1[fCount].empty() ) {
                    if ( !mVar1.insert( std::make_pair( mHeaderTrioComphetVar1[fCount], f ) ).second ) {
                        std::cout << "!!! ATTENTION !!! Field already presenty when trying to insert VAR1 VALUE in VAR1 sub-MAP: " << mHeaderTrioComphetVar1[fCount] << " --> " << f << " at position: " << fCount << '\n';
                    }
                }
            }

            if ( mHeaderTrioComphetVar2.find(fCount) != mHeaderTrioComphetVar2.end() ) {
                if ( !mHeaderTrioComphetVar2[fCount].empty() ) {
                    if ( !mVar2.insert( std::make_pair( mHeaderTrioComphetVar2[fCount], f ) ).second ) {
                        std::cout << "!!! ATTENTION !!! Field already presenty when trying to insert VAR2 VALUE in VAR2 sub-MAP: " << mHeaderTrioComphetVar2[fCount] << " --> " << f << '\n';
                    }
                }
            }
            if ( mHeaderTrioComphetVar12.find(fCount) != mHeaderTrioComphetVar12.end() ) {
                if ( mHeaderTrioComphetVar12.find(fCount) != mHeaderTrioComphetVar2.end() ) {
                    if ( !mVar12.insert( std::make_pair( mHeaderTrioComphetVar12[fCount], f ) ).second ) {
                        std::cout << "!!! ATTENTION !!! Field already presenty when trying to insert VAR12 VALUE in VAR12 sub-MAP: " << mHeaderTrioComphetVar12[fCount] << " --> " << f << '\n';
                    }
                }
            }

            // inserisco il nome GENE corretto nelle MAP VAR1 e 2
            if ( mHeaderTrioComphetVar1[fCount] == "gene_name" ) {
                if ( !mVar1.insert( std::make_pair( "gene_name_correct", f ) ).second ) {
                    std::cout << "!!! ATTENTION !!! Field already presenty when trying to insert VAR1 gene_name_correct in VAR1 sub-MAP: --> " << f << '\n';
                }
            }
            if ( mHeaderTrioComphetVar2[fCount] == "gene_name" ) {
                if ( !mVar2.insert( std::make_pair( "gene_name_correct", f ) ).second ) {
                    std::cout << "!!! ATTENTION !!! Field already presenty when trying to insert VAR1 gene_name_correct in VAR2 sub-MAP: --> " << f << '\n';
                }
            }


            fCount++;
        }

        // poi devo mettere in mVar12 il nome delle due VAR a cui fa riferitmento
        if ( !mVar12.insert( std::make_pair( "var1_name", mVar1["var_name"] ) ).second ) {
            std::cout << "!!! ATTENTION !!! Field already presenty when trying to insert VAR1 NAME " << mVar1["Variant ID"] << " in VAR12 sub-MAP" << '\n';
        }
        if ( !mVar12.insert( std::make_pair( "var2_name", mVar2["var_name"] ) ).second ) {
            std::cout << "!!! ATTENTION !!! Field already presenty when trying to insert VAR2 NAME " << mVar2["Variant ID"] << " in VAR12 sub-MAP" << '\n';
        }

        // SAMPLE NAME
        if ( mVar1["pz_name"] != mVar2["pz_name"] ) {
            std::cout << "!!! ATTENTION !!! sample name not matching for COMPHET TRIO VARS: " << mVar1["var_name"] << ": " << mVar1["pz_name"] << " and " << mVar2["var_name"] << ": " << mVar2["pz_name"] << " adding the FIRST ONE" << '\n';
        }
        if ( !mVar12.insert( std::make_pair( "pz_name", mVar1["pz_name"] ) ).second ) {
            std::cout << "!!! ATTENTION !!! Field already presenty when trying to insert VAR2 PZ_NAME " << mVar1["pz_name"] << " in VAR12 sub-MAP" << '\n';
        }

        // GENE NAME
        if ( mVar1["gene_name"] != mVar2["gene_name"] ) {
            std::cout << "!!! ATTENTION !!! gene name not matching for COMPHET TRIO VARS: " << mVar1["var_name"] << ": " << mVar1["gene_name"] << " and " << mVar2["var_name"] << ": " << mVar2["gene_name"] << " adding the FIRST ONE" << '\n';
        }
        // tolgo single quote
        if ( !mVar12.insert( std::make_pair( "gene_name_correct", mVar1["gene_name"] ) ).second ) {
            std::cout << "!!! ATTENTION !!! Field already presenty when trying to insert VAR2 PZ_NAME " << mVar1["gene_name_correct"] << " in VAR12 sub-MAP" << '\n';
        }
        replaceAll(mVar1["gene_name_correct"],"\'", "" );
        replaceAll(mVar2["gene_name_correct"],"\'", "" );
        replaceAll(mVar12["gene_name_correct"],"\'", "" );

        // aggiungo un campo FIELD COMPHET from TRIO FILE
        if ( !mVar1.insert( std::make_pair( "comphet_flag", "TRIO_file" ) ).second ) {
            std::cout << "!!! ATTENTION !!! Field already presenty when trying to insert VAR2 PZ_NAME " << mVar1["pz_name"] << " in VAR12 sub-MAP" << '\n';
        }
        if ( !mVar2.insert( std::make_pair( "comphet_flag", "TRIO_file" ) ).second ) {
            std::cout << "!!! ATTENTION !!! Field already presenty when trying to insert VAR2 PZ_NAME " << mVar1["pz_name"] << " in VAR12 sub-MAP" << '\n';
        }

        // infine inserisco la relativa MAP (ogni VAR) nella relativa V MAP
        vmVar1.push_back( mVar1 );
        mVar1.clear();
        vmVar2.push_back( mVar2 );
        mVar2.clear();
        vmVar12.push_back( mVar12 );
        mVar12.clear();

    }

    // devo poi convertire ogni MAP in un OBJECT della classe Variant
    //std::vector<Variant> vVarClass;
    //std::vector<std::string> vVarClassNames;
    for ( auto m : vmVar1 )
    {
        std::string var_name = m["var_name"];
        std::string pz_name = m["pz_name"];
        std::string gene_name = m["gene_name_correct"];

        // se NON in gene di interessa skippa al prossimo
        if (std::find (vGeneClassNames.begin(), vGeneClassNames.end(), gene_name) == vGeneClassNames.end()) {
            continue;
        }

        // se NON gia inserita
        if (std::find (vVarClassNames.begin(), vVarClassNames.end(), var_name) == vVarClassNames.end()) {
            Variant tmp;
            // apro loop su MAP della var
            for ( auto const &field : m )
            {
                if (std::find (variant_abbreviations.begin(), variant_abbreviations.end(), field.first) != variant_abbreviations.end()) {
                    tmp.addAttr(field.first, field.second);
                }
            }
            tmp.addAttr( "gene_name_correct", m["gene_name_correct"] );
            tmp.addPz( pz_name );
            vVarClass.push_back( tmp );
            vVarClassNames.push_back( var_name );
        }
        // se VAR gia in classe ci aggiungo il paziente
        else {
            std::ptrdiff_t posVar = find( vVarClassNames.begin(), vVarClassNames.end(), var_name) - vVarClassNames.begin();
            if ( vVarClass[ posVar ].getValue("gene_name_correct") == "NA" ) {
                vVarClass[ posVar ].addAttr("gene_name_correct", m["gene_name_correct"]);
            }
            vVarClass[ posVar ].addPz( pz_name );
        }

        /*
            DEVO anche capire se PZ nuovo e aggiungerlo alla Classe PAtient con le sue caratteristiche
            e idem per pzVars e idem per Gene
        */
        // se PZ NON ancora nella Classe Pazienti
        pzVars tmpPzVar;
        for ( std::string var_char : patient_var_abbreviations ) {
            // se la ABBR PzVar nella map
            if (m.find( var_char ) != m.end()) {
                tmpPzVar.addAttr( var_char, m[var_char] );
            }
        }
        // e aggiungo anche il nome CORRETTO Del gene
        // NB questo punto e quello che fa casino se ho COMPHET TRIO FILE o meno !!!
        tmpPzVar.addAttr( "gene_name_correct", m["gene_name_correct"] );
        if (std::find (vPatientClassNames.begin(), vPatientClassNames.end(), pz_name) == vPatientClassNames.end()) {
            Patient tmpPz;
            for ( auto var_char : patient_abbreviations ) {
                // se la ABBR Patient nella map
                if (m.find( var_char ) != m.end()) {
                    tmpPz.addAttr( var_char, m[var_char] );
                }
            } // chiudo loop su ABBR della classe paziente
            tmpPz.addVar( var_name );
            tmpPz.addPzVarVector( tmpPzVar );
            vPatientClass.push_back( tmpPz );
            vPatientClassNames.push_back( pz_name );
        } // chiudo IF Patient NEW
        // se PATIENT NON-NEW devo aggiungergli Variant e basta
        else {
            // se VAR MAI inserita in quel Patient la inserisce
            std::ptrdiff_t posPz = std::find( vPatientClassNames.begin(), vPatientClassNames.end(), pz_name ) - vPatientClassNames.begin();
            std::vector<std::string> vPzVars = vPatientClass[ posPz ].getVars();
            if ( std::find ( vPzVars.begin(), vPzVars.end(), var_name ) == vPzVars.end() ) {
                vPatientClass[ posPz ].addVar( var_name );
                vPatientClass[ posPz ].addPzVarVector( tmpPzVar );
            }
        } // chiudo ELSE (patient NON-NEW)

        if (std::find (vGeneClassNames.begin(), vGeneClassNames.end(), gene_name) == vGeneClassNames.end()) {
            Gene tmpGene;
            tmpGene.addVar( var_name );
            vGeneClass.push_back( tmpGene );
            vGeneClassNames.push_back( gene_name );
        } // chiudo IF Patient NEW
        // se PATIENT NON-NEW devo aggiungergli Variant e basta
        else {
            // se VAR MAI inserita in quel Patient la inserisce
            std::ptrdiff_t posGene = std::find( vGeneClassNames.begin(), vGeneClassNames.end(), gene_name ) - vGeneClassNames.begin();
            std::vector<std::string> vGeneVars = vGeneClass[ posGene ].getVars();
            if ( std::find ( vGeneVars.begin(), vGeneVars.end(), var_name ) == vGeneVars.end() ) {
                vGeneClass[ posGene ].addVar( var_name );
            }
        } // chiudo ELSE (patient NON-NEW)

    }

    for ( auto m : vmVar2 )
    {
        std::string var_name = m["var_name"];
        std::string pz_name = m["pz_name"];
        std::string gene_name = m["gene_name_correct"];

        if (std::find (vGeneClassNames.begin(), vGeneClassNames.end(), gene_name) == vGeneClassNames.end()) {
            continue;
        }

        // se NON gia inserita
        if (std::find (vVarClassNames.begin(), vVarClassNames.end(), var_name) == vVarClassNames.end()) {
            Variant tmp;
            // apro loop su MAP della var
            for ( auto const &field : m )
            {
                if (std::find (variant_abbreviations.begin(), variant_abbreviations.end(), field.first) != variant_abbreviations.end()) {
                    tmp.addAttr(field.first, field.second);
                }
            }
            tmp.addPz( pz_name );
            tmp.addAttr("gene_name_correct", m["gene_name_correct"]);
            vVarClass.push_back( tmp );
            vVarClassNames.push_back( var_name );
        }
        // se VAR gia in classe ci aggiungo il paziente
        else {
            std::ptrdiff_t posVar = find( vVarClassNames.begin(), vVarClassNames.end(), var_name) - vVarClassNames.begin();
            if ( vVarClass[ posVar ].getValue("gene_name_correct") == "NA" ) {
                vVarClass[ posVar ].addAttr("gene_name_correct", m["gene_name_correct"]);
            }
            vVarClass[ posVar ].addPz( pz_name );
        }

        /*
            DEVO anche capire se PZ nuovo e aggiungerlo alla Classe PAtient con le sue caratteristiche
            e idem per pzVars e idem per Gene
        */
        // se PZ NON ancora nella Classe Pazienti
        pzVars tmpPzVar;
        for ( std::string var_char : patient_var_abbreviations ) {
            // se la ABBR PzVar nella map
            if (m.find( var_char ) != m.end()) {
                tmpPzVar.addAttr( var_char, m[var_char] );
            }
        }
        // e aggiungo anche il nome CORRETTO Del gene
        // NB questo punto e quello che fa casino se ho COMPHET TRIO FILE o meno !!!
        tmpPzVar.addAttr( "gene_name_correct", m["gene_name_correct"] );
        if (std::find (vPatientClassNames.begin(), vPatientClassNames.end(), pz_name) == vPatientClassNames.end()) {
            Patient tmpPz;
            for ( auto var_char : patient_abbreviations ) {
                // se la ABBR Patient nella map
                if (m.find( var_char ) != m.end()) {
                    tmpPz.addAttr( var_char, m[var_char] );
                }
            } // chiudo loop su ABBR della classe paziente
            tmpPz.addVar( var_name );
            tmpPz.addPzVarVector( tmpPzVar );
            vPatientClass.push_back( tmpPz );
            vPatientClassNames.push_back( pz_name );
        } // chiudo IF Patient NEW
        // se PATIENT NON-NEW devo aggiungergli Variant e basta
        else {
            // se VAR MAI inserita in quel Patient la inserisce
            std::ptrdiff_t posPz = std::find( vPatientClassNames.begin(), vPatientClassNames.end(), pz_name ) - vPatientClassNames.begin();
            std::vector<std::string> vPzVars = vPatientClass[ posPz ].getVars();
            if ( std::find ( vPzVars.begin(), vPzVars.end(), var_name ) == vPzVars.end() ) {
                vPatientClass[ posPz ].addVar( var_name );
                vPatientClass[ posPz ].addPzVarVector( tmpPzVar );
            }
        } // chiudo ELSE (patient NON-NEW)

        if (std::find (vGeneClassNames.begin(), vGeneClassNames.end(), gene_name) == vGeneClassNames.end()) {
            Gene tmpGene;
            tmpGene.addVar( var_name );
            vGeneClass.push_back( tmpGene );
            vGeneClassNames.push_back( gene_name );
        } // chiudo IF Patient NEW
        // se PATIENT NON-NEW devo aggiungergli Variant e basta
        else {
            // se VAR MAI inserita in quel Patient la inserisce
            std::ptrdiff_t posGene = std::find( vGeneClassNames.begin(), vGeneClassNames.end(), gene_name ) - vGeneClassNames.begin();
            std::vector<std::string> vGeneVars = vGeneClass[ posGene ].getVars();
            if ( std::find ( vGeneVars.begin(), vGeneVars.end(), var_name ) == vGeneVars.end() ) {
                vGeneClass[ posGene ].addVar( var_name );
            }
        } // chiudo ELSE (patient NON-NEW)

    }




    /*
        a questo punto devo andare a creare una MAP< std::string, std::vector<std::string>>
        per ricordarmi quale paziente ha quali VAR in COMPHET e stessa cosa per il GENE
    */
    //std::map < std::string, std::map < std::string, std::vector<std::string >>> mPatientsComphetVars;

    int nVar = 0;
    for ( auto m : vmVar12 )
    {
        std::vector<std::string> vVars = { m["var1_name"], m["var2_name"] };
        std::string pzName = m["pz_name"];
        std::string geneName = m["gene_name_correct"];

        // se PAZIENTE NON ancora nella MAP per i pazienti
        if ( mPatientsComphetVars.find(pzName) == mPatientsComphetVars.end() ) {
            // creo la sotto-MAP del GENE da inserire poi nel PAZIENTE
            std::map < std::string, std::vector<std::string >> mGeneComphetVars;
            if ( !mGeneComphetVars.insert( std::make_pair( geneName, vVars ) ).second ) {
                std::cout << "!!! ATTENTION !!! Field already present when trying to insert in PZ_NAME sub-MAP " << pzName << " GENE: " << geneName << " in PATIENTS COMPHET MAP" << '\n';
            }
            if ( !mPatientsComphetVars.insert( std::make_pair( pzName, mGeneComphetVars ) ).second ) {
                std::cout << "!!! ATTENTION !!! Field already present when trying to insert PZ_NAME " << pzName << " GENE MAP: " << geneName << " in PATIENTS COMPHET MAP" << '\n';
            }
            nVar += 2;
        }
        // se invece PAZIENTE gia presente
        else {
            // se GENE non ancora presente fa la MAP e la inserisce
            if ( mPatientsComphetVars[pzName].find(geneName) == mPatientsComphetVars[pzName].end() ) {
                // creo la sotto-MAP del GENE da inserire poi nel PAZIENTE
                if ( !mPatientsComphetVars[pzName].insert( std::make_pair( geneName, vVars ) ).second ) {
                    std::cout << "!!! ATTENTION !!! Field already present when trying to insert in PZ_NAME sub-MAP " << pzName << " GENE: " << geneName << " in PATIENTS COMPHET MAP" << '\n';
                }
                nVar += 2;
            }
            // se invece GENE gia presente per il paziente devo solo appenddere le VAR se non ci sono gia
            else {
                if (std::find (mPatientsComphetVars[ pzName ][ geneName ].begin(), mPatientsComphetVars[ pzName ][ geneName ].end(), m["var1_name"]) == mPatientsComphetVars[ pzName ][ geneName ].end()) {
                    mPatientsComphetVars[ pzName ][ geneName ].push_back( m["var1_name"] );
                }
                nVar ++;
                if (std::find (mPatientsComphetVars[ pzName ][ geneName ].begin(), mPatientsComphetVars[ pzName ][ geneName ].end(), m["var2_name"]) == mPatientsComphetVars[ pzName ][ geneName ].end()) {
                    mPatientsComphetVars[ pzName ][ geneName ].push_back( m["var2_name"] );
                }
                nVar ++;
            }
        }
    }

    std::cout << "----> in file TRIO COMPHET found: " << nVar << " inserted" << '\n';
}
