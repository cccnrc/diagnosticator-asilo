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


#include <dirent.h>
#include <sys/stat.h>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
/*
    questo sara lo script per comparare il file VAR del user con
    il central DB delle VAR (nella cartella specificata)
    e aggiungere a quest ultimo le mancanti
    con in INFO il campo user=<username1:proj1,username2:proj2,...,usernameN:projI>
    se il file CORE per il chromosoma in questione NON esiste allora lo CREA con le var del file user
    ha inoltre bisogno di cartella "core_changed" in cui inserire (APPEND) le var cambiate e added
    per ogni chromosma con le caratt e utente:proj
*/
/*
gcc -O3 -lstdc++ --std=c++11 -g -lm  -L/opt/X11/lib -lboost_iostreams -lz \
    core.db.0.1.cpp\
    -o core.db.0.1 \
    && time ./core.db.0.1 \
        results/var_results/ \
        core_db \
        core_changed \
        enrico0 \
        skz0 &>core.out
*/
//vettorializzo la colonna dalla stringa
std::vector<std::string> vectorize_column( std::string sCol, char separator )
{
    std::vector<std::string> v;
    std::istringstream ss(sCol); // string stream delle posizioni
    while (getline ( ss, sCol, separator) ) v.push_back(sCol); //riempo array con la pos convertita a int
    return v;
}

/*--------------- convert GZ FILE in std::vector<std::vector<std::string>> -------------------*/
std::vector<std::vector<std::string>> get_VV_GZ_file (std::string fileName)
{
    std::vector<std::vector<std::string>> vv;

    std::ifstream file(fileName, std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf; //iniziallizzo filtering_streambuf inbuf
    inbuf.push(boost::iostreams::gzip_decompressor()); //ci metto dentro decompressore GZIP (se file GZIP)
    inbuf.push(file); //ci metto dentro file
    //boost::iostreams::copy(inbuf, std::cout); //copio in stdout

    //Convert streambuf to istream
    std::istream instream(&inbuf);
    //Iterate lines
    std::string line;
    while(std::getline(instream, line)) {
        std::vector<std::string> vLine = vectorize_column( line, '\t' );
        vv.push_back(vLine);
    }
    return vv;
}
/* ------------------------------------------------------------------------------------ */

// recreate VAR name
std::string recreate_var_name( std::string chr, std::string pos, std::string ref, std::string alt )
{
    std::stringstream streamVarName;
    streamVarName << chr << "-" << pos << "-" << ref << "-" << alt;
    return streamVarName.str();
}

std::vector<std::string> extract_var_from_VV_core( std::vector<std::vector<std::string>> vv )
{
    std::vector<std::string> v;
    for ( std::vector<std::string> line : vv )
    {
        std::string varName = recreate_var_name( line[0], line[1], line[3], line[4] );
        if (std::find (v.begin(), v.end(), varName) == v.end()) {
            v.push_back( varName );
        }
    }
    return v;
}

// add the DIR path to a file name (remove / if last character in DIR name)
//definisco funzione per aggungere davanti directory a nome dei file risultati
std::string recreate_chr_core_filename( std::string prefix, std::string chr, std::string suffix )
{
    std::stringstream streamVarName;
    streamVarName << prefix << chr << suffix;
    return streamVarName.str();
}

// CHECK file existence and accessibility
inline bool exists_test (const std::string& name) {
    struct stat buffer;
    return (stat (name.c_str(), &buffer) == 0);
}

// add the DIR path to a file name (remove / if last character in DIR name)
//definisco funzione per aggungere davanti directory a nome dei file risultati
std::string add_dir_path( std::string nameDir, std::string nameFile )
{
    if ( nameDir.back() == '/' ) {
        nameDir = nameDir.substr(0, nameDir.size()-1);
    }
    std::stringstream streamVarName;
    streamVarName << nameDir << "/" << nameFile;
    return streamVarName.str();
}

std::string create_nome_file_changed_chromosome( std::string chromosome )
{
    std::stringstream streamVarName;
    streamVarName << "changed." << chromosome << ".txt";
    return streamVarName.str();
}

std::string create_nome_file_added_chromosome( std::string chromosome )
{
    std::stringstream streamVarName;
    streamVarName << "added." << chromosome << ".txt";
    return streamVarName.str();
}

std::string create_userProj_name( std::string user, std::string proj )
{
    std::stringstream streamVarName;
    streamVarName << user << ":" << proj;
    return streamVarName.str();
}


bool get_files_dir ( std::string nomeDir, std::vector<std::string> &v, std::string source )
{
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir (nomeDir.c_str())) != NULL) {
        std::cout << "" << '\n';
        std::cout << "-------------------- starting " << source <<  " ---------------------" << '\n';
        /* print all the files and directories within directory */
        while ((ent = readdir(dir)) != NULL) {
            // converto l array doi CHAR in string
            std::string nameFile(ent->d_name);
            // per togliere i file accessori "."
            if ( nameFile.length() > 3 ) {
                v.push_back(nameFile);
                std::cout << " ----> " << source << ": reading file " << nameFile << " from " << nomeDir << '\n';
            }
      }
        closedir (dir);
        std::cout << "-------------------- " << source <<  " finished ---------------------" << '\n';
        std::cout << "" << '\n';
        return true;
    } else {
        /* could not open directory */
        perror ("");
    }
    return false;
}




/*--------------- convert FILE in std::vector<std::vector<std::string>> -------------------*/
void extract_VV(
  std::vector<std::vector<std::string>> &vvString
)
{

    std::string line;
    // apro loop sulle righe
    while(std::getline(std::cin, line))
    {
        std::vector<std::string> vLine = vectorize_column( line, '\t' );
        vvString.push_back(vLine);
    //chiudo LOOP sulle linee
    }
}

//questa cambia STD::CIN con il risultato della lettura del file
//poi chiama funzione su quella
std::vector<std::vector<std::string>> get_file_user( std::string nomeFile ) {

    std::vector<std::vector<std::string>> vvString;

    std::ifstream in(nomeFile.c_str());
    std::streambuf *cinbuf = std::cin.rdbuf(); //save old buf
    std::cin.rdbuf(in.rdbuf()); //redirect std::cin to in.txt!

    extract_VV( vvString ); //call function

    std::cin.rdbuf(cinbuf);   //reset to standard input again

    return vvString;
}
/*-----------------------------------------------------*/

//funzione per WRITE FILE di vec of vec
void write_file_vv_GZ( std::string nomeFile, std::vector< std::vector<std::string>> &vv ) {

    //std::ofstream out(nomeFile.c_str(), std::fstream::app);
    std::ofstream file(nomeFile.c_str(), std::ios_base::out | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::output> out;
    out.push(boost::iostreams::gzip_compressor());
    out.push(file);
    // apro lo stram per scriverci il file
    std::stringstream sstr;

    for (auto i : vv) {
      int iCount = 0;
      for (auto vi : i) {
        sstr << vi;
        if ( iCount < i.size()-1 ) {sstr << "\t";}
        iCount++;
      }
      sstr << '\n';
    }

    boost::iostreams::copy(sstr, out);
}

//funzione per WRITE FILE di vec of vec
void write_cin_file_vv( std::string &nomeFile, std::vector< std::vector<std::string>> &vvGnomadfounds ) {

    std::ofstream out(nomeFile.c_str(), std::fstream::app);
    //std::ofstream out(nomeFile.c_str());
    std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
    std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!

    for (auto i : vvGnomadfounds) {
      int iCount = 0;
      for (auto vi : i) {
        std::cout << vi;
        if ( iCount < i.size()-1 ) {std::cout << "\t";}
        iCount++;
      }
      std::cout << '\n';
    }

    std::cout.rdbuf(coutbuf); //reset to standard output again
}

//funzione per WRITE FILE di vec of vec
void write_cin_file_v( std::string &nomeFile, std::vector<std::string> &v ) {

    std::ofstream out(nomeFile.c_str(), std::fstream::app);
    //std::ofstream out(nomeFile.c_str());
    std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
    std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!

      for (std::string vi : v) {
        std::cout << vi << '\n';
      }

    std::cout.rdbuf(coutbuf); //reset to standard output again
}

// questa per creare un VV dalla std::map<std::string, std::map<std::string, std::map<std::string, std::string>>>
std::vector<std::vector<std::string>> vectorize_VV_MMM (
    std::map<std::string, std::map<std::string, std::map<std::string, std::string>>> &m,
    std::string userName
)
{

    // HEADER: "#var", "field", "from", "to", "user:proj"
    std::vector<std::vector<std::string>> vv;

    for ( auto const &i : m )
    {
        std::vector<std::string> v;
        for ( auto const &x : i.second )
        {
            v.push_back( i.first );
            v.push_back( x.first );
            for ( auto const &y : x.second )
            {
                v.push_back( y.first );
                v.push_back( y.second );
            }
            v.push_back( userName );
            vv.push_back(v);
            v.clear();
        }
    }
    return vv;
}


// trasforma vettore in stringa con delimiter specificato
std::string stringize_vector (const std::vector<std::string> & vec, std::string delimiter) {
    return std::accumulate(std::next(vec.begin()), vec.end(),
        vec[0],
        [&delimiter](std::string& a, std::string b) {
        return a + delimiter + b;
    });
}

// questo splitta in due la stringa SOLO sulla prima occurence del char separator
std::vector<std::string> split_on_first_ocurrence( std::string string, char sep )
{
    std::vector<std::string> v;
    if (string.find(sep) != std::string::npos) {
        const auto equals_idx = string.find_first_of(sep);
        if (std::string::npos != equals_idx) {
            v.push_back( string.substr(0, equals_idx) );
            v.push_back( string.substr(equals_idx + 1) );
        }
    }
    else {
        std::cout << "!!! ERROR !!! separator: " << sep << " absent in string: " << string << " when try to split it" << '\n';
    }
    return v;
}

// converto il campo INFO in una MAP <str, str>
std::map<std::string, std::string> mapize_info_field (std::string varUserInfo, std::string varName)
{
    std::vector<std::string> vFieldToExclude = {
        "het_cases",
        "hom_cases",
        "mother_dp_bin",
        "mother_gt",
        "father_dp_bin",
        "father_gt",
        "pz_all",
        "pz_comp_het",
        "pz_hom",
        "pz_single_het",
        "cases_af",
        "controls_af"
    };
    std::vector<std::string> vVarUserInfo = vectorize_column( varUserInfo, ';' );
    if ( vVarUserInfo.back().empty() ) {
        vVarUserInfo.pop_back();
    }
    std::map<std::string, std::string> mVarUserInfo;
    for (std::string field : vVarUserInfo )
    {
        std::vector<std::string> vField = split_on_first_ocurrence( field, '=' );
        // se in uno dei field da scartare lo ignora
        if (std::find (vFieldToExclude.begin(), vFieldToExclude.end(), vField[0]) != vFieldToExclude.end()) {
            continue;
        }
        if ( !mVarUserInfo.insert( std::make_pair( vField[0], vField[1] ) ).second ) {
            std::cout << " !!! ERROR !!! key: " << field[0] << " already present in USER MAP for var: " << varName << '\n';
        }
    }
    return mVarUserInfo;
}


void replaceAll(std::string& str, const std::string& from, const std::string& to) {
    if(from.empty())
        return;
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
    }
}


// questa stringa la std::map<std::string, std::string> con i due delimiter specificati (inter e intra column)
std::string stringize_map_str_str( std::map<std::string, std::string> m, std::string sep1, std::string sep2 )
{
    std::stringstream stream;
    int xCount = 0;
    for ( auto const &x : m )
    {
        std::string k = x.first;
        std::string v = x.second;
        replaceAll(v, ";", "_");
        stream << k << sep1 << v;
        if ( xCount < m.size()-1 ) {
            stream << sep2;
        }
        xCount++;
    }
    return stream.str();
}




int main(int argc, char const *argv[]) {

    // controllo argomenti passati allo script
    if ( argc == 6 ) {
      std::cout << "num arg correct" << '\n';
      std::cout << "> USER: " << argv[4] << '\n';
      std::cout << "> PROJECT: " << argv[5] << '\n';
      std::cout << "--> reading input user DBs from: " << argv[1] << '\n';
      std::cout << "--> reading input core DBs from: " << argv[2] << '\n';
      std::cout << "--> writing changes core DBs in: " << argv[3] << '\n';
    }
    else {
        std::cout << "!!! ERROR !!! Please call it with: ./core.db.0.1 <user-DB-dir> <core-DB-dir> <username>" << '\n';
        return 1;
    }

    std::string userDir = argv[1];
    std::string coreDir = argv[2];
    std::string coreChangeDir = argv[3];
    std::string userOnlyName = argv[4];
    std::string projName = argv[5];

    // gli USER saranno segnati come username:projname
    std::string userName = create_userProj_name( userOnlyName, projName );


    // funzione per CONTROLLARE che le DIR esistano e siano DIR
    // funzione per prendere TUTTI i filename nelle directory specificate
    std::vector<std::string> vFileUser;
    if ( !get_files_dir ( userDir, vFileUser, "user" ) ) {
        return 1;
    }
    std::vector<std::string> vFileCoreDB;
    if ( !get_files_dir ( coreDir, vFileCoreDB, "CORE" ) ) {
        return 1;
    }

    for ( std::string fileUser : vFileUser )
    {
        std::string pathFileUser = add_dir_path( userDir, fileUser);
        if ( exists_test (pathFileUser) ) {
            std::vector<std::vector<std::string>> vvFileUser = get_file_user(pathFileUser);
            std::string chromosome = vvFileUser[0][0];
            std::cout << "> CHROMOSOME: " << chromosome << '\n';
            std::cout << " --> " << vvFileUser.size() << " variants in USER file for chr: " << chromosome << " (" << pathFileUser << ")" << '\n';
            // ricostruisco nome del CHR corrispondente nel CORE DB
            std::string chrFileCore = recreate_chr_core_filename( "core.var.", chromosome, ".vcf.gz" );
            std::string pathFileCore = add_dir_path( coreDir, chrFileCore);

            // apro V per tenere conto di QUALI var ho aggiunto
            std::vector<std::string> vVarAdded;
            // apro MMM per tenere conto di:
            // key1: variante changed - value1: map sotto
            // key2: field changed - value2: map sotto
            // key3: from - value3: to
            std::map<std::string, std::map<std::string, std::map<std::string, std::string>>> mmmVarChanged;
            std::vector<std::vector<std::string>> vvVarChanged;

            // se esiste file CORE corrispondente
            if ( exists_test (pathFileCore) ) {
                std::vector<std::vector<std::string>> vvFileCore = get_VV_GZ_file( pathFileCore );
                std::vector<std::string> vCoreVar = extract_var_from_VV_core( vvFileCore );
                // devo aprire anche un V per tenere conto delle VAR trovate del file USER e NON nel file CORE per poter capire quali invece del file CORE
                // NON saranno trovate in file USER (visto che il loop e su USER e non su CORE mi serve per fare gli adattamenti anche a queste)
                // che all inizio sara uguale a quello delle VAR del file CORE al quale togliero una ad una quelle trovate, cosi da rimanere con SOLO quelle NON trovate in file USER
                std::vector<std::string> vNotFoundCore = vCoreVar;

                std::cout << " --> " << vvFileCore.size() << " variants in CORE file for chr: " << chromosome << " (" << pathFileCore << ")" << '\n';

                // COMPARA le varianti nei due VV (user vs. core)
                int numVarAdded = 0; int numVarChanged = 0; int numVarUnchanged = 0;
                for ( std::vector<std::string> v : vvFileUser )
                {
                    std::string varName = recreate_var_name( v[0], v[1], v[3], v[4] );

                    // apro V per tenere conto dei FIELD changed della VAR
                    // se NON ce ne sara nessuno NON la aggiungo al VV delle vvVarChanged
                    // se no invece la aggiungo
                    std::vector<std::string> vVarChanged;

                    // se NON esiste in VARIANTI del CORE aggiunge user attuale e aggiunge la linea a questo VV
                    if (std::find (vCoreVar.begin(), vCoreVar.end(), varName) == vCoreVar.end()) {
                        // prendo info della VAR USER
                        std::string varUserInfo = v.back();
                        std::map<std::string, std::string> mVarUserInfo = mapize_info_field( varUserInfo, varName );
                        if ( mVarUserInfo.find( "user" ) != mVarUserInfo.end() ) {
                            // vettorializza gli users e ci aggiunge user attuale se assente
                            std::vector<std::string> vUsers = vectorize_column( mVarUserInfo["user"], ',' );
                            // se NON gia presente tra gli user glielo aggiunge e lo mette come nuovo value del campo USER della map CORE
                            if (std::find (vUsers.begin(), vUsers.end(), userName) == vUsers.end()) {
                                vUsers.push_back( userName );
                                std::string usersStringized = stringize_vector( vUsers, "," );
                                mVarUserInfo["user"] = usersStringized;
                            }
                        } // chiudo IF field user gia nel file CORE
                        else {
                            if ( !mVarUserInfo.insert( std::make_pair( "user", userName ) ).second ) {
                                std::cout << "!!! ERROR !!! Trying to insert a USER " << userName << " in field already in USER map, for var: " << varName << '\n';
                            }
                        }
                        v.pop_back();
                        v.push_back( stringize_map_str_str( mVarUserInfo, "=", ";" ) );
                        // e inserisco quindi la map stringized come INFO della var che poi metto nel file CORE
                        vvFileCore.push_back(v);
                        vVarAdded.push_back(varName);
                        numVarAdded++;
                    } // chiudo IF var NON in file CORE

                    // se invece esiste gia la VAR nel file CORE deve comparare alcuni campi di INFO per vedere che non sia cambiato nulla
                    else {

                        // tolgo quindi la VAR dal vettore che serve a tenermi ricordate quali nel file CORE assenti nel file USER per cui devo fare le modifche
                        vNotFoundCore.erase(std::remove(vNotFoundCore.begin(), vNotFoundCore.end(), varName), vNotFoundCore.end());

                        // apro la MM per tenere conto delle KEY che cambiano (la SECONDA delle mmmVarChanged)
                        std::map<std::string, std::map<std::string, std::string>> mmChanged;

                        // prendo info della VAR USER
                        std::string varUserInfo = v.back();
                        // prendo info della VAR CORE
                        std::ptrdiff_t posVarCore = find( vCoreVar.begin(), vCoreVar.end(), varName) - vCoreVar.begin();
                        std::string varCoreInfo = vvFileCore[posVarCore].back();
                        // e tolgo il campo INFO dal VV CORE (perche lo sostituiro con quello modificato)
                        vvFileCore[posVarCore].pop_back();
                        // prendo la MAP<str,str> della var (TRANNE I CAMPI da ESCLUDERE)
                        std::map<std::string, std::string> mVarUserInfo = mapize_info_field( varUserInfo, varName );
                        std::map<std::string, std::string> mVarCoreInfo = mapize_info_field( varCoreInfo, varName );

                        // se nella map del file CORE per quella VAR NON si trova il FIELD user allora lo aggiunge con user attuale
                        if ( mVarCoreInfo.find( "user" ) == mVarCoreInfo.end() ) {
                            if ( !mVarCoreInfo.insert( std::make_pair( "user", userName ) ).second ) {
                                std::cout << "!!! ERROR !!! Trying to insert a USER " << userName << " in field already in USER map, for var: " << varName << '\n';
                            }
                        }
                        // se invece lo trova nel CORE MA NON nel USER per quella VAR (non lo aggiornerebbe quindi nel loop sotto)
                        else {
                            // agguiungo il suo valore alla map USER INFO
                            if ( mVarUserInfo.find( "user" ) == mVarUserInfo.end() ) {
                                if ( !mVarUserInfo.insert( std::make_pair( "user", userName ) ).second ) {
                                    std::cout << "!!! ERROR !!! Trying to insert a USER " << userName << " in field already in USER map, for var: " << varName << '\n';
                                }
                            }
                        } // chiudo IF field USER nel file CORE ma NON nel file USER

                        // apro quindi loop su map delle INFO del file USER
                        for (auto const &m : mVarUserInfo)
                        {

                            // se parliamo del campo USER
                            if (m.first == "user") {
                                // se gia presente il field nel file CORE
                                if ( mVarCoreInfo.find( "user" ) != mVarCoreInfo.end() ) {
                                    // vettorializza gli users e ci aggiunge user attuale se assente
                                    std::vector<std::string> vUsers = vectorize_column( mVarCoreInfo["user"], ',' );
                                    // se NON gia presente tra gli user glielo aggiunge e lo mette come nuovo value del campo USER della map CORE
                                    if (std::find (vUsers.begin(), vUsers.end(), userName) == vUsers.end()) {
                                        vUsers.push_back( userName );
                                        std::string usersStringized = stringize_vector( vUsers, "," );
                                        mVarCoreInfo["user"] = usersStringized;
                                    }
                                } // chiudo IF field user gia nel file CORE
                                // se invece field user assente inserisce il campo con value user attuale
                                else {
                                    if ( !mVarCoreInfo.insert( std::make_pair( "user", userName ) ).second ) {
                                        std::cout << "!!! ERROR !!! Trying to insert a USER " << userName << " in field already in CORE map, for var: " << varName << '\n';
                                    }
                                } // chiude ELSE field user assente in map CORE
                                continue; // perche poi deve skippare quanto sotto
                            } // chiude IF KEY user

                            // se parliamo di altro campo:
                            // se presente la KEY in CORE file
                            if ( mVarCoreInfo.find( m.first ) != mVarCoreInfo.end() ) {
                                // controlla che anche il VALUE sia identico, altrimenti lo cambia con il nuovo
                                if ( m.second != mVarCoreInfo[m.first] ) {
                                    // apro una map per ricordare il valore che cambio da PRE a POST (che sara la map piu interna della mmmVarChanged)
                                    std::map<std::string, std::string> mChanged;
                                    //std::cout << varName << "> key: " << m.first << ", user: " << m.second << ", CORE: " << mVarCoreInfo[m.first] << " changed" << '\n';
                                    if ( !mChanged.insert( std::make_pair( mVarCoreInfo[m.first], m.second ) ).second ) {
                                        std::cout << "!!! ERROR !!! Trying to insert a VALUE PRE " << mVarCoreInfo[m.first] << " already changed in CORE map CHANGED: " << m.first << " for var: " << varName << '\n';
                                    }
                                    mVarCoreInfo[m.first] = m.second;
                                    // e poi inserisco la map con il change in apposita sovra map: FIELD changed <pre, post>
                                    if ( !mmChanged.insert( std::make_pair( m.first, mChanged ) ).second ) {
                                        std::cout << "!!! ERROR !!! Trying to insert a key already in CORE map CHANGED: " << m.first << " for var: " << varName << '\n';
                                    }
                                }
                            } // chiudo IF KEY presente nel CORE file
                            // se invece la KEY manca nel file CORE la aggiungo alla sua MAP
                            else {
                                //std::cout << varName << "> key: " << m.first << ", value: " << m.second << ", ABSENT in CORE. Added." << '\n';
                                if ( !mVarCoreInfo.insert( std::make_pair( m.first, m.second ) ).second ) {
                                    std::cout << "!!! ERROR !!! Trying to insert a key already in CORE map: " << m.first << " for var: " << varName << '\n';
                                }
                            } // chiudo ELSE KEY di field absent nel file CORE
                        } // chiudo LOOP sui campi field fi INFO (della var USER) mapized

                        // infine, dopo aver fatto le modifiche, stirngizo la map delle INFO del file CORE e poi lo aggiungo al VV da scrivere
                        vvFileCore[posVarCore].push_back( stringize_map_str_str( mVarCoreInfo, "=", ";" ) );

                        // se per la VAR ha aggiunto almeno un campo nei CHANGED
                        // la aggiunge alla SOVRA MMM per ricordare a CHE VAR (KEY) corrispondono i cambi
                        if ( mmChanged.size() > 0 ) {
                            if ( !mmmVarChanged.insert( std::make_pair( varName, mmChanged ) ).second ) {
                                std::cout << "!!! ERROR !!! Trying to insert a VAR CHANGED in CORE map MMM: " << varName << '\n';
                            }
                            numVarChanged++;
                        }
                        // se NON ha cambiato nulla aumenta il conteggio delle unchanged
                        else {
                            numVarUnchanged++;
                        }
                    } // chiudo ELSE var gia esistente nel file CORE

                } // chiudo LOOP su righe del file USER

                // apro quindi LOOP sulle VAR non trovate di CORE in USER file per farci le modifiche necessarie
                int numAddedUsercore = 0; int numVarUnchangedCore = 0;
                for ( std::string varCore : vNotFoundCore )
                {
                    // riprendo POS e info della VAR CORE
                    std::ptrdiff_t posVarCore = find( vCoreVar.begin(), vCoreVar.end(), varCore) - vCoreVar.begin();
                    std::string varCoreInfo = vvFileCore[posVarCore].back();
                    std::map<std::string, std::string> mVarCoreInfo = mapize_info_field( varCoreInfo, varCore );
                    // se nella map del file CORE per quella VAR NON si trova il FIELD user allora lo aggiunge con user attuale
                    if ( mVarCoreInfo.find( "user" ) == mVarCoreInfo.end() ) {
                        if ( !mVarCoreInfo.insert( std::make_pair( "user", userName ) ).second ) {
                            std::cout << "!!! ERROR !!! Trying to insert a USER " << userName << " in field already in USER map, for var: " << varCore << '\n';
                        }
                        numAddedUsercore++;
                    }
                    // se INVECE presente lo aggiorna con user attuale
                    else {
                        // vettorializza gli users e ci aggiunge user attuale se assente
                        std::vector<std::string> vUsers = vectorize_column( mVarCoreInfo["user"], ',' );
                        // se NON gia presente tra gli user glielo aggiunge e lo mette come nuovo value del campo USER della map CORE
                        if (std::find (vUsers.begin(), vUsers.end(), userName) == vUsers.end()) {
                            vUsers.push_back( userName );
                            std::string usersStringized = stringize_vector( vUsers, "," );
                            mVarCoreInfo["user"] = usersStringized;
                            numAddedUsercore++;
                        }
                        // se invece user gia presente NON fa nulla e aumento il contatore
                        else {
                            numVarUnchangedCore++;
                        }
                    } // chiudo ELSE field USER gia presente in VAR del file CORE
                    vvFileCore[posVarCore].pop_back();
                    vvFileCore[posVarCore].push_back( stringize_map_str_str( mVarCoreInfo, "=", ";" ) );
                } // chiudo LOOP su VAR NOT FOUND di CORE in file USER

                std::cout << " -------> added " << numVarAdded <<  " variants to CORE file for chromosome: " << chromosome << " (" << pathFileCore << ") from USER file: " << pathFileUser << '\n';
                std::cout << " ---------------> changed " << numVarChanged <<  " variants to CORE file for chromosome: " << chromosome << " (" << pathFileCore << ") from USER file: " << pathFileUser << '\n';
                std::cout << " ---------------> CORE added user " << userName << ": " << numAddedUsercore <<  " variants of CORE file for chromosome: " << chromosome << " (" << pathFileCore << ") from USER file: " << pathFileUser << '\n';
                std::cout << " -----------------------> unchanged " << numVarUnchanged <<  " variants to CORE file for chromosome: " << chromosome << " (" << pathFileCore << ") from USER file: " << pathFileUser << '\n';
                std::cout << " -----------------------> CORE unchanged " << numVarUnchangedCore <<  " variants of CORE file for chromosome: " << chromosome << " (" << pathFileCore << ") from USER file: " << pathFileUser << '\n';
                std::cout << " --------------------------------> writing a TOTAL of: " << vvFileCore.size() <<  " variants to CORE file for chromosome: " << chromosome << " (" << pathFileCore << ")" << '\n';

                // controllo infine che i contatori siano corretti
                if ( numVarAdded + numVarChanged + numAddedUsercore + numVarUnchanged + numVarUnchangedCore != vvFileCore.size() ) {
                    std::cout << " ==================================> !!! ATTENTION !!! TOTAL NUMBERS NOT EQUAL to VAR writte in CORE updated file !!! CHECK IT !!! <================================== " << '\n';
                }

                // scrivo quindi il FILE CORE NUOVO del chromosoma
                write_file_vv_GZ( pathFileCore, vvFileCore );

                // e scrivo anche il file delle VAR AGGIUNTE e cambiate per quel chromosoma
                std::string pathFileChangedChromosome = add_dir_path( coreChangeDir, create_nome_file_changed_chromosome( chromosome ) );
                std::string pathFileAddedChromosome = add_dir_path( coreChangeDir, create_nome_file_added_chromosome( chromosome ) );
                std::vector<std::vector<std::string>> vMMM = vectorize_VV_MMM( mmmVarChanged, userName );
                write_cin_file_vv( pathFileChangedChromosome, vMMM );
                // metto a fianco della var added il nome del utenete e proj che la ha aggiunta
                std::vector<std::vector<std::string>> vvVarAdded;
                for ( std::string var : vVarAdded )
                {
                    std::vector<std::string> v;
                    v.push_back(var);
                    v.push_back(userName);
                    vvVarAdded.push_back(v);
                    v.clear();
                }
                write_cin_file_vv( pathFileAddedChromosome, vvVarAdded );
            } //chiudo IF EXISTS file CORE corrispondente del chromosome

            // se invece non esiste proprio il file core corrispondente per lo stesso chromosoma allora
            // deve sistemare le caratt della var user e inserirla come nuovo file core per il chromosoma
            else {
                std::cout << " -----> !!! ATTENTION !!! CORE DB file for chr: " << chromosome << " does NOT exists! Creating it!" << '\n';
                std::vector<std::vector<std::string>> vvFileCore;
                // COMPARA le varianti nei due VV (user vs. core)
                int numVarAdded = 0;
                for ( std::vector<std::string> v : vvFileUser )
                {
                    std::string varName = recreate_var_name( v[0], v[1], v[3], v[4] );
                    // se NON esiste in VARIANTI del CORE aggiunge user attuale e aggiunge la linea a questo VV
                    // prendo info della VAR USER
                    std::string varUserInfo = v.back();
                    std::map<std::string, std::string> mVarUserInfo = mapize_info_field( varUserInfo, varName );
                    if ( mVarUserInfo.find( "user" ) != mVarUserInfo.end() ) {
                        // vettorializza gli users e ci aggiunge user attuale se assente
                        std::vector<std::string> vUsers = vectorize_column( mVarUserInfo["user"], ',' );
                        // se NON gia presente tra gli user glielo aggiunge e lo mette come nuovo value del campo USER della map CORE
                        if (std::find (vUsers.begin(), vUsers.end(), userName) == vUsers.end()) {
                            vUsers.push_back( userName );
                            std::string usersStringized = stringize_vector( vUsers, "," );
                            mVarUserInfo["user"] = usersStringized;
                        }
                    } // chiudo IF field user gia nel file CORE
                    else {
                        if ( !mVarUserInfo.insert( std::make_pair( "user", userName ) ).second ) {
                            std::cout << "!!! ERROR !!! Trying to insert a USER " << userName << " in field already in USER map, for var: " << varName << '\n';
                        }
                    }
                    v.pop_back();
                    v.push_back( stringize_map_str_str( mVarUserInfo, "=", ";" ) );
                    // e inserisco quindi la map stringized come INFO della var che poi metto nel file CORE
                    vvFileCore.push_back(v);
                    vVarAdded.push_back(varName);
                    numVarAdded++;
                } // chiudo LOOP sulle righe del file del paziente del chromomosoma
                // scrivo quindi il FILE CORE NUOVO del chromosoma
                write_file_vv_GZ( pathFileCore, vvFileCore );
                // metto a fianco della var added il nome del utenete e proj che la ha aggiunta
                std::string pathFileAddedChromosome = add_dir_path( coreChangeDir, create_nome_file_added_chromosome( chromosome ) );
                std::vector<std::vector<std::string>> vvVarAdded;
                for ( std::string var : vVarAdded )
                {
                    std::vector<std::string> v;
                    v.push_back(var);
                    v.push_back(userName);
                    vvVarAdded.push_back(v);
                    v.clear();
                }
                write_cin_file_vv( pathFileAddedChromosome, vvVarAdded );
                std::cout << " -------> added " << numVarAdded <<  " variants to NEWLY CREATED CORE file for chromosome: " << chromosome << " (" << pathFileCore << ") from USER file: " << pathFileUser << '\n';
            } // chiudo ELSE file core corrispondente al chromosoma NON esistente
        } // chiudo IF EXISTS file del PAZIENTE
        std::cout << "" << '\n';
    } // chiudo loop su FILEs nella directory del PAZIENTE

    return 0;
}
