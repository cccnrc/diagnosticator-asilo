#include "basicFunctions.h"

/*
  questo file contiene invece le funzioni basilari di scrittura lettura ecc che mi servono per analisi
*/




void remove_duplicates_from_V (std::vector<std::string> &v)
{
    std::sort(v.begin(), v.end());
    auto it = std::unique(v.begin(), v.end());
    v.erase(it, v.end());
}



std::string stringize_float( float f )
{
    std::ostringstream ss;
    ss << f;
    std::string s(ss.str());
    return s;
}

std::vector<std::string> convert_map_keys_to_V ( std::map<std::string, std::vector<std::string>> m )
{
    std::vector<std::string> v;
    for (auto const &x : m)
    {
        v.push_back( x.first );
    }
    return v;
}


// questa invece scrive le linee VV corrispondente alla KEY della MAP fornita
// in una serie di file che sono nominati <nome-file><separator><KEY><suffix>
void write_map_VV
(
    std::map<std::string, std::vector<std::vector<std::string>>> mVV,
    std::string nomeFile,
    std::string separator,
    std::string suffix
)
{
    // apro quindi loop sulla map<std::string, VV>
    for ( auto const &chr : mVV )
    {
        std::string chr_name = chr.first;
        // creo il nome file definitivo
        std::string chr_filename = append_str_to_str( append_str_to_str( append_str_to_str( nomeFile, separator ), chr_name ), suffix );
        // estraggo la VV delle linee del suddetto file
        std::vector<std::vector<std::string>> vvChr = chr.second;
        // e lo scrivo
        write_cin_file_vv( chr_filename, vvChr );
    }
}


//definisco funzione per appender stringa ad altra stringa
std::string append_str_to_str ( std::string a, std::string b )
{
  std::stringstream streamVarName;
  streamVarName << a << b;
  return streamVarName.str();
}


// funzione per creare la DIR NB: se DIR esiste già va in ERROR
// se NON riesce a creare la directory result esce dal programma con errore
bool create_dir( std::string nomeDirResults )
{
    if (mkdir( nomeDirResults.c_str(), 0777 ) == -1 ) {
        std::cerr << "Error :  " << strerror(errno) << std::endl;
        return false;
    }
    else {
        std::cout << "------> " << nomeDirResults << " directory created" << '\n';
    }
    return true;
}

//definisco funzione per aggungere davanti directory a nome dei file risultati
std::string add_results_dir_path( std::string nameDir, std::string nameFile )
{
  std::stringstream streamVarName;
  streamVarName << nameDir << "/" << nameFile;
  return streamVarName.str();
}


// trasforma vettore in stringa con delimiter specificato
std::string stringize_vector (const std::vector<std::string> & vec, std::string delimiter) {
    return std::accumulate(std::next(vec.begin()), vec.end(),
        vec[0],
        [&delimiter](std::string& a, std::string b) {
        return a + delimiter + b;
    });
}

// questa stringa la std::map<std::string, std::string> con i due delimiter specificati (inter e intra column)
std::string stringize_map_str_str( std::map<std::string, std::string> m, std::string sep1, std::string sep2 )
{
    std::stringstream stream;
    for ( auto const &x : m )
    {
        std::string k = x.first;
        std::string v = x.second;
        replaceAll(v, ";", "_");
        stream << k << sep1 << v << sep2;
    }
    return stream.str();
}

// questa stringa la std::map<std::string, std::vector<std::string>> con i tre delimiter specificati (inter-objects, k:v e intra-v)
std::string stringize_map_str_V( std::map<std::string, std::vector<std::string>> m, std::string sep1, std::string sep2, std::string sep3 )
{
    std::stringstream stream;
    int xCount = 0;
    for ( auto const &x : m )
    {
        std::string k = x.first;
        // inserisce K - sep2
        stream << k << sep2;
        int yCount = 0;
        for ( std::string y : x.second )
        {
            // inserisce quindi ogni singolo valore del std::vector<std::string>
            replaceAll(y, ";", "_");
            stream << y;
            // se NON ultimo mette il sep3 (intra-vector)
            if (yCount < x.second.size()-1 ) {
                stream << sep3;
            }
            // se invece ultimo (e NON ulitmo elemento della map) mette il sep1 (inter-elements)
            else {
                if (xCount < m.size() - 1) {
                    stream << sep1;
                }
            }
            yCount++;
        } // chiufo LOOP su elementi del V value
        xCount++;
    } // chiudo LOOP su elementi della map
    return stream.str();
}


// questa stringa la std::map< std::map<std::string, std::vector<std::string>>> con i 5 delimiter specificati
// in questo modo: KEY1_1 (sep4) - KEY2_1 (sep2) - values in vector (sep3), (sep1) KEY2_2 ..., (sep5) KEY1_2
std::string stringize_map_map_str_V
(
    std::map< std::string, std::map<std::string, std::vector<std::string>>> mm,
    std::string sep1,
    std::string sep2,
    std::string sep3,
    std::string sep4,
    std::string sep5

)
{
    std::stringstream stream;
    int mCount = 0;
    for ( auto const &m : mm )
    {
        std::string z = m.first;
        stream << z << sep4;
        int xCount = 0;
        for ( auto const &x : m.second )
        {
            std::string k = x.first;
            // inserisce K - sep2
            stream << k << sep2;
            int yCount = 0;
            for ( std::string y : x.second )
            {
                // inserisce quindi ogni singolo valore del std::vector<std::string>
                replaceAll(y, ";", "_");
                stream << y;
                // se NON ultimo mette il sep3 (intra-vector)
                if (yCount < x.second.size()-1 ) {
                    stream << sep3;
                }
                // se invece ultimo (e NON ulitmo elemento della map) mette il sep1 (inter-elements)
                else {
                    if (xCount < m.second.size() - 1) {
                        stream << sep1;
                    }
                    else {
                        if (mCount < mm.size() - 1) {
                            stream << sep5;
                        }
                    }
                }
                yCount++;
            } // chiufo LOOP su elementi del V value
            xCount++;
        } // chiudo LOOP su elementi della map INT
        mCount++;
    } // chiudo LOOP su elementi della map EXT
    return stream.str();
}


// questa mi ritorna un TRUE se la std::string è convertibile in FLOAT
bool floatable( std::string &f_string ) {
  try {
    std::stof( f_string );
    return true;
  } catch(...) { return false; }
}

// funzione per comparare FLOAT (serve un indicatore di tolleranza se no il == 0 può dare risultati strani (https://noobtuts.com/cpp/compare-float-values))
bool cmpf(float A, float B, float epsilon) {
    return (fabs(A - B) < epsilon);
}

// funzione per comparare FLOAT (serve un indicatore di tolleranza se no il == 0 può dare risultati strani (https://noobtuts.com/cpp/compare-float-values))
bool cmpf_greater(float A, float B, float epsilon) {
    return (fabs(A - B) > epsilon);
}


// aggiunge a vettore se elemento non gia presente e mi restituisce 1 se new e 0 se gia presente
int add_absent_to_vector( std::vector<std::string> &v, std::string &value ) {
  int new_value = 0;
  if (std::find (v.begin(), v.end(), value) == v.end()) { v.push_back(value); new_value = 1; }
  return new_value;
}

//vettorializzo la colonna dalla stringa
void vectorize_column( std::string sCol, std::vector<std::string> &v, char separator ) {
  std::istringstream ss(sCol); // string stream delle posizioni
  while (getline ( ss, sCol, separator) ) v.push_back(sCol); //riempo array con la pos convertita a int
}

//vettorializzo la colonna dalla stringa
std::vector<std::string> vectorize_column_comma_quote( std::string sCol, char separator ) {
    std::vector<std::string> csvColumn;
    const char *mystart=sCol.c_str();    // prepare to parse the line - start is position of begin of field
    bool instring{false};
    for (const char* p=mystart; *p; p++) {  // iterate through the string
        if (*p=='"')                        // toggle flag if we're btw double quote
            instring = !instring;
        else if (*p==',' && !instring) {    // if comma OUTSIDE double quote
            csvColumn.push_back(std::string(mystart,p-mystart));  // keep the field
            mystart=p+1;                    // and start parsing next one
        }
    }
    csvColumn.push_back(std::string(mystart));   // last field delimited by end of line instead of comma
    return(csvColumn);
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

    //std::ofstream out(nomeFile.c_str(), std::fstream::app);
    std::ofstream out(nomeFile.c_str());
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


/*--------------- READ HEADER -------------------*/
void extract_header(
  std::vector<std::string> &vString
) {

    std::string line;
    int lineCount = 0;
    // apro loop sulle righe
    while(std::getline(std::cin, line)) {

      if (lineCount > 0) break;
      //if (line[0] == '#') continue;
      char fSep = ',';
      std::vector<std::string> vRiga;
      vectorize_column(line, vRiga, fSep);

      for (auto vi : vRiga) vString.push_back(vi);
      lineCount++;

    //chiudo LOOP sulle linee
    }
}

//questa cambia STD::CIN con il risultato della lettura del file
//poi chiama funzione su quella
std::vector<std::string> getHeader( std::string &nomeFile ) {

    std::vector<std::string> vString;

    std::ifstream in(nomeFile.c_str());
    std::streambuf *cinbuf = std::cin.rdbuf(); //save old buf
    std::cin.rdbuf(in.rdbuf()); //redirect std::cin to in.txt!

    extract_header( vString ); //call function

    std::cin.rdbuf(cinbuf);   //reset to standard input again

    return vString;
}
/*-----------------------------------------------------*/


/*--------------- convert FILE in std::vector<std::string> -------------------*/
void extract_vector(
  std::vector<std::string> &vString
) {

    std::string line;
    // apro loop sulle righe
    while(std::getline(std::cin, line)) {
      vString.push_back(line);
    //chiudo LOOP sulle linee
    }
}

//questa cambia STD::CIN con il risultato della lettura del file
//poi chiama funzione su quella
std::vector<std::string> getVectorFile( std::string &nomeFile ) {

    std::vector<std::string> vString;

    std::ifstream in(nomeFile.c_str());
    std::streambuf *cinbuf = std::cin.rdbuf(); //save old buf
    std::cin.rdbuf(in.rdbuf()); //redirect std::cin to in.txt!

    extract_vector( vString ); //call function

    std::cin.rdbuf(cinbuf);   //reset to standard input again

    return vString;
}
/*-----------------------------------------------------*/


/*--------------- convert FILE in std::vector<std::vector<std::vector<std::string>>> -------------------*/
// funzione per caricare il file come VVV specificando i due separatori da usare: es. inter-column (\t) e intra-column (,)
void extract_vector( std::vector<std::vector<std::vector<std::string>>> &vvvString, char sep1, char sep2 )
{
    std::string line;
    // apro loop sulle righe
    while(std::getline(std::cin, line)) {

        std::vector<std::vector<std::string>> vvLine;

        std::vector<std::string> vRiga;
        vectorize_column(line, vRiga, sep1);

        // apro loop in ogni colonna del file e trasformo i valori separati da virgola in vettori
        for (std::string v : vRiga)
        {
            std::vector<std::string> vColonna;
            vectorize_column(v, vColonna, sep2);
            vvLine.push_back(vColonna);
        }

        vvvString.push_back(vvLine);
        vvLine.clear();
    //chiudo LOOP sulle linee
    }
}

  //questa cambia STD::CIN con il risultato della lettura del file
  //poi chiama funzione su quella
std::vector<std::vector<std::vector<std::string>>> get_vvFile( std::string &nomeFile, char sep1, char sep2 ) {

    std::vector<std::vector<std::vector<std::string>>> vvvString;

    std::ifstream in(nomeFile.c_str());
    std::streambuf *cinbuf = std::cin.rdbuf(); //save old buf
    std::cin.rdbuf(in.rdbuf()); //redirect std::cin to in.txt!

    extract_vector( vvvString, sep1, sep2 ); //call function

    std::cin.rdbuf(cinbuf);   //reset to standard input again

    return vvvString;
}
/*-----------------------------------------------------*/

// questa prende un file messo in VVV e ne estrae la colonna desiderata come: std::vector<std::string>
std::vector<std::string> take_vvv_column( std::vector<std::vector<std::vector<std::string>>> &vvv, int numCol )
{
    // riporto il num colonna a posizione vettoriale
    numCol = numCol - 1;
    std::vector<std::string> v;
    for ( auto x : vvv )
    {
        v.push_back( x[numCol][0] );
    }
    return v;
}

std::map<std::string, std::vector<std::string>> create_vvv_map( std::vector<std::vector<std::vector<std::string>>> &vvv, int col_key, int col_values )
{
    std::map<std::string, std::vector<std::string>> m;
    col_key = col_key - 1;
    col_values = col_values - 1;
    for ( auto x : vvv )
    {
        if ( !m.insert( std::make_pair( x[col_key][0], x[col_values] ) ).second ) {
            std::cout << " !!! ERROR !!! key: " << x[col_key][0] << " repeated from original FILE (while converting in std::map)" << '\n';
        }
    }
    return m;
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

// questa rimpiazza TUTTE le occurences
void ReplaceStringInPlace(std::string& subject, std::string search, std::string replace)
{
    size_t pos = 0;
    while ((pos = subject.find(search, pos)) != std::string::npos) {
         subject.replace(pos, search.length(), replace);
         pos += replace.length();
    }
}

// NB questa rimpiazza SOLO la PRIMA occurence
bool replace_substring (std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}



//vettorializzo la colonna dalla stringa ---- v2
std::vector<std::string> vectorize_column_v2( std::string sCol, char separator ) {
    std::vector<std::string> v;
  std::istringstream ss(sCol); // string stream delle posizioni
  while (getline ( ss, sCol, separator) ) v.push_back(sCol); //riempo array con la pos convertita a int
  if ( v.back().empty() ) {
      v.pop_back();
  }
  return v;
}



/*--------------- convert FILE in std::vector<std::vector<std::string>> -------------------*/
void extract_VV(
  std::vector<std::vector<std::string>> &vvString,
  char sep
) {

    std::string line;
    // apro loop sulle righe
    while(std::getline(std::cin, line)) {
        std::vector<std::string> vRiga = vectorize_column_v2( line, sep );
        vvString.push_back( vRiga );
    //chiudo LOOP sulle linee
    }
}

//questa cambia STD::CIN con il risultato della lettura del file
//poi chiama funzione su quella
std::vector<std::vector<std::string>> getVVectorFile( std::string nomeFile, char sep ) {

    std::vector<std::vector<std::string>> vvString;

    std::ifstream in(nomeFile.c_str());
    std::streambuf *cinbuf = std::cin.rdbuf(); //save old buf
    std::cin.rdbuf(in.rdbuf()); //redirect std::cin to in.txt!

    extract_VV( vvString, sep ); //call function

    std::cin.rdbuf(cinbuf);   //reset to standard input again

    return vvString;
}
/*-----------------------------------------------------*/




// Returns true if x is in range [low..high], else false
bool inRange(unsigned low, unsigned high, unsigned x)
{
    return  ((x-low) <= (high-low));
}
