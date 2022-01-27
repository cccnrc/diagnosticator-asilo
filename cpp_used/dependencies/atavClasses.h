#ifndef ATAVCLASSES_H
#define ATAVCLASSES_H

#include <vector>
#include <string>
#include <map>
#include "basicFunctions.h"

class Variant {
  std::map<std::string, std::string> attributes;
  // un vettore per tenere riordato tutti i pazienti con la variante
  std::vector<std::string> vPatients;
  // un vettore per tenere le varie categorie in cui finisce di ACMG
  std::vector<std::string> vAcmg;
  // questo per le CLASSI invece
  std::vector<std::string> vAcmgClass;
  // e la map per tenere riassunto numerico delle subclassi
  std::map<std::string, int> mAcmgSubclasses;
  /* ---
   questi 2 vettori e la map sotto mi servono solo per geni AR
   infatti se la var sarà finita in qualsiasi categoria e il gene è AD i basta prendere il vettore vPatients con tutti i pazienti
   --- */
  // un vettore per sapere quali pazienti hanno la var in HOM e quali in SINGLE-HET (mi servira solo se gene AR)
  std::vector<std::string> vPzHom;
  std::vector<std::string> vPzSingleHet;
  // se GENE AR una map<str,vector<str>> per tenere ricordato quale paziente ha quali varianti in comphet (non mi basta vettore perchè devo sapere con quali ALTRE Var si trova in comphet) (mi servira solo se gene AR)
  std::map< std::string, std::map<std::string, std::vector<std::string>>> mPzComphet;

  public:
    /* --- 1: CHARACTERISTICS --- */
    void addAttr(std::string attr, std::string value) {
      if ( !attributes.insert( std::make_pair( attr, value ) ).second ) std::cout << "!!! ATTENTION !!! insert failed for: " << attr << ": " << value << " (attribute already present)" << '\n';
    }
    //void addAttr(std::string attr, std::string value){ attributes[attr] = value; }
    std::string getValue(std::string attr)
    {
        try {
          return attributes.at(attr);;
      } catch(...) { return "NA"; }
    }
    std::map<std::string, std::string> getAll() { return attributes; }
    void removeAttr( std::string key )
    {
        // se presente la KEY la elimina
        if (attributes.find(key) != attributes.end()) {
            attributes.erase(key);
        }
        else {
            std::cout << "!!! ATTENTION !!! attribute: " << key << " not found in VAR " << attributes.at("var_name") << " std::map attributes" << '\n';
        }
    }


    /* --- 2: ALL VARIANTS --- */
    // funzione per aggiungere i pazienti totali ( e di DEFAULT anche in SINGLE-HET
    // dal quale toglierò i pazienti dopo se HOM o COMPHET, lasciando alla fine solo i SINGLE-HET )
    void addPz( std::string pz_name ) {
      if (std::find (vPatients.begin(), vPatients.end(), pz_name) == vPatients.end()) { vPatients.push_back(pz_name); }
      if (std::find (vPzSingleHet.begin(), vPzSingleHet.end(), pz_name) == vPzSingleHet.end()) { vPzSingleHet.push_back(pz_name); }
    }
    // funzione per ritornare vettore di tutti i pazienti
    std::vector<std::string> & getPzAll() { return vPatients; };

    /* --- 3: SINGLE-HET --- */
    // funzione per ritornare vettore SINGLE-HET
    std::vector<std::string> & getPzSingleHet() { return vPzSingleHet; };
    // funzione per rimuovere i pazienti da vettore SINGLE-HET
    void removePzSingleHet( std::string pz_name ) {
      vPzSingleHet.erase(std::remove(vPzSingleHet.begin(), vPzSingleHet.end(), pz_name), vPzSingleHet.end());
    }
    // check rapido paziente in SINGLE_HET
    bool check_pzSingleHet ( std::string pz_name )
    {
        if (std::find (vPzSingleHet.begin(), vPzSingleHet.end(), pz_name) != vPzSingleHet.end()) {
            return true;
        }
        return false;
    }

    /* --- 4: HOM --- */
    // funzione per aggiungere i pazienti in cui HOM (e di concerto rimuoverli quindi da SINGLE-HET)
    void addPzHom( std::string pz_name ) {
      if (std::find (vPzHom.begin(), vPzHom.end(), pz_name) == vPzHom.end()) {
        vPzHom.push_back(pz_name);
        removePzSingleHet( pz_name );
      }
    }
    // controllo veloce se pz_name tra gli HOM della variante
    bool check_pzHom( std::string pz_name )
    {
        if (std::find (vPzHom.begin(), vPzHom.end(), pz_name) != vPzHom.end()) {
            return true;
        }
        return false;
    }
    std::vector<std::string> & getPzHom() { return vPzHom; };

    /* --- 5: COMPHET --- */
    // funzione per aggiungere un paziente (1st key) e gene (2nd key) e le altre var in comphet con la var in questione nello stesso
    void addPzComphet( std::string pz_name, std::string gene_name, std::vector<std::string> vVarsComphet ) {
        // ricreo la sub-map del GENE
        std::map<std::string, std::vector<std::string>> mGene;
        // se PZ NON ancora in MAP mPzComphet allora ci mette la MAP <gene, vars>
        if ( mPzComphet.find(pz_name) == mPzComphet.end() ) {
            mGene.insert(std::make_pair(gene_name, vVarsComphet));
            mPzComphet.insert(std::make_pair(pz_name, mGene));
        }
        else {
            // se GENE non ancora in mPzComphet[pz_name] lo aggiungo con le sue VAR
            if ( mPzComphet[ pz_name].find(gene_name) == mPzComphet[ pz_name ].end() ) {
                mPzComphet[pz_name].insert(std::make_pair(gene_name, vVarsComphet));
            }
        }
      // e lo rimuovo quindi da SINGLE-HET
      removePzSingleHet( pz_name );
    }
    // check rapido se un paziente nei COMPHET
    bool check_pzComphet( std::string pz_name )
    {
        if ( mPzComphet.find( pz_name ) != mPzComphet.end() ) {
            return true;
        }
        return false;
    }
    // funzione per avere la map con paziente:vars in COMPHET
    std::map< std::string , std::map<std::string, std::vector<std::string>>> & getPzComphet() { return mPzComphet; }

    /*
        quindi alla fine delle 3 sopra avrò i pazienti ben divisi in 3 diversi V/map: SINGLE HET (V) - HOM (V) - COMPHET (std::map<std::string, std::vector<std::string>>)
    */

    /* --- 5: ACMG --- */
    // 5a: sottocategorie
    void addAcmgCategory( std::string acmg_category ) {
      if (std::find (vAcmg.begin(), vAcmg.end(), acmg_category) == vAcmg.end()) {
        vAcmg.push_back(acmg_category);
      }
    }
    std::vector<std::string> & getAcmg() { return vAcmg; };
    // 5b: classi
    void addAcmgClass( std::string acmg_class ) {
      if (std::find (vAcmgClass.begin(), vAcmgClass.end(), acmg_class) == vAcmgClass.end()) {
        vAcmgClass.push_back(acmg_class);
      }
    }
    std::vector<std::string> & getAcmgClass() { return vAcmgClass; };
    // 5c: subclasses resumen
    void addAcmgSubclassesMap( std::map<std::string, int> &mapSubclasses ) { mAcmgSubclasses = mapSubclasses; }

	// empty ACMG classes vector
    void emptyAcmgClass() {
	  vAcmgClass.clear();
    }
	
};

// sottoclasse di VAR-PZ che mi serve per tenere le caratteristiche
// che sono tipiche della variante in relazione al paziente
// e non isolate della variante
class pzVars {
  std::map<std::string, std::string> attributes;
  public :
    void addAttr(std::string attr, std::string value) {
        if ( !attributes.insert( std::make_pair( attr, value ) ).second ) {
            // per gene_name_correct lo disabilito perche per alcune var o mette gia in TRIO COMPHET e lo ripeterebbe quando legge file DENOVO
            if ( attr != "gene_name_correct") {
                std::cout << "!!! ATTENTION !!! insert failed PZ_VAR for: " << attr << ": " << value << " (attribute already present)" << '\n';
            }
        }
    }
//    void addAttr(std::string attr, std::string value){ attributes[attr] = value; }
    std::string getValue(std::string attr){ return attributes.at(attr); }
    std::map<std::string, std::string> & getAll() { return attributes; }
    void removeAttr( std::string key )
    {
        // se presente la KEY la elimina
        if (attributes.find(key) != attributes.end()) {
            attributes.erase(key);
        }
        else {
            std::cout << "!!! ATTENTION !!! attribute: " << key << " not found in PT sub-VAR " << attributes.at("var_name") << " std::map attributes" << '\n';
        }
    }
};

class Patient {
  std::map<std::string, std::string> attributes;
  std::vector<std::string> vVar;
  // questo vettore con oggetti variante sottocaratt del paziente
  std::vector<pzVars> vPzVarsClass;
  public:
    void addAttr(std::string attr, std::string value) {
      if ( !attributes.insert( std::make_pair( attr, value ) ).second ) std::cout << "!!! ATTENTION !!! insert failed for: " << attr << ": " << value << " (attribute already present)" << '\n';
    }
//  void addAttr(std::string attr, std::string value){ attributes[attr] = value; }
    std::string getValue(std::string attr){ return attributes.at(attr); }
    std::map<std::string, std::string> & getAll() { return attributes; }
    void removeAttr( std::string key )
    {
        // se presente la KEY la elimina
        if (attributes.find(key) != attributes.end()) {
            attributes.erase(key);
        }
        else {
            std::cout << "!!! ATTENTION !!! attribute: " << key << " not found in PT " << attributes.at("pz_name") << " std::map attributes" << '\n';
        }
    }

    // per aggiunere la variante
    void addVar( std::string var ) { vVar.push_back(var); }
    std::vector<std::string> getVars() { return vVar; };
    // per aggiungere la var con sottocaratt tipiche del paziente ad apposito vettore
    void addPzVarVector( pzVars var ) { vPzVarsClass.push_back(var); }
    std::vector<pzVars> & getPzVarVector() { return vPzVarsClass; };
};

class Gene {
  std::map<std::string, std::string> attributes;
  std::vector<std::string> vVar;
  std::vector<std::string> vInhMode;
  std::vector<std::string> vGenelist;
  public:
    void addAttr(std::string attr, std::string value) {
      // su nome GENE CORRECT lo disabilito perche me lo mette gia quando fa trio analysis e lo ripeterebbe su altro file
      if ( !attributes.insert( std::make_pair( attr, value ) ).second ) {
          std::cout << "!!! ATTENTION !!! insert failed for: " << attr << ": " << value << " (attribute already present)" << '\n';
      }
    }
//  void addAttr(std::string attr, std::string value){ attributes[attr] = value; }
    //std::string getValue(std::string attr){ return attributes.at(attr); }
    std::string getValue(std::string attr)
    {
        try {
          return attributes.at(attr);;
      } catch(...) { return "NA"; }
    }
    std::map<std::string, std::string> &getAll() { return attributes; }
    void removeAttr( std::string key )
    {
        // se presente la KEY la elimina
        if (attributes.find(key) != attributes.end()) {
            attributes.erase(key);
        }
        else {
            std::cout << "!!! ATTENTION !!! attribute: " << key << " not found in GENE " << attributes.at("gene_name") << " std::map attributes" << '\n';
        }
    }

    // per aggiunere la variante
    void addVar( std::string var ) { vVar.push_back(var); }
    std::vector<std::string> & getVars() { return vVar; };
    void addInhModes( std::vector<std::string> vInhModes ) {
        vInhMode = vInhModes;
    }
    std::vector<std::string> & getInhModes() {
        return vInhMode;
    }
    void addGenelists( std::vector<std::string> vGenelists ) { vGenelist = vGenelists; }
    std::vector<std::string> & getGenelists() { return vGenelist; }
};

#endif
