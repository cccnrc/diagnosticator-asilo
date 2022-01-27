#include "filterFunctions.h"

/*
  questo file contiene i moduli di filtraggio di ogni categoria
  in modo da restituire un int: pass (0/1)
*/

// questa devo pre-dichiararla quassu perche mi serve chiamarla in questo stesso file per ACMG (e l ho messa in fondo al file la funzione)
std::map< std::string, std::vector<std::vector<std::string>> > get_ucsc_hotspot( std::string nomeFileInput );

bool contains_number(const std::string &c)
{
    return (c.find_first_of("0123456789") != std::string::npos);
}

// ritorna TRUE per i geni che hanno ALMENO un AD-XLD-U nei loro vInhMode (ovvero non solo RECESSIVE)
bool check_dominant_inh_mode ( const std::vector<std::string> &vInhMode ) {
  if (std::find (vInhMode.begin(), vInhMode.end(), "AD") != vInhMode.end()) {
      return true;
  }
  else if (std::find (vInhMode.begin(), vInhMode.end(), "XLD") != vInhMode.end()) {
      return true;
  }
  else if (std::find (vInhMode.begin(), vInhMode.end(), "U") != vInhMode.end()) {
      return true;
  }
  else {
      return false;
  }
}

// questo invece mi da TRUE se il gene è SOLO RECESSIVE
bool check_recessive_inh_mode ( const std::vector<std::string> &vInhMode ) {
    if (std::find (vInhMode.begin(), vInhMode.end(), "AR") != vInhMode.end()) {
        return true;
    }
    else {
        return false;
    }
}


// questo mi ritorna TRUE se presenti entrambi contemporaneamente
bool check_domRec_inh_mode( const std::vector<std::string> &vInhMode ) {
    if ( check_dominant_inh_mode( vInhMode ) && check_recessive_inh_mode( vInhMode ) ) {
        return true;
    }
    return false;
}

// questo invece mi da TRUE se il gene è SOLO XL/XLR RECESSIVE
bool check_xlinked_recessive_inh_mode( const std::vector<std::string> &vInhMode ) {
    if (std::find (vInhMode.begin(), vInhMode.end(), "XL") != vInhMode.end()) {
        return true;
    }
    else if (std::find (vInhMode.begin(), vInhMode.end(), "XLR") != vInhMode.end()) {
        return true;
    }
    else {
        return false;
    }
}



/* ----------
  -----------------------------------------------------------------------------
  ---------- singole funzioni di filtraggio delle car delle varianti ----------
  -----------------------------------------------------------------------------
--------- */

// filtraggio per vedere se annotate benigne su clinvar
int benign_filter(
  std::string &clinvar_clinsig,
  std::string &clinvar_disease,
  std::string &hgmd_class,
  std::string &hgmd_disease
) {
  int result;
  if ( clinvar_clinsig.find("enign") != std::string::npos && clinvar_disease.find("?Site") == std::string::npos) {
    if ( hgmd_class != "DM" && hgmd_class != "DM!" && hgmd_disease.find("?Site") == std::string::npos ) {
      result = 1;
    }
    else result = 0;
  }
  else result = 0;
  return result;
}

bool clinvar_patho_filter( std::string &clinvar_clinsig, std::string &clinvar_disease )
{
    // se clinvar trova la sottofrase athogenic e NON trova "onflict" e NON trova ?Site in Disease
    if ( clinvar_clinsig != "NA" ) {
        if ( clinvar_clinsig.find("athogenic") != std::string::npos && clinvar_clinsig.find("onflict") == std::string::npos && clinvar_disease.find("?Site") == std::string::npos) {
            return true;
        }
    }
    return false;
}

// filtraggio per vedere se annotate benigne su clinvar
int patho_filter(
  std::string &clinvar_clinsig,
  std::string &clinvar_disease,
  std::string &hgmd_class,
  std::string &hgmd_disease
) {
  int result = 0;
  int clinvar_result = 0; int hgmd_result = 0;
  // se clinvar trova la sottofrase athogenic e NON trova "onflict" e NON trova ?Site in Disease
  if ( clinvar_patho_filter(clinvar_clinsig, clinvar_disease) ) {
      clinvar_result = 1;
  }
  else {
      clinvar_result = 2;
  }
  /*
  if ( clinvar_clinsig != "NA" ) {
    if ( clinvar_clinsig.find("athogenic") != std::string::npos && clinvar_clinsig.find("onflict") == std::string::npos && clinvar_clinsig.find("onflict") == std::string::npos && clinvar_disease.find("?Site") == std::string::npos) { clinvar_result = 1; }
  } else { clinvar_result = 2; }
  */
  // se HGMD DM! or DM (ma NON DM?) e NON ?Site in Disease
  if ( hgmd_class != "NA" ) {
    if ( ( hgmd_class.find("DM!") != std::string::npos || hgmd_class.find("DM") != std::string::npos ) && hgmd_class.find("DM?") == std::string::npos && hgmd_disease.find("?Site") == std::string::npos ) {
        hgmd_result = 1;
    }
  }
  else {
      hgmd_result = 2;
  }
  // se almeno uno dei due PATHO e altro NA-PATHO
  if ( (clinvar_result == 1 && hgmd_result != 2) || (clinvar_result != 2 && hgmd_result == 1) ) {
      result = 1;
  }
  return result;
}

// filtraggio per AF 0 (compresi NA e 0)
int absent_filter(
  std::string &af
){
  int result;
  if ( af != "NA" ) {
    if ( cmpf(std::stof(af), 0.0f )) { result = 1; }
    else result = 0;
  } else result = 1;
  return result;
}

// filtraggio per AF RARE (compresi < soglia e 0)
int rare_filter(
  std::string &af,
  float &threshold
){
  int result;
  if ( af != "NA" ) {
    if ( std::stof(af) < threshold ) { result = 1; }
    else result = 0;
  } else result = 1;
  return result;
}

// filtraggio per HGMD 2bp
int hgmd_2bp_filter (
  std::string &hgmd_m2,
  std::string &hgmd_m1,
  std::string &hgmd_site,
  std::string &hgmd_class,
  std::string &hgmd_disease,
  std::string &hgmd_p1,
  std::string &hgmd_p2
) {
  int result = 0;
  std::vector<std::string> vHgmdFlank2bp = { hgmd_m2, hgmd_m1, hgmd_p1, hgmd_p2 };
  for ( auto pos : vHgmdFlank2bp ) {
    // se presente una entry per pos
    if (pos != "NA") {
      // separo in singola entry per la posizione
      std::vector<std::string> vSingleHgmdPos; char hSep = '|'; vectorize_column( pos, vSingleHgmdPos, hSep );
      // questo counter mi serve perche il PRIMO dei valori ha anche numero come inizio, mentre gli altri no
      int pCounter = 0;
      for ( auto single_pos : vSingleHgmdPos ) {
        std::vector<std::string> vSingleHgmdPosValue; char sSep = ';'; vectorize_column( single_pos, vSingleHgmdPosValue, sSep );
        // alcuni hanno solo PMID e basta e devo saltarli
        if ( vSingleHgmdPosValue.size() > 2 ) {
            // prendo il relativo campo Class e vedo se patho (se ci trova DM/DM! ma NON DM?)
            if ( pCounter == 0 ) {
                if ( ( vSingleHgmdPosValue[2].find("DM!") != std::string::npos || vSingleHgmdPosValue[2].find("DM") != std::string::npos ) && vSingleHgmdPosValue[2].find("DM?") == std::string::npos ) {
                    result = 1;
                    return result;
                }
            }
            else {
                if ( ( vSingleHgmdPosValue[1].find("DM!") != std::string::npos || vSingleHgmdPosValue[1].find("DM") != std::string::npos ) && vSingleHgmdPosValue[1].find("DM?") == std::string::npos ) {
                    result = 1;
                    return result;
                }
            }
        }
        pCounter++;
      } // chiudo loop su single entry
    } // chiudo IF NON NA
  } // chiudo loop sulle posizioni HGMD
  // se in hgmd_site ci trovo un numero maggiore di 0
  if ( hgmd_site != "NA" ) {
    hgmd_site = vectorize_column_v2( hgmd_site, ';' )[0];
    if ( std::stoi(hgmd_site) > 0 ) {
      // se in class ci trova DM/DM! ma NON DM?
      if ( ( hgmd_class.find("DM!") != std::string::npos || hgmd_class.find("DM") != std::string::npos ) && hgmd_class.find("DM?") == std::string::npos ) {
        // se trova anche ?Site in Disease significa che hgmd_class patho fa riferimento a questo
        if ( hgmd_disease.find("?Site") ) { result = 1; return result; }
        // se trova >1 in site significa che ne abbiamo altri patho oltre a lui (questo sarebbe pero gia finito in patho)
        if ( std::stoi(hgmd_site) > 1 ) { result = 1; return result; }
      // chiudo IF hgmd_class patho
      }
    // chiudo IF > 0
    }
  // chiudo SITE NON NA
  }
  return result;
}

// HGMD 9bp
int hgmd_9bp_filter (
  std::string &hgmd_9bp
){
  int result = 0;
  if ( hgmd_9bp != "NA" ) { if ( std::stoi(hgmd_9bp) > 0 ) { result = 1; return result; } }
  return result;
}

// con questa vedo se la var finisce in uno SPLICING_SITE
bool splice_filter(
  std::string &effect
) {
    // divido la colonna effect
    std::vector<std::string> vEffect; char eSep = '+'; vectorize_column( effect, vEffect, eSep );
    // apro loop sui singoli effetti e vedo se li trova nel vettore LOF e ne caso chiude funzione e ritorna 1
    for ( std::string eff : vEffect )
    {
        // cerco se parola splice nella std::string effect
        if ( eff.find("splic") != std::string::npos ) {
            return true;
        }
    }
    return false;
}

// this to see if the splice affects the gene
bool splice_patho_filter(
    std::string &ada_score,
    std::string &rf_score
) {
    if ( ada_score != "NA" ) {
        if ( floatable( ada_score ) ) {
            if ( std::stof(ada_score) > 0.59 ) {
                return true;
            }
        }
    }
    if ( rf_score != "NA" ) {
        if ( floatable( rf_score ) ) {
            if ( std::stof(rf_score) > 0.59 ) {
                return true;
            }
        }
    }
    return false;
}

// guardo quali varianti sono LOF
int lof_filter (
  std::string &var_effect,
  std::string &var_trap
){
  int result = 0;
  std::vector<std::string> vLof = {
      "exon_loss_variant",
      "frameshift_variant",
      "disruptive_inframe_deletion",
      "rare_amino_acid_variant",
      "stop_gained",
      "stop_lost",
      "start_lost",
      "gene_fusion",
      "bidirectional_gene_fusion",
      "splice_acceptor_variant",
      "splice_donor_variant"
  };
  // divido la colonna effect
  std::vector<std::string> vEffect; char eSep = '+'; vectorize_column( var_effect, vEffect, eSep );
  // apro loop sui singoli effetti e vedo se li trova nel vettore LOF e ne caso chiude funzione e ritorna 1
  for ( std::string eff : vEffect ) {
    if (std::find (vLof.begin(), vLof.end(), eff) != vLof.end()) {
        result = 1;
        return result;
    }
  }
  // se arriva fino qui vede se ha trap di almeno 0.9
  if ( floatable( var_trap ) ) {
      if ( std::stof(var_trap) > 0.89f ) {
          result = 1;
          return result;
      }
  }
  return result;
}


// guardo quali varianti sono LOF
bool pm4_effect_filter (
  std::string var_effect
){
  std::vector<std::string> vPM4effect = {
      "exon_loss_variant",
      "frameshift_variant",
      "stop_gained",
      "stop_lost",
      "start_lost"
  };
  // divido la colonna effect
  std::vector<std::string> vEffect; char eSep = '+'; vectorize_column( var_effect, vEffect, eSep );
  // apro loop sui singoli effetti e vedo se li trova nel vettore LOF e ne caso chiude funzione e ritorna 1
  for ( std::string eff : vEffect ) {
    if (std::find (vPM4effect.begin(), vPM4effect.end(), eff) != vPM4effect.end()) {
        return true;
    }
  }
  return false;
}


// guardo quali varianti sono TRAP > 0.45
int trap_filter (
  std::string &var_trap
){
  int result = 0;
  // se arriva fino qui vede se ha trap di almeno 0.9
  if ( var_trap != "NA" ) {
      if ( std::stof(var_trap) > 0.44f ) {
          result = 1;
          return result;
      }
  }
  return result;
}

// con questa vedo se la var sia sinonima o intronica
int synonymous_filter(
  std::string &effect
) {
  int result = 0;
  // non vettorializzo effect perche voglio scartare fuori quelle che sono SOLO SYN/INTR
  if ( effect == "synonymous_variant" || effect == "intron_variant" ) {
      result = 1;
      return result;
  }
  return result;
}


// guardo per PREDICTED
int pred_filter(
  std::string &humdiv_cat,
  std::string &humvar_cat,
  std::string &revel,
  std::string &subrvis_exon_centile,
  std::string &subrvis_domain_centile,
  std::string &oe_evs,
  std::string &oe_exac,
  std::string &cadd_phred,
  std::string &sift_prediction
) {
    // qui devo anche tenere conto degli NA e vedere che se unico available sia P allora considerarla tale
  int na_count = 0;
  if ( humdiv_cat == "NA" ) {
      na_count += 1;
  }
  if ( humvar_cat == "NA" ) {
      na_count += 1;
  }
  if ( sift_prediction == "NA" ) {
      na_count += 1;
  }
  if ( subrvis_exon_centile == "NA" ) {
      na_count += 1;
  }
  if ( subrvis_domain_centile == "NA" ) {
      na_count += 1;
  }
  if ( revel == "NA" ) {
      na_count += 1;
  }
  if ( oe_evs == "NA" ) {
      na_count += 1;
  }
  if ( oe_exac == "NA" ) {
      na_count += 1;
  }
  if ( cadd_phred == "NA" ) {
      na_count += 1;
  }

  int result = 0;

  // se una delle due PROBABLY o se entrambe POSSIBLY
  if ( humdiv_cat == "D" || humvar_cat == "D" || humdiv_cat == "P" || humvar_cat == "P" || humdiv_cat == "probably" || humvar_cat == "probably" || (humdiv_cat == "possibly" && humvar_cat == "possibly") ) {
      result += 1;
  }
  if ( sift_prediction == "deleterious_low_confidence" || sift_prediction == "deleterious" ) {
      result += 1;
  }
  // guarda se revel > 0.5
  if ( floatable( revel ) ) {
      if ( std::stof(revel) > 0.5f ) {
          result += 1;
    }
  }
  // guarda se subrvis_exon_centile < 35
  if ( floatable( subrvis_exon_centile ) ) {
      if ( std::stof(subrvis_exon_centile) < 35 ) {
          result += 1;
    }
  }
  // guarda se subrvis_domain_centile < 35
  if ( floatable( subrvis_domain_centile ) ) {
      if ( std::stof(subrvis_domain_centile) < 35 ) {
          result += 1;
      }
  }
  // guarda se OE[EVS] < 10
  if ( floatable( oe_evs ) ) {
      if ( std::stof(oe_evs) < 10 ) {
          result += 1;
      }
  }
  // guarda se OE[EXAC] < 10
  if ( floatable( oe_exac ) ) {
      if ( std::stof(oe_exac) < 10 ) {
          result += 1;
      }
  }
  // guarda se cadd_phred > 15
  if ( floatable( cadd_phred ) ) {
      if ( std::stof(cadd_phred) > 15.0f ) {
          result += 1;
      }
    // se CADD >22 la da comunque come patho anche se unico dato
    if ( std::stof(cadd_phred) > 22.0f ) {
        return 1;
    }
  }
  // if at least 2 predicted deleterious return 1, else 0
  if ( na_count < 5 ) {
      if ( result > 1 ) {
          return 1;
      }
  }
  else {
      if ( result > 0 ) {
          return 1;
      }
  }
  return 0;
}

// questa guarda se la var è PRED come BENIGN
bool pred_benign_filter (
    std::string &humdiv_cat,
    std::string &humvar_cat,
    std::string &revel,
    std::string &subrvis_exon_centile,
    std::string &subrvis_domain_centile,
    std::string &oe_evs,
    std::string &oe_exac,
    std::string &cadd_phred,
    std::string &sift_prediction
) {
    if ( humdiv_cat == "benign" || humvar_cat == "benign" || humdiv_cat == "B" || humvar_cat == "B" || humdiv_cat == "B" || humvar_cat == "B" ) {
        return true;
    }
    if ( floatable( cadd_phred ) ) {
        if ( std::stof( cadd_phred ) < 5 ) {
            return true;
        }
    }
    if ( floatable( revel ) ) {
        if ( std::stof( revel ) < 0.2 ) {
            return true;
        }
    }
    if ( sift_prediction == "tolerated"  ) {
        return true;
    }
    // threshold P < 35
    if ( floatable( subrvis_exon_centile ) ) {
        if ( std::stof(subrvis_exon_centile) > 50 ) {
            return true;
      }
    }
    // threshold P < 35
    if ( floatable( subrvis_domain_centile ) ) {
        if ( std::stof(subrvis_domain_centile) > 50 ) {
            return true;
        }
    }

    // threshold P OE[EVS] < 10
    if ( floatable( oe_evs ) ) {
        if ( std::stof(oe_evs) > 20 ) {
            return true;
        }
    }
    // threshold P OE[EXAC] < 10
    if ( floatable( oe_exac ) ) {
        if ( std::stof(oe_exac) > 20 ) {
            return true;
        }
    }
    return false;
}


/* ----------
  -----------------------------------------------------------------------------
  ---------- funzione iterativa per fargli passare TUTTE le varianti ----------
  -----------------------------------------------------------------------------
--------- */
void filter_variants (
    std::vector<Variant> &vVarClass,
    std::vector<std::string> &vPopulationSelected,
    std::vector<Gene> &vGeneClass,
    std::vector<std::string> &vGeneClassNames,
    float rare_threshold,
    int casesMaxAC
) {

  // non posso usare (auto var : vVarClass) perche se no non mi da errore ma non inserisce visto che fa una copia
  // devo mandare puntatore & a cella di memoria (https://stackoverflow.com/questions/56517006/c-auto-iterator-vector-of-class-elements-accession)
  for (auto &var : vVarClass) {

    /* --- CLINVAR BENIGN --- */
    // devo convertirli a std::string se no mi da errore pr const accession
    std::string var_clinsig = var.getValue("clinvar_clinsig");
    std::string var_cvDisease = var.getValue("clinvar_disease");
    std::string hgmd_class = var.getValue("hgmd_class");
    std::string hgmd_disease = var.getValue("hgmd_disease");
    // poi prendo risultato del filtraggio e lo inserisco nell elemento nella classe
    int benign_result = benign_filter( var_clinsig, var_cvDisease, hgmd_class, hgmd_disease );
    var.addAttr( "benign_cv", std::to_string(benign_result) );

    // se VAR benigna skippa a prossima iterazione
    if ( benign_result == 1 ) {
      continue;
    }

    // riprendo GENE e suoi inh_modes
    std::string var_gene_name = var.getValue("gene_name_correct");
    if (var.getValue("gene_name_correct") == "NA") {
        std::cout << "NO GENE for: " << var.getValue("var_name") << '\n';
    }

    // riprendo corrispettivo GENE in vettore classe Gene
    std::ptrdiff_t posGene = find( vGeneClassNames.begin(), vGeneClassNames.end(), var_gene_name) - vGeneClassNames.begin();
    Gene &var_gene_class = vGeneClass[ posGene ];
    // se gene REC cambio il valore di gene_inh_mode
    std::string gene_inh_mode = "dom";
    if ( check_recessive_inh_mode ( var_gene_class.getInhModes() ) || check_xlinked_recessive_inh_mode ( var_gene_class.getInhModes() ) ) {
        gene_inh_mode = "rec";
    }

    /* --- check for CASE AC MAX threshold --- */
    std::string hom_cases = var.getValue("hom_cases");
    std::string het_cases = var.getValue("het_cases");
    std::string ac_cases = var.getValue("cases_ac");
    int cases_ac = 0;
    int numNA = 0;
    if ( ac_cases != "NA" ) {
        cases_ac += std::stoi( ac_cases );
    }
    else {
        numNA += 1;
        if ( hom_cases != "NA" ) {
            cases_ac += std::stoi( hom_cases ) * 2;
        }
        else {
            numNA += 1;
        }
        if ( het_cases != "NA" ) {
            cases_ac += std::stoi( het_cases );
        }
        else {
            numNA += 1;
        }
    }


    int absent_overall = 1;
    int rare_overall = 1;
    // se gia oltre questa soglia allora imposta i valori iniziali di absent_overall e rare_overall a 0 cosi fallisce subito il filtraggio dopo
    // idem se TUTTI NA (quindi verosimilmente dato non inserito)
    if ( cases_ac > casesMaxAC || numNA > 2 ) {
        absent_overall = 0;
        rare_overall = 0;
    }

    /* --- DATABASES AF --- */
    // apro loop sulle popolazioni selezionate per GnomaAD e inserisco nel vettore risultati il relativo filter result
    for (auto pop : vPopulationSelected) {
      std::string pop_af = var.getValue( pop );
      int absent_result = absent_filter( pop_af );
      int rare_result = rare_filter( pop_af, rare_threshold );
      // se non trova NESSUNO 0 (tutti hanno passato filtro allora mette overall a 1)
      // se non passa RARE rompe il loop (così passa a VAR successiva senza farli tutti)
      // con absent non posso farlo perche invece ci sono var che possono essere RARE ma NON ABSENT (non viceversa)
      // cosi alla fine del loop di tutte le popolazioni se ancora 1 significa che ha passato filtri in tutte le popolazioni
      // se invece 0 ne ha fallita almeno una
      if (absent_result == 0) {
          absent_overall = 0;
      }
      if (rare_result == 0) {
          rare_overall = 0;
          // SOLO SE DOM rompe il loop (sui REC non voglio filtraggio di MAF)
          if ( gene_inh_mode == "dom" ) {
              break;
          }
      }
    // chiudo loop su popolazioni selezionate
    }
    // finito il loop delle popolazioni selezionate per ogni variante inserisce il risultato corrispondente nella classe
    if ( absent_overall == 1 && rare_overall == 1 ) {
        var.addAttr("af_absent","1");
        var.addAttr("af_rare","NA");
    }
    //se NON ASSENTE MA RARA mette 1 alla relativa categoria (per non avere 1 sia in ASSENTI che RARE ma solo in RARE)
    else {
      var.addAttr("af_absent","0");
      if ( rare_overall == 1 ) {
          var.addAttr("af_rare","1");
      }
      else {
          var.addAttr( "af_rare", "0");
      }
    }

    /* --- CLINVAR PATHO --- */
    int patho_result = patho_filter(var_clinsig, var_cvDisease, hgmd_class, hgmd_disease);
    var.addAttr( "patho", std::to_string(patho_result) );

    /* --- LOF --- */
    std::string effect = var.getValue("effect");
    std::string trap_score = var.getValue("trap_score");
    int lof_result = lof_filter( effect, trap_score );

    // if NON-lof effect directly calculate if splice effective on the gene transcript
    std::string ada_score = var.getValue("ada_score");
    std::string rf_score = var.getValue("rf_score");
    int splice_patho_result = 0;
    if ( lof_result == 0 ) {
        if ( splice_filter( effect ) ) {
            if ( splice_patho_filter( ada_score, rf_score ) ) {
                splice_patho_result = 1;
            }
        }
    }

    /* --- TRAP --- */
    // solo se NON LOF calcola score
    int trap_result;
    if ( lof_result == 0 ) {
        trap_result = trap_filter( trap_score );
        var.addAttr( "trap", std::to_string(trap_result) );
    }
    else {
        var.addAttr( "trap", "NA" );
        trap_result = 0;
    }

    /* --- SYN/INTRON --- */
    // solo se NON PATHO-LOF-TRAP
    int synonymous_result;
    if ( patho_result == 0 && lof_result == 0 && trap_result == 0 ) { synonymous_result = synonymous_filter( effect ); var.addAttr( "synonymous", std::to_string(synonymous_result) ); }
    else { var.addAttr( "synonymous", "NA" ); synonymous_result = 0; }

    /* --- PREDICTED --- */
    // se NON patho - LOF - TRAP e NON SYN calcolo PRED e FLANK
    if ( patho_result == 0 && splice_patho_result == 0 && lof_result == 0 && trap_result == 0 && synonymous_result == 0 ) {
      std::string humdiv_cat = var.getValue("humdiv_cat");
      std::string humvar_cat = var.getValue("humvar_cat");
      std::string revel = var.getValue("revel");
      std::string subrvis_exon_centile = var.getValue("subrvis_exon_centile");
      std::string subrvis_domain_centile = var.getValue("subrvis_domain_centile");
      std::string oe_evs = var.getValue("oe_evs");
      std::string oe_exac = var.getValue("oe_exac");
      std::string cadd_phred = var.getValue("cadd_phred");
      std::string sift_prediction = var.getValue("sift_prediction");
      int pred_result = pred_filter( humdiv_cat, humvar_cat, revel, subrvis_exon_centile, subrvis_domain_centile, oe_evs, oe_exac, cadd_phred, sift_prediction );
      var.addAttr( "pred", std::to_string(pred_result) );

      /* --- HGMD 2bp - 9 bp --- */
      std::string hgmd_m2 = var.getValue("hgmd_m2");
      std::string hgmd_m1 = var.getValue("hgmd_m1");
      std::string hgmd_site = var.getValue("hgmd_site");
      std::string hgmd_p1 = var.getValue("hgmd_p1");
      std::string hgmd_p2 = var.getValue("hgmd_p2");
      std::string hgmd_9bp = var.getValue("hgmd_9bp");
      // fa FLANK solo se anche PRED
      if ( pred_result == 1) {
        int hgmd_2bp_result = hgmd_2bp_filter ( hgmd_m2, hgmd_m1, hgmd_site, hgmd_class, hgmd_disease, hgmd_p1, hgmd_p2 );
        int hgmd_9bp_result = hgmd_9bp_filter ( hgmd_9bp );
        var.addAttr( "hgmd_flank_2bp", std::to_string( hgmd_2bp_result ) );
        var.addAttr( "hgmd_flank_9bp", std::to_string( hgmd_9bp_result ) );
      } else { var.addAttr( "hgmd_flank_2bp", "NA" ); var.addAttr( "hgmd_flank_9bp", "NA" ); }

    // se invece la var e almeno 1 tra PATHO - LOF - TRAP - SYN
  } else { var.addAttr( "pred", "NA" ); var.addAttr( "hgmd_flank_2bp", "NA" ); var.addAttr( "hgmd_flank_9bp", "NA" ); }

  /* ------- assegnazione CATEGORIA ------- */
  // se var ASSENTE
  if ( var.getValue("af_absent") == "1" ) {
    if ( patho_result == 1 ) var.addAttr( "dom", "patho0" );
    else if ( lof_result == 1 || splice_patho_result == 1 ) var.addAttr( "dom", "lof0" );
    else if ( trap_result == 1 ) var.addAttr( "dom", "trap0" );
    else if ( var.getValue("hgmd_flank_2bp") == "1" ) var.addAttr( "dom", "aa0" );
    else if ( var.getValue("hgmd_flank_9bp") == "1" ) var.addAttr( "dom", "flank0" );
    else if ( var.getValue("pred") == "1" ) var.addAttr( "dom", "pred0" );
    else { var.addAttr( "dom", "NA" ); }
  // se invece solo RARA
  } else if ( var.getValue("af_rare") == "1" ) {
    if ( patho_result == 1 ) var.addAttr( "dom", "patho5" );
    else if ( lof_result == 1 || splice_patho_result == 1 ) var.addAttr( "dom", "lof5" );
    else if ( trap_result == 1 ) var.addAttr( "dom", "trap5" );
    else if ( var.getValue("hgmd_flank_2bp") == "1" ) var.addAttr( "dom", "aa5" );
    else if ( var.getValue("hgmd_flank_9bp") == "1" ) var.addAttr( "dom", "flank5" );
    else if ( var.getValue("pred") == "1" ) var.addAttr( "dom", "pred5" );
    else { var.addAttr( "dom", "NA" ); }
  }
  else {
      if ( patho_result == 1 ) var.addAttr( "dom", "pathoR" );
      else if ( lof_result == 1 || splice_patho_result == 1  ) var.addAttr( "dom", "lofR" );
      else if ( trap_result == 1 ) var.addAttr( "dom", "trapR" );
      else if ( var.getValue("hgmd_flank_2bp") == "1" ) var.addAttr( "dom", "aaR" );
      else if ( var.getValue("hgmd_flank_9bp") == "1" ) var.addAttr( "dom", "flankR" );
      else if ( var.getValue("pred") == "1" ) var.addAttr( "dom", "predR" );
      else { var.addAttr( "dom", "NA" ); }
  }

  // chiudo LOOP sulle singole var della classe varianti
  }
}




/* ----------
  -----------------------------------------------------------------------------
  --------------- funzioni di filtraggio per le ACMG guidelines ---------------
  -----------------------------------------------------------------------------
--------- */
// guardo per PREDICTED
bool denovo_filter( std::string var_denovo ) {
  std::vector<std::string> vDenovo = {
      "COMPOUND DELETION",
      "DE NOVO",
      "DE NOVO; IN REF",
      "NEWLY HEMIZYGOUS",
      "NEWLY HEMIZYGOUS OR POSSIBLY DE NOVO",
      "NEWLY HEMIZYGOUS OR POSSIBLY DE NOVO; IN REF",
      "NEWLY HEMIZYGOUS; IN REF",
      "NEWLY HOMOZYGOUS",
      "POSSIBLY DE NOVO",
      "POSSIBLY DE NOVO; IN REF",
      "POSSIBLY NEWLY HEMIZYGOUS",
      "POSSIBLY NEWLY HEMIZYGOUS; IN REF"
  };
  // if denovo variant
  if (std::find (vDenovo.begin(), vDenovo.end(), var_denovo) != vDenovo.end()) { return true; }
  else { return false; }
}

bool ps3_filter(
    std::string clinvar_clinsig,
    std::string clinvar_disease,
    std::string clinvar_clinrevstars,
    std::string hgmd_class,
    std::string hgmd_disease
) {
  // se clinvar trova la sottofrase athogenic e NON trova "onflict" e NON trova ?Site in Disease e con almeno 1 stella di raccomandazione
  if ( clinvar_clinsig != "NA" && clinvar_clinrevstars != "NA" ) {
    if ( clinvar_clinsig.find("athogenic") != std::string::npos && clinvar_clinsig.find("onflict") == std::string::npos && clinvar_clinsig.find("onflict") == std::string::npos && clinvar_disease.find("?Site") == std::string::npos) {
      if ( std::stoi( clinvar_clinrevstars ) > 1 || ( ( hgmd_class == "DM" || hgmd_class == "DM!" )  && hgmd_disease.find("?Site") == std::string::npos ) ) {
          return true;
      }
      else { return false; }
    }
    else { return false; }
  }
  else { return false; }
}

bool pp5_filter(
    std::string clinvar_clinsig,
    std::string clinvar_disease,
    std::string clinvar_clinrevstars
) {
  // se clinvar trova la sottofrase athogenic e NON trova "onflict" e NON trova ?Site in Disease e con almeno 1 stella di raccomandazione
  if ( clinvar_clinsig != "NA" && clinvar_clinrevstars != "NA" ) {
    if ( clinvar_clinsig.find("athogenic") != std::string::npos && clinvar_clinsig.find("onflict") == std::string::npos && clinvar_clinsig.find("onflict") == std::string::npos && clinvar_disease.find("?Site") == std::string::npos) {
      if ( std::stoi( clinvar_clinrevstars ) < 2 ) {
          return true;
      }
      else { return false; }
    }
    else { return false; }
  }
  else { return false; }
}


// funzione per filtrare ASSENZA da controlli
bool pm2_filter(
    std::string gnomad_ex_controls_an,
    std::string gnomad_gen_controls_an,
    std::string gnomad_ex_controls_af,
    std::string gnomad_gen_controls_af,
    std::string gnomad_ex_global_af,
    std::string gnomad_gen_global_af,
    std::string gnomad_ex_global_nhet,
    std::string gnomad_ex_global_nhomalt,
    std::string gnomad_ex_controls_nhomalt,
    std::string gnomad_gen_global_nhet,
    std::string gnomad_gen_global_nhomalt,
    std::string gnomad_gen_controls_nhomalt,
    std::string &gene_inh_mode
) {
  if ( gnomad_ex_controls_an == "NA" && gnomad_gen_controls_an == "NA" && gnomad_ex_controls_af == "NA" && gnomad_gen_controls_af == "NA" && gnomad_ex_global_af == "NA" && gnomad_gen_global_af == "NA" ) {
      return true;
  }
  else {
    if ( gene_inh_mode == "dom" || gene_inh_mode == "NA" ) {
      // questo riquadro sotto da usare solo per i VCF
      // NB da VCF il gnomad_ex_global_af in realta si tratta del popmax
      if ( ( gnomad_ex_global_nhet == "NA" && gnomad_ex_global_nhomalt == "NA" ) && gnomad_ex_global_af != "NA" ) {
          if ( floatable( gnomad_ex_global_af ) ) {
              if ( std::stof(gnomad_ex_global_af) < 1E-4f ) {
                  return true;
              }
          }
      }
      // questo riquadro sotto da usare solo per i VCF
      // NB da VCF il gnomad_ex_global_af in realta si tratta del popmax
      // quindi devo fare in modo di farglielo capire automaticamente
      // e lo faccio usando il fatto che sia HET che NHOMALT sarebbero NA e non lo sarebbe gnomad_gen_global_af
      if ( ( gnomad_gen_global_nhet == "NA" && gnomad_gen_global_nhomalt == "NA" ) && gnomad_gen_global_af != "NA" ) {
          if ( floatable( gnomad_gen_global_af ) ) {
              if ( std::stof(gnomad_gen_global_af) < 1E-4f ) {
                  return true;
              }
          }
      }
      // if from ATAV
      else {
          if ( gnomad_ex_controls_nhomalt != "NA" ) {
              if (std::stoi(gnomad_ex_controls_nhomalt) == 0) {
                  if ( gnomad_ex_global_nhet != "NA" ) {
                      if (std::stoi(gnomad_ex_global_nhet) < 6) {
                        return true;
                    }
                }
              }
          }
          if ( gnomad_gen_controls_nhomalt != "NA" ) {
              if (std::stoi(gnomad_gen_controls_nhomalt) == 0) {
                  if ( gnomad_gen_global_nhet != "NA" ) {
                      if (std::stoi(gnomad_gen_global_nhet) < 6) {
                        return true;
                    }
                }
              }
          }
      }
    // se recessive deve controllare che siano RARI
  } else if ( gene_inh_mode == "rec" ) {
      int passed = 0;
      if ( gnomad_ex_controls_af != "NA" && floatable(gnomad_ex_controls_af) ) {
          if ( std::stof(gnomad_ex_controls_af) < 1E-4f ) {
              passed += 1;
          }
      }
      if ( gnomad_gen_controls_af != "NA" && floatable(gnomad_gen_controls_af) ) {
          if (std::stof(gnomad_gen_controls_af) < 1E-4f ) {
              passed += 1;
          }
      }
      if ( gnomad_ex_global_af != "NA" && floatable(gnomad_ex_global_af) ) {
          if ( std::stof(gnomad_ex_global_af) < 1E-4f ) {
              passed += 1;
          }
      }
      if ( gnomad_gen_global_af != "NA" && floatable(gnomad_gen_global_af) ) {
          if ( std::stof(gnomad_gen_global_af) < 1E-4f ) {
              passed += 1;
          }
      }
      if ( passed > 0 ) {
          return true;
      }
      else {
          return false;
      }
      return false;
    }
  // chiudo IF non-NA entrambi
  }
  return false;
}

bool frequent_filter(
  Variant &var,
  std::vector<std::string> &vPopulationSelected,
  float soglia
) {
  int result = 0;
  for ( std::string pop : vPopulationSelected ) {
    std::string value = var.getValue(pop);
    if ( value != "NA" ) {
      if ( floatable(value) ) {
        // se >soglia
        if ( std::stof( value ) > soglia  ) {
            return false;
        }
      }
    }
  }
  return true;
}

bool gerp_filter( std::string gerp_score ) {
  if ( gerp_score == "NA" ) {
      return false;
  }
  else {
    if (floatable(gerp_score)) {
      if (std::stof(gerp_score) > 2) {
          return true;
      }
      else {
          return false;
      }
    } else {
        return false;
    }
  }
}

bool bs2_filter (
    std::string gnomad_ex_controls_nhet,
    std::string gnomad_ex_controls_nhomalt,
    std::string gnomad_ex_global_af,
    std::string gnomad_gen_controls_nhet,
    std::string gnomad_gen_controls_nhomalt,
    std::string gnomad_gen_global_af,
    std::string &gene_inh_mode
){
    if ( gene_inh_mode == "dom" || gene_inh_mode == "NA" ) {
      // questo riquadro sotto da usare solo per i VCF
      // NB da VCF il gnomad_ex_global_af in realta si tratta del popmax
      if ( ( gnomad_ex_controls_nhet == "NA" && gnomad_ex_controls_nhomalt == "NA" ) && gnomad_ex_global_af != "NA" ) {
      }
      // se invece si tratta di ATAV
      else {
          if ( gnomad_ex_controls_nhet != "NA" ) {
              if ( std::stoi(gnomad_ex_controls_nhet) > 2 ) {
                  return true;
              }
          }
      }
      if ( ( gnomad_gen_controls_nhet == "NA" && gnomad_gen_controls_nhomalt == "NA" ) && gnomad_gen_global_af != "NA" ) {
      }
      else {
          if ( gnomad_gen_controls_nhet != "NA" ) {
              if ( std::stoi(gnomad_gen_controls_nhet) > 2 ) {
                  return true;
              }
          }
      }
  }
  // se gene RECESSIVE invece
  else {
      if ( ( gnomad_ex_controls_nhet == "NA" && gnomad_ex_controls_nhomalt == "NA" ) && gnomad_ex_global_af != "NA" ) {
      }
      // se invece si tratta di ATAV
      else {
          if ( gnomad_ex_controls_nhomalt != "NA" ) {
              if ( std::stoi(gnomad_ex_controls_nhomalt) > 5 ) {
                  return true;
              }
          }
      }

      if ( ( gnomad_gen_controls_nhet == "NA" && gnomad_gen_controls_nhomalt == "NA" ) && gnomad_gen_global_af != "NA" ) {
      }
      // se invece si tratta di ATAV
      else {
          if ( gnomad_gen_controls_nhomalt != "NA" ) {
              if ( std::stoi(gnomad_gen_controls_nhomalt) > 5 ) {
                  return true;
              }
          }
      }
  }
      return false;
}


// reputable benign in ClinVar
bool bs3_filter( std::string clinvar_clinsig, std::string clinvar_disease, std::string clinvar_clinrevstars ) {
  // se clinvar trova la sottofrase athogenic e NON trova "onflict" e NON trova ?Site in Disease e con almeno 1 stella di raccomandazione
  if ( clinvar_clinsig != "NA" && clinvar_clinrevstars != "NA" ) {
    if ( clinvar_clinsig.find("enign") != std::string::npos && clinvar_clinsig.find("onflict") == std::string::npos && clinvar_disease.find("?Site") == std::string::npos) {
      if ( std::stoi( clinvar_clinrevstars ) > 0 ) { return true; }
      else { return false; }
    } else { return false; }
  } else { return false; }
}

// questo per rimetter hgvs_p come annotato da ATAV in stesso formato che ha mio file da ClinVar
std::string adapt_atav_HGVSp ( std::string hgvs_p ) {
      replace_substring( hgvs_p, "p.", "" );
      replace_substring( hgvs_p, "*", "Ter" );
      std::string::size_type pos1 = hgvs_p.find_first_of("0123456789");
      std::string p_sub_ref = hgvs_p.substr(0, pos1);
      std::string::size_type pos2;
      std::string p_sub_alt;
      // se finisce con numero ritorna NA
      if ( hgvs_p.find_first_not_of("0123456789", pos1) == std::string::npos ) {
          pos2 = hgvs_p.length();
          p_sub_alt = "";
      }
      else {
          pos2 = hgvs_p.find_first_not_of("0123456789", pos1);
          p_sub_alt = hgvs_p.substr(pos2);
      }
      std::string p_sub_pos = hgvs_p.substr(pos1, pos2-pos1);
      if ( p_sub_ref == p_sub_alt ) {
          p_sub_alt = "%3D";
      }
      std::stringstream streamSub;
      streamSub << p_sub_ref << p_sub_pos << p_sub_alt;
      std::string new_sub = streamSub.str();
      return new_sub;
}

// questo per prendere AA ref e POS in modo da vedere se oltre alla nota ce ne sono altre
std::string adapt_atav_HGVSp_return_refPos ( std::string hgvs_p ) {
      replace_substring( hgvs_p, "p.", "" );
      replace_substring( hgvs_p, "*", "Ter" );
      std::string::size_type pos1 = hgvs_p.find_first_of("0123456789");
      std::string p_sub_ref = hgvs_p.substr(0, pos1);
      std::string::size_type pos2;
      std::string p_sub_alt;
      // se finisce con numero ritorna NA
      if ( hgvs_p.find_first_not_of("0123456789", pos1) == std::string::npos ) {
          pos2 = hgvs_p.length();
          p_sub_alt = "";
      }
      else {
          pos2 = hgvs_p.find_first_not_of("0123456789", pos1);
          p_sub_alt = hgvs_p.substr(pos2);
      }
      std::string p_sub_pos = hgvs_p.substr(pos1, pos2-pos1);
      if ( p_sub_ref == p_sub_alt ) {
          p_sub_alt = "%3D";
      }
      std::stringstream streamSub;
      streamSub << p_sub_ref << p_sub_pos;
      std::string new_sub = streamSub.str();
      return new_sub;
}


// questo per rimetter hgvs_p come annotato da ATAV in stesso formato che ha mio file da ClinVar
std::string extract_pos_atav_HGVSp ( std::string hgvs_p ) {
      replace_substring( hgvs_p, "p.", "" );
      std::string::size_type pos1 = hgvs_p.find_first_of("0123456789");
      std::string::size_type pos2;
      // se finisce con numero ritorna NA
      if ( hgvs_p.find_first_not_of("0123456789", pos1) == std::string::npos ) {
          pos2 = hgvs_p.length();
      }
      else {
          pos2 = hgvs_p.find_first_not_of("0123456789", pos1);
      }
      std::string p_sub_pos = hgvs_p.substr(pos1, pos2-pos1);
      return p_sub_pos;
}

// questa per inserire le sottoclassi ACMG in oggetto VAR della classe Variant
void add_acmg_subclasses(
    Variant &var,
    std::vector<Gene> &vGeneClass,
    std::vector<std::string> &vGeneClassNames,
    std::vector<std::string> &vLOFgenes,
    std::vector<std::string> &vMISSgenes,
    std::vector<std::vector<std::vector<std::string>>> &vvvHGVSPgenes,
    std::vector<std::string> &vHGVSPgenes,
    std::map<std::string, std::vector<std::string>> &mHGVSPpatho,
    std::map<std::string, std::vector<std::string>> &mHGVSPpos,
    int &parents_aff,
    std::map< std::string, std::vector<std::vector<std::string>> > &mHtospotUcsc,
    std::vector<std::string> &vHotspotCv,
    std::vector<std::string> vGeniMapHotspot,
    std::map<std::string, float> &mBS1genes,
    std::vector<std::string> &vBA1exceptions
) {

  // prendo il nome del gene
  std::string gene_name = var.getValue("gene_name_correct");

  // riprendo corrispettivo GENE in vettore classe Gene
  std::ptrdiff_t posGene = find( vGeneClassNames.begin(), vGeneClassNames.end(), gene_name) - vGeneClassNames.begin();
  Gene &var_gene_class = vGeneClass[ posGene ];

  /* ----------- 1: PVS ------------- */
  /* ----- 1a: PVS1 ----- */
  /* --- LOF in a LOF gene --- */
  // se gene in LOF fa il check per LOF-TRAP function
  int pvs1_variant = 0;
  if (std::find (vLOFgenes.begin(), vLOFgenes.end(), gene_name) != vLOFgenes.end()) {
    std::string effect = var.getValue("effect");
    std::string trap_score = var.getValue("trap_score");
    std::string ada_score = var.getValue("ada_score");
    std::string rf_score = var.getValue("rf_score");
    int lof_result = lof_filter( effect, trap_score );
    // se LOF in gene LOF mette in categoria acmg PVS1_lof
    if ( lof_result == 1 ) {
        var.addAcmgCategory("PVS1_lof");
		pvs1_variant=1;
    }
    // se non LOF ma trap>0.45 mette in acmg PVS1_trap
    else {
      //int trap_result = trap_filter( trap_score );
      int trap_result = 0;
      if (trap_result == 1 ) {
          var.addAcmgCategory("PVS1_trap");
		  pvs1_variant=1;
      }
      else {
          if ( splice_filter( effect ) ) {
              if ( splice_patho_filter( ada_score, rf_score ) ) {
                  var.addAcmgCategory("PVS1_splice");
				  pvs1_variant=1;
              }
          }
      }
    }
  // chiudo IF gene LOF
  }

  /* ----------- 2: PS ------------- */
  /* ----- 2a: PS1 & PM5 ----- */
  // se NON ClinVar nota benign
  int ps1_variant = 0;
  if (std::find (vBA1exceptions.begin(), vBA1exceptions.end(), var.getValue("var_name")) != vBA1exceptions.end()) {
    var.addAcmgCategory("PS1");
    ps1_variant = 1;
  }
  std::string var_clinsig = var.getValue("clinvar_clinsig"); std::string var_clindisease = var.getValue("clinvar_disease");
  if ( var.getValue("benign_cv") == "0" && var.getValue("effect") != "synonymous_variant" ) {
      // se in un gene tra quelli di cui abbiamo HGVS.p pathogeniche in ClinVar
      if (std::find (vHGVSPgenes.begin(), vHGVSPgenes.end(), gene_name) != vHGVSPgenes.end()) {
          // se non-NA HGVS.p
          if (var.getValue("hgvs_p") != "NA" && contains_number(var.getValue("hgvs_p"))) {
              // lo rimetto in stesso formato di mio file
              std::string var_hgvs_p = adapt_atav_HGVSp(var.getValue("hgvs_p"));
              // se la variante stessa NON gia trovata in clinvar cerca se altre var hanno = sostituzione AA
              if (! clinvar_patho_filter( var_clinsig, var_clindisease )) {
                  if (std::find (mHGVSPpatho[gene_name].begin(), mHGVSPpatho[gene_name].end(), var_hgvs_p ) != mHGVSPpatho[gene_name].end()) {
                    if ( ps1_variant == 0 ) {
                      var.addAcmgCategory("PS1");
                    }
                  }
                  // se NON trovata stessa sostituzione AA vedo se POS
                  else {
                      std::string var_hgvs_p_pos = extract_pos_atav_HGVSp(var.getValue("hgvs_p"));
                      if (std::find (mHGVSPpos[gene_name].begin(), mHGVSPpos[gene_name].end(), var_hgvs_p_pos ) != mHGVSPpos[gene_name].end()) {
                          var.addAcmgCategory("PM5");
                      }
                  }
              }
              // se gia nota in clinvar guarda comunque se ce ne sono altre su stesso AA prima come refPos poi come pos e basta
              else {
                  std::string var_hgvsp_posRef = adapt_atav_HGVSp_return_refPos(var.getValue("hgvs_p"));
                  // open loop on HGVSp and look if find another with same AA and POS substitution
                  for ( std::string hgvsp_patho : mHGVSPpatho[gene_name] )
                  {
                      if ( hgvsp_patho.find( var_hgvsp_posRef ) != std::string::npos ) {
                          // if NOT exactly the same variant
                          if ( hgvsp_patho != var_hgvs_p ) {
                              var.addAcmgCategory("PM5");
                          }
                      }
                  }
                  if (std::find (mHGVSPpatho[gene_name].begin(), mHGVSPpatho[gene_name].end(), var_hgvs_p ) == mHGVSPpatho[gene_name].end()) {
                      std::string var_hgvs_p_pos = extract_pos_atav_HGVSp(var.getValue("hgvs_p"));
                      if (std::find (mHGVSPpos[gene_name].begin(), mHGVSPpos[gene_name].end(), var_hgvs_p_pos ) != mHGVSPpos[gene_name].end()) {
                          var.addAcmgCategory("PM5");
                      }
                  }
              }
          }
      }
  }

  /* ----- 2b: PS2 ----- this need to be per patient !!! */
  /* --- De novo (both maternity and paternity confirmed) in a patient with the disease and no family history --- */
  // if denovo variant e PARENTS NOT AFFECTED
  if (denovo_filter(var.getValue("denovo")) == 1 && parents_aff == 0) { var.addAcmgCategory("PS2"); }

  /* ----- 2c: PS3 ----- */
  /* --- Well-established in vitro or in vivo functional studies supportive of a damaging effect on the gene or gene product --- */
  int ps3_variant = 0;
  if ( ps3_filter( var.getValue("clinvar_clinsig"), var.getValue("clinvar_disease"), var.getValue("clinvar_clinrevstars"), var.getValue("hgmd_class"), var.getValue("hgmd_disease") ) ) {
      var.addAcmgCategory("PS3");
      ps3_variant = 1;
  }

  /* ----- 2c: PS4 ----- */
  /* --- The prevalence of the variant in affected individuals is significantly increased compared with the prevalence in controls. --- */
  // ... users defined set of RISK VARIANTS ...

  /* ----------- 3: PM ------------- */
  /* ----- 3a: PM1 ----- */
  /* --- Located in a mutational hot spot and/or critical and well-established functional domain (e.g., active site of an enzyme) without benign variation --- */
  // solo SE gene tra quelli di cui conosco HOTSPOT faccio filtraggio per questo fattore (per ora la faccio con EXON ClinVar REV)
  int pm1_variant = 0;
  if (std::find (vGeniMapHotspot.begin(), vGeniMapHotspot.end(), gene_name) != vGeniMapHotspot.end()) {
      // riprendo la POS della var dal NOME
      int pos = std::stoi( vectorize_column_v2( var.getValue("var_name"), '-' )[1] );
      for ( std::vector<std::string> vExon : mHtospotUcsc[ gene_name ] )
      {
          // guardo se la POS della VAR cade nell intervallo dell esone
          if ( inRange( std::stoi(vExon[1]), std::stoi(vExon[2]), pos ) ) {
              var.addAcmgCategory("PM1");
              pm1_variant = 1;
              //std::cout << "in HOTSPOT " << gene_name << ", chr " << vExon[0] << " : " << var.getValue("var_name") << " [" << vExon[1] << " - " << vExon[2] << "]" << '\n';
          }
      }
  }
  // if not found previously try with ClinVar 60bp flanking region
  if ( pm1_variant == 0 ) {
      if (std::find (vHotspotCv.begin(), vHotspotCv.end(), var.getValue("var_name")) != vHotspotCv.end()) {
          var.addAcmgCategory("PM1");
          pm1_variant = 1;
      }
  }


  /* ----- 3b: PM2 ----- */
  // I need this to avoid filtering BS1 and BA1 if already a PM2 (mutually exclusive, see below)
  int pm2_variant = 0;
  /* --- Absent from controls (or at extremely low frequency if recessive) in Exome Sequencing Project, 1000 Genomes Project, or Exome Aggregation Consortium --- */
  std::string gene_inh_mode = "dom";
  // se gene REC cambio il valore di gene_inh_mode
  if ( check_recessive_inh_mode ( var_gene_class.getInhModes() ) || check_xlinked_recessive_inh_mode ( var_gene_class.getInhModes() ) || check_domRec_inh_mode( var_gene_class.getInhModes() ) ) {
      gene_inh_mode = "rec";
  }
  if ( pm2_filter(
            var.getValue("gnomad_ex_controls_an"),
            var.getValue("gnomad_gen_controls_an"),
            var.getValue("gnomad_ex_controls_af"),
            var.getValue("gnomad_gen_controls_af"),
            var.getValue("gnomad_ex_global_af"),
            var.getValue("gnomad_gen_global_af"),
            var.getValue("gnomad_ex_global_nhet"),
            var.getValue("gnomad_ex_global_nhomalt"),
            var.getValue("gnomad_ex_controls_nhomalt"),
            var.getValue("gnomad_gen_global_nhet"),
            var.getValue("gnomad_gen_global_nhomalt"),
            var.getValue("gnomad_gen_controls_nhomalt"),
            gene_inh_mode
        ) ) {
      var.addAcmgCategory("PM2");
      pm2_variant = 1;
  }

  /* ----- 3c: PM3 ----- this must be per patient !!! */
  /* --- For recessive disorders, detected in trans with a pathogenic variant --- */
  /* -------> INTEGRATA DOPO (perche deve fare tutte le VAR prima e fare analisi COMPHET). See function: acmg_pm3_screening() <-------- */

  /* ----- 3d: PM4 ----- */
  if ( pm4_effect_filter( var.getValue("effect") ) && var.getValue("repeated") == "false" ) {
      var.addAcmgCategory("PM4");
  }
  /* --- Protein length changes as a result of in-frame deletions/insertions in a non-repeat region or stop-loss variants --- */


  /* ----- 3e: PM5 ----- */
  /* --- Novel missense change at an amino acid residue where a different missense change determined to be pathogenic has been seen before --- */
  // integrata con PS1 !!!

  /* ----- 3f: PM6 ----- this must be per patient !!! */
  /* --- Assumed de novo, but without confirmation of paternity and maternity --- */

  /* ----------- 4: PP ------------- */
  /* ----- 4a: PP1 -----  this must be per patient !!!*/
  /* --- Cosegregation with disease in multiple affected family members in a gene definitively known to cause the disease --- */

  /* ----- 4b: PP2 ----- */
  /* --- Missense variant in a gene that has a low rate of benign missense variation and in which missense variants are a common mechanism of disease --- */
  // se in gene MISS
  int pp2_variant = 0;
  if (std::find (vMISSgenes.begin(), vMISSgenes.end(), gene_name) != vMISSgenes.end()) {
    std::string effect = var.getValue("effect");
    std::string trap_score = var.getValue("trap_score");
    int lof_result = lof_filter( effect, trap_score );
    int trap_result = trap_filter( trap_score );
    // se NON-LOF in gene LOF guarda anche che non sia synonymous-intronic
    if ( lof_result == 0 && trap_result == 0 ) {
        if ( synonymous_filter( effect ) == 0 ) {
            if ( effect != "conservative_inframe_deletion" ) {
                var.addAcmgCategory("PP2");
                pp2_variant = 1;
            }
        }
    }
  // chiudo IF gene MISS
  }

  /* ----- 4c: PP3 ----- */
  /* --- Multiple lines of computational evidence support a deleterious effect on the gene or gene product (conservation, evolutionary, splicing impact, etc.) --- */
  int pp3_var = 0;
  std::string humdiv_cat = var.getValue("humdiv_cat");
  std::string humvar_cat = var.getValue("humvar_cat");
  std::string revel = var.getValue("revel");
  std::string subrvis_exon_centile = var.getValue("subrvis_exon_centile");
  std::string subrvis_domain_centile = var.getValue("subrvis_domain_centile");
  std::string oe_evs = var.getValue("oe_evs");
  std::string oe_exac = var.getValue("oe_exac");
  std::string cadd_phred = var.getValue("cadd_phred");
  std::string sift_prediction = var.getValue("sift_prediction");
  if ( pred_filter(
            humdiv_cat,
            humvar_cat,
            revel,
            subrvis_exon_centile,
            subrvis_domain_centile,
            oe_evs,
            oe_exac,
            cadd_phred,
            sift_prediction
        ) == 1 ) {
      var.addAcmgCategory("PP3");
      pp3_var = 1;
  }

  /* ----- 4d: PP4 ----- this must be per patient !!! */
  /* --- Patient’s phenotype or family history is highly specific for a disease with a single genetic etiology --- */

  /* ----- 4e: PP5 ----- */
  /* --- Reputable source recently reports variant as pathogenic, but the evidence is not available to the laboratory to perform an independent evaluation --- */
  // if the var FAILS the PS3 (ClinVar 2 star at least) but is annotated as patho in ClinVar/HGMD
  int pp5_variant = 0;
  if ( ! ps3_filter( var.getValue("clinvar_clinsig"), var.getValue("clinvar_disease"), var.getValue("clinvar_clinrevstars"), var.getValue("hgmd_class"), var.getValue("hgmd_disease")  ) ) {
      if ( pp5_filter( var.getValue("clinvar_clinsig"), var.getValue("clinvar_disease"), var.getValue("clinvar_clinrevstars") ) ) {
          if ( var.getValue("benign_cv") == "0" ) {
          var.addAcmgCategory("PP5");
          pp5_variant = 1;
      }
    }
  }
  if ( pp5_variant == 0 ) {
      if ( (var.getValue("hgmd_class") == "DM" || var.getValue("hgmd_class") == "DM!"  )&& var.getValue("hgmd_disease").find("?Site") == std::string::npos  ) {
          var.addAcmgCategory("PP5");
          pp5_variant = 1;
      }
  }


  /* ----------- 5: BA ------------- */
  /* ----- 5a: BA1 ----- */
  int ba1_variant = 0;
  /* --- Allele frequency is >5% in Exome Sequencing Project, 1000 Genomes Project, or Exome Aggregation Consortium --- */
  std::vector<std::string> vPopulationSelected = {"evs_maf", "exac_global_af", "gnomad_ex_global_af", "gnomad_gen_global_af" };
  // BA1 is mutually exclusive with PM2 (https://www.clinicalgenome.org/site/assets/files/3677/clingen_variant-curation_sopv1.pdf)
  // only if not already a PM2 (mutually exclusive)
  // if is in BA1 exception do not calculate Benign values
  if (std::find (vBA1exceptions.begin(), vBA1exceptions.end(), var.getValue("var_name")) == vBA1exceptions.end()) {
      if ( ! frequent_filter(var, vPopulationSelected, 0.05f) ) {
          var.addAcmgCategory("BA1");
          ba1_variant = 1;
      }
  }

  /* ----------- 6: BS ------------- */
  /* ----- 6a: BS1 ----- */
  /* --- Allele frequency is greater than expected for disorder --- */
  // se <5% ma superiore a soglia MAF per la malattia
  if (std::find (vBA1exceptions.begin(), vBA1exceptions.end(), var.getValue("var_name")) == vBA1exceptions.end()) {
      if ( pm2_variant == 0 && ba1_variant == 0 ) {
          std::map<std::string, float>::iterator position = mBS1genes.find(gene_name);
          if (position != mBS1genes.end()) {
              float bs1_threshold = position->second;
              // il max sia comunque 1e-3 anche se diversamente specificato nel file map
              if ( bs1_threshold < 1E-3f ) {
                  if ( ! frequent_filter( var, vPopulationSelected, bs1_threshold ) ) {
                      var.addAcmgCategory("BS1");
                  }
              }
          }
          else {
            if ( ! frequent_filter(var, vPopulationSelected, 1E-3f) ) {
                var.addAcmgCategory("BS1");
            }
          }
      }
  }

  /* ----- 6b: BS2 ----- */
  /* --- Observed in a healthy adult individual for a recessive (homozygous), dominant (heterozygous), or X-linked (hemizygous) disorder, with full penetrance expected at an early age --- */
  /*
  if (std::find (vBA1exceptions.begin(), vBA1exceptions.end(), var.getValue("var_name")) == vBA1exceptions.end()) {
      if ( pm2_variant == 0 ) {
          if (bs2_filter (
              var.getValue("gnomad_ex_global_nhet"),
              var.getValue("gnomad_ex_global_nhomalt"),
              var.getValue("gnomad_ex_global_af"),
              var.getValue("gnomad_gen_global_nhet"),
              var.getValue("gnomad_gen_global_nhomalt"),
              var.getValue("gnomad_gen_global_af"),
              gene_inh_mode
          )) {
              var.addAcmgCategory("BS2");
          }
      }
  }
  */

  /* ----- 6c: BS3 ----- */
  /* --- Well-established in vitro or in vivo functional studies show no damaging effect on protein function or splicing --- */
  if (std::find (vBA1exceptions.begin(), vBA1exceptions.end(), var.getValue("var_name")) == vBA1exceptions.end()) {
      if ( var.getValue("benign_cv") == "1" ) {
          if ( bs3_filter( var.getValue("clinvar_clinsig"), var.getValue("clinvar_disease"), var.getValue("clinvar_clinrevstars") ) ) { var.addAcmgCategory("BS3"); }
      }
  }

  /* ----- 6d: BS4 ----- */
  /* --- Lack of segregation in affected members of a family --- */

  /* ----------- 7: BP ------------- */
  /* ----- 7a: BP1 ----- */
  /* --- Missense variant in a gene for which primarily truncating variants are known to cause disease --- */
  // if in a LOF known gene but is a missense variant (mutually exclusive with PP2)
  if (std::find (vBA1exceptions.begin(), vBA1exceptions.end(), var.getValue("var_name")) == vBA1exceptions.end()) {
      if ( pp2_variant  == 0 ) {
          if (std::find (vLOFgenes.begin(), vLOFgenes.end(), gene_name) != vLOFgenes.end()) {
              std::string effect = var.getValue("effect");
              std::string trap_score = var.getValue("trap_score");
              int lof_result = lof_filter( effect, trap_score );
              int trap_result = trap_filter( trap_score );
              // se NON-LOF in gene LOF guarda anche che non sia synonymous-intronic ma missense
              if ( lof_result == 0 && trap_result == 0 ) {
                  if ( !splice_filter( effect ) ) {
                      if ( synonymous_filter( effect ) == 0 ) {
                          var.addAcmgCategory("BP1");
                      }
                  }
              }
              // chiudo IF gene MISS
          }
      }
  }

  /* ----- 7b: BP2 ----- */
  /* --- Observed in trans with a pathogenic variant for a fully penetrant dominant gene/disorder or observed in cis with a pathogenic variant in any inheritance pattern --- */

  /* ----- 7c: BP3 ----- */
  /* --- In-frame deletions/insertions in a repetitive region without a known function --- */
  // ... find repeat regions of genes ...

  /* ----- 7d: BP4 ----- */
  /* --- Multiple lines of computational evidence suggest no impact on gene or gene product (conservation, evolutionary, splicing impact, etc.) --- */
  if (std::find (vBA1exceptions.begin(), vBA1exceptions.end(), var.getValue("var_name")) == vBA1exceptions.end()) {
      if ( pp3_var == 0 ) {
		  if ( pvs1_variant == 0 ) {
			  if ( pred_benign_filter( humdiv_cat, humvar_cat, revel, subrvis_exon_centile, subrvis_domain_centile, oe_evs, oe_exac, cadd_phred, sift_prediction ) ) {
				  var.addAcmgCategory("BP4");
			  }
		  }
      }
  }

  /* ----- 7e: BP5 ----- */
  /* --- Variant found in a case with an alternate molecular basis for disease --- */

  /* ----- 7f: BP6 ----- */
  /* --- Reputable source recently reports variant as benign, but the evidence is not available to the laboratory to perform an independent evaluation. --- */
  if (std::find (vBA1exceptions.begin(), vBA1exceptions.end(), var.getValue("var_name")) == vBA1exceptions.end()) {
      if ( var.getValue("benign_cv") == "1" ) {
          if ( ! bs3_filter( var.getValue("clinvar_clinsig"), var.getValue("clinvar_disease"), var.getValue("clinvar_clinrevstars") ) ) { var.addAcmgCategory("BP6"); }
      }
  }

  /* ----- 7g: BP7 ----- */
  /* --- A synonymous (silent) variant for which splicing prediction algorithms predict no impact to the splice consensus sequence nor the creation of a new splice site AND the nucleotide is not highly conserved --- */
  // ... to ADD some splicing algorithm score to predict it BENIGN
  if (std::find (vBA1exceptions.begin(), vBA1exceptions.end(), var.getValue("var_name")) == vBA1exceptions.end()) {
      std::string effect = var.getValue("effect");
      std::string ada_score = var.getValue("ada_score");
      std::string rf_score = var.getValue("rf_score");
      if ( ( synonymous_filter(effect) == 1 || splice_filter(effect) ) )  {
          if ( ada_score != "NA" || rf_score != "NA" ) {
              if ( !splice_patho_filter( ada_score, rf_score ) ) {
                  var.addAcmgCategory("BP7");
              }
          }
      }
  }
}

// questa funzione serve per ottenere la std::string della classe ACMG corrispondente
std::string get_acmg_class( Variant &var ) {
  // faccio un riassunto numerico delle sottocategorie
  int pvsCount = 0; int psCount = 0; int pmCount = 0; int ppCount = 0; int baCount = 0; int bsCount = 0; int bpCount = 0;
  // se almeno una categoria ACMG
  if (var.getAcmg().size() > 0 ) {
    // copio il vettore delle categorie ACMG
    std::vector<std::string> vAcmg = var.getAcmg();
    for (std::string sub : vAcmg) {
      if ( sub.find("PVS") != std::string::npos ) { pvsCount++; }
      if ( sub.find("PS") != std::string::npos ) { psCount++; }
      if ( sub.find("PM") != std::string::npos ) { pmCount++; }
      if ( sub.find("PP") != std::string::npos ) { ppCount++; }
      if ( sub.find("BA") != std::string::npos ) { baCount++; }
      if ( sub.find("BS") != std::string::npos ) { bsCount++; }
      if ( sub.find("BP") != std::string::npos ) { bpCount++; }
    }
  // chiudo IF ACMG vAcmg.size() > 1
  }
  // metto in map da aggiungere a oggetto var della classe Variant
  std::map<std::string, int> mAcmgSubclasses = {
    {"PVS", pvsCount},
    {"PS", psCount},
    {"PM", pmCount},
    {"PP", ppCount},
    {"BA", baCount},
    {"BS", bsCount},
    {"BP", bpCount}
  };
  // aggiungo la map di riassunto alla classe Variant
  var.addAcmgSubclassesMap(mAcmgSubclasses);

  // SVUOTO prima di mettere il nuovo valore
  // questo perche se no quando faccio il re-screening PM3 cerca di reinserirlo di nuovo
  if ( var.getValue("ACMG_bayesian") != "NA" ) {
      var.removeAttr("ACMG_bayesian");
  }
  /* --- BAYESIAN ACMG (articolo Hila) --- */
  // calcolo la Odds Pathogenicity (OP) della var basata sulle singole categorie
  float bay_exponent; float odds_patho; float post_probability;
  float odds_patho_PVS = 350.0;
  float prior_probablity = 0.10;
  if ( mAcmgSubclasses["BA"] == 0 ) {
    bay_exponent = ( (ppCount/8) + (pmCount/4) + (psCount/2) + (pvsCount/1) - (bpCount/8) - (bsCount/2) );
    odds_patho = pow(odds_patho_PVS, bay_exponent);
    post_probability = ( (odds_patho*prior_probablity) / ((odds_patho-1)*prior_probablity+1) );
//    std::stringstream stream;
//    stream << std::fixed << std::setprecision(2) << post_probability;
//    std::string s = stream.str();
//    std::cout << s << '\n';
  }
  else {
      post_probability = 10E-10;
  }
  // convert in string
  var.addAttr("ACMG_bayesian" , stringize_float(post_probability) );

  /* --- CLASS assignation --- */
  // https://www.acmg.net/docs/Standards_Guidelines_for_the_Interpretation_of_Sequence_Variants.pdf
  // Benign
  if ( mAcmgSubclasses["BA"] > 0 || mAcmgSubclasses["BS"] > 1 ) {
      var.addAcmgClass("B");
  }
  // Likely Benign
  else if (
    (mAcmgSubclasses["BS"] > 0 && mAcmgSubclasses["BP"] > 0) ||
    mAcmgSubclasses["BP"] > 1
  ) { var.addAcmgClass("LB"); }
  // Pathogenic
  if (bsCount < 2 ) {
      if ( pvsCount > 0 && psCount > 0 ) { var.addAcmgClass("P"); }
      else if ( pvsCount > 0 && pmCount > 1 ) { var.addAcmgClass("P"); }
      else if ( pvsCount > 0 && ppCount > 1 ) { var.addAcmgClass("P"); }
      else if ( pvsCount > 0 && pmCount > 0 && ppCount > 0 ) { var.addAcmgClass("P"); }
      else if ( pvsCount > 1 ) { var.addAcmgClass("P"); }
      else if ( pvsCount > 0 && pmCount > 2 ) { var.addAcmgClass("P"); }
      else if ( pvsCount > 0 && pmCount > 1 && ppCount > 2 ) { var.addAcmgClass("P"); }
      else if ( pvsCount > 0 && pmCount > 0 && ppCount > 3 ) { var.addAcmgClass("P"); }
      // Likely Pathogenic
      else if (
        (mAcmgSubclasses["PVS"] > 0 && mAcmgSubclasses["PM"] > 0) ||
        (mAcmgSubclasses["PS"] > 0 && ( mAcmgSubclasses["PM"] == 1 || mAcmgSubclasses["PM"] == 2)) ||
        (mAcmgSubclasses["PS"] > 0 && mAcmgSubclasses["PP"] > 1) ||
        (mAcmgSubclasses["PM"] > 2) ||
        (mAcmgSubclasses["PM"] > 1 && mAcmgSubclasses["PP"] > 1) ||
        (mAcmgSubclasses["PM"] >0 && mAcmgSubclasses["PP"] > 3)
      ) { var.addAcmgClass("LP"); }
  }

  std::string result;
  // se in NESSUNA CLASSE o se in PIU di 1 CLASSE = VUS
  if ( var.getAcmgClass().empty() || var.getAcmgClass().size() > 1 ) { result = "US"; }
  // altrimenti assegna ad ACMG il valore della CLASSE univoca trovata
  else { result = var.getAcmgClass()[0]; }

  return result;
}


/* ---------- ACMG
  -----------------------------------------------------------------------------
  ---------- funzione iterativa per fargli passare TUTTE le varianti ----------
  ---------- a farne associazione per le ACMG guidelines categories -----------
  -----------------------------------------------------------------------------
--------- */
void filter_variants_acmg (
  std::vector<Variant> &vVarClass,
  std::vector<Gene> &vGeneClass,
  std::vector<std::string> &vGeneClassNames,
  int &parents_aff
) {

    //caricamento dei file necessari (per non farli caricare ogni volta per ogni variante)
    // prendo in forma di vettore i file con i geni LOF e con i geni MISS
    // NB: se lo mando direttamente da solo da cartella asilo2 questa va modificata togliendo: app/dependencies/analisi/ da inizio
    std::string fileLOFgenes = "/usr/src/app/dependencies/databases/genes_lof_miss/clinvar_pvs1_genes.txt";
    std::string fileMISSgenes = "/usr/src/app/dependencies/databases/genes_lof_miss/clinvar_pp2_genes.txt";
    std::vector<std::string> vLOFgenes = getVectorFile( fileLOFgenes );
    std::vector<std::string> vMISSgenes = getVectorFile( fileMISSgenes );
    std::cout << "-------> found " << vLOFgenes.size() << " genes PVS1" << '\n';
    std::cout << "-------> found " << vMISSgenes.size() << " genes PP2" << '\n';
    // riprendo la MAP da file HOTSPOT region (che per ora sono gli esoni e le flank 60bp di ClinVar)
    std::string fileHOTSPOTucsc = "/usr/src/app/dependencies/databases/genes_lof_miss/ucsc.exons.patho_REV.hotspot.txt";
    std::map< std::string, std::vector<std::vector<std::string>> > mHtospotUcsc  = get_ucsc_hotspot( fileHOTSPOTucsc );
    std::cout << "-------> found " << mHtospotUcsc.size() << " HOTSPOTs (esoni)" << '\n';
    // questo lo ottengo per il CSV con /home/enrico/columbia/clinvar/db_jun_2019/launch_get_PM1.sh
    std::string fileHOTSPOTcv = "/usr/src/app/dependencies/databases/genes_lof_miss/baylopLMM_PLP_notFound_merged.pm1";
    std::vector<std::string> vHotspotCv  = getVectorFile( fileHOTSPOTcv );
    std::cout << "-------> found " << vHotspotCv.size() << " HOTSPOTs (ClinVar 60bp)" << '\n';

    // isolo la lista dei GENI dalla MAP HOTSPOT (gene_name = KEY)
    std::vector<std::string> vGeniMapHotspot;
    for( auto const &m : mHtospotUcsc ) {
        // NO DUPLICATI
        if (std::find (vGeniMapHotspot.begin(), vGeniMapHotspot.end(), m.first) == vGeniMapHotspot.end()) {
            vGeniMapHotspot.push_back(m.first);
        }
    }

    //prendo il file con i geni e le relative HGVS.p di geni PATHO ClinVar
    std::string fileHGVSPgenes = "/usr/src/app/dependencies/databases/genes_lof_miss/genefile_hgvsp_patho_miss_REV.txt";
    std::vector<std::vector<std::vector<std::string>>> vvvHGVSPgenes = get_vvFile( fileHGVSPgenes, '\t', ',' );
    std::cout << "-------> found " << vvvHGVSPgenes.size() << " HGVS.p mutations" << '\n';

    //prendo il file con le varianti nella ClinGen BA1 exception list
    std::string fileBA1exception = "/usr/src/app/dependencies/databases/genes_lof_miss/clingen_BA1_exception.list";
    std::vector<std::string> vBA1exceptions = getVectorFile( fileBA1exception );
    std::cout << "-------> found " << vBA1exceptions.size() << " BA1 ClinGen exceptions" << '\n';

    // questo creato con /media/enrico/vep_disk/columbia/emerge/extract_clinvar_max_patho_af.py
    std::string fileBS1genes = "/usr/src/app/dependencies/databases/genes_lof_miss/clinvar_20200203_bs1.tsv";
    std::vector<std::vector<std::string>> vvBS1genes = getVVectorFile( fileBS1genes, '\t' );
    // transformo il VV in una map (gene : threshold BS1)
    std::map<std::string, float> mBS1genes;
    for ( std::vector<std::string> v : vvBS1genes )
    {
        mBS1genes.insert(std::make_pair(v[0],  std::stof(v[1]) ));
    }
    std::cout << "-------> found " << mBS1genes.size() << " BS1 gene thresholds" << '\n';

    std::vector<std::string> vHGVSPgenes = take_vvv_column( vvvHGVSPgenes, 1 );
    std::map<std::string, std::vector<std::string>> mHGVSPpatho = create_vvv_map( vvvHGVSPgenes, 1, 2 );
    std::map<std::string, std::vector<std::string>> mHGVSPpos = create_vvv_map( vvvHGVSPgenes, 1, 3 );

    // apro loop sulle varianti (tutte)
    for (auto &var : vVarClass)
    {
        // inserisco le subclassi ACMG nel oggetto Var della classe Variant
        add_acmg_subclasses( var, vGeneClass, vGeneClassNames, vLOFgenes, vMISSgenes, vvvHGVSPgenes, vHGVSPgenes, mHGVSPpatho, mHGVSPpos, parents_aff, mHtospotUcsc, vHotspotCv, vGeniMapHotspot, mBS1genes, vBA1exceptions );

        // ottengo la classe ACMG e la metto come attribute ACMG della classe Variant
        var.addAttr("ACMG", get_acmg_class( var ));
        // chiudo loop sulle singole VAR della classe
    }
}




/* ----------
  -----------------------------------------------------------------------------
  ---------- funzione iterativa per ASSEGNARE PRIORITA FINALE alla var --------
  -----------------------------------------------------------------------------
--------- */
void assign_priority ( std::vector<Variant> &vVarClass ) {

  // assegno la rispettiva priority alle mie classi in questa map
  std::map<std::string, int> mPriorityMine = {
    {"patho0",18},
    {"lof0",17},
    {"trap0",12},
    {"aa0",11},
    {"flank0",10},
    {"pred0",9},
    {"patho5",16},
    {"lof5",15},
    {"trap5",8},
    {"aa5",7},
    {"flank5",6},
    {"pred5",5},
    {"pathoR",14},
    {"lofR",13},
    {"trapR",4},
    {"aaR",3},
    {"flankR",2},
    {"predR",1}
  };
  // e poi assegno priority a gruppo ACMG classification
  std::map<std::string, int> mPriorityACMG = {
    {"P",9},
    {"LP",7},
    {"US",5},
    {"LB",3},
    {"B",1}
  };

  // apro loop sulle varianti
  for (Variant &var : vVarClass ) {
    if ( dom_pass(var) && var.getValue("ACMG") != "B" && var.getValue("ACMG") != "LB") {
      // prendo la prioita della mia classe e di classe ACMG
      int prior_mine = mPriorityMine[var.getValue("dom")];
      int prior_acmg = mPriorityACMG[var.getValue("ACMG")];
      int prior_tot = prior_mine*prior_acmg;
      var.addAttr("priority", std::to_string(prior_tot));
    // chiudo IF DOM and ACMG assigned, se NON metto priority a 0
} else { var.addAttr("priority","0"); }
  // chiudo loop su singolo oggetto VAR
  }
}





/* ----------
  -----------------------------------------------------------------------------
  ---------- funzione iterativa per filtraggio varianti-pz RECESSIVE ----------
  -----------------------------------------------------------------------------
--------- */

/* ---
  NB: SOLO se gene AR devo fare tutto questo sotto per capire HOM e COMPHET
  infatti se la var sarà finita in almeno una categoria e il gene è AD
  mi basterà poi prendere il vettore con tutti i pazienti .getPzAll()
--- */
/*
    quindi alla fine di questa funzione mi trovo con la VAR che ha una std::map<std::string, std::vector<std::string>> mPzComphet
    in cui sono inseriti i pazienti in cui ha resistito in COMPHET (K) e le relative var con cui associata (V) (.getPzComphet())
*/
void filter_variants_recessive (
  std::vector<Patient> &vPatientClass,
  std::vector<Gene> &vGeneClass,
  std::vector<std::string> &vGeneClassNames,
  std::vector<Variant> &vVarClass,
  std::vector<std::string> &vVarClassNames,
  int comphet_analysis
) {

  // apro loop sui pazienti
  for ( Patient &pz : vPatientClass ) {

    std::string pz_name = pz.getValue("pz_name");

    // contatore per sapere quante recessive si sono salvate per ogni gene del paziente
    std::map< std::string, std::vector<std::string>> mGeneRecessiveResults;

    // apro loop sulla classe pzVars delle var del paziente
    for ( pzVars &varPz : pz.getPzVarVector() ) {

      // prendo nome, zigosi, denovo e gene della variante in questione
      std::string var_name = varPz.getValue("var_name");
      std::string var_zigosi = varPz.getValue("gt");
      std::string var_denovo = varPz.getValue("denovo");
      std::string var_gene_name = varPz.getValue("gene_name_correct");

      // riprendo corrispettiva VAR in vettore classe Variant
      std::ptrdiff_t posVar = find( vVarClassNames.begin(), vVarClassNames.end(), var_name) - vVarClassNames.begin();
      Variant &var_class = vVarClass[ posVar ];

      // riprendo corrispettivo GENE in vettore classe Gene
      std::ptrdiff_t posGene = find( vGeneClassNames.begin(), vGeneClassNames.end(), var_gene_name) - vGeneClassNames.begin();
      Gene &var_gene_class = vGeneClass[ posGene ];

      // capisco se inh mode del gene in questione AR o XLR e nel caso compilo i modelli rec_hom e rec_comphet
      if ( check_recessive_inh_mode ( var_gene_class.getInhModes() ) || check_xlinked_recessive_inh_mode ( var_gene_class.getInhModes() ) ) {

        // SOLO SE in almeno un modello di interesse
        if ( dom_pass(var_class) ) {
          if ( var_zigosi == "hom" ) {
            // riprendo la categoria in cui finito
            std::string category = var_class.getValue("dom");
            // SOLO se NON già presente (es per paziente trovato HOM in precedenza) inserisco il nome della categoria in rec_hom
            if ( var_class.getValue("rec_hom") == "NA" && category != "NA" ) {
                var_class.addAttr( "rec_hom", category );
            }
            var_class.addPzHom( pz_name );
          }

          // se invece in HET devo prima capire quante varianti in HET si sono salvate per ogni gene in ogni paziente
          else {

              // SOLO se richiesta ANALISI COMPHET (NO FILE COMPHET TRIO)
              if ( comphet_analysis == 1 ) {
                // se gene key non ancora in MAP
                if (mGeneRecessiveResults.find(var_gene_name) == mGeneRecessiveResults.end()) {
                  // creo un vettore con var_name per inserirlo nella map con nome gene e li inserisco
                  std::vector<std::string> vVarName = {var_name};
                  mGeneRecessiveResults.insert(std::make_pair(var_gene_name, vVarName));
                // se gene key gia in map
                }
                else {
                    if (std::find (mGeneRecessiveResults[var_gene_name].begin(), mGeneRecessiveResults[var_gene_name].end(), var_name) == mGeneRecessiveResults[var_gene_name].end()) {
                      // appendo la var al vettore .second corrispondente al .first gene
                      mGeneRecessiveResults[var_gene_name].push_back(var_name);
                    }
                }
            } // chiudo IF richiesta analisi COMPHET
            // se invece ho FILE COMPHET TRIO allora devo anche controllare che la VAR siano tra quelle con FLAG TRIO_file in attributes
            else {
                if ( var_class.getValue("comphet_flag") == "TRIO_file" ) {
                    // se gene key non ancora in MAP
                    if (mGeneRecessiveResults.find(var_gene_name) == mGeneRecessiveResults.end()) {
                        // creo un vettore con var_name per inserirlo nella map con nome gene e li inserisco
                        std::vector<std::string> vVarName = {var_name};
                        mGeneRecessiveResults.insert(std::make_pair(var_gene_name, vVarName));
                        // se gene key gia in map
                    }
                    else {
                        if (std::find (mGeneRecessiveResults[var_gene_name].begin(), mGeneRecessiveResults[var_gene_name].end(), var_name) == mGeneRecessiveResults[var_gene_name].end()) {
                            // appendo la var al vettore .second corrispondente al .first gene
                            mGeneRecessiveResults[var_gene_name].push_back(var_name);
                        }
                    }
                } // chiudo if FLAG TRIO_file
            } // chiudo ELSE (cioe se TRIO COMPHET File presente)

          } // chiudo ELSE (quindi if het)
        } // chiudo IF dom_pass
      } // chiudo IF AR-XLR
      // se gene DOM devo comunque vedere se sia HOM e metterlo in relativa categoria
      else {
        if ( var_zigosi == "hom" ) {
          // SOLO se NON già presente (es per paziente trovato HOM in precedenza) inserisco il nome della categoria in rec_hom
          var_class.addPzHom( pz_name );
        }
      }

    // chiudo loop su varPz
    }

    /*
        quindi se analisi comphet NON richiesta NON ci sara nessuno in questa MAP mGeneRecessiveResults e tutta la parte COMPHET successiva NON succedera
    */
    // apro quindi loop su mGeneRecessiveResults
    for ( auto const &gene : mGeneRecessiveResults ) {
      // se almeno 2 var salvate per il gene
      if ( gene.second.size() > 1 ) {
        // apro loop su queste var per aggiungerci le relative caratt agli obj classe var
        for ( std::string var_comphet : gene.second ) {
          // riprendo corrispettiva VAR in COMPHET dal vettore classe Variant per aggiungerci il modello rec_comphet
          // aggiungendo anche alla map mPzComphet della var itself il pz (.first) e l elenco var in comphet (.second)
          std::ptrdiff_t posVarCompHet = find( vVarClassNames.begin(), vVarClassNames.end(), var_comphet) - vVarClassNames.begin();
          Variant &var_class_comphet = vVarClass[ posVarCompHet ];
          // se NON già presente (es. se finita qui per un paziente precedente)
          if ( var_class_comphet.getValue("rec") == "NA" ) {
              var_class_comphet.addAttr("rec","comphet");
          }
          var_class_comphet.addPzComphet(pz_name, gene.first, gene.second);
        // chiudo loop var in comphet del gene
        }
      // chiudo IF > 1
      }
    // chiudo loop su mGeneRecessiveResults
    }

  // chiudo loop sui pazienti
  }
}

/*
    -------------------------------
    ------ RIFINITURE FINALI ------
    -------------------------------
*/

// questa controlla che ci sia almeno un paziente in ogni categoria (SINGLE/COMP HET and HOM)
// e nel caso la aggiunge come attributo alla relativa Variant (altrimenti NA)
void check_and_add_patients( Variant &var )
{
    // NON dovrebbe essercene nessuna Var che NON sia in almeno un paziente
    if ( var.getPzAll().size() > 0 ) {
        var.addAttr("pz_all",stringize_vector(var.getPzAll(), ","));
    }
    else {
        std::cout << "!!! ATTENTION !!! NO patient found in Var: " << var.getValue("var_name") << '\n';
        var.addAttr("pz_all","NA");
    }
    if ( var.getPzSingleHet().size() > 0 ) {
        var.addAttr("pz_single_het",stringize_vector(var.getPzSingleHet(), ","));
    }
    else {
        var.addAttr("pz_single_het","NA");
    }
    if ( var.getPzHom().size() > 0 ) {
        var.addAttr("pz_hom",stringize_vector(var.getPzHom(), ","));
    }
    else {
        var.addAttr("pz_hom","NA");
    }
    if ( var.getPzComphet().size() > 0 ) {
        var.addAttr("pz_comp_het",stringize_map_map_str_V( var.getPzComphet(), "|", ":", ",", ">", "!"));
    }
    else {
        var.addAttr("pz_comp_het","NA");
    }
}












/*
    questa invece la funzione per andare a mettere l'annotazione PM3
    che deve essere fatta DOPO lo screening delle altre var e il collezionamento COMPHET
    essendo relativa a essere in TRANS con var PATHO (solo per geni REC)
*/
void acmg_pm3_screening( std::vector<Variant> &vVarClass, std::vector<std::string> &vVarClassNames )
{
    // mi serve anche un V per tenere conto di quali VAR gia screenate per non ripeter (e dare errore inserimento)
    std::vector<std::string> vPM3screened;
    // qui voglio controllare i geni finiti in COMPHET per vedere se ci sono VAR P/LP in modo che le altre facciano un PM
    for ( Variant &var : vVarClass )
    {
        // se NON gia screenata
        if (std::find (vPM3screened.begin(), vPM3screened.end(), var.getValue("var_name")) != vPM3screened.end()) {
            continue;
        }
        vPM3screened.push_back( var.getValue("var_name"));

        // mi serve riprendere la posizione per vedere che siano in TRANS (almeno 50bp)
        std::vector<std::string> vNome; vectorize_column(var.getValue("var_name"), vNome, '-');
        int pos = std::stoi( vNome[1] );

        // inserisco i pazienti nei relativi campi come attribute della Var
        check_and_add_patients( var );
        // se almeno un paziente in COMPHET
        if ( var.getValue("pz_comp_het") != "NA" ) {
            //std::cout << var.getValue("pz_comp_het") << '\n';
            // apro loop su MAP< std::string (pz_name) ,std::vector<std::string> (vars_in_comphet) >
            for ( auto const &gene : var.getPzComphet() )
            {
                for ( auto const &pz_vars : gene.second )
                {
                    // rinomino pz_name e vettore VAR in comphet in quel paziente
                    std::string pz_name = pz_vars.first;
                    std::vector<std::string> vPzVars = pz_vars.second;
                    // apro loop sulle var in comphet di quel paziente
                    for ( std::string var_name : vPzVars )
                    {
                        // riprendo corrispettiva VAR
                        std::ptrdiff_t posVarCompHet = find( vVarClassNames.begin(), vVarClassNames.end(), var_name) - vVarClassNames.begin();
                        Variant &var_class_comphet = vVarClass[ posVarCompHet ];
                        if ( var_class_comphet.getValue("ACMG") == "P" || var_class_comphet.getValue("ACMG") == "LP" ) {

                            // mi serve riprendere la posizione per vedere che siano in TRANS (almeno 50bp)
                            std::vector<std::string> vNomeComphet; vectorize_column(var_class_comphet.getValue("var_name"), vNomeComphet, '-');
                            int posComphet = std::stoi( vNomeComphet[1] );
                            // controllo che siano in TRANS (50 bp di distanza) e non stessa INDEL x es
                            if ( std::abs( pos - posComphet ) > 50 ) {
                                // aggiungo il PM3 e poi devo ricalcolare e nel caso cambiare classe ACMG della VAR
                                var.addAcmgCategory("PM3");
								var.emptyAcmgClass();
								var.removeAttr("ACMG");
								var.addAttr("ACMG", get_acmg_class( var ));
								/*
								if ( var.getValue("var_name") == "1-5924454-GGTGGCAATCC-G" ) {
									std::cout << var.getValue("var_name") << '\n';
									std::cout << var.getValue("ACMG") << '\n';
									std::vector<std::string> vAcmg = var.getAcmg();
									for (std::string sub : vAcmg) {
											std::cout << sub << '\n';
									}
								}
								*/
                            } // chiudo if distanza almeno 50bp
                        } // chiudo IF VAR in comphet ACMG P/LP
                    } // chiudo loop sulle VAR in comphet del paziente
                }
            } // chiudo loop sulla MAP comphet della VAR
        } // chiudo IF VAR COMPHET NON-NA
    } // chiudo loop su VAR della classe Variant
    vPM3screened.clear(); // svuoto vettore delle gia corrette che non mi serve piu
}


// questa vede se la var finita come P/LP e se priority < 11 la rimette come US
void acmg_priority_screening( std::vector<Variant> &vVarClass, std::vector<std::string> &vVarClassNames )
{
    for ( Variant &var : vVarClass )
    {
      if ( var.getValue("ACMG") == "P" || var.getValue("ACMG") == "LP" ) {
        if ( var.getValue("priority") != "NA" ) {
          if ( std::stoi(var.getValue("priority")) < 0 ) {
            var.removeAttr("ACMG");
            var.addAttr("ACMG", "US");
          }
        }
      }
    }
}







/*
    questa invece la funzione per andare a riprendere gli esoni HOTSPOT
    come:  std::map< std::string, std::vector<std::vector<std::string>> >
    con una map con tutti i GENI come KEY
    e come VALUE un VV dove ogni sotto-V ha come valori: chr, start, stop
*/
std::map< std::string, std::vector<std::vector<std::string>> > get_ucsc_hotspot (
    std::string nomeFileInput
) {

    auto vvFileInput =  getVVectorFile( nomeFileInput, '\t' );

    // inizializzo la map in cui mettero il file input come <GENE , VV del <chr, start, stop> >
    std::map< std::string, std::vector<std::vector<std::string>> > mFileInput;
    for ( auto v : vvFileInput )
    {
        std::string geneName = v[0];
        std::string chrPos = v[2];
        // isolo i valori di chr e V start-end della posizione
        std::string chr = vectorize_column_v2( chrPos, ':' )[0];
        // ricreo il V chr, start, stop da inserire come elemento del VV degli esoni fi ogni gene nella MAP
        std::vector<std::string> vExon = { chr, vectorize_column_v2( vectorize_column_v2( chrPos, ':' )[1], '-' )[0], vectorize_column_v2( vectorize_column_v2( chrPos, ':' )[1], '-' )[1] };
        // se gene NON ancora in map lo inserisco con VV come value
        if ( mFileInput.find(geneName) == mFileInput.end() ) {
            std::vector<std::vector<std::string>> vvExons = { vExon };
            if ( !mFileInput.insert( std::make_pair( geneName, vvExons ) ).second ) {
                std::cout << "!!! ATTENTION !!! insert failed for: " << geneName << " in MAP mFileInput (attribute already present)" << '\n';
            }
        }
        // se invece geneName gia presente deve appendere il vExon al corrispondente vvExons
        else {
            mFileInput[ geneName ].push_back( vExon );
        }

        /*
            quindi a questo punto mi trovo con una map con tutti i GENI come KEY
            e come VALUE un VV dove ogni sotto-V ha come valori: chr, start, stop

        for ( auto const &m : mFileInput )
        {
            std::cout << "--> " << m.first << '\n';
            for ( auto vv : m.second )
            {
                for ( auto v : vv )
                {
                    std::cout << '\t' << v << ' ';
                }
                std::cout << "" << '\n';
            }
            std::cout << "" << '\n';
        }
        */
//        std::cout << geneName << ": " << vPosStartEnd[0] << "-" << vPosStartEnd[1] << '\n';
    }

    return mFileInput;
}
