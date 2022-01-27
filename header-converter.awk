BEGIN{
  FS = "\t";

  ### IGNORE GT VALUES ARRAY
  IGNORE_GT = ".,./.,.|.,0/0,0|0"
  split( IGNORE_GT, IGNORE_GT_ARRAY, /[,]/ );
  for ( g=1; g<=length( IGNORE_GT_ARRAY ); g++ )
  {
    IGNORE_GT_K_ARRAY[ IGNORE_GT_ARRAY[g] ] = g;
  }

  ### SAMPLE ARRAY
  split("", SAMPLE_ARRAY);
  split("", SAMPLE_K_ARRAY);

  ### VEP_HEADER ARRAY
  split("", VEP_HEADER_ARRAY);
  split("", VEP_HEADER_K_ARRAY);

  ### to store VEP array for each variant
  split( "", VEP_VARIANTS_ARRAY );

  ### FORMAT COLUMN
  FORMAT_COLUMN = 0;

  ### INFO COLUMN
  INFO_COLUMN = 0;

  ### START a VAR_COUNTER outside
  VAR_COUNTER = 0;

  ### store last sample column number
  SAMPLE_LAST_COL = 0;

  ### define CLINVAR ClinRevStat to ClinRevStars conversion
  CLINVAR_CLINREVSTA2CLINREVSTARS_ARRAY["practice_guideline"] = 4;
  CLINVAR_CLINREVSTA2CLINREVSTARS_ARRAY["reviewed_by_expert_panel"] = 3;
  CLINVAR_CLINREVSTA2CLINREVSTARS_ARRAY["criteria_provided&_multiple_submitters&_no_conflicts"] = 2;
  CLINVAR_CLINREVSTA2CLINREVSTARS_ARRAY["criteria_provided&_conflicting interpretations"] = 1;
  CLINVAR_CLINREVSTA2CLINREVSTARS_ARRAY["criteria_provided&_single_submitter"] = 1;
  CLINVAR_CLINREVSTA2CLINREVSTARS_ARRAY["no_assertion_for_the_individual_variant"] = 0;
  CLINVAR_CLINREVSTA2CLINREVSTARS_ARRAY["no_assertion_criteria_provided"] = 0;
  CLINVAR_CLINREVSTA2CLINREVSTARS_ARRAY["no_assertion_provided"] = 0;

  ### store HEADER conversion as array
  # VEP_HEADER_CONVERSION_ARRAY["gnomAD_exomes_AC"] = "gnomad_ex_global_an";
  VEP_HEADER_CONVERSION_ARRAY["gnomAD_exomes_POPMAX_AF"] = "gnomAD Exome global_AF";
  VEP_HEADER_CONVERSION_ARRAY["gnomAD_exomes_AN"] = "gnomAD Exome global_AN";
  VEP_HEADER_CONVERSION_ARRAY["gnomAD_exomes_AFR_AF"] = "gnomAD Exome afr_AF";
  VEP_HEADER_CONVERSION_ARRAY["gnomAD_exomes_AMR_AF"] = "gnomAD Exome amr_AF";
  VEP_HEADER_CONVERSION_ARRAY["gnomAD_exomes_ASJ_AF"] = "gnomAD Exome asj_AF";
  VEP_HEADER_CONVERSION_ARRAY["gnomAD_exomes_EAS_AF"] = "gnomAD Exome eas_AF";
  VEP_HEADER_CONVERSION_ARRAY["gnomAD_exomes_FIN_AF"] = "gnomAD Exome fin_AF";
  VEP_HEADER_CONVERSION_ARRAY["gnomAD_exomes_NFE_AF"] = "gnomAD Exome nfe_AF";
  VEP_HEADER_CONVERSION_ARRAY["gnomAD_exomes_SAS_AF"] = "gnomAD Exome sas_AF";
  # VEP_HEADER_CONVERSION_ARRAY["gnomAD_exomes_POPMAX_nhomalt"] = "";
  # VEP_HEADER_CONVERSION_ARRAY["gnomAD_exomes_controls_AC"] = "";
  VEP_HEADER_CONVERSION_ARRAY["gnomAD_exomes_controls_AF"] = "gnomAD Exome controls_AF";
  VEP_HEADER_CONVERSION_ARRAY["gnomAD_exomes_controls_AN"] = "gnomAD Exome controls_AN";

  VEP_HEADER_CONVERSION_ARRAY["gnomAD_genomes_POPMAX_AF"] = "gnomAD Genome global_AF";
  VEP_HEADER_CONVERSION_ARRAY["gnomAD_genomes_AC"] = "gnomAD Genome global_AN";
  # VEP_HEADER_CONVERSION_ARRAY["gnomAD_genomes_POPMAX_nhomalt"] = "";
  VEP_HEADER_CONVERSION_ARRAY["gnomAD_genomes_AFR_AF"] = "gnomAD Genome afr_AF";
  VEP_HEADER_CONVERSION_ARRAY["gnomAD_genomes_AMR_AF"] = "gnomAD Genome amr_AF";
  VEP_HEADER_CONVERSION_ARRAY["gnomAD_genomes_ASJ_AF"] = "gnomAD Genome asj_AF";
  VEP_HEADER_CONVERSION_ARRAY["gnomAD_genomes_EAS_AF"] = "gnomAD Genome eas_AF";
  VEP_HEADER_CONVERSION_ARRAY["gnomAD_genomes_FIN_AF"] = "gnomAD Genome fin_AF";
  VEP_HEADER_CONVERSION_ARRAY["gnomAD_genomes_NFE_AF"] = "gnomAD Genome nfe_AF";
  VEP_HEADER_CONVERSION_ARRAY["gnomAD_genomes_AN"] = "gnomAD Genome controls_AN";
  VEP_HEADER_CONVERSION_ARRAY["MAX_AF"] = "gnomAD Genome controls_AF";   ### this is a trick

  # VEP_HEADER_CONVERSION_ARRAY["ExAC_AC"] = "";
  VEP_HEADER_CONVERSION_ARRAY["ExAC_AF"] = "ExAC global af";
  VEP_HEADER_CONVERSION_ARRAY["ExAC_AFR_AF"] = "ExAC afr af";
  VEP_HEADER_CONVERSION_ARRAY["ExAC_AMR_AF"] = "ExAC amr af";
  VEP_HEADER_CONVERSION_ARRAY["ExAC_EAS_AF"] = "ExAC eas af";
  VEP_HEADER_CONVERSION_ARRAY["ExAC_FIN_AF"] = "ExAC fin af";
  VEP_HEADER_CONVERSION_ARRAY["ExAC_NFE_AF"] = "ExAC nfe af";
  VEP_HEADER_CONVERSION_ARRAY["ExAC_SAS_AF"] = "ExAC sas af";
  # VEP_HEADER_CONVERSION_ARRAY["ExAC_Adj_AF"] = "";
  # VEP_HEADER_CONVERSION_ARRAY["ExAC_AFR_AC"] = "";
  # VEP_HEADER_CONVERSION_ARRAY["ExAC_AMR_AC"] = "";
  # VEP_HEADER_CONVERSION_ARRAY["ExAC_Adj_AC"] = "";
  # VEP_HEADER_CONVERSION_ARRAY["ExAC_EAS_AC"] = "";
  # VEP_HEADER_CONVERSION_ARRAY["ExAC_FIN_AC"] = "";
  # VEP_HEADER_CONVERSION_ARRAY["ExAC_NFE_AC"] = "";
  # VEP_HEADER_CONVERSION_ARRAY["ExAC_SAS_AC"] = "";
  VEP_HEADER_CONVERSION_ARRAY["1000Gp3_AC"] = "1000G AC";
  VEP_HEADER_CONVERSION_ARRAY["1000Gp3_AF"] = "1000G AF";
  VEP_HEADER_CONVERSION_ARRAY["1000Gp3_AFR_AC"] = "1000G Afr AC";
  VEP_HEADER_CONVERSION_ARRAY["1000Gp3_AFR_AF"] = "1000G Afr AF";
  VEP_HEADER_CONVERSION_ARRAY["1000Gp3_AMR_AC"] = "1000G Amr AC";
  VEP_HEADER_CONVERSION_ARRAY["1000Gp3_AMR_AF"] = "1000G Amr AF";
  VEP_HEADER_CONVERSION_ARRAY["1000Gp3_EAS_AC"] = "1000G Eas AC";
  VEP_HEADER_CONVERSION_ARRAY["1000Gp3_EAS_AF"] = "1000G Eas AF";
  VEP_HEADER_CONVERSION_ARRAY["1000Gp3_EUR_AC"] = "1000G Eur AC";
  VEP_HEADER_CONVERSION_ARRAY["1000Gp3_EUR_AF"] = "1000G Eur AF";
  VEP_HEADER_CONVERSION_ARRAY["1000Gp3_SAS_AC"] = "1000G Sas AC";
  VEP_HEADER_CONVERSION_ARRAY["1000Gp3_SAS_AF"] = "1000G Sas AF";
  VEP_HEADER_CONVERSION_ARRAY["ClinVar"] = "ClinVar";
  VEP_HEADER_CONVERSION_ARRAY["ClinVar_CLNSIG"] = "ClinVar ClinSig";
  VEP_HEADER_CONVERSION_ARRAY["ClinVar_CLNREVSTAT"] = "ClinVar ClinRevStat";
  VEP_HEADER_CONVERSION_ARRAY["ClinVar_CLNREVSTARS"] = "ClinVar ClinRevStar";
  VEP_HEADER_CONVERSION_ARRAY["ClinVar_CLNDN"] = "ClinVar DiseaseName";
  VEP_HEADER_CONVERSION_ARRAY["ClinVar_CLNDISDB"] = "ClinVar DiseaseDB";
  VEP_HEADER_CONVERSION_ARRAY["UK10K_AF"] = "UK10k AF";
  VEP_HEADER_CONVERSION_ARRAY["ESP6500_AA_AF"] = "ESP6500 AA AF";
  VEP_HEADER_CONVERSION_ARRAY["ESP6500_EA_AF"] = "ESP6500 EA AC";
  VEP_HEADER_CONVERSION_ARRAY["REVEL_score"] = "REVEL";
  VEP_HEADER_CONVERSION_ARRAY["Polyphen2_HDIV_score"] = "Polyphen Humdiv Score (CCDS)";
  VEP_HEADER_CONVERSION_ARRAY["Polyphen2_HDIV_pred"] = "Polyphen Humdiv Prediction (CCDS)";
  VEP_HEADER_CONVERSION_ARRAY["Polyphen2_HVAR_pred"] = "Polyphen Humvar Score (CCDS)";
  VEP_HEADER_CONVERSION_ARRAY["Polyphen2_HVAR_score"] = "Polyphen Humvar Prediction (CCDS)";
  VEP_HEADER_CONVERSION_ARRAY["Consequence"] = "Effect";
  VEP_HEADER_CONVERSION_ARRAY["HGVSc"] = "HGVS_c";
  VEP_HEADER_CONVERSION_ARRAY["HGVSp"] = "HGVS_p";
  # VEP_HEADER_CONVERSION_ARRAY["SYMBOL"] = "Gene Name";
  VEP_HEADER_CONVERSION_ARRAY["SYMBOL"] = "UpToDate Gene Name";
  # VEP_HEADER_CONVERSION_ARRAY["GERP++_RS"] = "";
  VEP_HEADER_CONVERSION_ARRAY["SIFT"] = "SIFT Prediction";
  VEP_HEADER_CONVERSION_ARRAY["ada_score"] = "ADA score";
  VEP_HEADER_CONVERSION_ARRAY["rf_score"] = "RF score";
  VEP_HEADER_CONVERSION_ARRAY["CADD_PHRED"] = "CADD phred";
  VEP_HEADER_CONVERSION_ARRAY["CADD_RAW"] = "CADD raw";

  # VEP_HEADER_CONVERSION_ARRAY[""] = "";


  ### specify printing ARRAY
  PRINT_ARRAY_STRING="VAR_NAME,SAMPLE_NAME,SAMPLE_GT_CONVERTED,AD_REF_VALUE,AD_ALT_VALUE,DP_VALUE,GENE_NAME";
  split( PRINT_ARRAY_STRING, PRINT_ARRAY, /[,]/ );

  ### sample-specific header conversion
  SAMPLE_HEADER_CONVERSION_ARRAY["VAR_NAME"] = "Variant ID";
  SAMPLE_HEADER_CONVERSION_ARRAY["SAMPLE_NAME"] = "Sample Name";
  SAMPLE_HEADER_CONVERSION_ARRAY["SAMPLE_GT_CONVERTED"] = "GT";
  SAMPLE_HEADER_CONVERSION_ARRAY["AD_REF_VALUE"] = "AD REF";
  SAMPLE_HEADER_CONVERSION_ARRAY["AD_ALT_VALUE"] = "AD ALT";
  SAMPLE_HEADER_CONVERSION_ARRAY["DP_VALUE"] = "DP";
  SAMPLE_HEADER_CONVERSION_ARRAY["GENE_NAME"] = "Gene Name";

  ### print HEADER
  for ( ss = 1; ss <= length( PRINT_ARRAY ); ss++ )
  {
    SAMPLE_HEADER_VALUE = PRINT_ARRAY[ss];
    if ( ss == 1 ) {
      printf "%s", SAMPLE_HEADER_CONVERSION_ARRAY[SAMPLE_HEADER_VALUE];
    }
    else {
      printf ",%s", SAMPLE_HEADER_CONVERSION_ARRAY[SAMPLE_HEADER_VALUE];
    }
  }
  for ( vh in VEP_HEADER_CONVERSION_ARRAY ) {
    printf ",%s", VEP_HEADER_CONVERSION_ARRAY[vh];
  }
  printf "\n";
  ### END print HEADER


}


{

  ### if HEADER
  if ( substr($0,1,1) == "#" ) {
    # print;
    if ( substr($0,1,6) == "#CHROM") {
      ### extract sample names and relative column number
      for ( i=1; i<=NF; i++ )
      {
        ### extract FORMAT column number (the one before samples)
        if ( $i == "FORMAT" ) {
          FORMAT_COLUMN = i;
        }
        ### extract INFO column number (the one containing VEP annotation)
        else if ( $i == "INFO" ) {
          INFO_COLUMN = i;
        }
        ### if FORMAT column passed start memorize samples
        if ( FORMAT_COLUMN > 0 ) {
          if ( i > FORMAT_COLUMN ) {
            ### sample COLUMN = NAME
            SAMPLE_ARRAY[i] = $i;
            ### sample NAME = COLUMN
            SAMPLE_K_ARRAY[$i] = i;
          }
        }
      }
      SAMPLE_LAST_COL = FORMAT_COLUMN + length(SAMPLE_ARRAY);
      ### DEVELOPMENT
      if (0) {
        # printf "FORMAT column: %s\n", FORMAT_COLUMN;
        SAMPLE_LAST_COL = FORMAT_COLUMN + length(SAMPLE_ARRAY);
        for ( ss=1; ss<=SAMPLE_LAST_COL; ss++ )
        {
          printf "  - SAMPLE: %s --> column: %s\n", SAMPLE_ARRAY[ss], ss;
        }
      }
    }
    else if ( substr($0,1,14) == "##INFO=<ID=CSQ" ) {
      split( $0, VEP_HEADER_LINE_ARRAY, "Format: " );
      gsub( /">/, "", VEP_HEADER_LINE_ARRAY[2] );
      split( VEP_HEADER_LINE_ARRAY[2], VEP_HEADER_ARRAY_TMP, /[|]/ );
      for ( h=1; h<=length(VEP_HEADER_ARRAY_TMP); h++ )
      {
        VEP_HEADER_ARRAY[h] = VEP_HEADER_ARRAY_TMP[h];
        VEP_HEADER_K_ARRAY[VEP_HEADER_ARRAY_TMP[h]] = h;
        if ( VEP_HEADER_ARRAY[h] in VEP_HEADER_CONVERSION_ARRAY ) {
          # printf "  - VEP field %s: %s --> %s\n", h, VEP_HEADER_ARRAY[h], VEP_HEADER_CONVERSION_ARRAY[VEP_HEADER_ARRAY[h]];
        }
      }
    }
  }
  ### if NOT HEADER
  else {

    ### initialize VAR-specific values
    GT_POS = 0;
    AD_POS = 0;
    DP_POS = 0;
    split("",VARIANT_VEP_ARRAY);

    ### use to DEVELOPMENT to print only some
    if ( VAR_COUNTER > 0 ) {

      ### extract VAR NAME
      VAR_NAME = $1"-"$2"-"$4"-"$5
      # gsub( "chr", "", VAR_NAME );
      ### check for multiallelic
      if ( $5 ~ /,/ ) {
        printf "  - line: %s, identified MULTIALLELIC VAR: %s\n", NR, VAR_NAME;
      }
      ### split FORMAT column
      split( $FORMAT_COLUMN, FORMAT_ARRAY, /[:]/ );
      for ( f=1; f<=length(FORMAT_ARRAY); f++ )
      {
        if ( FORMAT_ARRAY[f] == "GT" ) {
          GT_POS = f;
        }
        else if ( FORMAT_ARRAY[f] == "AD" ) {
          AD_POS = f;
        }
        else if ( FORMAT_ARRAY[f] == "DP" ) {
          DP_POS = f;
        }
      }

      ### split INFO column and get VEP annotation
      split( $INFO_COLUMN, INFO_ARRAY, /[;=]/ );
      for ( z=1; z<=length(INFO_ARRAY); z++ )
      {
        if ( INFO_ARRAY[z] == "CSQ" ) {
          VEP_ANNOTATION = INFO_ARRAY[z+1];
          VEP_ANNOTATION_ARRAY_N = split( VEP_ANNOTATION, VEP_ANNOTATION_ARRAY, /[|]/ );

          ### 1p. this is to store, before LOOP and IF ELSE IF loops, the transcript position in dbNSFP
          #   NB: not necessarily they will agree, often dbNSFP is empty and VEP is not
          ### extract POSITION and VALUES in VEP of PolyPhen and Polyphen2_HVAR_score (PolyPhen is the VEP default and corresponds to HVAR)
          DBNSFP_TRANSCRIPT_POSITION = 1; ### default value
          POLYPHEN_VEP_HVAR_VALUE = VEP_ANNOTATION_ARRAY[ VEP_HEADER_K_ARRAY["PolyPhen"] ];
          if ( POLYPHEN_VEP_HVAR_VALUE != "" ) {
            split( POLYPHEN_VEP_HVAR_VALUE, POLYPHEN_VEP_HVAR_VALUE_ARRAY, /[()]/ );
            POLYPHEN_VEP_HVAR_VALUE_SCORE = POLYPHEN_VEP_HVAR_VALUE_ARRAY[2];
          }
          POLYPHEN_DBNSFP_HVAR_VALUE = VEP_ANNOTATION_ARRAY[ VEP_HEADER_K_ARRAY["Polyphen2_HVAR_score"] ];
          # printf " ==> polyphen VEP value: %s - polyphen dbNSFP value: %s\n", POLYPHEN_VEP_HVAR_VALUE, POLYPHEN_DBNSFP_HVAR_VALUE;
          split( POLYPHEN_DBNSFP_HVAR_VALUE, POLYPHEN_DBNSFP_HVAR_VALUE_ARRAY, /[&]/ );
          for ( ppv = 1; ppv <= length( POLYPHEN_DBNSFP_HVAR_VALUE_ARRAY ); ppv++ ) {
            if ( POLYPHEN_DBNSFP_HVAR_VALUE_ARRAY[ppv] == POLYPHEN_VEP_HVAR_VALUE_SCORE ) {
              DBNSFP_TRANSCRIPT_POSITION = ppv;
              # printf "  --> original dbNSFP Polyphen2_HVAR_score transcript position: %s\n", DBNSFP_TRANSCRIPT_POSITION;
              break;
            }
          }

            ### LOOP over VEP annotation fields
            for ( v=1; v<=VEP_ANNOTATION_ARRAY_N; v++ )
            {
              ### get the correspondent FIELD name
              VEP_FIELD = VEP_HEADER_ARRAY[v];
              ### store in VARIANT-specific VARIANT_VEP_ARRAY with K=VEP_FIELD : V=VEP_FIELD_VALUE
              if ( VEP_ANNOTATION_ARRAY[v] != "" ) {
                VARIANT_VEP_ARRAY[VEP_FIELD] = VEP_ANNOTATION_ARRAY[v];
              }
              else {
                VARIANT_VEP_ARRAY[VEP_FIELD] = "NA";
              }

              ### ARRANGE VEP values that need it
              ##  1. arrange HGVSc
              if ( VEP_FIELD == "HGVSc" ) {
                split("", HGVSC_VALUE_ARRAY);    ### re-initialize VEP_VALUE_ARRAY to ensure is empty
                HGVSC_VALUE = "";
                HGVSC_VALUE_CONVERTED = "";
                HGVSC_VALUE = VEP_ANNOTATION_ARRAY[v];
                if ( HGVSC_VALUE != "" ) {
                  split( HGVSC_VALUE, HGVSC_VALUE_ARRAY, /[:]/ );
                  HGVSC_VALUE_CONVERTED = HGVSC_VALUE_ARRAY[2];
                }
                else {
                  HGVSC_VALUE_CONVERTED = "NA";
                }
                # printf "HGVSc: %s\n", HGVSC_VALUE_CONVERTED;
                VARIANT_VEP_ARRAY[VEP_FIELD] = HGVSC_VALUE_CONVERTED;
              }
              ##  2. arrange HGVSp
              else if ( VEP_FIELD == "HGVSp" ) {
                split("", HGVSP_VALUE_ARRAY);    ### re-initialize VEP_VALUE_ARRAY to ensure is empty
                HGVSP_VALUE = "";
                HGVSP_VALUE_CONVERTED = "";
                HGVSP_VALUE = VEP_ANNOTATION_ARRAY[v];
                if ( HGVSP_VALUE != "" ) {
                  split( HGVSP_VALUE, HGVSP_VALUE_ARRAY, /[:]/ );
                  HGVSP_VALUE_CONVERTED = HGVSP_VALUE_ARRAY[2];
                  gsub( "%3D", "=", HGVSP_VALUE_CONVERTED );
                }
                else {
                  HGVSP_VALUE_CONVERTED = "NA";
                }
                # printf "  - HGVSp: %s\n", HGVSP_VALUE_CONVERTED;
                VARIANT_VEP_ARRAY[VEP_FIELD] = HGVSP_VALUE_CONVERTED;
              }
              ##  3. extract Polyphen HVAR
              else if ( VEP_FIELD == "PolyPhen" ) {
                split("", PPHVAR_VALUE_ARRAY);    ### re-initialize VEP_VALUE_ARRAY to ensure is empty
                PPHVAR_VALUE_PRED = "NA";
                PPHVAR_VALUE_SCORE = "NA";
                PPHVAR_VALUE = VEP_ANNOTATION_ARRAY[v];
                if ( PPHVAR_VALUE != "" ) {
                  split( PPHVAR_VALUE, PPHVAR_VALUE_ARRAY, /[()]/ );
                  PPHVAR_VALUE_PRED = PPHVAR_VALUE_ARRAY[1];
                  PPHVAR_VALUE_SCORE = PPHVAR_VALUE_ARRAY[2];
                }
                else {
                  PPHVAR_VALUE_PRED = "NA";
                  PPHVAR_VALUE_SCORE = "NA";
                }
                # printf "  --> splitted PolyPhen HVAR: %s - %s\n", PPHVAR_VALUE_PRED, PPHVAR_VALUE_SCORE;
                # printf "  - HGVSp: %s\n", HGVSP_VALUE_CONVERTED;
                if ( PPHVAR_VALUE_PRED == "." ) {
                  PPHVAR_VALUE_PRED = "NA";
                }
                else if ( PPHVAR_VALUE_PRED == "probably_damaging" ) {
                  PPHVAR_VALUE_PRED = "probably";
                }
                else if ( PPHVAR_VALUE_PRED == "possibly_damaging" ) {
                  PPHVAR_VALUE_PRED = "possibly";
                }
              }
              ##  4. extract Polyphen HDIV PRED
              else if ( VEP_FIELD == "Polyphen2_HDIV_pred" ) {
                split("", PPHDIV_VALUE_ARRAY);    ### re-initialize VEP_VALUE_ARRAY to ensure is empty
                PPHDIV_VALUE_PRED = "";
                PPHDIV_VALUE = VEP_ANNOTATION_ARRAY[v];
                # printf "  --> original dbNSFP Polyphen2_HDIV_pred: %s\n", PPHDIV_VALUE;
                if ( PPHDIV_VALUE != "" ) {
                  split( PPHDIV_VALUE, PPHDIV_VALUE_ARRAY, /[&]/ );
                  PPHDIV_VALUE_PRED = PPHDIV_VALUE_ARRAY[DBNSFP_TRANSCRIPT_POSITION];
                  if ( PPHDIV_VALUE_PRED == "." ) {
                    PPHDIV_VALUE_PRED = "NA";
                  }
                  else if ( PPHDIV_VALUE_PRED == "D" ) {
                    PPHDIV_VALUE_PRED = "probably";
                  }
                  else if ( PPHDIV_VALUE_PRED == "P" ) {
                    PPHDIV_VALUE_PRED = "possibly";
                  }
                  else if ( PPHDIV_VALUE_PRED == "B" ) {
                    PPHDIV_VALUE_PRED = "benign";
                  }
                }
                else {
                  PPHDIV_VALUE_PRED = "NA";
                }
              }
              ##  5. extract Polyphen HDIV SCORE
              else if ( VEP_FIELD == "Polyphen2_HDIV_score" ) {
                split("", PPHDIV_SCORE_VALUE_ARRAY);    ### re-initialize VEP_VALUE_ARRAY to ensure is empty
                PPHDIV_VALUE_SCORE = "";
                PPHDIV_VALUE_SCORE_ORIGINAL = VEP_ANNOTATION_ARRAY[v];
                # printf "  --> original dbNSFP Polyphen2_HDIV_score: %s\n", PPHDIV_VALUE_SCORE_ORIGINAL;
                if ( PPHDIV_VALUE_SCORE_ORIGINAL != "" ) {
                  split( PPHDIV_VALUE_SCORE_ORIGINAL, PPHDIV_SCORE_VALUE_ARRAY, /[&]/ );
                  PPHDIV_VALUE_SCORE = PPHDIV_SCORE_VALUE_ARRAY[DBNSFP_TRANSCRIPT_POSITION];
                  if ( PPHDIV_VALUE_SCORE == "." ) {
                    PPHDIV_VALUE_SCORE = "NA";
                  }
                }
                else {
                  PPHDIV_VALUE_SCORE = "NA";
                }
              }
              ##  6. extract SIFT
              else if ( VEP_FIELD == "SIFT" ) {
                split("", SIFT_VALUE_ARRAY);    ### re-initialize VEP_VALUE_ARRAY to ensure is empty
                SIFT_VALUE_PRED = "NA";
                SIFT_VALUE_SCORE = "NA";
                SIFT_VALUE_ORIGINAL = VEP_ANNOTATION_ARRAY[v];
                # printf "  --> original dbNSFP Polyphen2_HDIV_score: %s\n", PPHDIV_VALUE_SCORE_ORIGINAL;
                if ( SIFT_VALUE_ORIGINAL != "" ) {
                  split( SIFT_VALUE_ORIGINAL, SIFT_VALUE_ORIGINAL_ARRAY, /[()]/ );
                  SIFT_VALUE_PRED = SIFT_VALUE_ORIGINAL_ARRAY[1];
                  SIFT_VALUE_SCORE = SIFT_VALUE_ORIGINAL_ARRAY[2];
                  if ( SIFT_VALUE_PRED == "." ) {
                    SIFT_VALUE_PRED = "NA";
                  }
                }
                else {
                  SIFT_VALUE_PRED = "NA";
                }
              }
              ##  7. extract ClinVar_CLNREVSTARS
              else if ( VEP_FIELD == "ClinVar_CLNREVSTAT" ) {
                CLINVAR_CLNREVSTAT_VALUE_ORIGINAL = VEP_ANNOTATION_ARRAY[v];
                CLINVAR_CLNREVSTARS = "NA";
                if ( CLINVAR_CLNREVSTAT_VALUE_ORIGINAL != "" ) {
                  if ( CLINVAR_CLNREVSTAT_VALUE_ORIGINAL != "NA" ) {
                    if ( CLINVAR_CLNREVSTAT_VALUE_ORIGINAL in CLINVAR_CLINREVSTA2CLINREVSTARS_ARRAY ) {
                      CLINVAR_CLNREVSTARS = CLINVAR_CLINREVSTA2CLINREVSTARS_ARRAY[CLINVAR_CLNREVSTAT_VALUE_ORIGINAL];
                    }
                  }
                }
              }
              ### END ARRANGE VEP values that need it

            }   ### STOP LOOP over VEP annotation fields

            ### insert Polyphen HVAR
            VARIANT_VEP_ARRAY["Polyphen2_HVAR_pred"] = PPHVAR_VALUE_PRED;
            VARIANT_VEP_ARRAY["Polyphen2_HVAR_score"] = PPHVAR_VALUE_SCORE;
            ### insert Polyphen HDIV
            VARIANT_VEP_ARRAY["Polyphen2_HDIV_pred"] = PPHDIV_VALUE_PRED;
            VARIANT_VEP_ARRAY["Polyphen2_HDIV_score"] = PPHDIV_VALUE_SCORE;
            ### insert SIFT
            VARIANT_VEP_ARRAY["SIFT"] = SIFT_VALUE_PRED;
            ### insert ClinVar_CLNREVSTARS
            VARIANT_VEP_ARRAY["ClinVar_CLNREVSTARS"] = CLINVAR_CLNREVSTARS;
            ### duplicate gene name
            GENE_NAME = "'"VARIANT_VEP_ARRAY["SYMBOL"]"'";


            ### LOOP over VARIANT-specific VEP-values
            if (0) {
              for ( vv=1; vv<=length(VEP_HEADER_ARRAY); vv++ )
              {
                ### store VARIANT VEP array
                VEP_VARIANTS_ARRAY[ VAR_NAME ][ VEP_HEADER_ARRAY[vv] ] = VARIANT_VEP_ARRAY[VEP_HEADER_ARRAY[vv]];
                # printf "      --> field N %s, %s: %s\n", vv, VEP_HEADER_ARRAY[vv], VARIANT_VEP_ARRAY[VEP_HEADER_ARRAY[vv]];
              }   ### STOP LOOP over VARIANT-specific VEP-values
            }

        }
      }
      if (0) {
        printf "  --> VEP annotation: %s\n", VEP_ANNOTATION;
        printf "    --> VEP NF: %s\n", VEP_ANNOTATION_ARRAY_N;
      }

      ### loop over samples
      for ( sc=FORMAT_COLUMN+1; sc<=SAMPLE_LAST_COL; sc++ )
      {

        # printf "- CURRENT COLUMN: %s\n", sc;

        ### extract sample name
        SAMPLE_NAME = SAMPLE_ARRAY[sc];

        ### RE-initialize SAMPLE-SPECIFIC values
        AD_REF_POS = 0;
        AD_ALT_POS = 0;
        AD_REF_VALUE = 0;
        AD_ALT_VALUE = 0;
        SAMPLE_GT = "";
        SAMPLE_GT_CONVERTED = "";

        ### LOOP over SAMPLE_INFO
        split( $sc, SAMPLE_INFO_ARRAY, /[:]/);
        for ( x=1; x<=length(SAMPLE_INFO_ARRAY); x++ )
        {

          ### extract SAMPLE GT
          if ( x == GT_POS ) {
            SAMPLE_GT = SAMPLE_INFO_ARRAY[x];

            ### if SAMPLE NOT NEGATIVE for the VAR
            if ( SAMPLE_GT in IGNORE_GT_K_ARRAY == 0 ) {

              if ( SAMPLE_GT == "0/1" || SAMPLE_GT == "0|1" ) {
                SAMPLE_GT_CONVERTED = "het";
              }
              else if ( SAMPLE_GT == "1/1" || SAMPLE_GT == "1|1" ) {
                SAMPLE_GT_CONVERTED = "hom";
              }
              ### identify where the GT has 0 or 1/2/3 for subsequent AD extraction
              split( SAMPLE_GT, SAMPLE_GT_ARRAY, /[|/]/);
              GT_REF_VALUE = SAMPLE_GT_ARRAY[1];
              GT_ALT_VALUE = SAMPLE_GT_ARRAY[2];
              ### if NOT hom
              if ( SAMPLE_GT_CONVERTED != "hom" ) {
                ### if first value is 0 then store it as REF AND the value of the second as ALT
                if ( GT_REF_VALUE == 0 ) {
                  AD_REF_POS = 1;
                  if ( GT_ALT_VALUE != 0 ) {
                    AD_ALT_POS = GT_ALT_VALUE+1;
                  }
                }
                ### if 1st != 0 AND GT_ALT_VALUE == 0 store 1st value as ALT and 2nd as REF
                else {
                  if ( GT_ALT_VALUE == 0 ) {
                    AD_ALT_POS = GT_REF_VALUE+1;
                    AD_REF_POS = 2;
                  }
                }
              }
              ### if HOM
              else {
                AD_REF_POS = 1;
                AD_ALT_POS = 2;
              }
            }   ### CLOSE if SAMPLE NOT NEGATIVE for the VAR
          }
          ### extract AD
          else if ( x == AD_POS ) {
            ### if SAMPLE NOT NEGATIVE for the VAR
            if ( SAMPLE_GT in IGNORE_GT_K_ARRAY == 0 ) {
              split( SAMPLE_INFO_ARRAY[x], SAMPLE_AD_ARRAY, /[,]/ );
              for ( a=1; a<=length(SAMPLE_AD_ARRAY); a++ )
              {
                AD_REF_VALUE = SAMPLE_AD_ARRAY[AD_REF_POS];
                AD_ALT_VALUE = SAMPLE_AD_ARRAY[AD_ALT_POS];
              }
            }   ### CLOSE if SAMPLE NOT NEGATIVE for the VAR
          }
          ### extract DP
          else if ( x == DP_POS ) {
            ### if SAMPLE NOT NEGATIVE for the VAR
            if ( SAMPLE_GT in IGNORE_GT_K_ARRAY == 0 ) {
              DP_VALUE = SAMPLE_INFO_ARRAY[x];
            }
          }
        } ### STOP LOOP over SAMPLE_INFO

        ### if SAMPLE NOT NEGATIVE for the VAR
        if ( SAMPLE_GT in IGNORE_GT_K_ARRAY == 0 ) {
          # printf "  - SAMPLE: %s, VAR: %s, GT: %s, SAMPLE_GT_CONVERTED: %s, AD_REF_POS: %s, AD_ALT_POS: %s, AD REF: %s, AD ALT: %s\n", SAMPLE_NAME, VAR_NAME, SAMPLE_GT, SAMPLE_GT_CONVERTED, AD_REF_POS, AD_ALT_POS, AD_REF_VALUE, AD_ALT_VALUE;
        }   ### CLOSE if SAMPLE NOT NEGATIVE for the VAR

        ### print results as specified in PRINT_ARRAY only for NON-HOMREF
        if ( SAMPLE_GT_CONVERTED != "" ) {
          for ( p=1; p<=length(PRINT_ARRAY); p++ )
          {
            gsub( ",", "_", PRINT_ARRAY[p] );   ### remove any possible comma from generated CSV
            if ( p == 1 ) {
              printf "%s", SYMTAB[PRINT_ARRAY[p]];
            }
            else {
              printf ",%s", SYMTAB[PRINT_ARRAY[p]];
            }
            # printf " --> value: %s : %s\n", PRINT_ARRAY[p], SYMTAB[PRINT_ARRAY[p]];
          }
          for ( pv in VEP_HEADER_CONVERSION_ARRAY )
          {
            gsub( ",", "_", VARIANT_VEP_ARRAY[pv] );   ### remove any possible comma from generated CSV
            ### print VARIANT VEP array
            printf ",%s", VARIANT_VEP_ARRAY[pv];
            # printf "      --> field N %s, %s: %s\n", pv, VEP_HEADER_ARRAY[pv], VARIANT_VEP_ARRAY[VEP_HEADER_ARRAY[pv]];
          }   ### STOP LOOP over VARIANT-specific VEP-values
          printf "\n";
        }


      } ### CLOSE SAMPLEs LOOP

    }
    VAR_COUNTER++;
  }

}


END{}
