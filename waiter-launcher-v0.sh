#!/usr/bin/env bash

DIR='/INPUT'
AWK='/usr/src/app/header-converter.awk'

while true; do
  COUNT=$( ls -1 $DIR/*.asilo_input 2>/dev/null | wc -l )
  if [ $COUNT -gt 0 ]; then
    INPUT_FILE=$( ls -t $DIR/*.asilo_input 2>/dev/null | head -n1 )
    INPUT="${DIR}/$( awk -F'\t' '{ print $1 }' $INPUT_FILE )"
    GENELIST="$( awk -F'\t' '{ print $2 }' $INPUT_FILE )"
    FILTER_AF=$( awk -F'\t' '{ print $3 }' $INPUT_FILE )
    CSV="$( dirname $INPUT )/$( basename $INPUT .vcf ).csv"
    AWK_COMMAND="gawk -f $AWK $INPUT > $CSV"
    ASILO_COMMAND_RM="rm -rf ${DIR}/analisi_result"
    ASILO_COMMAND_EXE="./asilo.1.0
                            ${CSV}
                            NONE
                            ${DIR}/${GENELIST}
                            ${DIR}/analisi_result
                            prova.var
                            prova.pz
                            prova.pzvar
                            prova.gene
                            ${FILTER_AF}
                            0
                            0
                            gnomad_ex_global_af,gnomad_ex_afr_af,gnomad_gen_global_af,gnomad_gen_afr_af
                            5"
     echo -e " --------------------- INPUTS --------------------- "
     echo -e "  - INPUT:\t$INPUT_FILE"
     echo -e "    - VCF:\t$INPUT"
     echo -e "    - CSV:\t$CSV"
     echo -e "    - GENES:\t$GENELIST"
     echo -e "    - MAF:\t$FILTER_AF"
     echo -e " -------------------------------------------------- "
     echo -e "   1. converting VCF2CSV ..."
     eval $AWK_COMMAND
     wait
     sleep 2
     if [ -d "${DIR}/analisi_result" ]; then
       echo -e "   2a. removing previous analisi_result DIR ..."
       eval $ASILO_COMMAND_RM
       wait
       sleep 2
     fi
     echo -e "   3. annotating with ASILO ..."
     eval $ASILO_COMMAND_EXE
     wait
     sleep 2
     echo -e "   4. create ASILO-NEW flag ..."
     date '+%Y/%m/%d %H:%M:%S' > "${DIR}/ASILO.NEW"
     echo -e " ---------------------- END ---------------------- "
     mv $INPUT_FILE "$( dirname $INPUT_FILE )/$( basename $INPUT_FILE .input ).output"
     echo "    - done: $OUTPUT"
  else
    echo "  - NO INPUT file"
    sleep 10
  fi
done


























exit 0
