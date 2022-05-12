#!/bin/env bash

source /data/UHTS/2backup/tools/covid_pipeline/covpi.conf
source ${SRC_DIR}/src/utils.sh

source $1 
if [ "x${PROJECT}" == "x" ]; then
    echo "Missing analysis configuration file"
    exit 1
fi

HUG_SPSP_XLSX="$2"

export PERL5LIB=${PERL5LIB}:${SRC_DIR}/src/perl5lib

cd ${OUT_DIR}

###############################################################
# GET CONSENSI FROM HUG SPSP xlsx file
###############################################################

# consensi OK
dir=${PROJECT}_SPSP
echo_info "Prepare $dir"
rm -rf ${dir} ${dir}.tar.gz
mkdir -p $dir
${SRC_DIR}/src/parseXlsxSPSP.pl "../${HUG_SPSP_XLSX}" |  while read l; do
    if [ -f ${PROJECT}_GISAID_OK/${l} ]; then 
        cp -v ${PROJECT}_GISAID_OK/${l} ${dir}/${l}
    fi
done

gpg -q --import --fingerprint $PUBKEY
gpg -q --batch --yes --always-trust -o ${dir}/${HUG_SPSP_XLSX}.gpg -e -r $RECIPIENT ../${HUG_SPSP_XLSX}
tar cvf - $dir | gzip -9c > ${dir}.tar.gz
sha256sum ${dir}.tar.gz > ${dir}.tar.gz.sha256
cd ..

