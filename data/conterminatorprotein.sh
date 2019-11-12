#!/bin/sh -e
# Sequence search workflow script
fail() {
    echo "Error: $1"
    exit 1
}

notExists() {
	[ ! -f "$1" ]
}

#pre processing
[ -z "$MMSEQS" ] && echo "Please set the environment variable \$MMSEQS to your MMSEQS binary." && exit 1;
# check amount of input variables
[ "$#" -ne 4 ] && echo "Please provide <sequence.fasta> <mappingFile> <result> <tmp>" && exit 1;
# check if files exists
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[ ! -f "$2" ] &&  echo "$2 not found!" && exit 1;
[   -f "$3" ] &&  echo "$3 exists already!" && exit 1;
[ ! -d "$4" ] &&  echo "tmp directory $4 not found!" && mkdir -p "$4";

TMP_PATH="$4"


if notExists "$TMP_PATH/sequencedb"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createdb "$1" "$TMP_PATH/sequencedb" ${CREATEDB} \
        || fail "createdb step died"
fi

if notExists "$TMP_PATH/sequencedb_mapping"; then
    # shellcheck disable=SC2086
    "$MMSEQS" createtaxdb "$TMP_PATH/sequencedb" "${TMP_PATH}/createtaxdb" --tax-mapping-file "${TAXMAPPINGFILE}" ${ONLYVERBOSITY} \
        || fail "createtaxdb step died"
fi

if notExists "$TMP_PATH/clu.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" linclust "$TMP_PATH/sequencedb" "$TMP_PATH/clu" "$TMP_PATH/linclust" ${LINCLUST_PAR} \
        || fail "kmermatcher step died"
fi

if notExists "$TMP_PATH/aln.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" align "$TMP_PATH/sequencedb" "$TMP_PATH/sequencedb" "$TMP_PATH/clu" "$TMP_PATH/aln" ${ALN_PAR} \
        || fail "kmermatcher step died"
fi

if notExists "$TMP_PATH/conterm_aln.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" crosstaxonfilter "$TMP_PATH/sequencedb" "$TMP_PATH/aln" "$TMP_PATH/conterm_aln" ${CROSSTAXA_PAR} \
        || fail "kmermatcher step died"
fi

if notExists  "$TMP_PATH/conterm_aln_stats.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" createstats "$TMP_PATH/sequencedb" "$TMP_PATH/conterm_aln" "$TMP_PATH/conterm_aln_stats" ${CREATESTATS_PAR} \
        || fail "createtsv step died"
fi

if notExists "${2}_stats"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" prefixid "$TMP_PATH/conterm_aln_stats" "${3}_stats" --threads 1 --tsv \
        || fail "prefixid step 1  died"
fi

if notExists "$TMP_PATH/conterm_aln_all.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" createallreport "$TMP_PATH/sequencedb" "$TMP_PATH/conterm_aln" "$TMP_PATH/conterm_aln_all" ${CREATESTATS_PAR} \
        || fail "createtsv step died"
fi

if notExists "${2}_all"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" prefixid "$TMP_PATH/conterm_aln_all" "${3}_all" --threads 1 --tsv \
        || fail "prefixid step 2 died"
fi

if [ -n "$REMOVE_TMP" ]; then
  echo "Remove temporary files"
    echo "Remove temporary files"
  $MMSEQS rmdb "$TMP_PATH/conterm_aln_all"
  $MMSEQS rmdb "$TMP_PATH/conterm_aln_stats"
  $MMSEQS rmdb "$TMP_PATH/sequencedb"
  $MMSEQS rmdb "$TMP_PATH/sequencedb_h"
  $MMSEQS rmdb "$TMP_PATH/conterm_aln"
  $MMSEQS rmdb "$TMP_PATH/aln"
  $MMSEQS rmdb "$TMP_PATH/clu"
fi

