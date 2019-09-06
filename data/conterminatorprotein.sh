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
[ "$#" -ne 3 ] && echo "Please provide <sequenceDB> <outDB> <tmp>" && exit 1;
# check if files exists
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[   -f "$2" ] &&  echo "$2 exists already!" && exit 1;
[ ! -d "$3" ] &&  echo "tmp directory $3 not found!" && mkdir -p "$3";


DB="$1"
TMP_PATH="$3"


if notExists "$TMP_PATH/clu.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" linclust "$DB" "$TMP_PATH/clu" "$TMP_PATH/linclust" ${LINCLUST_PAR} \
        || fail "kmermatcher step died"
fi

if notExists "$TMP_PATH/aln.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" align "$DB" "$DB" "$TMP_PATH/clu" "$TMP_PATH/aln"  ${ALN_PAR} \
        || fail "kmermatcher step died"
fi

if notExists "$TMP_PATH/conterm_aln.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" crosstaxonfilter "$DB" "$TMP_PATH/aln" "$TMP_PATH/conterm_aln" ${CROSSTAXA_PAR} \
        || fail "kmermatcher step died"
fi

if notExists  "$TMP_PATH/conterm_aln_stats.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" createstats "${DB}" "$TMP_PATH/conterm_aln" "$TMP_PATH/conterm_aln_stats" ${CREATESTATS_PAR} \
        || fail "createtsv step died"
fi

if notExists "${2}_stats"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" prefixid "$TMP_PATH/conterm_aln_stats" "${2}_stats" --threads 1 --tsv \
        || fail "prefixid step 1  died"
fi

if notExists "$TMP_PATH/conterm_aln_all.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" createallreport "${DB}" "$TMP_PATH/conterm_aln" "$TMP_PATH/conterm_aln_all" ${CREATESTATS_PAR} \
        || fail "createtsv step died"
fi

if notExists "${2}_all"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" prefixid "$TMP_PATH/conterm_aln_all" "${2}_all" --threads 1 --tsv \
        || fail "prefixid step 2 died"
fi

if [ -n "$REMOVE_TMP" ]; then
  echo "Remove temporary files"
  rm -f "$3/aln_distance"    "$4/aln_distance.index"
fi

