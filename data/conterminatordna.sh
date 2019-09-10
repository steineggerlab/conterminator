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


if notExists "$TMP_PATH/db_rev_split"; then
    # shellcheck disable=SC2086
    "$MMSEQS" splitsequence "$DB" "$TMP_PATH/db_rev_split"  ${SPLITSEQ_PAR} \
        || fail "splitsequence step died"
fi

if notExists "$TMP_PATH/pref"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" kmermatcher "$TMP_PATH/db_rev_split"  "$TMP_PATH/pref"  ${KMERMATCHER_PAR} \
        || fail "kmermatcher step died"
fi

if notExists "$TMP_PATH/aln.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" rescorediagonal "$TMP_PATH/db_rev_split" "$TMP_PATH/db_rev_split" "$TMP_PATH/pref" "$TMP_PATH/aln" ${RESCORE_DIAGONAL_PAR} \
        || fail "rescorediagonal step died"
fi

if notExists "$TMP_PATH/aln_offset.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" offsetalignment "${DB}" "$TMP_PATH/db_rev_split" "${DB}" "$TMP_PATH/db_rev_split"  "$TMP_PATH/aln" "$TMP_PATH/aln_offset" ${OFFSETALIGNMENT_PAR} \
        || fail "rescorediagonal step died"
fi

if notExists "$TMP_PATH/contam_aln.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" extractalignments "${DB}"  "$TMP_PATH/aln_offset" "$TMP_PATH/contam_aln" ${EXTRACTALIGNMENTS_PAR} \
        || fail "extractalignment step died"
fi

if notExists "$TMP_PATH/contam_region.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" extractalignedregion "${DB}" "${DB}" "$TMP_PATH/contam_aln" "$TMP_PATH/contam_region" ${THREADS_PAR} \
        || fail "extractalignedregion step died"
fi

cp "$TMP_PATH/contam_region.index" "$TMP_PATH/contam_region.old.index"
awk '{print NR"\t"$2"\t"$3}' "$TMP_PATH/contam_region.index" > "$TMP_PATH/contam_region.new.index"
mv "$TMP_PATH/contam_region.new.index" "$TMP_PATH/contam_region.index"

if notExists "$TMP_PATH/contam_region_rev.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" extractframes "$TMP_PATH/contam_region" "$TMP_PATH/contam_region_rev" ${EXTRACT_FRAMES_PAR}  \
        || fail "Extractframes died"
fi

#TODO
if notExists "$TMP_PATH/contam_region_pref.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" prefilter "$TMP_PATH/db_rev_split" "$TMP_PATH/contam_region_rev" "$TMP_PATH/contam_region_pref" ${PREFILTER_PAR} \
        || fail "createdb step died"
fi

if notExists "$TMP_PATH/contam_region_aln.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" rescorediagonal "$TMP_PATH/db_rev_split" "$TMP_PATH/contam_region_rev" "$TMP_PATH/contam_region_pref" "$TMP_PATH/contam_region_aln" ${RESCORE_DIAGONAL_PAR} \
        || fail "rescorediagonal2 step died"
fi

if notExists "${TMP_PATH}/contam_region_aln_swap.dbtype"; then
     # shellcheck disable=SC2086
    "$MMSEQS" swapresults "$TMP_PATH/db_rev_split" "$TMP_PATH/contam_region_rev" "${TMP_PATH}/contam_region_aln" "${TMP_PATH}/contam_region_aln_swap" ${SWAP_PAR} \
        || fail "Swapresults pref died"
fi

if notExists "$TMP_PATH/contam_region_aln_swap_offset.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" offsetalignment "$TMP_PATH/contam_region" "$TMP_PATH/contam_region_rev" "${DB}" "$TMP_PATH/db_rev_split"  "$TMP_PATH/contam_region_aln_swap" "$TMP_PATH/contam_region_aln_swap_offset" ${THREADS_PAR} \
        || fail "rescorediagonal step died"
fi

if notExists  "$TMP_PATH/contam_region_aln_swap_offset_stats.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" createstats "${DB}" "$TMP_PATH/contam_region_aln_swap_offset" "$TMP_PATH/contam_region_aln_swap_offset_stats" ${CREATESTATS_PAR} \
        || fail "createtsv step died"
fi

if notExists "${2}_stats"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" prefixid "$TMP_PATH/contam_region_aln_swap_offset_stats" "${2}_stats" --threads 1 --tsv \
        || fail "prefixid step 1  died"
fi

if notExists  "$TMP_PATH/contam_region_aln_swap_offset_all.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" createallreport "${DB}" "$TMP_PATH/contam_region_aln_swap_offset" "$TMP_PATH/contam_region_aln_swap_offset_all" ${CREATESTATS_PAR} \
        || fail "createtsv step died"
fi

if notExists "${2}_all"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" prefixid "$TMP_PATH/contam_region_aln_swap_offset_all" "${2}_all" --threads 1 --tsv \
        || fail "prefixid step 2 died"
fi

if [ -n "$REMOVE_TMP" ]; then
  echo "Remove temporary files"
  rm -f "$3/aln_distance"    "$4/aln_distance.index"
fi

