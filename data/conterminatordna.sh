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

if notExists "$TMP_PATH/db_rev_split"; then
    # shellcheck disable=SC2086
    "$MMSEQS" splitsequence "$TMP_PATH/sequencedb" "$TMP_PATH/db_rev_split"  ${SPLITSEQ_PAR} \
        || fail "splitsequence step died"
fi

if notExists "$TMP_PATH/pref"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" kmermatcher "$TMP_PATH/db_rev_split" "$TMP_PATH/pref"  ${KMERMATCHER_PAR} \
        || fail "kmermatcher step died"
fi

if notExists "$TMP_PATH/aln.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" rescorediagonal "$TMP_PATH/db_rev_split" "$TMP_PATH/db_rev_split" "$TMP_PATH/pref" "$TMP_PATH/aln" ${RESCORE_DIAGONAL1_PAR} \
        || fail "rescorediagonal step died"
fi

if notExists "$TMP_PATH/aln_offset.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" offsetalignment "$TMP_PATH/sequencedb" "$TMP_PATH/db_rev_split" "$TMP_PATH/sequencedb" "$TMP_PATH/db_rev_split"  "$TMP_PATH/aln" "$TMP_PATH/aln_offset" ${OFFSETALIGNMENT_PAR} \
        || fail "rescorediagonal step died"
fi

if notExists "$TMP_PATH/contam_aln.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" extractalignments "$TMP_PATH/sequencedb" "$TMP_PATH/aln_offset" "$TMP_PATH/contam_aln" ${EXTRACTALIGNMENTS_PAR} \
        || fail "extractalignment step died"
fi

if notExists "$TMP_PATH/contam_region.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" extractalignedregion "$TMP_PATH/sequencedb" "$TMP_PATH/sequencedb" "$TMP_PATH/contam_aln" "$TMP_PATH/contam_region" ${THREADS_PAR} \
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
    $RUNNER "$MMSEQS" rescorediagonal "$TMP_PATH/db_rev_split" "$TMP_PATH/contam_region_rev" "$TMP_PATH/contam_region_pref" "$TMP_PATH/contam_region_aln" ${RESCORE_DIAGONAL2_PAR} \
        || fail "rescorediagonal2 step died"
fi

if notExists "${TMP_PATH}/contam_region_aln_swap.dbtype"; then
     # shellcheck disable=SC2086
    "$MMSEQS" swapresults "$TMP_PATH/db_rev_split" "$TMP_PATH/contam_region_rev" "${TMP_PATH}/contam_region_aln" "${TMP_PATH}/contam_region_aln_swap" ${SWAP_PAR} \
        || fail "Swapresults pref died"
fi

if notExists "$TMP_PATH/contam_region_aln_swap_offset.dbtype"; then
    # shellcheck disable=SC2086
    "$MMSEQS" offsetalignment "$TMP_PATH/contam_region" "$TMP_PATH/contam_region_rev" "$TMP_PATH/sequencedb" "$TMP_PATH/db_rev_split"  "$TMP_PATH/contam_region_aln_swap" "$TMP_PATH/contam_region_aln_swap_offset" ${THREADS_PAR} \
        || fail "rescorediagonal step died"
fi

if notExists  "$TMP_PATH/contam_region_aln_swap_offset_all.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" createallreport "$TMP_PATH/sequencedb" "$TMP_PATH/contam_region_aln_swap_offset" "$TMP_PATH/contam_region_aln_swap_offset_all" ${CREATESTATS_PAR} \
        || fail "createtsv step died"
fi

if notExists "${3}_all"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" prefixid "$TMP_PATH/contam_region_aln_swap_offset_all" "${3}_all" --threads 1 --tsv \
        || fail "prefixid step 2 died"
fi

if notExists  "$TMP_PATH/contam_region_aln_swap_offset_predconterm.dbtype"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" predictcontamination "$TMP_PATH/contam_region_aln_swap_offset_all" "$TMP_PATH/contam_region_aln_swap_offset_predconterm" ${CREATESTATS_PAR} \
        || fail "createtsv step died"
fi

if notExists "${3}_stats"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" prefixid "$TMP_PATH/contam_region_aln_swap_offset_predconterm" "${3}_conterm_prediction" --threads 1 --tsv \
        || fail "prefixid step 1  died"
fi

if [ -n "$REMOVE_TMP" ]; then
  echo "Remove temporary files"
  $MMSEQS rmdb "$TMP_PATH/contam_region_aln_swap_offset_predconterm"
  $MMSEQS rmdb "$TMP_PATH/contam_region_aln_swap_offset_all"
  $MMSEQS rmdb "$TMP_PATH/contam_region_aln_swap_offset"
  $MMSEQS rmdb "$TMP_PATH/contam_region_aln_swap"
  $MMSEQS rmdb "$TMP_PATH/contam_region_pref"
  $MMSEQS rmdb "$TMP_PATH/contam_region_aln"
  $MMSEQS rmdb "$TMP_PATH/contam_region_rev"
  $MMSEQS rmdb "$TMP_PATH/contam_region"
  $MMSEQS rmdb "$TMP_PATH/db_rev_split"
  $MMSEQS rmdb "$TMP_PATH/contam_aln"
  $MMSEQS rmdb "$TMP_PATH/aln_offset"
  $MMSEQS rmdb "$TMP_PATH/sequencedb"
  $MMSEQS rmdb "$TMP_PATH/sequencedb_h"
  $MMSEQS rmdb "$TMP_PATH/pref"
  $MMSEQS rmdb "$TMP_PATH/aln"
fi

