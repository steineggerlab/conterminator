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


if notExists "$TMP_PATH/db_rev"; then
    # shellcheck disable=SC2086
    "$MMSEQS" extractframes "$DB" "$TMP_PATH/db_rev"  ${EXTRACT_FRAMES_PAR} \
        || fail "extractframes step died"
fi

if notExists "$TMP_PATH/db_rev_split"; then
    # shellcheck disable=SC2086
    "$MMSEQS" splitsequence "$TMP_PATH/db_rev" "$TMP_PATH/db_rev_split"  ${SPLITSEQ_PAR} \
        || fail "splitsequence step died"
fi

if notExists "$TMP_PATH/pref"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" kmermatcher "$TMP_PATH/db_rev_split"  "$TMP_PATH/pref"  ${KMERMATCHER_PAR} \
        || fail "kmermatcher step died"
fi

if notExists "$TMP_PATH/aln"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" rescorediagonal "$TMP_PATH/db_rev_split" "$TMP_PATH/db_rev_split" "$TMP_PATH/pref" "$TMP_PATH/aln" ${RESCORE_DIAGONAL_PAR} \
        || fail "rescorediagonal step died"
fi

if notExists "$TMP_PATH/aln_offset"; then
    # shellcheck disable=SC2086
    "$MMSEQS" offsetalignment "${DB}" "$TMP_PATH/db_rev_split" "${DB}" "$TMP_PATH/db_rev_split"  "$TMP_PATH/aln" "$TMP_PATH/aln_offset" ${OFFSETALIGNMENT_PAR} \
        || fail "rescorediagonal step died"
fi

if notExists "$TMP_PATH/aln_kingdom"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" multipletaxas "${DB}"  "$TMP_PATH/aln_offset" "$TMP_PATH/aln_kingdom" ${MULTIPLETAXA_PAR} \
        || fail "rescorediagonal step died"
fi

if notExists "$4/aln_offset_distance"; then
    # shellcheck disable=SC2086
    "$MMSEQS" distanceton  "${DB}" "${DB}" "$TMP_PATH/aln_kingdom" "$TMP_PATH/aln_kingdom_distance" ${DISTANCETON_PAR}  \
        || fail "distanceton step died"
fi

(mv -f "$TMP_PATH/aln_kingdom_distance" "$2" && mv -f "$TMP_PATH/aln_kingdom_distance.index" "$2.index") \
    || fail "Could not move result to $2"

if [ -n "$REMOVE_TMP" ]; then
  echo "Remove temporary files"
  rm -f "$3/aln_distance"    "$4/aln_distance.index"
fi

