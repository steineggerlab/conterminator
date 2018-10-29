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
[ "$#" -ne 4 ] && echo "Please provide <queryDB> <targetDB> <outDB> <tmp>" && exit 1;
# check if files exists
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[ ! -f "$2" ] &&  echo "$2 not found!" && exit 1;
[   -f "$3" ] &&  echo "$3 exists already!" && exit 1;
[ ! -d "$4" ] &&  echo "tmp directory $4 not found!" && mkdir -p "$4";


QUERY="$1"
TARGET="$2"
TMP_PATH="$4"

mkdir -p "$4/search"
if notExists "$4/aln"; then
    # shellcheck disable=SC2086
    $RUNNER "$MMSEQS" search "${QUERY}" "${TARGET}" "$4/aln" "$4/search" ${SEARCH1_PAR} \
        || fail "Search step died"
fi

if notExists "$4/aln_offset"; then
    # shellcheck disable=SC2086
    "$MMSEQS" distanceton "${QUERY}" "${TARGET}" "$4/aln"  "$4/aln_distance"  \
        || fail "Offset step died"
fi

(mv -f "$4/aln_distance" "$3" && mv -f "$4/aln_distance.index" "$3.index") \
    || fail "Could not move result to $3"

if [ -n "$REMOVE_TMP" ]; then
  echo "Remove temporary files"
  rm -f "$4/aln_distance"    "$4/aln_distance.index"
fi

