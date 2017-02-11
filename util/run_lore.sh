#!/bin/bash

# Purpose: run multilore non-interactively to save human time
#
# Depends on mktemp from debianutils
#
# Author: Donnie Berkholz <berkhold@science.oregonstate.edu>
# 12 April 2006

TMP="$(mktemp)"
PROG="multilore"

usage() {
	echo "${0##*/} [ -m modelnums ] pdb1 [ pdb2 ] [ pdb3 ] [ ... ]"
	echo "  Run $PROG non-interactively"
	echo "  on list of PDB files."
	echo "  Each PDB must contain multiple models."
	echo "   -m Change number of models in each group. Defaults to 1,1."
	exit 1
}

err_check() {
	if [[ $1 -ne 0 ]]; then
		return 1
	else
		return 0
	fi
}

exists() {
	local err
	for arg in $@; do
		if [[ ! -e $arg ]]; then
			echo "Input file '$arg' doesn't exist"
			err=1
		fi
	done
	err_check $err
}

available() {
	local err
	for arg in $@; do
		if [[ -e $arg ]]; then
			echo "Output file '$arg' already exists"
			err=1
		fi
	done
	err_check $err
}

process_args() {
	while getopts m: OPT; do
		case $OPT in
			m) MODELS=$OPTARG ;;
			*) usage ;;
		esac
	done
	shift $(( $OPTIND - 1 ))

	if [[ ! -n "$MODELS" ]]; then
		MODELS="1,1"
	fi

	INFILES=$@
	OUTFILES=${INFILES//pdb/lore}
}

if [[ $# -lt 1 ]]; then
	usage
fi

process_args $@

exists $INFILES || err="1"
available $OUTFILES || err="1"

if [[ "$err" -eq 1 ]]; then
	exit 1
fi

for INFILE in $INFILES; do
	OUTFILE=${INFILE//pdb/lore}

	cat << EOF >> $TMP
$INFILE
$MODELS
$OUTFILE
EOF

	$PROG < $TMP

	rm $TMP
done
