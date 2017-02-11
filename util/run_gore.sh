#!/usr/bin/env bash

# Purpose: run multicore/gore non-interactively to save human time
#
# Depends on mktemp from debianutils
#
# Author: Donnie Berkholz <berkhold@science.oregonstate.edu>
# 12 April 2006

# Max deviation of atoms to use in overlaying structures, in angstroms
MAX_DEV="1.0"

PROG_1="multicore"
PROG_2="multigore"

usage() {
	echo "${0##*/} [ -m MODELS ] pdb1 [ pdb2 ] [ pdb3 ] [ ... ]"
	echo "  Run $PROG_1 and $PROG_2 non-interactively"
	echo "  on list of PDB files."
	echo "  Each PDB must contain two models unless you pass -m MODELS."
	echo "   -m MODELS  Pass in a two-field, comma-delimited pair of numbers"
	echo "      that indicates how many models are in each group."
	exit 1
}

err_check() {
	if [[ $1 -ne 0 ]]; then
		return 1
	else
		return 0
	fi
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

	OUTFLAGS=${INFILES//pdb/flag}
	OUTFILES="${INFILES//pdb/flag} ${INFILES//pdb/gore}"
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
	OUTFLAG=${INFILE//pdb/flag}
	OUTGORE=${INFILE//pdb/gore}
	OUTMATRIX=${INFILE//pdb/rot}

	if type mktemp > /dev/null 2>&1; then
		TMP="$(mktemp)"
	else
		TMP=$RANDOM
	fi

	cat << EOF >> $TMP
$INFILE
$OUTFLAG
$MAX_DEV
EOF

	$PROG_1 < $TMP

	mv multicore.log ${INFILE/.pdb}-multicore.log

	rm $TMP fort.*

	if type mktemp > /dev/null 2>&1; then
		TMP="$(mktemp)"
	else
		TMP=$RANDOM
	fi

	# multigore outputs an average structure if either group of models is
	# greater than 1
	if [[ `echo $MODELS | cut -d, -f1` -gt 1 ]]; then
		OUTAVG_1=${INFILE/.pdb/-avg1.pdb}
		if [[ `echo $MODELS | cut -d, -f2` -gt 1 ]]; then
			OUTAVG_2=${INFILE/.pdb/-avg2.pdb}
		fi
	elif [[ `echo $MODELS | cut -d, -f2` -gt 1 ]]; then
		OUTAVG_1=${INFILE/.pdb/-avg2.pdb}
	fi

	cat << EOF >> $TMP
$INFILE
$MODELS
$OUTGORE
$OUTFLAG
$OUTAVG_1
$OUTAVG_2
EOF

	$PROG_2 < $TMP

	mv fort.80 ${OUTMATRIX}

	rm $TMP fort.*

done

