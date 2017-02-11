#!/bin/bash

# Purpose: concatenate every possible combination of PDB files.
#
# Author: Donnie Berkholz <berkhold@science.oregonstate.edu>
# 11 April 2006

INFILES=$@

usage() {
	echo "${0##*/} pdb1 pdb2 [ pdb3 ] [ ... ]"
	echo "  Concatenate multiple PDB files"
	echo "  and name the resulting PDB logically."
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
			echo "Input PDB '$arg' doesn't exist"
			err=1
		fi
	done
	err_check $err
}

available() {
	local err
	for arg in $@; do
		if [[ -e $arg ]]; then
			echo "PDB '$arg' already exists"
			err=1
		fi
	done
	err_check $err
}

if [[ $# -lt 2 ]]; then
	usage
fi

concat() {
	while [[ $# -ge 1 ]]; do
		IN_1=$1
		OUT_1=${1/-replicate}
		shift
		for IN_2 in $@; do
			OUT_2=${IN_2/-replicate}
			OUTFILE=${OUT_1%.pdb}-$OUT_2;
			if [[ -n "$REAL" ]]; then
				cat $IN_1 $IN_2 > $OUTFILE
			else
				exists $IN_1 $IN_2 || err=1
				available $OUTFILE || err=1
			fi
		done
	done
	err_check $err
}

concat $@ || exit 1
REAL="yes" concat $@
