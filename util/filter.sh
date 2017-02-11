#!/bin/bash

# Purpose: Filter an input file by any of its columns for numbers under a given
# value. One use of this is to filter PDBs for B values under e.g., 15.00.
#
# Author: Donnie Berkholz <berkhold@science.oregonstate.edu>
# Karplus lab, Oregon State University

maxval_default="15"
column_default="3"

usage() {
	echo "${0##*/} [ -m maxval ] [ -c column ] infile"
	echo "  infile should be any tab-delimited file;"
	echo "  for example, an output file from multigore."
	echo "  By default, ${0##*/} uses ${maxval_default} for maxval and looks in column ${column_default}."
}

# Check for wrong number of arguments
if [ $# -lt 1 -o $# -gt 5 ]; then
	usage
	exit 1
fi

while getopts m:c: options; do
	case $options in
		m) maxval="${OPTARG}"
			;;
		c) column="${OPTARG}"
			;;
	esac
done

infile="${!OPTIND}"

maxval="${maxval:-${maxval_default}}"
column="${column:-${column_default}}"

# Here's the magic
awk '{ while ( getline ) if ( $'"${column}"' < '"$maxval"' ) print; }' < $infile
