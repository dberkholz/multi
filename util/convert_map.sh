#!/bin/bash

# Purpose: Convert electron-density maps between various formats
# Depends on mapman, part of Gerard Kleywegt's RAVE package,
# and mktemp, part of debianutils.
#
# Author: Donnie Berkholz <berkhold@science.oregonstate.edu>
# 10 April 2006

# To get valid map types out of mapman, run:
# read map1 x ?

# Increase max memory allowed for mapman -- 20 megabytes
export MAPSIZE="$( echo $(( 20 * 1000 * 1000 )) )"

# Handle more than 2 maps at a time
export NUMMAPS="5"

# Create a temporary file for mapman script
TMP="$(mktemp)"

usage() {
	echo "${0##*/} -i informat -o outformat inmap outmap"
	echo "  Convert a map inmap from informat to outformat,"
	echo "  and put the output into outmap."
	exit 1
}

check_mapman() {
	which mapman >& /dev/null
	if [[ $? -ne 0 ]]; then
		echo "mapman not found!"
		echo "Please install the RAVE package"
		echo "from http://alpha2.bmc.uu.se/usf/rave.html"
		exit 1
	fi
}

process_args() {
	while getopts i:o: OPT; do
		case $OPT in
			i) INFORMAT=$OPTARG ;;
			o) OUTFORMAT=$OPTARG ;;
			*) usage ;;
		esac
	done
	shift $(( $OPTIND - 1 ))
	INFILE=$1
	OUTFILE=$2
}

check_files() {
	if [[ ! -e "$INFILE" ]]; then
		echo "Input map '$INFILE' doesn't exist"
		exit 1
	elif [[ -e $OUTFILE ]]; then
		echo "Output map '$OUTFILE' already exists"
		exit 1
	fi
}

if [[ $# -ne 6 ]]; then
	usage
fi

check_mapman

process_args $@

check_files

cat << EOF >> $TMP
read map1 $INFILE $INFORMAT
write map1 $OUTFILE $OUTFORMAT
quit
EOF

mapman < $TMP

rm $TMP
