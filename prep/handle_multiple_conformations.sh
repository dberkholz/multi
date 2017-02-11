#!/bin/bash

# This script turns multiple conformations in a single PDB file into an ensemble
# of PDB files, each containing a single conformation. This is for use with
# multigore, so the multiple conformations are taken into account and basically
# averaged out of the calculation. It also has the option of concatenating the
# created PDB files into a single file for feeding to multi* tools.
#
# Dependencies: >=bash-2.04
#
# Author: Donnie Berkholz <berkhold@science.oregonstate.edu>
# Karplus lab, Oregon State University

unset concatenate files specials

PREPPATH="${0%/*}/"
source ${PREPPATH}/color.sh
source ${PREPPATH}/fix_model.sh
source ${PREPPATH}/bash_check.sh

color_setup

# void usage(void)
usage() {
	echo
	cecho "  ${0##*/} [ -c concatenated-file ] pdb" b
	echo
	cecho "  Split a single PDB file with multiple conformations"
	cecho "  into multiple PDB files, each with a single conformation."
	cecho "  Output files will have the same name as the PDB, with an appended letter"
	cecho "  that's based on the letter of the conformation."
	echo
	cecho "  -c" green b
	cecho "    Concatenate the output files into one pdb"
	cecho "    for easy import into multigore and related tools."
	cecho "    This does consider the MODEL, TER and END pdb options."
	echo
	exit 1
}

bash_version

# Check for wrong number of arguments
if [[ $# -lt 1 ]] || [[ $# -gt 3 ]]; then
	usage
fi

# Handle command-line options
while getopts c: options; do
	case ${options} in
		c)	concatenate=${OPTARG}
			;;
		*)	usage
			;;
	esac
done
# Reset post-option stuff to positional parameters
shift $((OPTIND - 1))

pdb=${1}

pdbbase=$(basename ${pdb} .pdb)

cecho "Checking for multiple conformations in ${1} ... "

# First, find the number of conformations
count=0
while read line; do
	found=0
	special=${line:16:1}
	case ${special} in
		[A-Z]) 
				for letter in ${specials}; do
					if [[ "${special}" = "${letter}" ]]; then
						found=1
					fi
				done
				if [[ ${found} -ne 1 ]]; then
					specials="${specials} ${special}"
					(( count++ )) # C-like construct
				fi
			;;
	esac
done < $pdb
if [[ ! -n "${specials}" ]]; then
	specials="space"
	(( count++ ))
fi
cecho "  Conformations found: ${count}" b

# Set up all the files, so we get atoms before any multiple conformations
# written to every file
echo -n "  Will write to these files:"
for special in ${specials}; do
	files="${files} ${pdbbase}${special}.pdb"
	echo -n "  ${pdbbase}${special}.pdb"
done
echo

cecho "  Checking for existing copies..." b
for file in ${files}; do
	if [[ -e ${file} ]]; then
		cecho "    ${file} exists! Moving to ${file}.bak"
		mv ${file} ${file}.bak
	fi
done

# Scan for multiple conformations and write to files
cecho "  Creating PDB files..." b
while read line; do
	special=${line:16:1}
	case ${special} in
		[A-Z])	echo "${line:0:16} ${line:17}" >> ${pdbbase}${special}.pdb
			;;
		*)		for file in ${files}; do
					echo "${line}" >> ${file}
				done
			;;
	esac
done < $pdb

# Concatenate files
# Assume MODEL, TER, ENDMDL are in the right spots in the original.
if [[ -n "${concatenate}" ]]; then
	cecho "  Concatenating PDB files to ${concatenate}..." b
	if [[ ! -e "${concatenate}" ]]; then
		model=0
		for file in ${files}; do
			fix_model ${file} ${concatenate}
			# Get rid of ENDs
			sed -i '/^END$/d' ${concatenate}

			echo "END" >> ${concatenate}
		done
	else
		cecho "${concatenate} already exists! Please remove it."
	fi
fi
