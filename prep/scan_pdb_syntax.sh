#!/bin/bash

# The purpose of this script is to check for a PDB in compliance with what we
# need for the multi* tools.
# Things like the right format: MODEL #, atoms, TER, ENDMDL, MODEL #... END.
#
# Dependencies: >=bash-2.04
#
# Author: Donnie Berkholz <berkhold@science.oregonstate.edu>
# Karplus lab, Oregon State University

unset lastline

PREPPATH="${0%/*}/"
source ${PREPPATH}/color.sh
source ${PREPPATH}/bash_check.sh

color_setup

# void usage(void)
usage() {
	echo
	cecho "  ${0##*/} [ -m ] pdb" b
	echo
	cecho "  Check for PDB file syntax compliant with the multi* tools."
	echo
	cecho "  -m" green b
	cecho "    Check for invalid multiple conformation formats:"
	cecho "    AABB instead of ABAB."
	exit 1
}

modelerror() {
	local defaultmsg="${1:-${errortype} not allowed outside MODEL/ENDMDL.}"
	cecho "  line ${lineno}: ${defaultmsg}" red b
	error=1
}

bash_version

# Check for wrong number of arguments
if [[ $# -ne 1 ]] && [[ $# -ne 2 ]]; then
	usage
fi

# Handle command-line options
while getopts m options; do
	case ${options} in
		m)	CHUNK_SCAN=1
			;;
		*)	usage
			;;
	esac
done
# Reset post-option stuff to positional parameters
shift $((OPTIND - 1))

pdb=${1}

pdbbase=$(basename ${pdb} .pdb)

cecho "Scanning ${pdb##*/} ..." b
modeldepth=0
lineno=0
ter=0
while read line; do
	(( lineno++ ))

	# Depth into MODEL/ENDMDL blocks, should only ever be 0 or 1
	if [[ "${line/#MODEL}" != "${line}" ]]; then
		(( modeldepth++ ))
	elif [[ "${line/#ENDMDL}" != "${line}" ]]; then
		(( modeldepth-- ))
	fi

	if [[ "${modeldepth}" -ge 2 ]]; then
		modelerror "No ENDMDL to close MODEL/ENDMDL block."
		(( modeldepth-- )) # Recover from error
	elif [[ "${modeldepth}" -lt 0 ]]; then
		modelerror "ENDMDL with no initial MODEL."
		(( modeldepth++ )) # Recover from error
	fi

	if [[ "${line/#ENDMDL}" != "${line}" ]]; then
		if [[ "${ter}" -eq 0 ]]; then
			modelerror "No TER at end of MODEL/ENDMDL block."
		fi
		ter=0
	fi

	# Check for valid MODEL line -- number in columns 11-14, rather than
	# immediately following MODEL as would be logical
	if [[ "${line/#MODEL}" != "${line}" ]]; then
		if [[ "${line:5:5}" != "     " ]]; then
			modelerror "MODEL should have spaces in columns 6-10."
		# Hack: strings equal zero in a numerical-only comparison
		elif [[ ${line:10:4} -eq 0 ]]; then
			modelerror "MODEL should have number in columns 11-14."
		fi
	fi

	# No coordinate-related things allowed outside MODEL
	if [[ "${modeldepth}" -eq 0 ]]; then
		if [[ "${line/#ATOM}" != "${line}" ]]; then
			errortype=ATOM modelerror
		fi
		if [[ "${line/#ANISOU}" != "${line}" ]]; then
			errortype=ANISOU modelerror
		fi
		if [[ "${line/#TER}" != "${line}" ]]; then
			errortype=TER modelerror
		fi
	# No END allowed inside MODEL
	elif [[ "${modeldepth}" -ge 1 ]]; then
		if [[ "${line/%END}" != "${line}" ]]; then
			modelerror "END not allowed inside MODEL/ENDMDL block."
		fi
	fi

	# We need at least one TER within every MODEL/ENDMDL block, on the line
	# before ENDMDL
	if [[ "${line/#TER}" != "${line}" ]]; then
		ter=1
	fi

	if [[ -n ${CHUNK_SCAN} ]] && [[ "${line/#ATOM}" != "${line}" ]]; then
		if [[ "${lastline:16:1}" != " " ]] && [[ "${line:16:1}" != " " ]] \
			&& [[ -n "${lastline:16:1}" ]] && [[ -n "${line:16:1}" ]]; then
			if [[ "${lastline:16:1}" = "${line:16:1}" ]]; then
				modelerror "${line:17:3}${line:22:4}: Two ${lastline:16:1} conformations in a row."
			fi
		fi
	fi

	lastline="${line}"
done < ${pdb}

# If we're missing the END tag, check for other problems...
if [[ "${lastline}" != "END" ]]; then
	if [[ "${lastline}" = "ENDMDL" ]]; then
		modelerror "END tag missing."
	elif [[ "${lastline}" = "TER" ]]; then
		modelerror "ENDMDL and END tags missing."
	else
		modelerror "TER, ENDMDL and END tags missing."
	fi
fi

# Die if we got any errors
if [[ "${error}" -eq 1 ]]; then
	cecho "  Errors found during scan." red b
	exit 1
fi

# We made it through, so everything's OK
cecho "  PDB syntax clean!" b
