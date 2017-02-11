#!/bin/bash

# Purpose: Delete entries not present in all PDB's. This allows multi* to
# make an effective analysis across atoms common to all of the PDB files, and
# multi* will not work with non-replicate atoms.
#
# Dependencies: >=bash-2.04
#
# http://moongroup.com/pipermail/shell.scripting/2002-October/005754.html
# may be useful for understanding use of dynamic variables.
#
# Author: Donnie Berkholz <berkhold@science.oregonstate.edu>
# Karplus lab, Oregon State University

if [[ -n "${DEBUG}" ]]; then
	DEBUGATOMSEQUAL=1
	DEBUGFIRSTARRAY=1
	DEBUGREADARRAYS=1
	DEBUGCOMMONARRAY=1
	DEBUGSCAN=1
fi

unset concatenate files specials IN_ARRAY_STOP
LAST_NONMULTI="0"

PREPPATH="${0%/*}/"
source ${PREPPATH}/color.sh
source ${PREPPATH}/bash_check.sh

color_setup

# void usage(void)
usage() {
	echo
	cecho "  ${0##*/} pdb pdb [ pdb ... ]" b
	echo
	cecho "  Remove entries not common to all PDBs supplied as arguments."
	cecho "  This allows multi* to make an effective analysis across atoms"
	cecho "  common to all of the PDB files, and multi* will not work with"
	cecho "  non-replicate atoms."
	echo
	cecho "  Output will be written to files the same as the arguments, but"
	cecho "  will end with -replicate.pdb."
	echo
	exit 1
}

res_num() {
	echo ${1:22:4}
}

res_name() {
	echo ${1:17:3}
}

res_id() {
	echo $(res_name "${1}") $(res_num "${1}")
}

atom_name() {
	echo ${1:12:4}
}

atom_id() {
	echo $(res_id "${1}") $(atom_name "${1}")
}

residues_equal() {
	if [[ "${1:22:4}" = "${2:22:4}" ]]; then
		return 0
	fi

	return 1
}

atoms_equal() {
	if [[ -n "${DEBUGATOMSEQUAL}" ]]; then
		echo "1:0:26 = ${1:0:26}"
		echo "2:0:26 = ${2:0:26}"
	fi
	if residues_equal "${1}" "${2}"; then
		if [[ "${1:0:6}" = "${2:0:6}" ]] && [[ "${1:11:5}" = "${2:11:5}" ]] && [[ "${1:17:5}" = "${2:17:5}" ]]; then
			return 0
		fi
	fi

	return 1
}

is_atom() {
	if [[ ${1:0:4} = "ATOM" ]]; then
		return 0
	fi

	return 1
}

is_multi() {
	if [[ ${1:16:1} = [A-Z] ]]; then
		return 0
	fi

	return 1
}

# void in_array(char * arrayname[], char * elementname)
in_array() {
	local i
	local array="${1}"
	local array_length=$(eval echo \${#${array}[@]})
	for (( i = ${IN_ARRAY_STOP:-0}; i < ${array_length}; i++ )); do
		local myvar=${array}[${i}]
		local RES_MYVAR=$(res_num "${!myvar}")
		local RES_2=$(res_num "${!2}")

		if ! is_multi "${!myvar}"; then
			LAST_NONMULTI=${i}
		fi

		# If we aren't in the right residue, don't bother comparing
		if [[ ${RES_MYVAR} -lt ${RES_2} ]]; then
			continue
		fi

		# Stop looking if we get past the desired residue number in the primary
		# array
		# This is where AABB-style multiple conformations screw us.
		if [[ ${RES_MYVAR} -gt ${RES_2} ]]; then
			if ! is_multi "${!myvar}" && ! is_multi "${!2}"; then
				echo ${LAST_NONMULTI}
				return 2
			fi
		fi

		if atoms_equal "${!myvar}" "${!2}"; then
			echo ${LAST_NONMULTI}
			return 0
		fi
	done

	echo ${LAST_NONMULTI}
	return 1
}

# Handle command-line options
while getopts o: options; do
	case ${options} in
		o)	outfile=${OPTARG}
			;;
		*)	usage
			;;
	esac
done
# Reset post-option stuff to positional parameters
shift $((OPTIND - 1))

bash_version

# Check for wrong number of arguments
if [[ $# -lt 1 ]]; then
	usage
fi

infiles=$@
for i in $@; do
	if [[ "${i%.pdb}" = "${i}" ]]; then
		cecho "${i} doesn't end in .pdb"
		exit 1
	fi
done

# How to deal with identical parts of arbitrary number of files?
# We're looking for the first 26 columns to be identical.
# That's ATOM/ANISOU, atom #, atom type, residue type, residue #
#
# Basically what we want is an intersect algorithm.
#
# Something like this:
# 	1. Read in all files, one array per file
# 	2. Take first line in A. Is line in B? If so, is line in C? If in all,
# 		save line.
# 	3. Iterate through all lines in A. No need to iterate through other files,
# 		since all saved lines must be common to every PDB and therefore in A.


# Read files into arrays
for file in $@; do
	# Get rid of hyphens and periods in name, since bash doesn't allow these in
	# variables
	filename=${file##*/}
	arrayname=array${filename//-/}
	arrayname=${arrayname//.pdb/}
	# Collect a list of all the files
	arraynames="${arraynames} ${arrayname}"

	lineno=0 # Lines match array elements, not PDB lines
	declare -a ${arrayname}

	# Read in file
	while read line; do
		# To assign an element of a dynamic array, we need to:
		# 	1. Assign a new variable to that element
		# 	2. Run eval to dereference the new var, then assign the line to it
		myvar=${arrayname}[${lineno}]
		eval ${myvar}='"${line}"'

		(( lineno++ ))
	done < ${file}
done

# Print all elements of all files, in sequential order by file
if [[ -n "${DEBUGREADARRAYS}" ]]; then
	cecho "Entering DEBUGREADARRAYS"
	for array in ${arraynames}; do
		# To get the array length, we need to do a double dereference using
		# `eval echo`:
		# 	1. echo dereferences the unescaped variable
		# 	2. eval dereferences the other variable.
		# This is necessary because bash can't handle dereferencing a variable
		# inside a variable.
		array_length=$(eval echo \${#${array}[@]})
		for (( i = 0; i < ${array_length}; i++ )); do
			# To get at the actual contents of an array element, we need to:
			# 	1. Assign a new variable to that element
			# 	2. Indirectly dereference the new variable
			myvar=${array}[${i}]
			cecho " ${myvar} = ${!myvar}"
		done
	done
	cecho "Exiting DEBUGREADARRAYS"
fi

# Grab the first file as a standard to compare against
declare -a firstarray
for array in ${arraynames}; do
	if [[ -n "${firstarray}" ]]; then
		break
	fi

	array_length=$(eval echo \${#${array}[@]})
	for (( i = 0; i < ${array_length}; i++ )); do
		myvar=${array}[${i}]
		firstarray[${i}]="${!myvar}"
	done
done

if [[ -n "${DEBUGFIRSTARRAY}" ]]; then
	cecho "Entering DEBUGFIRSTARRAY"
	for (( i = 0; i < ${#firstarray[@]}; i++ )); do
		cecho " firstarray[${i}] = ${firstarray[${i}]}"
	done
	cecho "Exiting DEBUGFIRSTARRAY"
fi

# This just sounds like somewhere a recursive function might be useful, but I'm
# not sure how the details would work.

# OK, now we've got our "standard" array to compare others against.
# Was that even a good idea?

declare -a commonarray
element=0
match=0
# Skip the array we scanned into firstarray, and use positional parameters
# for arrays
set ${arraynames}
# Save first array name for later
firstarrayname=${1}
# Drop $1
shift

# "firstline" means a line in the first array, not the actual first line.
# Be careful to preserve multiple conformations in one vs. single in another
# 	FIXME: How to do this?
# 			Relevant case here: multiple in firstline, single in others
# 			Can resolve some of this by rescanning the last match line in others
# 			eval ${start}=${i} instead of eval ${start}="$((i + 1 ))"
# 			This will only solve the ABABAB case, not AAABBB.
#
# 			The problem is PDB's (can?) have all one conformation, then all
# 			another, e.g., 4 lines of A, then 4 lines of B. How to deal?
# 			Could create some mini array read-ahead for all of ${letter}
# 			Have to compensate for approaching ${array_length} etc.
#
# 			It's not clear to me that we even have to care about multiples in
# 			this loop -- common elements will be added regardless. We should
# 			deal with them later on, when we're actually scanning through arrays
# 			to write to the files.
#
# 			But a run with a multi-PDB first does produce incorrect output on
# 			both the ABAB and AABB-style multi's.
#
# 			Maybe we shouldn't deal with AABB-style multi at all and instead,
# 			catch them in scan_pdb_syntax.sh. That sounds like a good option,
# 			but it may mean lots of manual work if programs consistently produce
# 			AABB-style multi's.
#
# 			It fails on ABABAB -> AAAB _because_ the B is the final atom in
# 			the residue, so there must be some way of dealing with residue
# 			number mismatches between firstarray and the others.
cecho "Creating list of elements common to all files ... "
# Take a line in the first array.
for firstline in "${firstarray[@]}"; do
	# Scan through all arrays
	for array in ${@}; do
		array_length=$(eval echo \${#${array}[@]})
		# Scan through each line in an array
		for (( i = $(eval echo \${${array}start:-0}); i < ${array_length}; i++ )); do
			myvar=${array}[${i}]
			arrayline=${!myvar}
			# If arrayline = firstline, then we stop looking in this array
			# and move on to the next array.
			if atoms_equal "${firstline}" "${arrayline}"; then
				if [[ -n "${DEBUGCOMMONARRAY}" ]]; then
					cecho "  $(atom_id "${firstline}") found!"
				fi
				match=1
				start=${array}start
				# ${i} here instead of $((i+1)) actually makes things worse.
				# Instead of ABABAB->AAAB, it gives ABABAB->ABABAB (no change).
				# But perhaps that suggests we just aren't handling multi
				# correctly yet and isn't a true problem in and of itself.
				eval ${start}="$(( i + 1 ))"
				break
			fi
			match=0
			# If we get past the residue we're looking for, quit looking
			# altogether and move on to next firstline.
			# NOTE: This can be screwed up by things like randomly inserted
			# waters with high numbers -- deal with this in multiprep.sh.
			# Let's just skip this check while we're on a multiple conformation.
			# This skipping for multiples actually slows down the program by
			# about 15x, from ~1 min to ~15 min on my computer.
			# But we still get a blank model 3 with or without it.
			if is_atom "${arrayline}" && is_atom "${firstline}"; then
				if [[ $(res_num "${firstline}") -lt $(res_num "${arrayline}") ]]; then
					if [[ -n "${DEBUGCOMMONARRAY}" ]]; then
						cecho "  Got to residue ${arrayline:22:4} without matching all of residue ${firstline:22:4}."
						cecho "    This may be because of a missing residue or multiple conformations."
					fi
					break 2
				fi
			fi
		done
		if [[ ${match} -ne 1 ]]; then
			cecho "  $(atom_id "${firstline}") not found"
		fi
	done

	if [[ "${match}" -eq 1 ]]; then
		# Make sure it isn't already in the commonarray, to account for
		# possible multiple conformations
		IN_ARRAY_STOP=$(in_array commonarray firstline)
		IN_ARRAY_RET=$?
		if [[ ${IN_ARRAY_RET} -ne 0 ]]; then
			commonarray[${element}]="${firstline}"
			(( element++ ))
		fi
	fi
done

if [[ -n "${DEBUGCOMMONARRAY}" ]]; then
	cecho "Entering DEBUGCOMMONARRAY"
	for (( i = 0; i < ${#commonarray[@]}; i++ )); do
		cecho " commonarray[${i}] = ${commonarray[${i}]}"
	done
	cecho "Exiting DEBUGCOMMONARRAY"
fi

# Now we want to use commonarray's first 26 columns to determine which rows of
# the original arrays to write out.

# Be careful to preserve multiple conformations in one vs. single in another
# 	FIXME: How to do this? Two cases:
# 		1) multiple in firstline, single in others
# 			DO HERE
# 		2) single in firstline, multiple in others
# 			Don't stop after the first match. Scan until line !atom_eq lastline
#			DONE HERE, WORKS
# 		In both cases, we must ignore the multiple conformation column.
for infile in ${infiles}; do
	outfile=${infile##*/}
	outfile="${outfile%.pdb}-replicate.pdb"
	if [[ -n "${DEBUGSCAN}" ]]; then
		cecho "outfile=${outfile}"
	fi
	if [[ -e "${outfile}" ]]; then
		rm ${outfile}
	fi
	# Get rid of hyphens and periods in name, since bash doesn't allow these in
	# variables
	arrayinfile=${infile##*/}
	arrayname=array${arrayinfile//-/}
	arrayname=${arrayname//.pdb/}
	if [[ -n "${DEBUGSCAN}" ]]; then
		cecho "arrayname=$arrayname"
	fi

	## Let's run through the common array, then scan the current array for that.
	# We can leave a "stopping" element to start the next scan at.

	array_length=$(eval echo \${#${arrayname}[@]})

	unset start
	cecho "Scanning for common elements in ${infile} ... "
	LAST_NONMULTI=0
	for (( i = 0; i < ${#commonarray[@]}; i++ )); do
		found=0
		for (( j = ${start:-0}; j < ${array_length}; j++ )); do
			myvar=${arrayname}[${j}]

			# If we're on a higher common residue than the specific one we want
			# to find, and we're not a multiple conformation, then we're too far
			if [[ $(res_num "${!myvar}") -gt $(res_num "${commonarray[${i}]}") ]]; then
				if ! is_multi "${!myvar}" && ! is_multi "${commonarray[${i}]}"; then
					break
				fi
			fi

			# Do we have a match? If so, write out the line.
			if atoms_equal "${!myvar}" "${commonarray[${i}]}"; then
				if [[ -n "${DEBUGSCAN}" ]]; then
					cecho "Match at $(atom_id "${!myvar}")"
				fi

				found=1
				echo "${!myvar}" >> ${outfile}
			fi

			# No need to keep looking in this array if we already found all
			# valid atoms (i.e., if it isn't multiple).
			if [[ "${found}" -eq 1 ]]; then
				if ! is_multi "${!myvar}"; then
					break
				fi
			fi

			if ! is_multi "${!myvar}"; then
				LAST_NONMULTI=${j}
			fi

			if [[ "${j}" -ge $(( array_length-1 )) ]]; then
				cecho "$(atom_id "${!myvar}") not found!"
			fi
		done
		start="$(( LAST_NONMULTI+1 ))"
	done
done
