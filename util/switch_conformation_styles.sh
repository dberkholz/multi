#!/bin/bash

# This script interconverts ABAB- and AABB-style multiple conformations.
#
# Dependencies: >=bash-2.04
#
# Author: Donnie Berkholz <berkhold@science.oregonstate.edu>
# Karplus lab, Oregon State University

unset concatenate files specials

PREPPATH="${0%/*}/../prep/"
source ${PREPPATH}/color.sh
source ${PREPPATH}/bash_check.sh

color_setup

# void usage(void)
usage() {
	echo
	cecho "  ${0##*/} [ -r ] pdb" b
	echo
	cecho "  Interconvert between ABAB- and AABB-style multiple conformations."
	cecho "  By default, go AABB->ABAB."
	echo
	cecho "  Output file will have '-conv.pdb' at the end."
	echo
	cecho "  -r" green b
	cecho "    Convert in the reverse direction: ABAB->AABB."
	echo
	exit 1
}

# Checks for the first argument in a list of one or more arguments
# void has(char * item, char * list_item1, ...)
has() {
	[[ " ${*:2} " == *" $1 "* ]]
}

# Right now, we aren't renumbering the atom serials, we're just changing the
# order.
get_coord_num() {
	echo ${1:6:5}
}

get_multi_letter() {
	echo ${1:16:1}
}

is_atom() {
	if [[ ${1:0:4} = "ATOM" ]]; then
		return 0
	fi

	return 1
}

is_sigatm() {
	if [[ ${1:0:6} = "SIGATM" ]]; then
		return 0
	fi

	return 1
}

is_anisou() {
	if [[ ${1:0:6} = "ANISOU" ]]; then
		return 0
	fi

	return 1
}

is_siguij() {
	if [[ ${1:0:6} = "SIGUIJ" ]]; then
		return 0
	fi

	return 1
}

is_hetatm() {
	if [[ ${1:0:6} = "HETATM" ]]; then
		return 0
	fi

	return 1
}

is_coord() {
	if is_atom "${1}" || is_sigatm "${1}" || is_anisou "${1}" \
		|| is_siguij "${1}" || is_hetatm "${1}"; then

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

## Perhaps the superior way to implement this would be to have a very complete,
## more complex sort_array_within() that took an optional second-level parameter
## to sort. Then sort_array() would just act as a wrapper for it to simplify the
## calling. Maybe a recursive function that calls itself only if there's a
## second argument?

# void
# sort_array(char * in_arrayname, char * out_arrayname, char * sort_function)
sort_array() {
	local i myvar sort_function_out sort_by
	local in_array="${1}"
	local out_array="${2}"
	local sort_function="${3}"
	local array_length=$(eval echo \${#${in_array}[@]})
	# Set up the things to sort by
	for (( i = 0; i < ${array_length}; i++ )); do
		myvar=${in_array}[${i}]
		sort_function_out=$(${sort_function} "${!myvar}")
		if ! has "${sort_function_out}" ${sort_by}; then
			sort_by="${sort_by} ${sort_function_out}"
		fi
	done
#	sort_by=$(echo ${sort_by} | sort -u)

	local sort_item arrayvar
	local lineno=0
	# Run the sort
	for sort_item in ${sort_by}; do
		for (( i = 0; i < ${array_length}; i++ )); do
			local myvar=${in_array}[${i}]
			if [[ $(${sort_function} "${!myvar}") = ${sort_item} ]]; then
				arrayvar=${out_array}[${lineno}]
				eval ${arrayvar}='"${!myvar}"'
				(( lineno++ ))
			fi
		done
	done
}

# void
# sort_array_recursive(char * in_arrayname, char * out_arrayname,
# 					char * outer_sort_function, char * inner_sort_function)
sort_array_recursive() {
	local i myvar inner_sort_function_out inner_sort_by
	local in_array="${1}"
	local out_array="${2}"
	local outer_sort_function="${3}"
	local inner_sort_function="${4}"
	local array_length=$(eval echo \${#${in_array}[@]})
	# Set up the things to sort by
	for (( i = 0; i < ${array_length}; i++ )); do
		myvar=${in_array}[${i}]
		inner_sort_function_out=$(${inner_sort_function} "${!myvar}")
		if ! has "${inner_sort_function_out}" ${inner_sort_by}; then
			inner_sort_by="${inner_sort_by} ${inner_sort_function_out}"
		fi
	done
	cecho "Sort items: $inner_sort_by"
#	inner_sort_by=$(echo ${inner_sort_by} | sort -u)

	local sort_item arrayvar
	local lineno=0
	# Run the sort
	for sort_item in ${inner_sort_by}; do
		cecho "Sorting by $sort_item"
		for (( i = 0; i < ${array_length}; i++ )); do
			echo -n "Running $inner_sort_function in $sort_item ... "
			local myvar=${in_array}[${i}]
			if [[ $(${inner_sort_function} "${!myvar}") = ${sort_item} ]]; then
				cecho "$($inner_sort_function "${!myvar}") = ${sort_item}"
				arrayvar=${out_array}[${lineno}]
				eval ${arrayvar}='"${!myvar}"'
				(( lineno++ ))
			else
				cecho "$($inner_sort_function "${!myvar}") != ${sort_item}"
			fi
		done
	done
}

bash_version

# Check for wrong number of arguments
if [[ $# -lt 1 ]] || [[ $# -gt 2 ]]; then
	usage
fi

# Handle command-line options
while getopts r options; do
	case ${options} in
		r)	REVERSE=1
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
	cecho "  No multiple conformations found! Quitting ..."
	exit 1
fi
cecho "  Conformations found: ${count}" b

# Hold PDB contents
declare -a pdb_array
# Set up array to hold current set of multiple conformations
declare -a multi_array
# Array for sorting multi_array in
declare -a tmp_array
# Save size of array, so we don't need to unset elements beyond the end of
# the current array and generally make things easier
multi_lineno=0

# Set name of input file
file=${1}
outfile="${1%.pdb}-conv.pdb"

# Read file into array
lineno=0 # Lines match array elements, not PDB lines
while read line; do
	pdb_array[${lineno}]="${line}"
	(( lineno++ ))
done < ${file}

cecho "Read ${#pdb_array[@]} lines"

# Scan through each line in the array
for (( i = 0; i < ${#pdb_array[@]}; i++ )); do
	arrayline=${pdb_array[${i}]}
#	cecho "Current arrayline:"
#	cecho "    ${arrayline}"

	# If we're not a coordinate, or if we are a coordinate but aren't a
	# multiple conformation
	if ! is_coord "${arrayline}" || ! is_multi "${arrayline}"; then
#		cecho "Either !ATOM or !MULTI:"
#		cecho "    ${arrayline}"
		# Handle what happens when we were in multi, but the next thing is
		# single. This should be where we sort and write out the cached multi's.
		if [[ ${on_multi} -eq 1 ]]; then
			# Sort cached multi's -- allow for forward and reverse cases

			if [[ ${REVERSE} -eq 1 ]]; then
			# In ABAB->AABB case, sort by coord_num first, then by conformation
			# within each coord_num.
				cecho "Sorting ABAB->AABB"
				sort_array multi_array tmp_array get_coord_num
				# sort by conf
#				sort_array_within tmp_array multi_array get_multi_letter get_coord_num
			else
			# In AABB->ABAB case, sort by conformation first, then by coord_num
			# within each conformation.
				cecho "Sorting AABB->ABAB"
				sort_array multi_array tmp_array get_multi_letter
				# sort by coord
#				sort_array_within tmp_array multi_array get_coord_num get_multi_letter
			fi

			# Write out cached multi's
			for (( j=0; j<multi_lineno; j++ )); do
				echo "${multi_array[${j}]}" >> "${outfile}"
			done

			# Reset multi cache
			multi_lineno=0
			unset multi_array[@]
			unset tmp_array[@]
		fi

		# Write out single conformation
		echo "${arrayline}" >> "${outfile}"
		on_multi=0
	else
		cecho "Multiple conformation: $(get_coord_num "${arrayline}") $(get_multi_letter "${arrayline}")"
#		cecho "    ${arrayline}"
		# We are a multiple conformation
		multi_array[${multi_lineno}]="${arrayline}"
		(( multi_lineno++ ))
		on_multi=1
	fi
done
