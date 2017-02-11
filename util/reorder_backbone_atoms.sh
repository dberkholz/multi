#!/bin/bash

# This script reorders protein backbones to always be N-CA-C-O, as the PDB
# format requires. Its input file must have multiple conformations as AABB.
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
	cecho "  ${0##*/} pdb" b
	echo
	cecho "  Reorder protein backbones to always be N-CA-C-O,"
	cecho "  as the PDB format requires."
	echo
	cecho "  Input files must be run through"
	cecho "  switch_conformation_styles.sh -r."
	echo
	cecho "  Output file will have '-ordered.pdb' at the end."
	echo
	exit 1
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
## calling.

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
		sort_by="${sort_by} ${sort_function_out}"
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
# sort_array_within(char * in_arrayname, char * out_arrayname,
# 					char * sort_function, char * within_sort_function)
sort_array_within() {
	# call sort_array() somewhere in here on the subarray we generate from
	# within_sort_function
	local i myvar within_sort_function_out within_sort_by
	local in_array="${1}"
	local out_array="${2}"
	local sort_function="${3}"
	local within_sort_function="${4}"
	local array_length=$(eval echo \${#${in_array}[@]})
	# Set up the things to sort by
	for (( i = 0; i < ${array_length}; i++ )); do
		myvar=${in_array}[${i}]
		within_sort_function_out=$(${within_sort_function} "${!myvar}")
		within_sort_by="${within_sort_by} ${within_sort_function_out}"
	done
#	within_sort_by=$(echo ${within_sort_by} | sort -u)

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

bash_version

# Check for wrong number of arguments
if [[ $# -lt 1 ]] || [[ $# -gt 1 ]]; then
	usage
fi

pdb=${1}

pdbbase=$(basename ${pdb} .pdb)

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
outfile="${1%.pdb}-ordered.pdb"

# Read file into array
lineno=0 # Lines match array elements, not PDB lines
while read line; do
	pdb_array[${lineno}]="${line}"
	(( lineno++ ))
done < ${file}


# Scan through each line in the array
for (( i = 0; i < ${#pdb_array[@]}; i++ )); do
	arrayline=${pdb_array[${i}]}

	# If we're not a coordinate, or if we are a coordinate but aren't a
	# multiple conformation
	if ! is_coord "${arrayline}" || ! is_multi "${arrayline}"; then
		# Handle what happens when we were in multi, but the next thing is
		# single. This should be where we sort and write out the cached multi's.
		if [[ ${on_multi} -eq 1 ]]; then
			# Sort cached multi's -- allow for forward and reverse cases

			# In ABAB->AABB case, sort by coord_num first, then by conformation
			# within each coord_num.
			sort_array multi_array tmp_array get_coord_num
			# sort by conf
			sort_array_within tmp_array multi_array get_multi_letter get_coord_num

			# Write out cached multi's
			for (( i=0; i<=multi_lineno; i++ )); do
				echo "${multi_array[${i}]}" >> "${outfile}"
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
		# We are a multiple conformation
		multi_array[${multi_lineno}]="${arrayline}"
		(( multi_lineno++ ))
		on_multi=1
	fi
done
