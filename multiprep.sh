#!/bin/bash

# Purpose: Wrapper script for scan_pdb_syntax.sh, require_replicate_atoms.sh and
# handle_multiple_conformations.sh. This will:
#
# 	1) Check for errors in PDB syntax;
# 	2) Output *-replicate PDB files with all atoms not common to all PDBs
# 		deleted;
# 	3) Create a single concatenated PDB with a number of single-conformation
# 		structures from an individual PDB with multiple conformations.
#
# The output single PDB can then be fed into multicore, multilore and multigore
# for further analysis.
#
# Author: Donnie Berkholz <berkhold@science.oregonstate.edu>
# Karplus lab, Oregon State University


# Use case:
# We come in with 2 PDB files, from two different states. We want to compare
# them. We need to output a single PDB with only replicate atoms, and with
# multiple conformations split out. We need to notify the user of how many
# conformations of each original PDB were produced, so the user knows correct
# input information for multi*.
#
# Handling this: Take any number of input PDBs, scan, check for replicates.
# Output to same number of *-replicate.pdb files. With each file, then split
# out the multiple conformations. All of this is handled by the scripts.
#
# Next, we need to produce a single PDB from all of these replicate-only,
# single-conformation PDB files. This is not handled by the scripts.
#
# Beginning with two PDB files 1.pdb and 2.pdb, this is the process:
# 	1) Scan results in no output; may require manual fixing
# 	2) Replicate run:
#		Input: 1.pdb and 2.pdb
# 		Outpit: 1-replicate.pdb and 2-replicate.pdb
# 	3) Multiple run:
#		Input: 1-replicate.pdb and 2-replicate.pdb
#		Output: Any number of 1-replicate[A-Z].pdb and 2-replicate[A-Z].pdb
#			Also 1-replicate-all.pdb and 2-replicate-all.pdb
#
# We need to count the number of MODEL entries and output this. We then need
# to clip END off all MODEL entries, concatenate them and tag END onto the end
# of the file.
# 	Input: 1-replicate-all.pdb and 2-replicate-all.pdb
# 	Output: 1-2-multi.pdb and 1-2-multi.txt (conformation counts)
#
# This should be as easy as:
# 	1) grep for instances of ^MODEL in each PDB, count total
# 	2) sed /^END$/d to delete END but not ENDMDL
# 	3) cat the PDBs, echo END >> multi.pdb
#
# In addition, to be a good citizen we should put everything in a temporary
# directory instead of polluting the current one, and just put the end output
# in the current directory. This might be nice to have configurable.

OUTFILE="multi.pdb"
PREPPATH="${0%/*}/prep/"

source ${PREPPATH}/color.sh
source ${PREPPATH}/fix_model.sh
source ${PREPPATH}/bash_check.sh

color_setup

# void usage(void)
usage() {
	echo
	cecho "  ${0##*/} [ -t ] [ -m ] pdb pdb ..." b
	echo
	cecho "  Prepare PDB files for input into multi*. Takes any number of PDB"
	cecho "  files as arguments."
	echo
	cecho "  -t" green b
	cecho "    Save temporary work directory."
	echo
	cecho "  -m" green b
	cecho "    Split multiple conformations into separate models, each with a"
	cecho "    single conformation."
	exit 1
}

bash_version

# Check for wrong number of arguments
if [[ $# -lt 1 ]]; then
	usage
fi

# Handle command-line options
while getopts tm options; do
	case ${options} in
		t)	SAVE_TEMP=1
			;;
		m)	HANDLE_MULTI=1
			;;
		*)	usage
			;;
	esac
done
# Reset post-option stuff to positional parameters
shift $((OPTIND - 1))

# Error checking -- are inputs real files, etc.
for i in $@; do
	if [[ ! -e $i ]]; then
		cecho "$i doesn't exist!" red b
		usage
	fi
done

if [[ -e "${OUTFILE}" ]]; then
	cecho "${OUTFILE} already exists! Remove it to continue."
	exit 1
fi

# Set up temp directory
if [[ -e /bin/mktemp ]]; then
	TMPDIR="$(/bin/mktemp -d /tmp/multiprep.XXXXXX)"
else
	TMPDIR="/tmp/multiprep.$$"
fi

# Actual run of scripts
for i in $@; do
	${PREPPATH}/scan_pdb_syntax.sh $i
	if [[ $? -ne 0 ]]; then
		cecho "Error encountered in ${PREPPATH}/scan_pdb_syntax.sh"
		exit 1
	fi
done

ORIGDIR="$(pwd)"
pushd ${TMPDIR} &> /dev/null



# Clean out waters, to avoid problems with randomly inserted high-numbered
# waters in the middle of chains.
cecho "Cleaning out water, crystal, anisotropic B, remarks, etc ... "
for i in $@; do
	grep -v \
		-e HOH \
		-e CRYST \
		-e SCALE \
		-e ANISOU \
		-e REMARK \
		-e AUTHOR \
		-e COMPND \
		-e CONECT \
		-e FORMUL \
		-e FTNOTE \
		-e HEADER \
		-e HELIX \
		-e 'HET ' \
		-e JRNL \
		-e MASTER \
		-e ORIGX1 \
		-e ORIGX2 \
		-e ORIGX3 \
		-e REVDAT \
		-e SEQRES \
		-e SHEET \
		-e SOURCE \
		-e SPRSDE \
		-e SSBOND \
		${ORIGDIR}/${i} \
		> ${i##*/}
	LOCS="${LOCS} ${i##*/}"
done


cecho "Cleaning out atoms not common to all PDBs ... "
if [[ -n "${SAVE_TEMP}" ]]; then
	cecho "Output files will be at ${TMPDIR}/*-replicate.pdb"
else
	cecho "Not saving intermediate replicate output files."
fi
${PREPPATH}/require_replicate_atoms.sh ${LOCS}
if [[ $? -ne 0 ]]; then
	popd
	cecho "Error encountered in ${PREPPATH}/require_replicate_atoms.sh"
	exit 1
fi

if [[ -n "${HANDLE_MULTI}" ]]; then
	cecho "Splitting out multiple conformations into multiple models ... "
	for i in ${LOCS}; do
		${ORIGDIR}/${PREPPATH}/handle_multiple_conformations.sh -c ${i%.pdb}-all.pdb \
			 ${i%.pdb}-replicate.pdb
		if [[ $? -ne 0 ]]; then
			cecho "Error encountered in ${PREPPATH}/handle_multiple_conformations.sh"
			exit 1
		fi
	done
fi

# Concatenate all PDBs for suitable multi* input
echo
cecho "Creating ${OUTFILE} containing all models ... "

# The PDBs to concatenate change, depending on whether we split out multi's.
if [[ -n "${HANDLE_MULTI}" ]]; then
	SUFFIX="all"
else
	SUFFIX="replicate"
fi

for CONF in ${LOCS}; do
	# Count MODEL instances
	MODELTOT="$(grep '^MODEL' ${CONF%.pdb}-${SUFFIX}.pdb | wc -l)"
	cecho "Total models for ${CONF%.pdb}-${SUFFIX}.pdb: ${MODELTOT}"
	echo "${CONF%.pdb}-${SUFFIX}.pdb ${MODELTOT} models" >> ${ORIGDIR}/${CONF%.pdb}-${SUFFIX}.txt
	# Get rid of END
	sed -i '/^END$/d' ${CONF%.pdb}-${SUFFIX}.pdb
	# Add to multi.pdb
	cat ${CONF%.pdb}-${SUFFIX}.pdb >> ${OUTFILE}
done

echo "END" >> ${OUTFILE}

# Fix MODEL numbers
model=0
fix_model ${OUTFILE} ${ORIGDIR}/${OUTFILE}

# Cleanup
# Remove temp files, etc. unless (not implemented) told to save them.
echo
if [[ ! -n "${SAVE_TEMP}" ]]; then
	cecho "Cleaning temp directory ..."
	rm -rf ${TMPDIR}
else
	cecho "Temp directory is ${TMPDIR}."
fi

popd &> /dev/null
