# fix_model.sh
# For use by scripts in sourcing (.) only. Cannot run as a standalone script.

# Fix MODEL numbers
# void fix_model(char * input, char * output)
fix_model() {
	while read line; do
		# Check for MODEL number
		numbered=0
		# If we're in a MODEL line, and it is numbered, set flag
		if [[ "${line/#MODEL}" != "${line}" ]]; then
			if [[ "${line%[0-9]}" != "${line}" ]]; then
				numbered=1
			fi
			(( model++ ))
		fi
		# Set up MODEL number
		line=${line/#MODEL/MODEL ${model}}
		# Strip old MODEL number from end; shortest match to " #$"
		if [[ "${numbered}" -eq 1 ]]; then
			line=${line% [0-9]}
		fi
		echo "${line}" >> ${2}
	done < ${1}
}
