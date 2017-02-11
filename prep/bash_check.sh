# bash_check.sh
# For use by scripts in sourcing (.) only. Cannot run as a standalone script.

# We require >=bash-2.04 for the constructs used here. 2.04 was released in
# March 2000, so we should work on almost any system.
# void bash_version(void)
bash_version() {
	if [[ -n "${BASH_VERSION}" ]]; then
		if [[ ${BASH_VERSINFO[0]} -ge 3 ]]; then
			true
		elif [[ ${BASH_VERSINFO[1]} -eq 2 ]]; then
			if [[ ${BASH_VERSINFO[2]} -eq 0 ]]; then
				if [[ ${BASH_VERSINFO[2]} -le 3 ]]; then
					cecho "Bash version too old! This script requires >=2.04"
					exit 1
				fi
			fi
		else
			cecho "Bash version too old! This script requires >=2.04"
			exit 1
		fi
	else
		cecho "BASH_VERSION variable unset! Either the bash version is too old"
		cecho "(previous to 2.04) or the environment has been modified."
	fi
}
