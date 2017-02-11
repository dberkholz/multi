# color.sh
# For use by scripts in sourcing (.) only. Cannot run as a standalone script.

color_setup() {
	# Color setup, works only for interactive sessions
	# Can't quote ${PS1}
	if [[ -n ${PS1} ]]; then
		color="yes"
	else
		unset color
	fi

	bold="${color:+\033[1m}"
	unbold="${color:+\033[0m}"
	# shortcuts
	b="${bold}"
	ub="${unbold}"

	red="${color:+\E[40;31m}"
	green="${color:+\E[40;32m}"
	yellow="${color:+\E[40;33m}"
	blue="${color:+\E[40;34m}"
	magenta="${color:+\E[40;35m}"
	cyan="${color:+\E[40;36m}"
	white="${color:+\E[40;37m}"

	# Reset text attributes to normal without clearing screen.
	reset() {
		tput sgr0
	}

	# Color echo.
	# void cecho(char *message, char *color)
	cecho () {
		local default_msg="No message passed."
		message=${1:-$default_msg}
		color=${!2:-${!black}}
		color2=${!3}
	
		echo -ne "${color}${color2}"
		echo "$message"
		reset
	}
}
