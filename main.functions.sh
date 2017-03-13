#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 0.1.0
# Description:
# 	functions included in the main.sh script to mangae pipeline steps and inputs
# 
# ------------------------------------------------------------------------------



# CHECK FUNCTIONS ==============================================================

function check_settings() {
	echo -e '\nChecking settings...'

	# Check that there are no missing settings ---------------------------------

	# The experiment ID is required
	if [[ -z $experiment ]]; then
		msg="\n !ERROR!\n An experiment ID is required.\n"
		msg="$msg Please specify the 'experiment' field in the settings file.\n"
		echo -e $msg
		exit 1
	fi

	# The comma-separated conditions are required
	if [[ -z conditions ]]; then
		msg="\n !ERROR!\n At least one condition is required.\n"
		msg="$msg Please specify the 'conditions' field in the settings file.\n"
		echo -e $msg
		exit 1
	fi

	# Check that elective settings are present, otherwise fall back to default -
	

	# Elective settings
	vars=($numbproc $genome $aligner $rmX $rmY $cutsite $mapqthr $platform \
		$umi_length $binsize $binstep $csrange $pthr $emax $eperc)

	# Elective setting names
	fields=('numbproc' 'genome' 'aligner' 'rmX' 'rmY' 'cutsite' 'mapqthr' \
		'platform' 'umi_length' 'binsize' 'binstep' 'csrange' 'pthr' \
		'emax' 'eperc')

	# Elective setting default values
	defaults=(8 'hg19' 'bwa' false true 'AAGCTT' 30 'L' 8 1e6 1e5 40 0 1e-3 20)

	# Cycle counter
	field_id=0
	# Iterate through the elective settings
	while [[ $field_id -lt ${#fields[@]} ]]; do
		# If unset, fall back to default
		if [[ -z ${vars[$field_id]} ]]; then
			msg="\n !WARNING!\n"
			msg="$msg No '${fields[$field_id]}' field in the settings file.\n"
			msg="$msg Using default: ${defaults[$field_id]}\n"
			echo -e $msg
			eval "${fields[$field_id]}=${defaults[$field_id]}"
		fi

		# Increment counter
		field_id=$(($field_id + 1))
	done
}

# STEP FUNCTIONS ===============================================================

function execute_step() {
	# execute_step opt stepname funcname
	# 
	# Executes steps of the pipeline asking for confirmation.
	# If $opt == 1 then the step is automatically executed.
	
	if [ "$#" -lt 3 ]; then
		echo -e "Correct usage: execute_step $opt step_name funcname\n"
		exit 1
	fi
	
	opt=$1
	step=$2
	f=$3

	if [[ $opt == 1 ]]; then
		ans='y'
	else
		msg="\nRun $step?\nYes (y), Skip (s), Abort (a)"
		echo -e $msg
		read -e ans
	fi
	end=0
	while [[ 0 -eq $end ]]; do
		if [[ -z $ans ]]; then
			echo -e $msg
			read -e ans
		elif [[ 'a' == $ans ]]; then
			end=1
			echo "Aborted."
			exit 1
		elif [[ 's' == $ans ]]; then
			echo -e "Skipped $step.\n"
			end=1
		elif [[ 'y' == $ans ]]; then
			echo -e "\n"
			$f
			end=1
		else
			echo -e $msg
			read -e ans
		fi
	done
}

# INPUT FUNCTIONS ==============================================================

function input_int() {
	# input_int default var_name [val]
	# 
	# Gets an integer by input
	
	if [ "$#" -lt 2 ]; then
		echo -e "Correct usage: input_int default var_name [val]\n"
		exit 1
	fi

	default=$1
	k=$2
	if [ "$#" -lt 3 ]; then
		v=''
	else
		v=$3
	fi
	
	# Get mapq threshold by input
	if [[ -z $v ]]; then
		msg="\nPlease provide a $k [$default]"
		echo -e $msg
		read -e v
	fi

	# Default mapq threshold value
	if [[ -z $v ]]; then
		v="$default"
	fi

	# mapq threshold must be an integer
	[[ "$v" =~ ^([+-]?[0-9]+|[+-]?[.0-9]+[e]?[0-9]+)$ ]]; cond=$?
	while [[ $cond -ne true ]]; do
		msg="\n$k must be an integer.\nPlease provide a $k [$default]"
		echo -e $msg
		read -e v

		# Default mapq threshold value
		if [[ -z $v ]]; then
			v="$default"
		fi

		# Update condition
		[[ "$v" =~ ^([+-]?[0-9]+|[+-]?[.0-9]+[e]?[0-9]+)$ ]]; cond=$?
	done

	echo -e "$k: $v\n"
}

function input_fname() {
	# input_fname fname description [val]
	# 
	# Gets a filepath by input
	
	if [ "$#" -lt 2 ]; then
		echo -e "Correct usage: input_fname fname description [val]\n"
		exit 1
	fi

	fname=$1
	description=$2
	if [ "$#" -lt 3 ]; then
		v=''
	else
		v=$3
	fi
	
	if [[ -z $v ]]; then
		msg="\nDo you want to provide $description? [y/n]"
		echo -e $msg
		read -e ans

		end2=0
		while [[ 0 -eq $end2 ]]; do
			if [[ 'n' == $ans ]]; then
				echo -e "Proceeding without $description.\n"
				end2=1
			elif [[ 'y' == $ans ]]; then
				msg2="Path to $description:"
				echo -e $msg2
				read -e v

				if [[ ! -a $v ]]; then
					echo -e "File not found.\n"
					ans=''
				else
					end2=1
				fi
			else
				echo -e $msg
				read -e ans
			fi
		done
	else
		if [[ ! -a $v ]]; then
			echo -e "$fname not found. Proceeding without $description.\n"
			v=''
		else
			echo -e "$fname:\n$v\n"
		fi
	fi
}

################################################################################
