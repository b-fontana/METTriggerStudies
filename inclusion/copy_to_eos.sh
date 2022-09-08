#!/usr/bin/env bash

### Argument parsing help strings
help_description="prints this help message"
tag_description="select tags"
eos_description="/eos/ username (same used in lxplus)"
dryrun_description="dryrun: prints everything, runs nothing"
function print_usage_workflowClean {
    usage=" $(basename "$0") [-h] [-t -d]: removes folders
where:
    -h / --help  [ ${help_description} ]
    -t / --tag   [ ${tag_description} ]
    -f / --full  [ ${full_description} ]
    --eos  [ ${eos_description} ]
    -d / --dryrun [ ${dryrun_description} ]

    Run example: bash $(basename "$0") -t 000 -t 111 -t 222 -t 654 -f -d
"
    printf "${usage}"
}

### Argument parsing
DRYRUN=true # safety during script development
EOS_USER="bfontana"
declare -a TAGS;
while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
	-h|--help)
	    print_usage_workflowClean
	    exit 1
	    ;;	
	-t|--tag)
	    TAGS+=("$2")
	    shift # past argument
	    shift # past value
	    ;;
	--eos)
	    EOS_USER="$2"
	    shift # past argument
	    shift # past value
	    ;;
	-d|--dryrun)
	    DRYRUN=true
	    shift # past argument
	    ;;
	*)    # unknown option
	    echo "Wrong parameter ${1}."
	    exit 1
	    ;;
    esac
done

### Load functions
source "lib/funcs.sh"

### General parameters
BASE_LOCAL_DIR="/data_CMS/cms/${USER}/TriggerScaleFactors"
BASE_EOS_DIR="/eos/user/${EOS_USER:0:1}/${EOS_USER}/www/TriggerScaleFactors"
declare -a CHANNELS=( "etau" "mutau" "tautau" )

declare -a OLD_TAGS=( $(/bin/ls -1 "${BASE_PATH}") )
check_tags ${#TAGS[@]} ${TAGS[@]} ${OLD_TAGS[@]}

declare -a REMOVED_TAGS;
for tag1 in ${TAGS[@]}; do
	flag=0
	for tag2 in ${OLD_TAGS[@]}; do
		[[ ${tag1} == ${tag2} ]] && flag=1
	done
	if [ ${flag} -eq 0 ]; then
		REMOVED_TAGS+=( ${tag1} )
	fi
done

# print tags that do not exist and take them out of the list of tags to remove
NREM=${#REMOVED_TAGS[@]}
if [ ${NREM} -eq 1 ]; then
	echo "Tag ${REMOVED_TAGS[1]} was ignored since it is not present."
elif [ ${NREM} -ge 1 ]; then
	list_ignored_tags ${NREM} ${REMOVED_TAGS[@]}
	list_tags ${OLD_TAGS[@]}
fi
for tag in ${REMOVED_TAGS[@]}; do
	TAGS=( ${TAGS[@]/${tag}} ) # remove one element from the array
done

### Argument parsing: information for the user
echo "### Arguments"
echo "TAGS  = ${TAGS[*]}"
echo "DRYRUN = ${DRYRUN}"
echo "FULL  = ${FULL}"
echo "#############"

for tag in ${TAGS[@]}; do
	for chn in ${CHANNELS[@]}; do
		comm="cp ${BASE_LOCAL_DIR}/${tag}/Outputs/${chn}/Counts_Data_MET_EG_Mu/table.csv ${BASE_EOS_DIR}/New2/table.csv"
		[[ ${DRYRUN} -eq 1 ]] && echo ${comm} || ${comm}
	done
done
