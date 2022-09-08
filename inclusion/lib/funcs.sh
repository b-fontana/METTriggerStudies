#!/usr/bin/env bash

function list_tags() {
	if [ ${@} -eq 1 ]; then
		echo "The following tag is currently available:"
	else
		echo "The following tags are currently available:"
	fi
	for tag in "$@"; do
		echo " - ${tag}"
	done
}

function list_ignored_tags() {
	declare -a tags=( "${@:2}" )
	
	printf "Tags"
	for tag in ${tags[@]}; do
		if [ ${tag} == ${tags[ $((${1}-1)) ]} ]; then
			printf " and '${tag}' "
		elif [ ${tag} == ${tags[0]} ]; then
			printf " '${tag}'"
		else
			printf ", '${tag}'"
		fi
	done
	echo "were ignored since they are not present."
}
