#!/usr/bin/env bash

function setup() {
	cmsenv
	return 0
}

function action() {
    if setup "$@"; then
        echo -e "\x1b[0;49;35msf_inclusion successfully setup\x1b[0m"
        return "0"
    else
        local code="$?"
        echo -e "\x1b[0;49;31msf_inclusion failed with code ${code}\x1b[0m"
        return "$code"
    fi
}

action "$@"
