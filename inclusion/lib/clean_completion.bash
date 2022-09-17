#!/usr/bin/env bash

# Source: https://unix.stackexchange.com/questions/322047/bash-argument-autocomplete

# Define a base directory for relative paths.
BASE_PATH="/data_CMS/cms/alves/TriggerScaleFactors"

# Define the completion function
function _clean_completion() {

    # Define local variables to store adjacent pairs of arguments
    local prev_arg;
    local curr_arg;

    # If there are at least two arguments then we have a candidate
    # for path-completion, i.e. we need the option flag '-d' and
    # the path string that follows it.
    if [[ ${#COMP_WORDS[@]} -ge 2 ]]; then

        # Get the current and previous arguments from the command-line
        prev_arg="${COMP_WORDS[COMP_CWORD-1]}";
        curr_arg="${COMP_WORDS[COMP_CWORD]}";

        # We only want to do our custom path-completion if the previous
        # argument is the '-d' option flag
        if [[ "${prev_arg}" = "-t" ]]; then

            # We only want to do our custom path-completion if a base
            # directory is defined and the argument is a relative path
            if [[ -n "${BASE_PATH}" && ${curr_arg} != /* ]]; then

                # Generate the list of path-completions starting from BASE_PATH
                COMPREPLY=( $(compgen -d -o default -- "${BASE_PATH}/${curr_arg}") );

                # Don't append a space after the command-completion
                # This is so we can continue to apply completion to subdirectories
                compopt -o nospace;

                # Return immediately
                return 0;
            fi
        fi
    fi

    # If no '-d' flag is given or no base directory is defined then apply default command-completion
    COMPREPLY=( $(compgen -o default -- "${curr_arg}") );
    return 0;
}

# Activate the completion function
complete -F _clean_completion clean
