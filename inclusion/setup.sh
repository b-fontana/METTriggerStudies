#!/usr/bin/env bash

### Load functions
THIS_FILE="${BASH_SOURCE[0]}"
THIS_DIR="$( cd "$( dirname ${THIS_FILE} )" && pwd )"
source "${THIS_DIR}/lib/funcs.sh"
source "${THIS_DIR}/lib/clean_completion.bash"

export PATH="${PATH}:${THIS_DIR}/lib/"

chmod +x ${THIS_DIR}/clean.sh
alias incl_clean=${THIS_DIR}/clean.sh
