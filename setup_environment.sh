#!/bin/bash

function load_module() {
  local MODULE_NAME="$1"
  module load "${MODULE_NAME}" > /dev/null 2>&1
  if [[ "$?" -ne 0 ]]; then
    echo "failed to: module load ${MODULE_NAME}" 1>&2
    return 1
  fi
}

function load_modules() {
  module purge || return 1
  load_module cmake/3.15.4 || return 1
  load_module GCC/5.4.0-2.26 || return 1
  load_module gsl/2.5 || return 1
}

function main() {
  # For using conda
  source "${HOME}/.bashrc" || return 1

  # Uncomment to use environment modules
  # load_modules || return 1
}

main "$@"
