#!/bin/sh
# generic script to prepare the prmf conda environment for a program and then run the 
# program in that environment

# TODO implement getoptlong style interface
usage="prmf_runner.sh <env> <cmd> <args>"

# act like a login shell, load conda (and other) environment
# TODO what if the conda setup is not in ~/.bash_profile (but is in ~/.bashrc for example)
source ~/.bash_profile

# non-login bash shells such as the one expected to be created when invoking this script
# do not source conda's bash configuration files, source them if needed
# TODO sourcing ~/.bash_profile renders this unnecessary 
type -a conda | grep -q 'conda is a function'
if [ "$?" != "0" ]; then
  if [ -f $_CONDA_ROOT/etc/profile.d/conda.sh ]; then
    source $_CONDA_ROOT/etc/profile.d/conda.sh
  else
    echo 'cannot locate conda bash configuration file at $_CONDA_ROOT/etc/profile.d/conda.sh' 1>&2
    echo "value of \$_CONDA_ROOT: $_CONDA_ROOT" 1>&2
    exit 21
  fi
fi

env="$1"
conda activate "$env"
if [ "$?" != "0" ]; then
  # conda will print an error message
  exit 22
fi

cmd="$(which "$2")"
args=""
for arg in "$@"; do
  if [ "$arg" != "$1" ] && [ "$arg" != "$2" ]; then
    args="${args} ${arg}"
  fi
done

$cmd $args
