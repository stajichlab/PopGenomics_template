#!/usr/bin/bash -l
#SBATCH -p short -c 1 -N 1 -n 1 --mem 4gb --out logs/build_env.log

if [ ! -d ./env ]; then
  c=$(which mamba)
  if [ -z $c ]; then
    c=$(which conda)
  fi
  $c env create -p ./env -f environment.yml
fi
