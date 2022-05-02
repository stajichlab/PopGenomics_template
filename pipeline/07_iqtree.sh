#!/usr/bin/bash -l
#SBATCH -p intel -N 1 -n 8 --mem 16gb  --out logs/IQTREE2.log

module load iqtree/2.2.0
TREEDIR=strain_tree
source config.txt
if [ -z $TREEDIR ]; then
  echo "Need a TREEDIR defined in config.txt"
  exit
fi
for POPNAME in $(yq eval '.Populations | keys' $POPYAML | perl -p -e 's/^\s*\-\s*//')
do
  for TYPE in SNP INDEL
  do
    iqtree2 -s $TREEDIR/$PREFIX.SNP.mfa -m GTR+ASC -nt AUTO -bb 1000 -alrt 1000
   done
 done
