#!/usr/bin/bash
#SBATCH --mem=1gb --ntasks 1 --nodes 1
#SBATCH -p short
#SBATCH -J maketreePoppr --out logs/setup_tree_poppr.log

module load yq

CPU=2
if [ $SLURM_CPUS_ON_NODE ]; then
  CPU=$SLURM_CPUS_ON_NODE
fi

if [[ -f config.txt ]]; then
  source config.txt
else
  echo "Need a config.txt"
  exit
fi

if [[ -z $REFNAME ]]; then
  REFNAME=REF
fi

if [[ -z $POPYAML || ! -s $POPYAML ]]; then
  echo "Cannot find \$POPYAML variable - set in config.txt"
  exit
fi

module load bcftools
module load samtools
module load workspace/scratch


mkdir -p $TREEDIR
for POPNAME in $(yq eval '.Populations | keys' $POPYAML | grep All | perl -p -e 's/^\s*\-\s*//' )
do
  for TYPE in SNP
  do 
    root=$FINALVCF/$PREFIX.$POPNAME.$TYPE.combined_selected
    vcf=$root.vcf.gz
    tree=$TREEDIR/$PREFIX.$POPNAME.$TYPE.poppr.upgma.tre
    if [[ ! -s $tree || ${vcf} -nt $tree ]]; then
	    sbatch -p batch -N 1 -n 8 --mem 32gb --out logs/make_poppr_$POPNAME.upgma.%A.log --wrap "time Rscript ./scripts/poppr_tree.R --vcf $vcf --tree $tree --method upgma"
    fi
    tree=$TREEDIR/$PREFIX.$POPNAME.$TYPE.poppr.nj.tre
    if [[ ! -s $tree || ${vcf} -nt $tree ]]; then
	    sbatch -p batch -N 1 -n 8 --mem 32gb --out logs/make_poppr_$POPNAME.nj.%A.log --wrap "time Rscript ./scripts/poppr_tree.R  --vcf $vcf --tree $tree --method nj"
    fi
  done
done
