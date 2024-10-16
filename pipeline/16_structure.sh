#!/usr/bin/bash -l
#SBATCH -c 8 --mem 8gb --out logs/structure.log

module load plink
module load yq
module load faststructure

OUT=structure

mkdir -p $OUT


CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi
if [ -f config.txt ]; then
	source config.txt
else
	echo "need a config.txt"
fi

if [ -z $FINALVCF ]; then
	echo "Need to define FINALVCF in config.txt"
	exit
fi

for POPNAME in $(yq eval '.Populations | keys' $POPYAML | perl -p -e 's/^\s*\-\s*//')
do
    # create a filtered VCF containing only invariant sites
    for TYPE in SNP
    do
	VCF=$FINALVCF/$PREFIX.$POPNAME.$TYPE.combined_selected.vcf.gz
	plink --vcf $VCF --const-fid --allow-extra-chr  --vcf-idspace-to _ --keep-allele-order --make-bed --out $OUT/$PREFIX.$POPNAME.$TYPE
	K=${SLURM_ARRAY_TASK_ID}
	parallel -j 8 structure.py -K {} --seed 121 --input=$OUT/$PREFIX.$POPNAME.$TYPE --output=$OUT/$PREFIX.$POPNAME.$TYPE ::: $(seq 8)
	chooseK.py --input=$OUT/$PREFIX.$POPNAME.$TYPE
    done
done

Rscript scripts/plot_distruct.R
