#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 1 -c 96 --mem 8gb --out logs/18_pairiwise_SNPs.log

# this script generates number of SNPs each strain has with the ref

module load bcftools
module load workspace/scratch
module load yq
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

if [[ -z $POPYAML || ! -s $POPYAML ]]; then
	echo "Cannot find \$POPYAML variable - set in config.txt"
	exit
fi

OUTDIR=reports/pairwise_strain_compare
mkdir -p $OUTDIR
for POPNAME in $(yq eval '.Populations | keys' $POPYAML | perl -p -e 's/^\s*\-\s*//')
do
    for TYPE in SNP INDEL
    do
        OUT=$OUTDIR/$PREFIX.$POPNAME.$TYPE.pairwise_count.tsv
        IN=$FINALVCF/$PREFIX.$POPNAME.$TYPE.combined_selected.vcf.gz
        echo -e "STRAIN1\tSTRAIN2\tCOUNT" > $OUT
        parallel -j $CPU if [[ "{1}" != "{2}" ]]\; then bcftools view -s {1},{2} -Ob -o $SCRATCH/{1}-{2}.$TYPE.bcf $IN \; fi ::: $(bcftools query -l $IN | head -n 3) ::: $(bcftools query -l $IN | head -n 3)
        parallel -j $CPU if [[ "{1}" != "{2}" ]]\; then bcftools +fill-tags $SCRATCH/{1}-{2}.$TYPE.bcf -Ob -o $SCRATCH/{1}-{2}.$TYPE.filltags.bcf -- -t all \; fi ::: $(bcftools query -l $IN | head -n 3) ::: $(bcftools query -l $IN | head -n 3)
        parallel -j $CPU if [[ "{1}" != "{2}" ]]\; then bcftools view --exclude='AC=0' -f 'PASS,.' -Ob -o $SCRATCH/{1}-{2}.$TYPE.pair.bcf $SCRATCH/{1}-{2}.$TYPE.filltags.bcf \; fi ::: $(bcftools query -l $IN | head -n 3) ::: $(bcftools query -l $IN | head -n 3)
        parallel -j $CPU if [[ "{1}" != "{2}" ]]\; then bcftools query -f '%CHROM\\t%POS\\t%REF[\\t%TGT]\\n' -o $SCRATCH/{1}-{2}.$TYPE.tsv $SCRATCH/{1}-{2}.$TYPE.pair.bcf \; fi ::: $(bcftools query -l $IN | head -n 3) ::: $(bcftools query -l $IN | head -n 3)
	parallel -j $CPU if [[ "{1}" != "{2}" ]]\; then M=\$\(./scripts/count_pairwise_vcftab.py --input $SCRATCH/{1}-{2}.$TYPE.tsv\) \; echo -e "{}\\t{}\\t{}\\t\$M" >> \$OUT \; fi ::: $(bcftools query -l $IN | head -n 3) ::: $(bcftools query -l $IN | head -n 3)
    done
done
