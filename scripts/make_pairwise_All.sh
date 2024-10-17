#!/usr/bin/bash -l
#SBATCH -p short -N 1 -c 16 --mem 96gb --out logs/pairwise_SNPs.log

module load bcftools
module load workspace/scratch
source config.txt

IN=vcf/$PREFIX.All.SNP.combined_selected.vcf.gz
CPU=16
OFILE=pairwise_All_SNP_v6_count.tsv
echo -e "STRAIN\tCOUNT" > $OFILE
for strain in $(bcftools query -l $IN)
do
	bcftools view --threads $CPU -s $strain -Ob -o $SCRATCH/$strain.bcf $IN
	bcftools +fill-tags $SCRATCH/$strain.bcf -Ob -o $SCRATCH/$strain.tags.bcf -- -t all
	count=$(bcftools view -e "AF=0" $SCRATCH/$strain.tags.bcf | grep -c PASS)
	echo -e "$strain\t$count"
done >> $OFILE

IN=vcf/$PREFIX.All.INDEL.combined_selected.vcf.gz
OFILE=pairwise_All_INDEL_v6_count.tsv
echo -e "STRAIN\tCOUNT" > $OFILE
for strain in $(bcftools query -l $IN)
do
	bcftools view --threads $CPU -s $strain -Ob -o $SCRATCH/$strain.bcf $IN
	bcftools +fill-tags $SCRATCH/$strain.bcf -Ob -o $SCRATCH/$strain.tags.bcf -- -t all
	count=$(bcftools view -e "AF=0" $SCRATCH/$strain.tags.bcf | grep -c PASS)
	echo -e "$strain\t$count"
done >> $OFILE



