#!/usr/bin/bash
#SBATCH --mem=24gb --ntasks 24 --nodes 1
#SBATCH --time=2:00:00 -p short
#SBATCH -J maketree --out logs/make_tree.log

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

module load parallel
module load bcftools
module load samtools
module load iqtree
module load fasttree
module load workspace/scratch

print_fas() {
  printf ">%s\n%s\n" $1 $(bcftools query -s $1 -f '[%IUPACGT]' $2)
}

iqtreerun() {
	in=$1
	out=$in.treefile
	if [[ ! -f $out || $in -nt $out ]]; then
		sbatch -p intel -n 6 -N 1 --mem 16gb -J iqtree --wrap "module load IQ-TREE/2.1.1; iqtree2 -m GTR+ASC -s $in -st DNA -nt AUTO -bb 1000 -alrt 1000"
	fi
}

fasttreerun() {
        in=$1
	out=$(echo $in | perl -p -e 's/\.mfa/.fasttree.tre/')
        if [[ ! -f $out || $in -nt $out ]]; then
                sbatch -p short -n 32 -N 1 --mem 16gb -p short -J FastTree --wrap "module load fasttree; FastTreeMP -gtr -gamma -nt < $in > $out"
        fi
}

export -f print_fas fasttreerun iqtreerun
mkdir -p $TREEDIR
for POPNAME in $(yq eval '.Populations | keys' $POPYAML | perl -p -e 's/^\s*\-\s*//' )
do
  for TYPE in SNP
  do
    root=$FINALVCF/$PREFIX.$POPNAME.$TYPE.combined_selected
    FAS=$TREEDIR/$PREFIX.$POPNAME.$TYPE.mfa

    if [ -f $root.vcf ]; then
      bgzip $root.vcf
      tabix $root.vcf.gz
    fi

    vcf=$root.vcf.gz
    if [[ ! -f $FAS || ${vcf} -nt $FAS ]]; then
      vcftemp=$SCRATCH/$PREFIX.$POPNAME.$TYPE.vcf.gz
      bcftools filter -Oz -o $vcftemp --SnpGap 3 -e 'QUAL < 1000 || AF=1 || INFO/AF < 0.05 || F_MISSING > 0' $vcf
      bcftools index $vcftemp
      # no ref genome alleles
      printf ">%s\n%s\n" $REFNAME $(bcftools query -f '%REF' $vcftemp) > $FAS
      parallel -j $CPU print_fas ::: $(bcftools query -l ${vcf}) ::: $vcftemp >> $FAS
#      perl -ip -e 'if(/^>/){s/[\(\)#]/_/g; s/_+/_/g } else {s/[\*.]/-/g }' $FAS
    fi
  done
done
#parallel -j 2 fasttreerun ::: $(ls $TREEDIR/*.mfa)
#parallel -j 4 iqtreerun ::: $(ls $TREEDIR/*.mfa)
