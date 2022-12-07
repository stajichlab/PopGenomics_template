#!/usr/bin/bash -l
#SBATCH --mem=64G -N 1 -n 1 -c 2 --out logs/snpEff.log

# this module defines SNPEFFJAR and SNPEFFDIR
module load snpEff/4.3t
module load bcftools
module load yq

CPU=$SLURM_CPUS_ON_NODE
if [ -z $CPU ]; then
  CPU=1
fi
# Fix this for your genome
SNPEFFGENOME=Clusitianiae
GFFGENOME=candida_lusitaniae_1_fixed.gff3
TOPFOLDER=`pwd` # expecting to be run in top folder of the github checkout

MEM=64g

if [ -f config.txt ]; then
	source config.txt
fi
FASTAGENOMEFILE=$REFGENOME
GFFGENOMEFILE=$GFFGENOME
echo "GFF is $GFFGENOMEFILE"
if [ -z $SNPEFFJAR ]; then
 echo "need to defined \$SNPEFFJAR in module or config.txt"
 exit
fi
if [ -z $SNPEFFDIR ]; then
 echo "need to defined \$SNPEFFDIR in module or config.txt"
 exit
fi
# could make this a config

if [ -z $FINALVCF ]; then
	echo "need a FINALVCF variable in config.txt"
	exit
fi

mkdir -p $SNPEFFOUT
if [ ! -e $SNPEFFOUT/$snpEffConfig ]; then
	rsync -a $SNPEFFDIR/snpEff.config $SNPEFFOUT/$snpEffConfig
	echo "# Clusitaniae " >> $SNPEFFOUT/$snpEffConfig
	if [ -z $NAME 
  	echo "$SNPEFFGENOME.genome : $SNPEFFGENOME" >> $SNPEFFOUT/$snpEffConfig
	echo "$FASTAGENOMEFILE"
	chroms=$(grep '^>' $FASTAGENOMEFILE | perl -p -e 's/>//; s/\n/, /g' | perl -p -e 's/,\s+$/\n/')
	echo -e "\t$SNPEFFGENOME.chromosomes: $chroms" >> $SNPEFFOUT/$snpEffConfig

	# THIS WOULD NEED SPEIFIC FIX BY USER - IN A.fumigatus the MT contig is called mito_A_fumigatus_Af293
	echo -e "\t$SNPEFFGENOME.Supercontig_1.1.codonTable : Alternative_Yeast_Nuclear" >> $SNPEFFOUT/$snpEffConfig
	echo -e "\t$SNPEFFGENOME.Supercontig_1.2.codonTable : Alternative_Yeast_Nuclear" >> $SNPEFFOUT/$snpEffConfig
	echo -e "\t$SNPEFFGENOME.Supercontig_1.3.codonTable : Alternative_Yeast_Nuclear" >> $SNPEFFOUT/$snpEffConfig
	echo -e "\t$SNPEFFGENOME.Supercontig_1.4.codonTable : Alternative_Yeast_Nuclear" >> $SNPEFFOUT/$snpEffConfig
	echo -e "\t$SNPEFFGENOME.Supercontig_1.5.codonTable : Alternative_Yeast_Nuclear" >> $SNPEFFOUT/$snpEffConfig
	echo -e "\t$SNPEFFGENOME.Supercontig_1.6.codonTable : Alternative_Yeast_Nuclear" >> $SNPEFFOUT/$snpEffConfig
	echo -e "\t$SNPEFFGENOME.Supercontig_1.7.codonTable : Alternative_Yeast_Nuclear" >> $SNPEFFOUT/$snpEffConfig
	echo -e "\t$SNPEFFGENOME.Supercontig_1.8.codonTable : Alternative_Yeast_Nuclear" >> $SNPEFFOUT/$snpEffConfig
	echo -e "\t$SNPEFFGENOME.MT_CBS_6936.codonTable : Mold_Mitochondrial" >> $SNPEFFOUT/$snpEffConfig

	mkdir -p $SNPEFFOUT/data/$SNPEFFGENOME
	pigz -c $GFFGENOMEFILE > $SNPEFFOUT/data/$SNPEFFGENOME/genes.gff.gz
	rsync -aL $REFGENOME $SNPEFFOUT/data/$SNPEFFGENOME/sequences.fa
	java -Xmx$MEM -jar $SNPEFFJAR build -datadir $TOPFOLDER/$SNPEFFOUT/data -c $SNPEFFOUT/$snpEffConfig -gff3 -noCheckCds -noCheckProtein -nodownload -v $SNPEFFGENOME
fi

POPYAML=$(realpath $POPYAML)
REFGENOME=$(realpath $REFGENOME)
pushd $SNPEFFOUT

makeMatrix() {
  POPNAME=$1
  echo "POPNAME is $POPNAME"
  mkdir -p $POPNAME
  COMBVCF="$TOPFOLDER/$FINALVCF/$PREFIX.$POPNAME.SNP.combined_selected.vcf.gz $TOPFOLDER/$FINALVCF/$PREFIX.$POPNAME.INDEL.combined_selected.vcf.gz"
  echo "COMBVCF is '$COMBVCF'"
  for n in $COMBVCF
  do
    echo $n
    st=$(echo $n | perl -p -e 's/\.gz//')
    if [ ! -f $n ]; then
      bgzip $st
    fi
    if [ ! -f $n.tbi ]; then
      tabix $n
    fi
  done

  pushd $POPNAME
  INVCF=$PREFIX.$POPNAME.allvariants_combined_selected.vcf
  OUTVCF=$PREFIX.$POPNAME.snpEff.vcf
  OUTTAB=$PREFIX.$POPNAME.snpEff.tab
  OUTMATRIX=$PREFIX.$POPNAME.snpEff.matrix.tsv
  DOMAINVAR=$PREFIX.$POPNAME.snpEff.domain_variant.tsv
  if [[ ! -s $INVCF || $TOPFOLDER/$FINALVCF/$PREFIX.$POPNAME.SNP.combined_selected.vcf.gz -nt $INVCF ]]; then
    bcftools concat -a -d both -o $INVCF -O v $COMBVCF
  fi
  echo " java -Xmx$MEM -jar $SNPEFFJAR eff -c $TOPFOLDER/$SNPEFFOUT/$snpEffConfig -dataDir $TOPFOLDER/$SNPEFFOUT/data -v $SNPEFFGENOME $INVCF > $OUTVCF"
  java -Xmx$MEM -jar $SNPEFFJAR eff -nodownload -c $TOPFOLDER/$SNPEFFOUT/$snpEffConfig -dataDir $TOPFOLDER/$SNPEFFOUT/data -v $SNPEFFGENOME $INVCF > $OUTVCF

  bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT{0}[\t%TGT]\t%INFO/ANN\n' $OUTVCF > $OUTTAB

  # this requires python3 and vcf script
  # this assumes the interpro domains were downloaded from FungiDB and their format - you will need to generalize this
  # Nedthis only if genome came from FungiDB otherwise we need to have a diff script
  #$TOPFOLDER/scripts/map_snpEff2domains.py --vcf $OUTVCF --domains $DOMAINS --output $DOMAINVAR

  # this requires Python and the vcf library to be installed
  module load pyvcf
  $TOPFOLDER/scripts/snpEff_2_tab.py $OUTVCF $REFGENOME > $OUTMATRIX
  module unload pyvcf
  popd
}
source $(which env_parallel.bash)
export -f makeMatrix
env_parallel -j $CPU --env _ makeMatrix ::: $(yq eval '.Populations | keys' $POPYAML | perl -p -e 's/^\s*\-\s*//' )
