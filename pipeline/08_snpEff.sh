#!/usr/bin/bash
#SBATCH --mem=64G -p batch --nodes 1 --ntasks 8 --out logs/snpEff.log
module unload miniconda2
module load miniconda3
module load snpEff
module load bcftools/1.11
module load tabix
module load yq

# THIS IS AN EXAMPLE OF HOW TO MAKE SNPEFF - it is for A.fumigatus
SNPEFFGENOME=AfumigatusAf293_FungiDB_50
GFFGENOME=$SNPEFFGENOME.gff
CPU=$SLURM_CPUS_ON_NODE
if [ -z $CPU ]; then
  CPU=1
fi
MEM=64g
TOPFOLDER=`pwd` # expecting to be run in top folder of the github checkout
DOMAINS=$TOPFOLDER/genome/$(basename ${GFFGENOME} .gff)_InterproDomains.txt

# this module defines SNPEFFJAR and SNPEFFDIR
if [ -f config.txt ]; then
  source config.txt
fi
GFFGENOMEFILE=$GENOMEFOLDER/$GFFGENOME
FASTAGENOMEFILE=$GENOMEFOLDER/$GENOMEFASTA
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

if [[ -z $POPYAML || ! -s $POPYAML ]]; then
  echo "Need populations.yaml file and variable POPYAML set in the config.txt"
  exit
fi

mkdir -p $SNPEFFOUT
## NOTE YOU WILL NEED TO FIX THIS FOR YOUR CUSTOM GENOME
if [ ! -e $SNPEFFOUT/$snpEffConfig ]; then
  rsync -a $SNPEFFDIR/snpEff.config $SNPEFFOUT/$snpEffConfig
  echo "# AfumAf293.fungidb " >> $SNPEFFOUT/$snpEffConfig
  # CHANGE Aspergillus fumigatus Af293 FungiDB to your genome name and source - though this is really not important - $SNPEFFGENOME.genome is really what is used
  echo "$SNPEFFGENOME.genome : Aspergillus fumigatus Af293 FungiDB" >> $SNPEFFOUT/$snpEffConfig
  chroms=$(grep '##sequence-region' $GFFGENOMEFILE | awk '{print $2}' | perl -p -e 's/\n/, /' | perl -p -e 's/,\s+$/\n/')
  echo -e "\t$SNPEFFGENOME.chromosomes: $chroms" >> $SNPEFFOUT/$snpEffConfig
  # THIS WOULD NEED SPEIFIC FIX BY USER - IN A.fumigatus the MT contig is called mito_A_fumigatus_Af293
  echo -e "\t$SNPEFFGENOME.mito_A_fumigatus_Af293.codonTable : Mold_Mitochondrial" >> $SNPEFFOUT/$snpEffConfig
  mkdir -p $SNPEFFOUT/data/$SNPEFFGENOME
  gzip -c $GFFGENOMEFILE > $SNPEFFOUT/data/$SNPEFFGENOME/genes.gff.gz
  rsync -aL $REFGENOME $SNPEFFOUT/data/$SNPEFFGENOME/sequences.fa

  java -Xmx$MEM -jar $SNPEFFJAR build -datadir $TOPFOLDER/$SNPEFFOUT/data -c $TOPFOLDER/$SNPEFFOUT/$snpEffConfig -gff3 -v $SNPEFFGENOME
fi
# get full path to YAML file
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
  java -Xmx$MEM -jar $SNPEFFJAR eff -c $TOPFOLDER/$SNPEFFOUT/$snpEffConfig -dataDir $TOPFOLDER/$SNPEFFOUT/data -v $SNPEFFGENOME $INVCF > $OUTVCF

  bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT{0}[\t%TGT]\t%INFO/ANN\n' $OUTVCF > $OUTTAB

  # this requires python3 and vcf script
  # this assumes the interpro domains were downloaded from FungiDB and their format - you will need to generalize this
  $TOPFOLDER/scripts/map_snpEff2domains.py --vcf $OUTVCF --domains $DOMAINS --output $DOMAINVAR

  # this requires Python and the vcf library to be installed
  $TOPFOLDER/scripts/snpEff_2_tab.py $OUTVCF $REFGENOME > $OUTMATRIX
  popd
}
source $(which env_parallel.bash)
export -f makeMatrix
env_parallel -j $CPU --env _ makeMatrix ::: $(yq eval '.Populations | keys' $POPYAML | perl -p -e 's/^\s*\-\s*//' )
