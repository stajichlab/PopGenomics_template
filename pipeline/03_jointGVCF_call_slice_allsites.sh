#!/usr/bin/bash
#SBATCH --mem 24G --nodes 1 --ntasks 4 -J slice.GVCFGeno --out logs/GVCFGenoGATK4.allsites.slice_%a.%A.log  -a 1-23
hostname
MEM=24g
module unload R
module unload java
module load picard
module load gatk/4.6.0.0
module load bcftools
module load parallel
module load yq
module load workspace/scratch

source config.txt

N=${SLURM_ARRAY_TASK_ID}

if [ -z $N ]; then
    N=$1
    if [ -z $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi
if [ -f config.txt ]; then
	source config.txt
fi
TEMPDIR=$SCRATCH
if [ ! -f $REFGENOME ]; then
    module load samtools
    samtools faidx $REFGENOME
fi
NSTART=$(perl -e "printf('%d',1 + $GVCF_INTERVAL * ($N - 1))")
NEND=$(perl -e "printf('%d',$GVCF_INTERVAL * $N)")
MAX=$(wc -l $REFGENOME.fai | awk '{print $1}')
if [ "$NSTART" -gt "$MAX" ]; then
	echo "NSTART ($NSTART) > $MAX"
	exit
fi
if [ "$NEND" -gt "$MAX" ]; then
	NEND=$MAX
fi
echo "$NSTART -> $NEND"

CPU=$SLURM_CPUS_ON_NODE
if [ ! $CPU ]; then
    CPU=2
fi
if [[ $(ls $GVCFFOLDER | grep -c -P "\.g.vcf$") -gt "0" ]]; then
   parallel -j $CPU bgzip {} ::: $GVCFFOLDER/*.g.vcf
  parallel -j $CPU tabix -f {} ::: $GVCFFOLDER/*.g.vcf.gz
fi

if [[ -z $POPYAML || ! -s $POPYAML ]]; then
	echo "Cannot find \$POPYAML variable - set in config.txt"
	exit
fi
if [ -z $SLICEVCF ]; then
	SLICEVCF=vcf_slice
fi
SLICEVCF=$SLICEVCF.all_sites
mkdir -p $SLICEVCF
for POPNAME in $(yq eval '.Populations | keys' $POPYAML | perl -p -e 's/^\s*\-\s*//')
do
	FILES=$(yq eval '.Populations.'$POPNAME'[]' $POPYAML | perl -p -e "s/(\S+)/-V $GVCFFOLDER\/\$1.g.vcf.gz/g"  )
	INTERVALS=$(cut -f1 $REFGENOME.fai  | sed -n "${NSTART},${NEND}p" | perl -p -e 's/(\S+)\n/--intervals $1 /g')
	INTERVALS=$(cut -f1 $REFGENOME.fai  | sed -n "${NSTART},${NEND}p" | perl -p -e 's/(\S+)\n/--intervals $1 /g')

	mkdir -p $SLICEVCF/$POPNAME
	STEM=$SLICEVCF/$POPNAME/$PREFIX.$N
	GENOVCFOUT=$STEM.all.vcf
	GENOALLVCFOUT=$STEM.all-sites.selected.vcf.gz
	FILTERSNP=$STEM.SNP.filter.vcf
	SELECTSNP=$STEM.SNP.selected.vcf
	echo "$STEM is stem; GENOVCFOUT=$STEM.all.vcf POPNAME=$POPNAME slice=$SLICEVCF"
	mkdir -p $TEMPDIR
	if [ ! -f $GENOVCFOUT.gz ]; then
	    if [ ! -f $GENOVCFOUT ]; then
		DB=$TEMPDIR/${GVCFFOLDER}_slice_$N
		rm -rf $DB

		echo gatk  --java-options "-Xmx$MEM -Xms$MEM" GenomicsDBImport --consolidate --merge-input-intervals \
		      --genomicsdb-workspace-path $DB $FILES $INTERVALS --tmp-dir $TEMPDIR --reader-threads $CPU

		gatk  --java-options "-Xmx$MEM -Xms$MEM" GenomicsDBImport --consolidate --merge-input-intervals \
		      --genomicsdb-workspace-path $DB $FILES $INTERVALS --tmp-dir $TEMPDIR --reader-threads $CPU

		time gatk  --java-options "-Xmx$MEM -Xms$MEM" GenotypeGVCFs --reference $REFGENOME --output $GENOVCFOUT \
		     -V gendb://$DB --tmp-dir $TEMPDIR -all-sites  -G StandardAnnotation -G AS_StandardAnnotation $INTERVALS
		ls -l $DB
		rm -rf $DB
	    fi
	    if [ -f $GENOVCFOUT ]; then
	    	bgzip $GENOVCFOUT
	    	tabix $GENOVCFOUT.gz
	    fi
	fi
	if [[ ! -f $GENOALLVCFOUT || $GENOVCFOUT.gz -nt $GENOALLVCFOUT ]]; then

	    module load vcftools
	    vcftools --gzvcf $GENOVCFOUT.gz --max-maf 0 --recode --stdout |  bgzip -c > $SCRATCH/invariant.vcf.gz
	    #	bcftools view -i "MAF[0] <= 0" -Oz -o $SCRATCH/invariant.vcf.gz $GENOVCFOUT.gz
	    vcftools --gzvcf $GENOVCFOUT.gz --mac 1 --recode --stdout | bgzip -c > $SCRATCH/variant.vcf.gz
	    #	bcftools view -i "MAC == 1" -Oz -o $SCRATCH/invariant.vcf.gz $GENOVCFOUT.gz
	    module unload vcftools
	    tabix $SCRATCH/invariant.vcf.gz
	    tabix $SCRATCH/variant.vcf.gz
	    bcftools concat --allow-overlaps $SCRATCH/variant.vcf.gz $SCRATCH/invariant.vcf.gz \
		     -O z -o $GENOALLVCFOUT
	    tabix $GENOALLVCFOUT
	fi
	if [[ ! -f $SELECTSNP.gz || $GENOVCFOUT.gz -nt $SELECTSNP.gz ]]; then

	    TYPE=SNP
	    gatk SelectVariants \
		 -R $REFGENOME \
		 --variant $GENOVCFOUT.gz \
		 -O $STEM.$TYPE.vcf \
		 --restrict-alleles-to BIALLELIC \
		--select-type-to-include $TYPE --create-output-variant-index false

	    bgzip -f $STEM.$TYPE.vcf
	    tabix $STEM.$TYPE.vcf.gz
	    # create a filtered VCF containing only variant sites
	    gatk VariantFiltration --output $FILTERSNP --tmp-dir $TEMPDIR \
		 --variant $STEM.$TYPE.vcf.gz -R $REFGENOME \
		 --cluster-window-size 10  \
		 --filter-expression "QD < 2.0" --filter-name QualByDepth \
		 --filter-expression "MQ < 40.0" --filter-name MapQual \
		 --filter-expression "SOR > 3.0" --filter-name StrandOddsRatio \
		 --filter-expression "FS > 60.0" --filter-name FisherStrandBias \
		 --missing-values-evaluate-as-failing --create-output-variant-index false
	    bgzip $FILTERSNP
	    tabix $FILTERSNP.gz
	    gatk SelectVariants -R $REFGENOME \
		 --variant $FILTERSNP.gz \
		 --output $SELECTSNP \
		 --exclude-filtered --create-output-variant-index false
	    bgzip $SELECTSNP
	    tabix $SELECTSNP.gz
	fi
done
