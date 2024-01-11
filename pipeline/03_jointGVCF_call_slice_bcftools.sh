#!/usr/bin/bash -l
#SBATCH --mem 24G -N 1 -n 1 -c 4 -J slice.GVCFGeno --out logs/GVCFGenoGATK4_bcf.slice_%a.%A.log -a 1-9 --time 5-0:0:0
hostname
MEM=24g
BCFTOOLSCPU=4
module load java
module load picard
module load gatk/4
module load bcftools
module load parallel
module load yq
module load workspace/scratch

source config.txt


GVCF_INTERVAL=1
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

if [ ! -f $REFGENOME ]; then
    module load samtools/1.11
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
mkdir -p $SLICEVCF
for POPNAME in $(yq eval '.Populations | keys' $POPYAML | perl -p -e 's/^\s*\-\s*//')
do
	FILES=$(yq eval '.Populations.'$POPNAME'[]' $POPYAML | perl -p -e "s/(\S+)/-V $GVCFFOLDER\/\$1.g.vcf.gz/g"  )
	INTERVALS=$(cut -f1 $REFGENOME.fai  | sed -n "${NSTART},${NEND}p" | perl -p -e 's/(\S+)\n/--intervals $1 /g')

	mkdir -p $SLICEVCF/$POPNAME
	STEM=$SLICEVCF/$POPNAME/$PREFIX.$N
	GENOVCFOUT=$STEM.all.vcf
	FILTERSNP=$STEM.SNP.filter.vcf
	FILTERINDEL=$STEM.INDEL.filter.vcf
	SELECTSNP=$STEM.SNP.selected.vcf
	SELECTINDEL=$STEM.INDEL.selected.vcf
	echo "$STEM is stem; GENOVCFOUT=$STEM.all.vcf POPNAME=$POPNAME slice=$SLICEVCF"
	if [ ! -f $GENOVCFOUT.gz ]; then
	    if [ ! -f $GENOVCFOUT ]; then
		DB=$SCRATCH/${GVCFFOLDER}_slice_$N
		rm -rf $DB
		gatk  --java-options "-Xmx$MEM -Xms$MEM" GenomicsDBImport --consolidate --merge-input-intervals --genomicsdb-workspace-path $DB $FILES $INTERVALS --tmp-dir $SCRATCH --reader-threads $CPU
		#--reader-threads $CPU
		#gatk  --java-options "-Xmx$MEM -Xms$MEM" GenomicsDBImport --genomicsdb-workspace-path $DB $FILES $INTERVALS  --reader-threads $CPU
		time gatk GenotypeGVCFs --reference $REFGENOME --output $GENOVCFOUT -V gendb://$DB --tmp-dir $SCRATCH
		ls -l $SCRATCH
		rm -rf $DB
	    fi
	    if [ -f $GENOVCFOUT ]; then
	    	bgzip $GENOVCFOUT
	    	tabix $GENOVCFOUT.gz
	    fi
	fi
	TYPE=SNP
	echo "VCF = $STEM.$TYPE.vcf.gz"
	if [[ ! -f $STEM.$TYPE.vcf.gz ]]; then
	    gatk SelectVariants \
		-R $REFGENOME \
		--variant $GENOVCFOUT.gz \
		-O $STEM.$TYPE.vcf \
		--restrict-alleles-to BIALLELIC \
		--select-type-to-include $TYPE --create-output-variant-index false

	    bgzip $STEM.$TYPE.vcf
	    tabix $STEM.$TYPE.vcf.gz
	fi

	if [[ ! -f $SELECTSNP.gz || $STEM.$TYPE.vcf.gz -nt $SELECTSNP.gz ]]; then
		START=$(date +%s)
		bcftools view -e "QD < 2.0" $STEM.$TYPE.vcf.gz -Ou --threads  4 | 
			bcftools view -e "MQ < 40.0" -Ou --threads 4 |
			bcftools view -e "QUAL < 100" -Ou |
			bcftools view -e "SOR > 4.0" -Ou | bcftools view -e "MQRankSum < -12.5" -Ou | bcftools view -e "ReadPosRankSum < -8.0" -Ou |
			bcftools view -e "FS > 60.00" -Oz -o $SELECTSNP.gz
	    END=$(date +%s)
	    echo "Elapsed Time bcftools SELECTSNP $SELECTSNP: $(($end-$start)) seconds"
	    tabix $FILTERSNP.gz
	fi

	TYPE=INDEL
	if [ ! -f $STEM.$TYPE.vcf.gz ]; then
	    gatk SelectVariants \
	        -R $REFGENOME \
	        --variant $GENOVCFOUT.gz \
	        -O $STEM.$TYPE.vcf  --select-type-to-include MIXED --select-type-to-include MNP \
	        --select-type-to-include $TYPE --create-output-variant-index false
	    bgzip $STEM.$TYPE.vcf
	    tabix $STEM.$TYPE.vcf.gz
	fi

	if [[ ! -f $SELECTINDEL.gz || $STEM.$TYPE.vcf.gz -nt $SELECTINDEL.gz ]]; then
		bcftools view -e "QD < 2.0" -Ou $STEM.$TYPE.vcf.gz | bcftools view -e "SOR > 10.0" -Ou |
			bcftools view -e "FS > 200.0" -Ou | bcftools view -e "InbreedingCoeff < -0.8" -Oz -o $SELECTINDEL.gz
	    tabix $SELECTINDEL.gz
	fi
done
