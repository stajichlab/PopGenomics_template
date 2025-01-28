#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 64 --mem 256gb --out logs/mmseqs_classify_unmappedreads.%a.log -a 1-90

module load workspace/scratch
module load mmseqs2

if [ -f config.txt ]; then
  source config.txt
fi

CPU=2
if [ $SLURM_CPUS_ON_NODE ]; then
  CPU=$SLURM_CPUS_ON_NODE
fi
N=${SLURM_ARRAY_TASK_ID}
if [ -z $N ]; then
  N=$1
fi
if [ -z $N ]; then
  echo "cannot run without a number provided either cmdline or --array in sbatch"
  exit
fi

MAX=$(wc -l $SAMPFILE | awk '{print $1}')
if [ $N -gt $MAX ]; then
  echo "$N is too big, only $MAX lines in $SAMPFILE"
  exit
fi
INDIR=unmapped
DB2=/srv/projects/db/ncbi/mmseqs/uniref50
DB2NAME=$(basename $DB2)
IFS=,
OUTSEARCH=results/mmseqs2_unmapped
mkdir -p $OUTSEARCH
tail -n +2 $SAMPFILE | sed -n ${N}p | while read STRAIN FILEBASE
do
    mkdir -p $OUTSEARCH/$STRAIN
    BASEPATTERN=$(echo $FILEBASE | perl -p -e 's/\;/ /g; ')
    if [ ! -s $OUTSEARCH/$STRAIN/mmseq_${DB2NAME}_report ]; then
	IDX=$SCRATCH/${STRAIN}_reads.idx
	mmseqs createdb $INDIR/$STRAIN.fastq.gz $INDIR/${STRAIN}_single.fastq.gz $IDX --dbtype 2
	mmseqs taxonomy $IDX $DB2 $SCRATCH/${DB2NAME}_taxo $SCRATCH -s 2 --threads $CPU
	mmseqs taxonomyreport $DB2 $SCRATCH/${DB2NAME}_taxo $OUTSEARCH/$STRAIN/mmseq_$DB2NAME.krona_native.html --report-mode 1
	mmseqs taxonomyreport $DB2 $SCRATCH/${DB2NAME}_taxo $OUTSEARCH/$STRAIN/mmseq_${DB2NAME}_report

#	mmseqs easy-taxonomy $INDIR/$STRAIN.fastq.gz $INDIR/${STRAIN}_single.fastq.gz $DB2 $OUTSEARCH/$STRAIN/mmseq_$DB2NAME $SCRATCH --threads $CPU --lca-ranks kingdom,phylum,family  --tax-lineage 1
#	pigz $OUTSEARCH/$STRAIN/*.tsv $OUTSEARCH/$STRAIN/mmseq_uniref50_tophit_aln $OUTSEARCH/$STRAIN/mmseq_uniref50_tophit_report
    fi
done
