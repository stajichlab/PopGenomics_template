#!/bin/bash -l
#SBATCH -p short -N 1 -n 64 --mem 256gb --out logs/unmapped_asm_mmseqs_classify.%a.log

module load mmseqs2
module load KronaTools
module load workspace/scratch

UNMAPPEDASM=unmapped_asm
OUTSEARCH=unmapped_asm_taxonomy
mkdir -p $OUTSEARCH
DB=/srv/projects/db/ncbi/mmseqs/swissprot
DB2=/srv/projects/db/ncbi/mmseqs/uniref50

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
    if [ -z $N ]; then
	echo "cannot run without a number provided either cmdline or --array in sbatch"
	exit
    fi
fi

MAX=$(wc -l $SAMPFILE | awk '{print $1}')
if [ $N -gt $MAX ]; then
  echo "$N is too big, only $MAX lines in $SAMPFILE"
  exit
fi

IFS=,
tail -n +2 $SAMPFILE | sed -n ${N}p | while read STRAIN FILEBASE 
do
    
    CONTIG=$UNMAPPEDASM/$STRAIN/scaffolds.fasta
    if [ ! -f $CONTIG ]; then
	echo "Need to finish running assembly for $STRAIN"
	exit
    fi
    mkdir -p $OUTSEARCH/$STRAIN
    if [ ! -f  $OUTSEARCH/$STRAIN/mmseq_uniref50_report ]; then
	mmseqs easy-taxonomy $CONTIG $DB2 $OUTSEARCH/$STRAIN/mmseq_uniref50 $SCRATCH --threads $CPU --lca-ranks kingdom,phylum,family  --tax-lineage 1
    fi
    if [ ! -f $OUTSEARCH/$STRAIN/mmseq_sprot_report ]; then
	mmseqs easy-taxonomy $CONTIG $DB $OUTSEARCH/$STRAIN/mmseq_sprot $SCRATCH --threads $CPU --lca-ranks kingdom,phylum,family --tax-lineage 1
    fi
     mmseqs taxonomyreport $DB2 $OUTSEARCH/$STRAIN/mmseq_uniref50_lca.tsv $OUTSEARCH/$STRAIN/mmseq_uniref50.krona_native.html --report-mode 1
     mmseqs taxonomyreport $DB $OUTSEARCH/$STRAIN/mmseq_sprot_lca.tsv $OUTSEARCH/$STRAIN/mmseq_sprot.krona_native.html --report-mode 1

    ktImportTaxonomy -o $OUTSEARCH/$STRAIN/mmseq_uniref50.krona.html $OUTSEARCH/$STRAIN/mmseq_uniref50_report
    ktImportTaxonomy -o $OUTSEARCH/$STRAIN/mmseq_sprot.krona.html $OUTSEARCH/$STRAIN/mmseq_sprot_report
done
