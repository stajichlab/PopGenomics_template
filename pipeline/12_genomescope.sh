#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 24 --mem 64gb --out logs/genomescope.%a.log -a 1-28

module load workspace/scratch
module load samtools
module load jellyfish
module load R

GENOMESCOPE=genomescope
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
JELLYFISHSIZE=1000000000
IFS=,
KMER=21
READLEN=150
tail -n +2 $SAMPFILE | sed -n ${N}p | while read STRAIN FILEBASE
do
    # for this project has only a single file base  but this needs fixing otherse
    BASEPATTERN=$(echo $FILEBASE | perl -p -e 's/\;/ /g; ')
    FINALFILE=$ALNFOLDER/$STRAIN.$HTCEXT
    jellyfish count -C -m $KMER -s $JELLYFISHSIZE -t $CPU -o $SCRATCH/$STRAIN.jf <(samtools fastq -F 4 --threads 4 $FINALFILE)
    #jellyfish count -C -m $KMER -s $JELLYFISHSIZE -t $CPU -o $SCRATCH/$STRAIN.jf <(pigz -dc $FASTQFOLDER/$BASEPATTERN)
    jellyfish histo -t $CPU $SCRATCH/$STRAIN.jf > $GENOMESCOPE/$STRAIN.histo
    Rscript scripts/genomescope.R $GENOMESCOPE/$STRAIN.histo $KMER $READLEN $GENOMESCOPE/$STRAIN/
done
