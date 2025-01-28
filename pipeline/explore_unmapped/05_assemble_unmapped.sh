#!/usr/bin/bash -l
#SBATCH -p batch -N 1 -n 8 --mem 48gb --out logs/assemble_unmapped.%a.log

module load spades
module load fastp
module load workspace/scratch
MEM=48
UNMAPPEDASM=unmapped_asm
UNMAPPED=unmapped
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

IFS=,
tail -n +2 $SAMPFILE | sed -n ${N}p | while read STRAIN FILEBASE
do
  UMAP=$UNMAPPED/${STRAIN}.$FASTQEXT
  UMAPSINGLE=$UNMAPPED/${STRAIN}_single.$FASTQEXT
  if [[ ! -s $UMAP && ! -s $UMAPSINGLE ]]; then
    echo "Need unmapped FASTQ file, skipping $STRAIN ($FILEBASE)"
  else
    fastp -y --in1 $UMAP --interleaved_in --average_qual 4 --correction --compression 5 --out1 $SCRATCH/${STRAIN}_filter_1.fq.gz --out2 $SCRATCH/${STRAIN}_filter_2.fq.gz --merged_out $SCRATCH/${STRAIN}_merged.fq.gz --merge --json $UNMAPPED/$STRAIN.fastp.PE.json --html $UNMAPPED/$STRAIN.fastp.PE.html 
    fastp -y --in1 $UMAPSINGLE -y --compression 5 --out1 $SCRATCH/${STRAIN}_single.fq.gz  --json $UNMAPPED/$STRAIN.fastp.SE.json --html $UNMAPPED/$STRAIN.fastp.SE.html

    spades.py --pe-1 1 $SCRATCH/${STRAIN}_filter_1.fq.gz --pe-2 1 $SCRATCH/${STRAIN}_filter_2.fq.gz \
	    --pe-s 1 $SCRATCH/${STRAIN}_single.fq.gz --pe-s 2 $SCRATCH/${STRAIN}_merged.fq.gz -o $UNMAPPEDASM/$STRAIN -t $CPU -m $MEM --careful --tmp-dir $SCRATCH
  fi
done
