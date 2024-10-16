#!/bin/bash -l
#SBATCH -N 1 -n 16 --mem 32gb --out logs/bwa.%a.log --time 8:00:00
module load bwa
module load samtools
module load picard
module load gatk/4
module load java
module load workspace/scratch

TMPOUTDIR=$SCRATCH

if [ -f config.txt ]; then
  source config.txt
fi
if [ -z $REFGENOME ]; then
  echo "NEED A REFGENOME - set in config.txt and make sure 00_index.sh is run"
  exit
fi
TEMP=$SCRATCH
if [ ! -f $REFGENOME.dict ]; then
  echo "NEED a $REFGENOME.dict - make sure 00_index.sh is run"
fi
mkdir -p $TMPOUTDIR $ALNFOLDER $UNMAPPED

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
cat $SAMPFILE | sed -n ${N}p | while read STRAIN FILEBASE
do
  PREFIX=$STRAIN
  FINALFILE=$ALNFOLDER/$STRAIN.$HTCEXT
  echo "infile is $FILEBASE"
  echo "To process $PREFIX and $FINALFILE"
  if [ ! -s $FINALFILE ]; then
    BAMSTOMERGE=()
    for BASEPATTERN in $(echo $FILEBASE | perl -p -e 's/\;/,/g');
    do
      BASE=$(basename $BASEPATTERN | perl -p -e 's/(\S+)\[12\].+/$1/g; s/_R\?(_\d+)?\.(fastq|fq)(\.gz)?$//g;')
      # END THIS PART IS PROBABLY PROJECT SPECIFIC
      echo "STRAIN is $STRAIN BASE is $BASE BASEPATTERN is $BASEPATTERN"      
      TMPBAMFILE=$TEMP/$BASE.unsrt.bam
      SRTED=$TMPOUTDIR/$BASE.srt.bam
      DDFILE=$TMPOUTDIR/$BASE.DD.bam

      FINALFILE=$ALNFOLDER/$STRAIN.$HTCEXT
      READGROUP="@RG\tID:$BASE\tSM:$STRAIN\tLB:$BASE\tPL:illumina\tCN:$RGCENTER"
      echo "$TMPBAMFILE $READGROUP"
#      echo "bwa mem -t $CPU -R $READGROUP $REFGENOME $FASTQFOLDER/$BASEPATTERN | samtools sort --threads $CPU -O bam -o $SRTED -T $TEMP -"

      if [ ! -s $DDFILE ]; then
        if [ ! -s $SRTED ]; then
          if [ -e $PAIR1 ]; then
            if [ ! -f $SRTED ]; then
              # potential switch this to bwa-mem2 for extra speed
              bwa mem -t $CPU -R $READGROUP $REFGENOME $FASTQFOLDER/$BASEPATTERN | samtools sort --threads $CPU -O bam -o $SRTED -T $TEMP -
            fi
          else
            echo "Cannot find $BASEPATTERN, skipping $STRAIN"
            exit
          fi
        fi # SRTED file exists or was created by this block
        time java -jar $PICARD MarkDuplicates -I $SRTED -O $DDFILE \
          -METRICS_FILE logs/$STRAIN.dedup.metrics -CREATE_INDEX true -VALIDATION_STRINGENCY SILENT
        if [ -f $DDFILE ]; then
          rm -f $SRTED
        fi
      fi # DDFILE is created after this or already exists
      BAMSTOMERGE+=( $DDFILE )
    done
    samtools merge --write-index -O $HTCFORMAT --threads $CPU --reference $REFGENOME $FINALFILE "${BAMSTOMERGE[@]}"
    if [ -f $FINALFILE ]; then
      rm -f "${BAMSTOMERGE[@]}"
      rm -f $(echo "${BAMSTOMERGE[@]}" | sed 's/bam$/bai/')
    fi

    if [ -f $FINALFILE ]; then
      rm -f $DDFILE
      rm -f $(echo $DDFILE | sed 's/bam$/bai/')
    fi
  fi #FINALFILE created or already exists
  FQ=$(basename $FASTQEXT .gz)
  UMAP=$UNMAPPED/${STRAIN}.$FQ
  UMAPF=$UNMAPPED/${STRAIN}.F.$FQ
  UMAPR=$UNMAPPED/${STRAIN}.R.$FQ

  UMAPSINGLE=$UNMAPPED/${STRAIN}_single.$FQ
  #echo "$UMAP $UMAPSINGLE $FQ"

  if [ ! -f $UMAP.gz ]; then
    module load BBMap
    samtools fastq -f 4 --threads $CPU -N -1 $UMAPF -2 $UMAPR $FINALFILE > $UMAPSINGLE
    if [ -s $UMAP_SINGLE ]; then
    	pigz $UMAPSINGLE
    fi
    repair.sh in=$UMAPF in2=$UMAPR out=$UMAP.gz
    unlink $UMAPF
    unlink $UMAPR
  fi
done
