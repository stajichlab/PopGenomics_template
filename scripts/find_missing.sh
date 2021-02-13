#!/usr/bin/bash
#SBATCH -p short
joinByString() {
  local separator="$1"
  shift
  local first="$1"
  shift
  printf "%s" "$first" "${@/#/$separator}"
}
export -f joinByString
source config.txt

ALIGN=()
GVCF=()
for n in $(cut -d, -f1 $SAMPFILE | tail -n +2); do 
	if [ ! -f $ALNFOLDER/$n.$HTCEXT ]; then 
#		echo "($N) $n"; 
		ALIGN+=( $N )
	fi
	if [ ! -f $GVCFFOLDER/$n.g.vcf.gz ]; then
		GVCF+=( $N )
	fi
	N=$(expr $N + 1)
done

ALNLIST=$(joinByString , "${ALIGN[@]}")
echo "sbatch -a $ALNLIST pipeline/01_align.sh"
VCFLIST=$(joinByString , "${GVCF[@]}")
echo "sbatch -a $VCFLIST pipeline/02_call_gvcf.sh"
