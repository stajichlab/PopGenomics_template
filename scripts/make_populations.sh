#!/usr/bin/bash -l

echo -e "Populations:\n\n  All:" > populations.yaml
IFS=,
while read POPNAME FILES
do
	echo "    - $POPNAME"
done < samples.csv >> populations.yaml
