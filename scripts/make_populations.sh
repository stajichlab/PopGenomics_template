#!/usr/bin/bash -l

echo -e "Populations:\n\n  All:" > population_sets.yaml
IFS=,
while read POPNAME FILES
do
	echo "    - $POPNAME"
done < samples.csv >> population_sets.yaml
