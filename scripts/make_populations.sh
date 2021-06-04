#!/usr/bin/bash -l

echo -e "Populations:\n\n  All:" > population_sets.yaml
IFS=,
tail -n +2 samples.csv | while read POPNAME FILES
do
	echo "    - $POPNAME"
done >> population_sets.yaml
