#!/bin/bash
#$ -cwd
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

while read -r line; do
	g=$(echo $line | cut -d "," -f 1 )
	t=$(echo $line | cut -d "," -f 2 )
	$SCRIPT_DIR/RNAhybrid -c -s 3utr_human ${g} ${t}
done < $1