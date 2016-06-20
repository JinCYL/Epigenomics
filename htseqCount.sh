#!/bin/bash
module load python

input="matureLeaf.filenames.txt"
while IFS= read -r j
do
  echo $j
  # echo ~/jfw-lab/Projects/Duplicated_Pathways/DupNetMap/.$j
  htseq-count -f bam --stranded=no -r name ~/jfw-lab/Projects/Duplicated_Pathways/DupNetMap/$j ~/jfw-lab/GenomicResources/pseudogenomes/D5.primaryOnly.gtf > $j.txt
done < "$input"
