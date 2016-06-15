# Epigenomic Control of Whole-Genome Duplicated Gene Expression
Analysis of FAIRE-seq, MNase-seq, etc.

## Preliminary RNA-seq analysis
RNA-seq libraries were acquired for mature cotton leaves under Long day (LD: 7am-9pm) and short day (SD:7am-5pm) conditions and sampled at both dawn and dusk. From four sampling stages LD7, LD9, SD7 and SD5, which one should we select for nucleosone occupancy analysis? Is these a better choice rather than picking anyone?

    # access to mapping results, get bam files started with LD7, LD9, SD7 and SD5, save into "matureLeaf.filenames.txt"
    cd ~/jfw-lab/Projects/Duplicated_Pathways/DupNetMap/
    # get read count
    bash htseqCount.sh


<htseqCount.sh>

    #!/bin/bash
    module load python
    input="matureLeaf.filenames.txt"
    while IFS= read -r j
    do
    echo $j
    # echo ~/jfw-lab/Projects/Duplicated_Pathways/DupNetMap/.$j
    htseq-count -f bam --stranded=no -r name ~/jfw-lab/Projects/Duplicated_Pathways/DupNetMap/$j ~/jfw-lab/GenomicResources/pseudogenomes/D5.primaryOnly.gtf > $j.txt
    done < "$input"
