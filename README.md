# Epigenomic Control of Whole-Genome Duplicated Gene Expression
Analysis of FAIRE-seq, MNase-seq, etc.

## Preliminary RNA-seq analysis
RNA-seq libraries were acquired for mature cotton leaves under Long day (LD: 7am-9pm) and short day (SD:7am-5pm) conditions and sampled at both dawn and dusk. From four sampling stages LD7, LD9, SD7 and SD5, which one should we select for nucleosone occupancy analysis? Is these a better choice rather than picking anyone?

First, analyze total read counts extracted by [htseqCount.sh](https://github.com/huguanjing/Epigenomics/blob/master/htseqCount.sh), using R script [deseq2.R](https://github.com/huguanjing/Epigenomics/blob/master/htseqCount.sh).

    # access to mapping results, get bam files started with LD7, LD9, SD7 and SD5, save into "matureLeaf.filenames.txt"
    ls ~/jfw-lab/Projects/Duplicated_Pathways/DupNetMap/
    # get read count
    bash htseqCount.sh


