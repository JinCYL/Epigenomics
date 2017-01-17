# Epigenomic Control of Whole-Genome Duplicated Gene Expression
Analysis of FAIRE-seq, MNase-seq, etc.

## Preliminary RNA-seq analysis
RNA-seq libraries were acquired for mature cotton leaves under Long day (LD: 7am-9pm) and short day (SD:7am-5pm) conditions and sampled at both dawn and dusk. From four sampling stages LD7, LD9, SD7 and SD5, which one should we select for nucleosone occupancy analysis? Is these a better choice rather than blindly picking one?

Total read counts were extracted by [htseqCount.sh](https://github.com/huguanjing/Epigenomics/blob/master/htseqCount.sh), and differetial expression between conditions X genomes were analyzed using R script [deseq2.R](https://github.com/huguanjing/Epigenomics/blob/master/htseqCount.sh).

    # access to mapping results, get bam files started with LD7, LD9, SD7 and SD5, save into "matureLeaf.filenames.txt"
    ls ~/jfw-lab/Projects/Duplicated_Pathways/DupNetMap/
    # get read count
    bash htseqCount.sh

Basically, the expression differences between sunrise and sunset dominate other varations, including LD vs SD and among genomes. Only 6.3% of DE genes were found between LD7 and SD7, and 11.9% were found between LD9 and SD5. I ended up using SD5 for MNase-seq analysis, (shamefully) because it is just easier to collect at that time .... 

