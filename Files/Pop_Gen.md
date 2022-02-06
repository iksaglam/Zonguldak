# A tutorial for some basic population genetic analyses using RADseq data.

In this tutorial we will be using 60 samples of *Phonochorion artvinensis* to perform basic population genetic analysis from low-depth RAD sequencing data using several programs including, [samtools](http://www.htslib.org/), [bwa](http://bio-bwa.sourceforge.net/), [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD), [PCangsd](http://www.popgen.dk/software/index.php/PCAngsd) and [NGSadmix](http://www.popgen.dk/software/index.php/NgsAdmix). Because we do not have a reference genome for this species or from any other closely related species we will use a denovo set of RAD contigs discovered using the OVIT population of *Ph. uvarovi* as both an outgroup and a reference genome for alignment. More information on this study can be found [here](https://www.sciencedirect.com/science/article/pii/S1055790319302428), and the whole data set is publically available [here](https://github.com/iksaglam/Phonochorion_data). 

The data for this practical can be found here `/egitim/iksaglam/data/artvinensis` and example shell scripts for all analysis and R scripts for manipulating and plotting results can be found [here](https://github.com/iksaglam/Zonguldak/tree/main/Scripts). Please be aware that shell and R scripts provided here are for illustrative purposes only and if you wish to use them on your own data/project they should be modified accordingly. 

Our goal here is to show how to go from raw fastq files to alignment files (BAMs) and from there to basic summary statistics and analysis in population genetics.
Specifically we will be covering the following topics:

- Sequence alignment and basic filtering
- Genotyping and variant calling
- Population structure (PCA)
- Admixture
- Site frequency spectrum
- Nucleotive diversity and neutrality
- Population genetic differentiation


Finally, all analysis, frameworks and programs covered here are not restricted to RADseq data and can be thought of as an application of probabilistics approaches for processing NGS data in general.


## stitle
### sstitle

normal text


```
code block
```

*italic*



[github](www.github.com)



**bold**




divider below
---



- list item
- list item 2





table below

a | b | c
-- | -- | --
xxx | yyy | zzz


Inline code block `im a code block`


And to create a Markdown file, just name your file with extension `.md` or `.MD`
