# A tutorial for some basic population genetic analyses using RADseq data.

In this tutorial we will be using 60 samples of *Phonochorion artvinensis* to perform basic population genetic analysis from low-depth RAD sequencing data using several programs including, samtools, bwa, ANGSD, PCangsd and NGSadmix. Because we do not have a reference genome for this species or from any other closely related species we will use a denovo set of RAD contigs discovered using the OVIT population of *Ph. uvarovi* as both an outgroup and a reference genome for alignment. We will also use several R scripts for manipulating and plotting results but please be aware that scripts provided here are for illustrative purposes only and if you wish to use them on different data sets they should be modified. Finally, all analysis, frameworks and programs covered here are not restricted to RADseq data and can be thought of as an application of probabilistics approaches for processing NGS data in general.


Our goal here is to show how to go from raw fastq files, to alignment files (BAMs) and from ther to basic summary statistics and analysis in population genetics. Nanemly we will cover the follwong topics:

- Sequence alignment and basic filtering
- Genotyping and variant calling
- Population structure (PCA)
- Admixture
- Site frequency spectrum
- Nucleotive diversity and neutrality
- Population genetic differentiation



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
