# A tutorial for some basic population genetic analyses using RADseq data.

In this tutorial we will be using 60 samples from 5 populations of *Phonochorion artvinensis* (Figure 1) to perform basic population genetic analysis from low-depth RAD sequencing data using several programs including, [samtools](http://www.htslib.org/), [bwa](http://bio-bwa.sourceforge.net/), [picard-tools](https://broadinstitute.github.io/picard/), [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD), [PCangsd](http://www.popgen.dk/software/index.php/PCAngsd) and [NGSadmix](http://www.popgen.dk/software/index.php/NgsAdmix). Because we do not have a reference genome for this species or from any other closely related species we will use a denovo set of RAD contigs discovered using the OVIT population of *Ph. uvarovi* as both an outgroup and a reference genome for alignment. More information on this study can be found [here](https://www.sciencedirect.com/science/article/pii/S1055790319302428), and the whole data set is publically available [here](https://github.com/iksaglam/Phonochorion_data).
<br/><br/>

![Figure 1](https://github.com/iksaglam/Zonguldak/blob/main/Files/1-s2.0-S1055790319302428-gr1_lrg.jpg)
**Figure 1.** Distribution of *Phonochorion artvinensis* and *Ph. uvarovi* populations within the Eastern Black Sea Region of Turkey. *Ph. artvinensis*: 1=Erikli Köyü; 4=Kama Köyü; 6=Çamlık Köyü; 10=Çaymakçur Yaylası; 12=Napopen Yaylası. *Ph. uvarovi*: 2=Boğalı Köy, 3=Anzer Yaylası, 5=Ovit Geçidi, 7=Çağrankaya Yaylası, 8=Alt Çat, 9=Kavron Yaylası, 11=Sırt Yaylası.
<br/><br/>

The data for this practical can be found here `/egitim/iksaglam/data/artvinensis` and the denovo partial reference genome which we will be using for mapping and variant calling can be found here `/egitim/iksaglam/ref`. In addition, example shell scripts for all analysis and R scripts for manipulating and plotting results can be found [here](https://github.com/iksaglam/Zonguldak/tree/main/Scripts). Please be aware that shell and R scripts provided here are for illustrative purposes only and if you wish to use them on your own data/project they should be modified accordingly. 

Our goal here is to show how to go from raw fastq files to alignment files (BAMs) and from there to basic summary statistics and analysis in population genetics.
Specifically we will be covering the following topics:

- Sequence alignment, cleanup and indexing
- Genotyping and variant calling
- Population structure (PCA)
- Admixture
- Site frequency spectrum
- Nucleotive diversity and neutrality
- Population genetic differentiation


Finally, all analysis, frameworks and programs covered here are not restricted to RADseq data and can be thought of as an application of probabilistics approaches for processing NGS data in general.


## Sequence alignment, cleanup and indexing
As a first step let us copy the data into our own directory and list to make sure everything was copied correctly.

```
cd ~/my_directory/
cp /egitim/iksaglam/data/artvinensis/*.fastq.gz ./
ls *.fastq.gz
```

Let us also make a file listing all individuals and another one listing all populations for future use
```
ls *.fastq.gz | cut -d'_' -f1-2 > artv.list
cut -c3-5 artv.list | uniq > pop.list
```
### Aligning reads with BWA

For aligning reads we will use [BWA](http://bio-bwa.sourceforge.net/) (Burrows-Wheeler Aligner) which is based on the [Burrows-Wheeler Transformed](https://academic.oup.com/bioinformatics/article/25/14/1754/225615?login=true) index of the reference genome. It performs gapped alignment for single-end reads, supports paired-end mapping, generates mapping quality and can give multiple hits if called for.


#### Creating BWA-MEM index
Like all other alignment tools for genome wide short reads the first step is to index the reference genome. BWA indexes the genome with an FM Index based on the Burrows-Wheeler Transform enabling memory efficient mapping. 

The basic command for indexing the genome using BWA is:
```
ref=/egitim/iksaglam/ref/uvar_ref_contigs_300.fasta
bwa index -a bwtsw $ref
```

This will result in several indexes including the BWT index `/egitim/iksaglam/ref/uvar_ref_contigs_300.fasta.bwt` which the alignments will be based on. The creation of the index only needs to be performed once and does not have to be recreated for every alignment job.

#### Aligning reads with BWA-MEM

Now that we have our indexes created, we can start aligning our reads (i.e. individual fastq files) to the reference genome. Let us first create a new directory where we will store our alignment files and cd (change directory) into it.
```
mkdir alignments
cd alignments/
```
Next we will use the BWA-MEM algorithm to align one of our paired-end reads to the reference genome and output results as a [SAM](https://samtools.github.io/hts-specs/SAMv1.pdf) (Sequence Alignment/Map Format) file.

```
reads=/egitim/iksaglam/data/artvinensis
ref=/egitim/iksaglam/ref/uvar_ref_contigs_300.fasta
out=/egitim/iksaglam/alignments
bwa mem $ref $reads/A_CAMD04_R1.fastq.gz $reads/A_CAMD04_R2.fastq.gz > $out/A_CAMD04.sam
```
As before we can take a look at our SAM files using less.

```
less -S A_CAMD04_R1.sam
```

#### SAM/BAM conversion, sorting and PCR clone (duplicate) markup/removal:

To save space we ideally want to transform our SAM file into a binary BAM format and sort by coordinates. Next when working with reduced representation libraries like RADseq data it is advantegous to keep only properly paired individuals. Finally we want to remove/mark PCR duplicates using [picard-tools](https://broadinstitute.github.io/picard/), so they do not bias variant calling and genotyping and index our final BAM file.

```
samtools view -bS A_CAMD04.sam > A_CAMD04.bam
samtools sort A_CAMD04.bam A_CAMD04_sorted.bam
samtools view -b -f 0x2 A_CAMD04_sorted.bam > A_CAMD04_sorted_proper.bam
java -jar /kuacc/apps/picard/2.22.1/picard.jar MarkDuplicates INPUT=A_CAMD04_sorted_proper.bam OUTPUT=A_CAMD04_sorted_proper_rmdup.bam METRICS_FILE=A_CAMD04_metrics.txt VALIDATION_STRINGENCY=LENIENT  REMOVE_DUPLICATES=True
samtools index A_CAMD04_sorted_proper_rmdup.bam
```
For the purposes of this tutorial we are directly removing duplicates `REMOVE_DUPLICATES=True` from our resulting bam files (to primarily save space) but usually we only want to mark our duplicates and not remove them `REMOVE_DUPLICATES=False`.

To view our resuling bam files we can use the following command
```
samtools view A_CAMD04_sorted_proper_rmdup.bam | less -S
```

We can also get alignment statistics using samtools flagstat
```
samtools flagstat A_CAMD04_sorted_proper_rmdup.bam
```

Instead of aligning each paired read one by one normally we would want to do this in bulk. You can find an example shell script that can do this [here] (https://github.com/iksaglam/Zonguldak/blob/main/Scripts/align_pe_reads.sh). We can use this script to execute the above pipeline and get alignment files for all individuals `artv.list` using the following command.

```
sh align_pe_reads.sh artv.list
```

Now that we have 60 BAM files at low/medium depth we can create a bamlist for the data set as a whole and for each population to use in downstream analysis in a new directory called `analyses`.

```
mkdir analyses
cd analyses
ls /egitim/iksaglam/alignments/*.sorted_proper_rmdup.bam > artv.bamlist
for i in `cat pop.list`; do grep $i artv.bamlist > ${i}.bamlist; done
```

## Genotyping and variant calling

To peform basic population genetic analysis in [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD), we first need to assign genotype probabilities at each site for each individual and determine allele frequencies for each site. The specific option in [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD) for calculating genotype probabilities and determining polymorphic sites is `-doGeno` together with `-doMAF`. When working with low/medium depth data it is always advantegous to work with genotype porbabilities but we can also call genotypes using the `-doPost` option with an appopriate probability cutoff `-postCutoff`. Additionally (or alternatively) we can also create output files in beagle `-doGlf`, vcf `-dovcf` and PLINK `-doPlink` format to have the option to use other NGS analysis toolkits like [PLINK](https://www.cog-genomics.org/plink/2.0/), [GATK](https://gatk.broadinstitute.org/hc/en-us) or [STACKS](https://catchenlab.life.illinois.edu/stacks/) among others.

### Basic filtering
When determining genotype probabilities and polymorphic sites we ideally want to work with only high quality sites. We can do this by filtering out low quality reads, reads with low mapping quality and removing sites where half of the individuals have no data. This can be achieved by setting appropriate parameters to the below options.

Parameter | Meaning |
-- | -- |
-minMapQ 10 | minimum mapping quality of 10 |
-minQ 20 | minimum base quality of 20 |
-minInd 15 | use only sites with data from at least 15 individuals |

Note that there are many more filtering options in ANGSD. For a list of all avaliable filters see links below:

- quality and depth, see [here](http://www.popgen.dk/angsd/index.php/Filters)
- SNP quality, see [here](http://popgen.dk/angsd/index.php/SnpFilters)
- sites, see [here](http://popgen.dk/angsd/index.php/Sites)

### Calculating genotype probabilities and determining polymorphic sites

As input for [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD) we give the list of BAM files with option `-bam` and then specify the references sequence with `-ref` and the prefix for output files with `-out`. For assigning genotype probabilities we set `-doGeno 5` and define a model to calculate posterior probability of genotypes `-doPost 1`. We estimate allele frequencies by setting `-doMAF 1`, specify major and minor alleles by setting `-doMajorMinor 1` and specify the genotype likelihood model to use `-GL 1`. For most cases, we want to restrict the analysis to polymorphic sites (SNPs), as non-variable sites across all samples will not carry information regarding population structure or differentiation. We can perform SNP calling by setting `-minMaf 0.05` and `-SNP_pval 1e-12`. Finally if we want to call genotypes we can set `-postCutoff 0.85`. Explaination of all options are summarized below.

Option | Meaning |
-- | -- |
-doGeno 5 | 1: write major and minor; 4: write the called genotype as: AA,AC etc  |
-doPost 1 | Use HWE frequencies as prior |
-doMAF 1 | Frequency (fixed major and minor) |
-doMajorMinor 1 | Infer major and minor alleles from GL |
-GL 1 | use SAMtools GL model |
-minMaf 0.05 | Remove sites with MAF below 0.05 |
-SNP_pval 1e-12 | Remove sites with a pvalue larger than 1e-12  |
-postCutoff 0.85 | Call genotypes with a posterior probability of over 85%  |

Now let us create a directory for carrying out all our anal
Recalling also our choice for data filtering and that we would like to output files in various formats, our final command line would look something like this:

```
mkdir results_geno
ref=/egitim/iksaglam/ref/uvar_ref_contigs_300.fasta
nInd=$(wc -l artv.bamlist | awk '{print $1}')
mInd=$((${nInd}/2))

angsd -bam artv.bamlist -out results_geno/artv -GL 1 -doMajorMinor 1 -doMaf 1 -doGlf 2 -doGeno 5 -dovcf 1 -doPlink 2 -doPost 1 -postCutoff 0.85 -minMapQ 10 -minQ 20 -minInd $mInd -SNP_pval 1e-12 -minMaf 0.05 -nThreads 2
```
















our list of indiviuals 

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
