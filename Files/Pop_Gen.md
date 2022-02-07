# A tutorial for some basic population genetic analyses using RADseq data.

In this tutorial we will be using 60 samples from 5 populations of *Phonochorion artvinensis* (Figure 1) to perform basic population genetic analysis from low-depth RAD sequencing data using several programs including, [samtools](http://www.htslib.org/), [bwa](http://bio-bwa.sourceforge.net/), [picard-tools](https://broadinstitute.github.io/picard/), [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD), [PCangsd](http://www.popgen.dk/software/index.php/PCAngsd) and [NGSadmix](http://www.popgen.dk/software/index.php/NgsAdmix). Because we do not have a reference genome for this species or from any other closely related species we will use a denovo set of RAD contigs discovered using the OVIT population of *Ph. uvarovi* as both an outgroup and a reference genome for alignment. More information on this study can be found [here](https://www.sciencedirect.com/science/article/pii/S1055790319302428), and the whole data set is publically available [here](https://github.com/iksaglam/Phonochorion_data).
<br/><br/>

![Figure 1](https://github.com/iksaglam/Zonguldak/blob/main/Files/1-s2.0-S1055790319302428-gr1_lrg.jpg)
**Figure 1.** Distribution of *Phonochorion artvinensis* and *Ph. uvarovi* populations within the Eastern Black Sea Region of Turkey. *Ph. artvinensis*: 1=Erikli Köyü; 4=Kama Köyü; 6=Çamlık Köyü; 10=Çaymakçur Yaylası; 12=Napopen Yaylası. *Ph. uvarovi*: 2=Boğalı Köy, 3=Anzer Yaylası, 5=Ovit Geçidi, 7=Çağrankaya Yaylası, 8=Alt Çat, 9=Kavron Yaylası, 11=Sırt Yaylası.
<br/><br/>

The data for this practical can be found here `/egitim/iksaglam/data/artvinensis` and the denovo partial reference genome which we will be using for mapping and variant calling can be found here `/egitim/iksaglam/ref`. In addition, example shell scripts for all analysis and R scripts for manipulating and plotting results can be found [here](https://github.com/iksaglam/Zonguldak/tree/main/Scripts). Please be aware that shell and R scripts provided here are for illustrative purposes only and if you wish to use them on your own data/project they should be modified accordingly. 

Our goal here is to show how to go from raw fastq files to alignment files (BAMs) and from there to basic summary statistics and analysis in population genetics.

#### Topics to be covered:

- [Sequence alignment, cleanup and indexing](https://github.com/iksaglam/Zonguldak/blob/main/Files/Pop_Gen.md#sequence-alignment-cleanup-and-indexing)
- [Genotyping and variant calling](https://github.com/iksaglam/Zonguldak/blob/main/Files/Pop_Gen.md#genotyping-and-variant-calling)
- [Population structure (PCA)](https://github.com/iksaglam/Zonguldak/blob/main/Files/Pop_Gen.md#population-structure-pca)
- [Admixture](https://github.com/iksaglam/Zonguldak/blob/main/Files/Pop_Gen.md#admixture)
- [The Site frequency spectrum](https://github.com/iksaglam/Zonguldak/blob/main/Files/Pop_Gen.md#the-site-frequency-spectrum-sfs)
- [Nucleotive diversity and neutrality](https://github.com/iksaglam/Zonguldak/blob/main/Files/Pop_Gen.md#nucleotive-diversity-and-neutrality)
- [Population genetic differentiation](https://github.com/iksaglam/Zonguldak/blob/main/Files/Pop_Gen.md#population-genetic-differentiation)


Finally, all analysis, frameworks and programs covered here are not restricted to RADseq data and can be thought of as an application of probabilistics approaches for processing NGS data in general.


## Sequence alignment, cleanup and indexing
As a first step let us copy the data into our own directory and list to make sure everything was copied correctly.

```Bash
cd ~/my_directory/
cp /egitim/iksaglam/data/artvinensis/*.fastq.gz ./
ls *.fastq.gz
```

Let us also make a file listing all individuals and another one listing all populations for future use
```Bash
ls *.fastq.gz | cut -d'_' -f1-2 > artv.list
cut -c3-5 artv.list | uniq > pop.list
```
### Aligning reads with BWA

For aligning reads we will use [BWA](http://bio-bwa.sourceforge.net/) (Burrows-Wheeler Aligner) which is based on the [Burrows-Wheeler Transformed](https://academic.oup.com/bioinformatics/article/25/14/1754/225615?login=true) index of the reference genome. It performs gapped alignment for single-end reads, supports paired-end mapping, generates mapping quality and can give multiple hits if called for.


#### Creating BWA-MEM index
Like all other alignment tools for genome wide short reads the first step is to index the reference genome. BWA indexes the genome with an FM Index based on the Burrows-Wheeler Transform enabling memory efficient mapping. 

The basic command for indexing the genome using BWA is:
```Bash
ref=/egitim/iksaglam/ref/uvar_ref_contigs_300.fasta
bwa index -a bwtsw $ref
```

This will result in several indexes including the BWT index `/egitim/iksaglam/ref/uvar_ref_contigs_300.fasta.bwt` which the alignments will be based on. The creation of the index only needs to be performed once and does not have to be recreated for every alignment job.

#### Aligning reads with BWA-MEM

Now that we have our indexes created, we can start aligning our reads (i.e. individual fastq files) to the reference genome. Let us first create a new directory where we will store our alignment files and cd (change directory) into it.
```Bash
mkdir alignments
cd alignments/
```
Next we will use the BWA-MEM algorithm to align one of our paired-end reads to the reference genome and output results as a [SAM](https://samtools.github.io/hts-specs/SAMv1.pdf) (Sequence Alignment/Map Format) file.

```Bash
reads=/egitim/iksaglam/data/artvinensis
ref=/egitim/iksaglam/ref/uvar_ref_contigs_300.fasta
out=/egitim/iksaglam/alignments
bwa mem $ref $reads/A_CAMD04_R1.fastq.gz $reads/A_CAMD04_R2.fastq.gz > $out/A_CAMD04.sam
```
As before we can take a look at our SAM files using less.

```Bash
less -S A_CAMD04_R1.sam
```

#### SAM/BAM conversion, sorting and PCR clone (duplicate) markup/removal:

To save space we ideally want to transform our SAM file into a binary BAM format and sort by coordinates. Next when working with reduced representation libraries like RADseq data it is advantegous to keep only properly paired individuals. Finally we want to remove/mark PCR duplicates using [picard-tools](https://broadinstitute.github.io/picard/), so they do not bias variant calling and genotyping and index our final BAM file.

```Bash
samtools view -bS A_CAMD04.sam > A_CAMD04.bam
samtools sort A_CAMD04.bam A_CAMD04_sorted.bam
samtools view -b -f 0x2 A_CAMD04_sorted.bam > A_CAMD04_sorted_proper.bam
java -jar /kuacc/apps/picard/2.22.1/picard.jar MarkDuplicates INPUT=A_CAMD04_sorted_proper.bam OUTPUT=A_CAMD04_sorted_proper_rmdup.bam METRICS_FILE=A_CAMD04_metrics.txt VALIDATION_STRINGENCY=LENIENT  REMOVE_DUPLICATES=True
samtools index A_CAMD04_sorted_proper_rmdup.bam
```
For the purposes of this tutorial we are directly removing duplicates `REMOVE_DUPLICATES=True` from our resulting bam files (to primarily save space) but usually we only want to mark our duplicates and not remove them `REMOVE_DUPLICATES=False`.

To view our resulting bam files we can use the following command
```Bash
samtools view A_CAMD04_sorted_proper_rmdup.bam | less -S
```

We can also get alignment statistics using samtools flagstat
```Bash
samtools flagstat A_CAMD04_sorted_proper_rmdup.bam
```
#### Running in parallel or bulk:
Instead of aligning each paired read one by one normally we would want to do this in bulk. You can find an example shell script that can do this [here](https://github.com/iksaglam/Zonguldak/blob/main/Scripts/align_pe_reads.sh). We can use this script to execute the above pipeline and get alignment files for all individuals `artv.list` simultaneously using the following command.

```Bash
sbatch align_pe_reads.sh artv.list /egitim/iksaglam/ref/uvar_ref_contigs_300.fasta
```

Now that we have 60 BAM files at low/medium depth we can create a bamlist for the data set as a whole and for each population to use in downstream analysis in a new directory called `analyses`.

```Bash
mkdir analyses
cd analyses
ls /egitim/iksaglam/alignments/*.sorted_proper_rmdup.bam > artv.bamlist
for i in `cat pop.list`; do grep $i artv.bamlist > ${i}.bamlist; done
```
[Topics](https://github.com/iksaglam/Zonguldak/blob/main/Files/Pop_Gen.md#topics-to-be-covered)

## Genotyping and variant calling

To peform basic population genetic analysis in [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD), we first need to assign genotype probabilities at each site for each individual and determine allele frequencies for each site. The specific option in [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD) for calculating genotype probabilities and determining polymorphic sites is `-doGeno` together with `-doMAF`. When working with low/medium depth data it is always advantegous to work with genotype probabilities but we can also call genotypes using the `-doPost` option with an appopriate probability cutoff `-postCutoff`. Additionally (or alternatively) we can also create output files in beagle `-doGlf`, vcf `-dovcf` and PLINK `-doPlink` format to have the option to use other NGS analysis toolkits like [PLINK](https://www.cog-genomics.org/plink/2.0/), [GATK](https://gatk.broadinstitute.org/hc/en-us) or [STACKS](https://catchenlab.life.illinois.edu/stacks/) among others.

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

Recalling  our choice for data filtering and that we would like to output files in various formats, our final command line would look something like this:

```Bash
mkdir results_geno
ref=/egitim/iksaglam/ref/uvar_ref_contigs_300.fasta
nInd=$(wc -l artv.bamlist | awk '{print $1}')
mInd=$((${nInd}/2))

angsd -bam artv.bamlist -ref $ref -out results_geno/artv -GL 1 -doMajorMinor 1 -doMaf 1 -doGlf 2 -doGeno 5 -dovcf 1 -doPlink 2 -doPost 1 -postCutoff 0.85 -minMapQ 10 -minQ 20 -minInd $mInd -SNP_pval 1e-12 -minMaf 0.05 -nThreads 2
```
An example shell script for running the analysis on the cluster can be found [here](https://github.com/iksaglam/Zonguldak/blob/main/Scripts/get_genos.sh) and can be executed as follows:
```Bash
sbatch get_genos.sh artv /egitim/iksaglam/ref/uvar_ref_contigs_300.fasta
```

Take a look at some of the resulting files using less. Can you makes sense of them?
```Bash
zless artv.mafs.gz 
zless artv.beagle.gz  
zless artv.geno.gz  
zless artv.vcf.gz
less artv.tped
```
[Topics](https://github.com/iksaglam/Zonguldak/blob/main/Files/Pop_Gen.md#topics-to-be-covered)

## Population structure (PCA)

Now that we have calculated our genotype likelihoods and probabilities we can use these to perform a principal component analyses (PCA) summarizing genetic structure between populations, taking genotype uncertainty into account. To this end we will first estimate the covariance matrix between individuals based on genotype likelihoods using PCangsd and then decompose this matrix into principle components using R and plot the results. 

`PCangsd` takes as input genotype likelihoods in `Beagle` format so to calculate our covarinace matrix we will be using the `artv.beagle.gz` file created earlier.

```Bash
mkdir results_pca
python pcangsd.py -beagle results_pca/artv.beagle.gz -o results_pca/artv -threads 2
```

This will create an output called `art.cov` which contains our covariance matrix which we will now use to conduct an eigendecomposition and plot our first two PC axes using the Rscript [here](https://github.com/iksaglam/Zonguldak/blob/main/Scripts/plotPCA.R). We will also create a simple cluster file so that we can color label our populations in the plot.

```Bash
cut -c3-5 artv.list | sed 1iCLUSTER > CLUSTER
seq 30 | sed 1iFID > FID
printf '1\n%.0s' {1..30} | sed 1iIID > IID
paste -d' ' FID IID CLUSTER > artv.clst
```

```Bash
scripts=/egitim/iksaglam/scripts
Rscript $scripts/plotPCA.R -i results_pca/artv.cov -c 1-2 -a artv.clst -o results_pca/artv_pca_1_2.pdf
```

![PCA](https://github.com/iksaglam/Zonguldak/blob/main/Files/artv_pca_1_2.png)

Download the resulting pdf file to your local computer and take a look.

An example shell script for running the analysis on the cluster can be found [here](https://github.com/iksaglam/Zonguldak/blob/main/Scripts/get_PCA_Angsd.sh) and can be executed as follows:

```Bash
sbatch get_PCA_Angsd.sh artv
```
[Topics](https://github.com/iksaglam/Zonguldak/blob/main/Files/Pop_Gen.md#topics-to-be-covered)

## Admixture

To determine ancestry (i.e. admixture) between populations we will use `NGSadmix` together with R for plotting the results. Like `PCangsd` `NGSadmix` requires a genotype likelihood file in `beagle` format. We will conduct an admixture analysis for K=2 to K=5 clusters and also run 10 replicates for each K value. These replicate runs will allow us to use Evanno's method ([Evanno et al. 2005](https://onlinelibrary.wiley.com/doi/10.1111/j.1365-294X.2005.02553.x)) to determine the most likely K based on the log likelihood of each K.

To run the analysis we will input our previously generated genotype likelihood file `artv.beagle.gz` into `NGSadmix` and loop over each different K values 10 times.

```Bash
mkdir results_admix
x=2
while [ $x -le 5 ] 
do
  y=1
  while [ $y -le 10 ]
  do
  NGSadmix -likes results_geno/artv.beagle.gz -K $x -P 4 -seed $[RANDOM] -o results_admix/artv_admix${x}_run${y}
  y=$(( $y + 1 )) 
  done 
x=$(( $x + 1 ))
done
```

An example shell script for running the analysis on the cluster and in parallel can be found [here](https://github.com/iksaglam/Zonguldak/blob/main/Scripts/get_admix.sh) and can be executed as follows:

```Bash
sbatch get_admix.sh artv 5
```

Next, we will take the likelihood value from each run of `NGSadmix` and prepare a [Clumpak](http://clumpak.tau.ac.il/bestK.html) file to determine the most like K based on Evanno's method ([Evanno et al. 2005](https://onlinelibrary.wiley.com/doi/10.1111/j.1365-294X.2005.02553.x)).
```Bash
cd results_admix
(for log in `ls *.log`; do grep -Po 'like=\K[^ ]+' $log; done) > logfile
(for log in `ls *.log`; do grep -Po 'nPop=\K[^ ]' $log; done) > noK
paste noK logfile > admix_runs_LH.txt
```

We can now import our formatted logfile into [Clumpak](http://clumpak.tau.ac.il/bestK.html) and determibe the most likely K value for our populations.

Lastly we will visulize our admixture results for best K in the form of a barplot plot using the Rscript [here](https://github.com/iksaglam/Zonguldak/blob/main/Scripts/plotAdmix.R). For this we need to import the `.qopt` file from one of our runs for the most likely K into R. In this case best K was K=3. Additionally we will create an `info` file so we can label our population in the bar plot. We can create such a file using the simple commands below.
```Bash
cut -c3-5 artv.list | paste - artv.list > artv.info
```
Now that our info file is ready we can plot our results like below
```Bash
scripts=/egitim/iksaglam/scripts
Rscript $scripts/plotAdmix.R artv_admix3_run1.qopt artv.info
```
![Admix](https://github.com/iksaglam/Zonguldak/blob/main/Files/artv_admix3_run1.qopt.png)


Download the resulting pdf file onto your local computer and view!

[Topics](https://github.com/iksaglam/Zonguldak/blob/main/Files/Pop_Gen.md#topics-to-be-covered)

## The site frequency spectrum (SFS)

The SFS is one of the most important summary statistics in population genetics as it summarizes the distributuion of different allele frequencies along the genome. From this distribution we can calculate such statistics as Watterson's theta, Pi, Tajima's D as well as conduct demographic analysis to model past evolutionary or ecological forces affecting genetic diversity of populations (i.e. like bottlenecks or selection). The SFS can be folded or unfolded, and the latter case implies the use of an outgroup species to define the ancestral state (derived vs ancestral allele).

We will use `ANGSD` to estimate the SFS based on the methods given [here](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0037558). Like before we will use genotype likelihoods as input and from these `ANGSD` computes posterior probabilities of the Sample Allele Frequency (SAF) `-doSaf`. Next the generated `saf` files are used to estimate the SFS using the program `realSFS`. 

The SFS should be computed for each population separately. Since we have an out group (i.e. *Ph. uvarovi*) we can polorise our alleles (i.e. determine ancestral and derived states) so we want to estimate the unfolded SFS. Basic commands for estimating the SFS for a single population is given below.

### Calculate SAF file
```Bash
mkdir results_sfs
ref=/egitim/iksaglam/ref/uvar_ref_contigs_300.fasta
angsd -bam CAM.bamlist -ref $ref -anc $ref -out results_sfs/CAM -GL 1 -doSaf 1 -doCounts 1 -minMapQ 10 -minQ 20
```

Let us take a look at the out file.
```Bash
realSFS print results_sfs/CAM.saf.idx | less -S
```
These values represent the sample allele frequency likelihoods at each site. So the first value (after the chromosome and position columns) is the likelihood of having 0 copies of the derived allele, the second indicates the probability of having 1 copy and so on. Note that these values are in log format and scaled so that the maximum is 0.

### Estimate the SFS
```Bash
realSFS results_sfs/CAM.saf.idx -maxIter 100 > results_sfs/CAM.sfs
```

Take a look at the output file:
```bash
less -S results_sfs/CAM.sfs
```
The first value represent the expected number of sites with derived allele frequency equal to 0, the second column the expected number of sites with frequency equal to 1 and so on.

Lastly we can plot the SFS using a simple R script given [here](https://github.com/iksaglam/Zonguldak/blob/main/Scripts/plotSFS.R).
```Bash
scripts=/egitim/iksaglam/scripts
Rscript $scripts/plotSFS.R results_sfs/CAM.sfs CAM 0 results_sfs/CAM.sfs.pdf
```
![SFS](https://github.com/iksaglam/Zonguldak/blob/main/Files/CAM.sfs.png)

Download the resulting pdf file onto your local computer and view!

Ideally we would of course want to cycle through all populations in parallel. An example shell script for running the analysis on the cluster and in parallel can be found [here](https://github.com/iksaglam/Zonguldak/blob/main/Scripts/get_sfs.sh) and can be executed as follows:

```Bash
sbatch get_sfs.sh pop.list
```
[Topics](https://github.com/iksaglam/Zonguldak/blob/main/Files/Pop_Gen.md#topics-to-be-covered)

## Nucleotive diversity and neutrality

Once we have computed the allele frequency posterior probabilities and the SFS for each population we can use these to calculate nucleotide diversity and related statistics for each population in ANGSD without relying on called genotypes.

Again it is better to calculate such statistics for each population separately. Calculation of diversity statistics for one population can be achieved using the following pipeline.
```Bash
mkdir results_diversity
realSFS saf2theta results_sfs/CAM.saf.idx -sfs results_sfs/CAM.sfs -outname results_diversity/CAM
thetaStat do_stat results_diversity/CAM.thetas.idx
```

This will result in an output file named `CAM.thetas.idx.pestPG`. We can view this file using less.
```Bash
less -S CAM.thetas.idx.pestPG
```

The output will contain 5 different estimators of theta: Watterson, Pi, FuLi, fayH, L. and 5 different neutrality test statistics: Tajima's D, Fu&Li F's, Fu&Li's D, Fay's H, Zeng's E. The final column is the effetive number of sites with data for each RAD contig.

```
#(indexStart,indexStop)(firstPos_withData,lastPos_withData)(WinStart,WinStop)   Chr     WinCenter       tW      tP      tF      tH      tL      Tajima  fuf     fud     fayh    zeng    nSites
(0,344)(1,345)(0,345)   R000002 172     2.268923        1.798771        3.482656        0.477495        1.138133        -0.808122       -0.936711       -0.775660       0.630364        -0.571185       344
(0,298)(1,308)(0,308)   R000004 154     14.467825       17.316357       4.444719        13.537845       15.427101       0.899640        1.510947        1.377009        0.315925        0.091875        298
(0,355)(1,356)(0,356)   R000006 178     4.619073        5.465138        2.793336        7.232932        6.349035        0.781203        0.854562        0.676516        -0.441681       0.477114        355
(0,424)(1,451)(0,451)   R000007 225     4.602208        4.815590        3.410900        5.577429        5.196509        0.197677        0.450634        0.442741        -0.190996       0.164437        424
(0,287)(1,290)(0,290)   R000009 145     4.320633        5.315162        2.142283        4.696980        5.006071        0.975362        1.071780        0.851829        0.164373        0.200529        287
(0,372)(1,373)(0,373)   R000010 186     4.796857        4.255665        6.827616        2.108999        3.182332        -0.482878       -0.797435       -0.729749       0.517732        -0.430593       372
(0,363)(1,364)(0,364)   R000011 182     3.576046        3.913829        2.765189        2.245555        3.079692        0.392182        0.451592        0.368236        0.528321        -0.171232       363
(0,371)(1,391)(0,391)   R000012 195     4.542996        4.553729        4.873980        3.009494        3.781611        0.010060        -0.103836       -0.124303       0.391854        -0.213097       371
.
.
.
```

Alternatively if you have full genome data you might want to do a sliding window analysis by adding -win/-step arguments to the last command `thetaStat`. In this case the final column above (nSites) would be the effetive number of sites with data in the window.
```
thetaStat do_stat results_diversity/CAM.thetas.idx -win 50000 -step 10000
```

As always we would like to run the analysis in paralell for all populations. An example shell script for running the analysis on the cluster and in parallel can be found [here](https://github.com/iksaglam/Zonguldak/blob/main/Scripts/get_diversity.sh) and can be executed as follows:
```Bash
sbatch get_diversity.sh pop.list
```

If we wish we can also plot or make a table comparing our results between populations for selected statistics. To do this we would first need to extract the relative columns from each population and write them into a new file. For example let us say we want to compare and plot ThetaW, ThetaPi and Tajima's D estimates between populations. To create such a spreadsheet we could do something like this:
```Bash
for i in `cat ../pop.list`; do cut -f2,4,5,9,14 ${i}.thetas.idx.pestPG | sed 1d | sed "s/^/$i\t/" >> all_pops_diversity.tsv ; done
sed -i 1i'Pop\tRADloci\ttW\ttP\tTajD\tNsites’ all_pops_diversity.tsv
```
Next we can use a simple R script as given [here](https://github.com/iksaglam/Zonguldak/blob/main/Scripts/plotDiv.R) to plot our results and create a summary table.
```Bash
scripts=/egitim/iksaglam/scripts
Rscript $scripts/plotDiv.R all_pops_diversity.tsv
```
![Tajima's D](https://github.com/iksaglam/Zonguldak/blob/main/Files/artv_pops_TajimaD.png)

[Topics](https://github.com/iksaglam/Zonguldak/blob/main/Files/Pop_Gen.md#topics-to-be-covered)

## Population genetic differentiation
Here we are going to estimate allele frequency differentiation between populations using the Fst metrics. We can achieve this in ANGSD without relying on genotype calls  by directly working with the sample allele frequencies likelihoods we calculated before (i.e. `saf` files). However to estimate Fst values between populations we also need to estimate the joint SFS between any two populations (2D-SFS) which will serve as prior information for estimating Fst values. As an example we will calculate the Fst between the CAM and CAY populations.

First we need to estimate the 2D-SFS between these two populations. This can be achieved by using the `realSFS` program. An important issue when doing this is to be sure that we are comparing the exact same sites between populations. Luckily ANGSD does this automatically and considers only a set of overlapping sites.
```Bash
realSFS results_sfs/CAM.saf.idx results_sfs/CAY.saf.idx > results_sfs/CAM.CAY.2dsfs
```

Next we can use `saf` files of both populations together with the `2D-SFS` to calculate joint allele frequency probabilities at each site.
```Bash
mkdir results_fst
realSFS fst index results_sfs/CAM.saf.idx results_sfs/CAY.saf.idx -sfs results_sfs/CAM.CAY.2dsfs -fstout results_fst/CAM.CAY
```
This command will create and output file ending with `*.fst.idx` which stores per-site FST indexes. We can take a look at this file as follows.
```Bash
realSFS fst print results_fst/CAM.CAY.fst.idx | less -S
```
Here columns are: `RAD_loci, position, (a), (a+b)` values, where FST is defined as `a/(a+b)`. Note that FST on multiple SNPs is calculated as `sum(a)/sum(a+b)`.

Lastly we can calculate the global Fst between these two populations using the following command.
```Bash
realSFS fst stats results_fst/CAM.CAY.fst.idx 2> results_fst/CAM.CAY_global.fst
```

We can take a look at this final estimate by using less.
```Bash
less CAM.CAY_global.fst
```

```
 -> Assuming idxname:results_fst/CAM.CAY.fst.idx
 -> Assuming .fst.gz file: results_fst/CAM.CAY.fst.gz
 -> FST.Unweight[nObs:13835235]:0.116738 Fst.Weight:0.260466
```

Note, if you have full genome data you might want to do a windowed scan analysis at this step.
```Bash
realSFS fst stats2 results_fst/CAM.CAY.fst.idx -win 50000 -step 10000 -whichFST 1 > results_fst/CAM.CAY.fst.txt
```

As before we can perform all pairwise comparisons in bulk using appropriate shell scripts found [here](https://github.com/iksaglam/Zonguldak/blob/main/Scripts/get_2Dsfs.sh) for calculating the 2D-SFS and [here](https://github.com/iksaglam/Zonguldak/blob/main/Scripts/get_fst.sh) for calculating Fst statistics.

2D-SFS between all populations
```Bash
sbatch get_2Dsfs.sh pop.list
```
Fst between all populations
```Bash
sbatch get_fst.sh pop.list
```

Lastly we can plot our results in the form of a heat map illustrating pairwise genetic differentiation (Fst) between populations. To do this first let us make a table summarizing global fst values between each population.
```
cd results_fst
(for fst in `ls *.fst`; do grep -Po 'Fst.Weight:\K[^ ]+' $fst; done) | sed 1iFst > fstfile
ls *.fst | cut -d'_' -f1 | tr '.' '\t' | sed 1i'Pop1\tPop2' > pops
paste pops fstfile > all_pops_fst.tsv
```
We can then input this table into an Rscript given [here](https://github.com/iksaglam/Zonguldak/blob/main/Scripts/plotFst.R) and plot a heat map.
```Bash
scripts=/egitim/iksaglam/scripts
Rscript $scripts/plotFst.R all_pops_fst.tsv
```
![Fst](https://github.com/iksaglam/Zonguldak/blob/main/Files/artv_pops_Fst.png)

[Top](https://github.com/iksaglam/Zonguldak/blob/main/Files/Pop_Gen.md#a-tutorial-for-some-basic-population-genetic-analyses-using-radseq-data)
