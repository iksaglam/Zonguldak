# Winter course in analysis of RADseq data for molecular ecology and population genetics


## Working without a reference genome: RAD loci discovery

Our goal in this tutorial is to discover individual RAD loci from raw RADseq data (fastq files). These RAD loci can then be assembled into longer 300-800 bp contigs using forward and reverse reads. This set of reference RAD-contigs than can be used as a de-novo reference genome for all subsequent downstream analysis. More information on this procedure can be found [here](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.13732) and a similar approach has been implemented [here](https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.15253) as well.


#### In this tutorial we will cover the following topics:
- Qaulity filtering of fastq files
- Building a hash table for all reads with ID and count information
- Building a pairwise alignment map between all reads
- Determining criteria for locus discovery

#### Data:
To discover a set of unique RAD loci we will use forward reads from 6 fastq files from the Ovit population of *Ph. uvarovi*. The data can be found here `/egitim/iksaglam/data/uvarovi`. As a first step let us now copy this data into our own directory and list to make sure everything was copied correctly.

```
cd ~/my_directory/
cp /egitim/iksaglam/data/uvarovi/*.fastq ./
ls *.fastq
```

Let us also make a file listing all individuals for future use

```
ls *.fastq | cut -d'_' -f1-2 > U_OVT.list
```


#### Qaulity filtering of fastq files:

Now let us truncate reads to 80bp and remove low quality reads (i.e. < Q20). Here we will be using a custom script but any quality filtering software can be used.

```
scripts=/egitim/iksaglam/scripts
for i in `cat U_OVT.list`
do
perl $scripts/QualityFilter.pl ${i}_R1.fastq > ${i}_L80P80.fastq
done
```

We can quickly take a look at on of these files using less.

```
less U_OVTD04_L80P80.fastq
```


<p><details><summary>Click to view first few lines of output</summary>

```
@HS1:444:C6AWNACXX:4:1101:1448:2123 2:N:0:TGACCA
CTACGTGCAGATCTTCAGGTGCTGGGGGAGCTGTGAGAAATGGGCGACCAGCATCGCAATATATTCGAATAGGGACAGTT
+
IIIIIEHIIIIIGGIIIIIFDAHGIIII@FEDG;7@==D@ACEE@B63>??<C>ABBD?@CCCDEEDA?ABCCDD?BD@C
@HS1:444:C6AWNACXX:4:1101:2428:2201 1:N:0:TGACCA
ATTACGCCGAGGTATTTGTATGAATTTACTACTTCTAGCCTCTGATTTCCACATGTGAAAACGTCCTTTGAAGAGATTCT
+
IIJJIIIGGFGGDGHGIGHHIICAHIGHIIJIIIGHC9CGGGGGGHHHE;?E@;??CBE>CCC=ABBCAC>CCC?CCCAC
@HS1:444:C6AWNACXX:4:1101:4611:2167 1:N:0:TGACCA
TGCTTTTAAAAAGTATGTTAGAAGAAGAATGTGACTCACGGTAGAATTATAATCACATCACGAATCCTAGTTAATTATTC
+
EGIJIJJIGHIIJCGIIIGIGIIIJIIJJJJIJJJJJJJJJFHHHIJIIIIIIGHJJGJHGHCEEFEEEECEECADDEFE
@HS1:444:C6AWNACXX:4:1101:5404:2136 1:N:0:TGACCA
AAGAGGAAAAATCAAATGCATGTGAAGGTCGTTTGAGCCTCTAGAGCTATCCGCCTTTTGGATGTCCAATACCACTACGA
```
</details><p>


#### Building a hash table for all reads:
Next we will use another custom script to build a hash table in fasta format containing ID and count information for all reads within each fastq file.  

```
scripts=/egitim/iksaglam/scripts
for i in `cat U_OVT.list`
do
perl $scripts/HashSeqs.pl ${i}_L80P80.fastq ${i} > ${i}_L80P80.hash
done
```
Let us again take a quick look at on of these files using less.

```
less U_OVTD04_L80P80.hash
```

<p><details><summary>Click to view first few lines of output</summary>

```
>U_OVTD04;1;4674
GCCATGACGGCCGAGCTCTGGAAGCGCAGGTCGGTCTTGAAGTCCTGAGCGATCTCCCTGACCAACCGCTGGAAGGGCAG
>U_OVTD04;2;4509
GATTTAAATTATATTTTAGGCCATAACTCAAAAACTATTGCAGATATTAACGTAATTCTTTCACGGAATTACACAGAAGA
>U_OVTD04;3;4185
AAGCCAGCGAGGCCTACTTGGTGGGACTTTTCGAAGACACCAACCTGTGCGCCATTCATGCCAAGCGAGTCACCATCATG
>U_OVTD04;4;4181
CGACACGCCCCATTCTGCTGTCACATAGATTGGAGCGCGCTGATGCCTCATTGCCGATTAATTACTACAGGCAAATTTAT
>U_OVTD04;5;3858
TCATCTGCATATAAGTGAAATGAACTGTATTTAAATACCGTGGATATGTCATTGATATATATGGAGAAACATAAAGGTCC
>U_OVTD04;6;3114
CTTACATCCACTATTCTCCCAATACGGCAGAAAATCAAGTGAAACTCTTCAGCGAAGATATAGATAGGCTGAGTAACTGG
>U_OVTD04;7;3095
ATCATCACTTCAACATTCTAGATGCCTGGAACAATCAATGGCAGGAACAAATTTGCCTGACACATCGCCACTTAATTGTA
>U_OVTD04;8;3029
TTGACTGCTGTTCGCCATGATGGTTGTCTGGACTTTAGTCGCTGCGGTAAAGGCTCTGGGAGCTCAGAGTGTATTGGCAG
>U_OVTD04;9;2416
CTGGTATTAAATAGATCAATATCGATAGAGTTTGAGAGGCGCATCATTCTATCAACTGGACTGTTGAAGCCATATTTAGT
>U_OVTD04;10;2129
TCATCTGAGACAACCCACAGAGAGATATCAACTGCTTAGAGAGTAGATCTAACTCTCGGACTCAGAGGGAGTAAGTTGTA
```
</details><p>  


#### Building a pairwise alignment map between all reads:
Now we are ready to build our alignment map. For the purposes of this tutorial we will be using the aligner novaalign. Firstly we need to concatenate all files into a single hash file (i.e. fasta file) and index this file for mapping.

```
cat *.hash > U_OVT.fasta
novoindex U_OVT.index U_OVT.fasta
```

Next we will map our hash file containing all reads (i.e. U_OVT.fasta) back on to its own index. This way we will obtain pairwise alignment statictics between each read in our hash file, essentially giving us an alignment map for all pairwise read comparisons with enough similarity to align. To make sure that novoalign does not report only the best alignments we will be using the "exhaustive" mapping option (-rE) with -t 180 (%67 percent identity). Here we will set -rE to 50 (reports the top 50 alignments) for quick results but feel free to go as high as you want.
  
```
novoalign  -r E 50 -t 180 -d U_OVT.index -f U_OVT.fasta > U_OVT.novo
```
  
  
  
  
![image](www.e0dad810-dc1a-4674-8759-8af6ff6bc7bc.png)




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


