## Working without a reference genome: RAD loci discovery

Our goal in this tutorial is to understand the basics of RAD loci discovery from raw RADseq data (fastq files). Once individual RAD loci have been discovered thet can then be assembled into longer 300-800 bp contigs using forward and reverse reads and  the resulting contigs can be used as a de-novo set of reference RAD-contigs (partial reference genome) for all subsequent downstream analysis (i.e. mapping, variant calling, population genetic analysis, see [here](https://github.com/iksaglam/Zonguldak/blob/main/Files/Pop_Gen.md)). More information on this procedure can be found [here](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.13732) and a similar approach has been implemented [here](https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.15253) as well. A good paper for the advantages of using such an approach can be found [here](https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13324).


#### In this tutorial we will cover the following topics:
- Qaulity filtering of fastq files
- Building a hash table for all reads with ID and count information
- Building a pairwise alignment map between all reads
- Locus discovery

#### Data:
To discover a set of unique RAD loci we will use forward reads from 6 fastq files from the Ovit population of *Ph. uvarovi*. The data can be found here `/egitim/iksaglam/data/uvarovi`. As a first step let us now copy this data into our own directory and list to make sure everything was copied correctly.

```Bash
cd ~/my_directory/
cp /egitim/iksaglam/data/uvarovi/*.fastq ./
ls *.fastq
```

Let us also make a file listing all individuals for future use

```Bash
ls *.fastq | cut -d'_' -f1-2 > U_OVT.list
```


#### Qaulity filtering of fastq files:

Now let us truncate reads to 80bp and remove low quality reads (i.e. < Q20). We will be using a custom perl script given [here](https://github.com/iksaglam/Zonguldak/blob/main/Scripts/QualityFilter.pl), but any quality filtering software can be used as well.

```Bash
scripts=/egitim/iksaglam/scripts
for i in `cat U_OVT.list`
do
perl $scripts/QualityFilter.pl ${i}_R1.fastq > ${i}_L80P80.fastq
done
```

We can quickly take a look at on of these files using less.

```Bash
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
Next we will use another custom perl script given [here](https://github.com/iksaglam/Zonguldak/blob/main/Scripts/HashSeqs.pl) to build a hash table in fasta format containing ID and count information for all reads within each fastq file.  

```Bash
scripts=/egitim/iksaglam/scripts
for i in `cat U_OVT.list`
do
perl $scripts/HashSeqs.pl ${i}_L80P80.fastq ${i} > ${i}_L80P80.hash
done
```
Let us again take a quick look at on of these files using less.

```Bash
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
Now we are ready to build our alignment map. For the purposes of this tutorial we will be using the aligner novaalign. Firstly we need to concatenate all files into a single hash file in fasta format and index this file for mapping.

```Bash
cat *.hash > U_OVT.fasta
novoindex U_OVT.index U_OVT.fasta
```

Next we will map our hash file containing all reads (i.e. U_OVT.fasta) back on to its own index. This way we will obtain pairwise alignment statictics between each read in our hash file, essentially giving us an alignment map for all pairwise read comparisons with enough similarity to align. To make sure that novoalign does not report only the best alignments we will be using the "exhaustive" mapping option `-rE` with `-t 180` (%67 percent identity). Here we will set `-rE` to `50` (reports the top 50 alignments) for quick results but feel free to go as high as you want.
  
```Bash
novoalign  -rE 50 -t 180 -d U_OVT.index -f U_OVT.fasta > U_OVT.novo
```

Let us take a look at the first few lines of our alignment map

<p><details><summary>Click to view first few lines of output</summary>

```
# novoalign (V2.08.03 - Build Oct 16 2012 @ 12:37:58 - A short read aligner with qualities.
# (C) 2008,2009,2010,2011 NovoCraft Technologies Sdn Bhd.
# License file not found.
# Licensed for evaluation, educational, and not-for-profit use only.
#  novoalign -r E 50 -t 180 -d U_OVT.index -f U_OVT.fasta 
# Starting at Mon Jan 31 16:16:21 2022
# Interpreting input files as FASTA.
# Index Build Version: 2.8
# Hash length: 11
# Step size: 1
>U_OVTD04;1;4674        S       GCCATGACGGCCGAGCTCTGGAAGCGCAGGTCGGTCTTGAAGTCCTGAGCGATCTCCCTGACCAACCGCTGGAAGGGCAG        .       R       0       1       >U_OVTE11;1;10984       1       F       .       .       .
>U_OVTD04;1;4674        S       GCCATGACGGCCGAGCTCTGGAAGCGCAGGTCGGTCTTGAAGTCCTGAGCGATCTCCCTGACCAACCGCTGGAAGGGCAG        .       R       0       1       >U_OVTE08;1;9513        1       F       .       .       .
>U_OVTD04;1;4674        S       GCCATGACGGCCGAGCTCTGGAAGCGCAGGTCGGTCTTGAAGTCCTGAGCGATCTCCCTGACCAACCGCTGGAAGGGCAG        .       R       0       1       >U_OVTD05;1;10052       1       F       .       .       .
>U_OVTD04;1;4674        S       GCCATGACGGCCGAGCTCTGGAAGCGCAGGTCGGTCTTGAAGTCCTGAGCGATCTCCCTGACCAACCGCTGGAAGGGCAG        .       R       0       1       >U_OVTD08;1;14778       1       F       .       .       .
>U_OVTD04;1;4674        S       GCCATGACGGCCGAGCTCTGGAAGCGCAGGTCGGTCTTGAAGTCCTGAGCGATCTCCCTGACCAACCGCTGGAAGGGCAG        .       R       0       1       >U_OVTD04;1;4674        1       F       .       .       .
>U_OVTD04;1;4674        S       GCCATGACGGCCGAGCTCTGGAAGCGCAGGTCGGTCTTGAAGTCCTGAGCGATCTCCCTGACCAACCGCTGGAAGGGCAG        .       R       0       1       >U_OVTE09;1;20320       1       F       .       .       .
>U_OVTD04;1;4674        S       GCCATGACGGCCGAGCTCTGGAAGCGCAGGTCGGTCTTGAAGTCCTGAGCGATCTCCCTGACCAACCGCTGGAAGGGCAG        .       R       30      0       >U_OVTD04;135299;2      1       F       .       .       .       75T>G
>U_OVTD04;1;4674        S       GCCATGACGGCCGAGCTCTGGAAGCGCAGGTCGGTCTTGAAGTCCTGAGCGATCTCCCTGACCAACCGCTGGAAGGGCAG        .       R       30      0       >U_OVTE08;80772;3       1       F       .       .       .       73G>A
>U_OVTD04;1;4674        S       GCCATGACGGCCGAGCTCTGGAAGCGCAGGTCGGTCTTGAAGTCCTGAGCGATCTCCCTGACCAACCGCTGGAAGGGCAG        .       R       30      0       >U_OVTE08;3589;12       1       F       .       .       .       77C>G
>U_OVTD04;1;4674        S       GCCATGACGGCCGAGCTCTGGAAGCGCAGGTCGGTCTTGAAGTCCTGAGCGATCTCCCTGACCAACCGCTGGAAGGGCAG        .       R       30      0       >U_OVTE09;22664;7       1       F       .       .       .       20T>G
>U_OVTD04;1;4674        S       GCCATGACGGCCGAGCTCTGGAAGCGCAGGTCGGTCTTGAAGTCCTGAGCGATCTCCCTGACCAACCGCTGGAAGGGCAG        .       R       30      0       >U_OVTD08;40935;6       1       F       .       .       .       38C>T
>U_OVTD04;1;4674        S       GCCATGACGGCCGAGCTCTGGAAGCGCAGGTCGGTCTTGAAGTCCTGAGCGATCTCCCTGACCAACCGCTGGAAGGGCAG        .       R       30      0       >U_OVTD04;134146;2      1       F       .       .       .       57G>C
>U_OVTD04;1;4674        S       GCCATGACGGCCGAGCTCTGGAAGCGCAGGTCGGTCTTGAAGTCCTGAGCGATCTCCCTGACCAACCGCTGGAAGGGCAG        .       R       30      0       >U_OVTD05;30610;6       1       F       .       .       .       1T>G
>U_OVTD04;1;4674        S       GCCATGACGGCCGAGCTCTGGAAGCGCAGGTCGGTCTTGAAGTCCTGAGCGATCTCCCTGACCAACCGCTGGAAGGGCAG        .       R       30      0       >U_OVTE11;108248;2      1       F       .       .       .       2T>C
>U_OVTD04;1;4674        S       GCCATGACGGCCGAGCTCTGGAAGCGCAGGTCGGTCTTGAAGTCCTGAGCGATCTCCCTGACCAACCGCTGGAAGGGCAG        .       R       30      0       >U_OVTE09;4653;12       1       F       .       .       .       7T>A
>U_OVTD04;1;4674        S       GCCATGACGGCCGAGCTCTGGAAGCGCAGGTCGGTCTTGAAGTCCTGAGCGATCTCCCTGACCAACCGCTGGAAGGGCAG        .       R       30      0       >U_OVTE11;43254;5       1       F       .       .       .       22T>A
>U_OVTD04;1;4674        S       GCCATGACGGCCGAGCTCTGGAAGCGCAGGTCGGTCTTGAAGTCCTGAGCGATCTCCCTGACCAACCGCTGGAAGGGCAG        .       R       30      0       >U_OVTD05;8641;9        1       F       .       .       .       77T>G
>U_OVTD04;1;4674        S       GCCATGACGGCCGAGCTCTGGAAGCGCAGGTCGGTCTTGAAGTCCTGAGCGATCTCCCTGACCAACCGCTGGAAGGGCAG        .       R       30      0       >U_OVTE09;130083;2      1       F       .       .       .       9A>G
>U_OVTD04;1;4674        S       GCCATGACGGCCGAGCTCTGGAAGCGCAGGTCGGTCTTGAAGTCCTGAGCGATCTCCCTGACCAACCGCTGGAAGGGCAG        .       R       30      0       >U_OVTE09;3742;13       1       F       .       .       .       48C>A
>U_OVTD04;1;4674        S       GCCATGACGGCCGAGCTCTGGAAGCGCAGGTCGGTCTTGAAGTCCTGAGCGATCTCCCTGACCAACCGCTGGAAGGGCAG        .       R       30      0       >U_OVTD05;96175;3       1       F       .       .       .       61G>A
>U_OVTD04;1;4674        S       GCCATGACGGCCGAGCTCTGGAAGCGCAGGTCGGTCTTGAAGTCCTGAGCGATCTCCCTGACCAACCGCTGGAAGGGCAG        .       R       30      0       >U_OVTE11;12356;8       1       F       .       .       .       29T>G
>U_OVTD04;1;4674        S       GCCATGACGGCCGAGCTCTGGAAGCGCAGGTCGGTCTTGAAGTCCTGAGCGATCTCCCTGACCAACCGCTGGAAGGGCAG        .       R       30      0       >U_OVTE08;54664;4       1       F       .       .       .       80T>G
>U_OVTD04;1;4674        S       GCCATGACGGCCGAGCTCTGGAAGCGCAGGTCGGTCTTGAAGTCCTGAGCGATCTCCCTGACCAACCGCTGGAAGGGCAG        .       R       30      0       >U_OVTE09;81380;3       1       F       .       .       .       10T>G
>U_OVTD04;1;4674        S       GCCATGACGGCCGAGCTCTGGAAGCGCAGGTCGGTCTTGAAGTCCTGAGCGATCTCCCTGACCAACCGCTGGAAGGGCAG        .       R       30      0       >U_OVTD08;28511;7       1       F       .       .       .       20T>G
>U_OVTD04;1;4674        S       GCCATGACGGCCGAGCTCTGGAAGCGCAGGTCGGTCTTGAAGTCCTGAGCGATCTCCCTGACCAACCGCTGGAAGGGCAG        .       R       30      0       >U_OVTD08;3367;15       1       F       .       .       .       26T>G
>U_OVTD04;1;4674        S       GCCATGACGGCCGAGCTCTGGAAGCGCAGGTCGGTCTTGAAGTCCTGAGCGATCTCCCTGACCAACCGCTGGAAGGGCAG        .       R       30      0       >U_OVTE09;45064;5       1       F       .       .       .       39T>G
>U_OVTD04;1;4674        S       GCCATGACGGCCGAGCTCTGGAAGCGCAGGTCGGTCTTGAAGTCCTGAGCGATCTCCCTGACCAACCGCTGGAAGGGCAG        .       R       30      0       >U_OVTD08;2035;18       1       F       .       .       .       77T>G
>U_OVTD04;1;4674        S       GCCATGACGGCCGAGCTCTGGAAGCGCAGGTCGGTCTTGAAGTCCTGAGCGATCTCCCTGACCAACCGCTGGAAGGGCAG        .       R       30      0       >U_OVTE09;34851;6       1       F       .       .       .       50A>C
>U_OVTD04;1;4674        S       GCCATGACGGCCGAGCTCTGGAAGCGCAGGTCGGTCTTGAAGTCCTGAGCGATCTCCCTGACCAACCGCTGGAAGGGCAG        .       R       30      0       >U_OVTE08;85807;3       1       F       .       .       .       9A>G
>U_OVTD04;1;4674        S       GCCATGACGGCCGAGCTCTGGAAGCGCAGGTCGGTCTTGAAGTCCTGAGCGATCTCCCTGACCAACCGCTGGAAGGGCAG        .       R       30      0       >U_OVTE11;35151;5       1       F       .       .       .       76T>G
>U_OVTD04;1;4674        S       GCCATGACGGCCGAGCTCTGGAAGCGCAGGTCGGTCTTGAAGTCCTGAGCGATCTCCCTGACCAACCGCTGGAAGGGCAG        .       R       30      0       >U_OVTD05;5820;10       1       F       .       .       .       29T>G
```
</details><p>  

#### Locus discovery:
Now that we have an alignment map we are ready to discover unique (i.e. individual) RAD loci. To do this we will use a custom perl script given [here](https://github.com/iksaglam/Zonguldak/blob/main/Scripts/IdentifyLoci3.pl). The first lines of this script contains important criteria which we will have to set and which will influence our results (i.e. the number loci we end up with).

```Perl
#!/usr/bin/perl

#######################################################################################
$max_alignment_score = 30;
$divergence_factor = 2;
$min_count = 0; #zero turns off to use min_flag instead
$min_flag = 0; #zero turns off to use min_count instead
$max_internal_alignments = 1;
$min_internal_alignments = 0;
$max_external_alignments = 10; # twice the df of the total number of samples (i.e. individuals) in the alignment
$min_external_alignments = 3; # approximately half of the df of the total number of samples (i.e. individuals) in the alignment
$max_total_alignments = 11; # max_ext_align + max_int_align
$min_total_alignments = 3; # min_ext_align + min_int_align
$min_alleles = 1;
$max_alleles = 4;
$min_samples_per_allele = 2;
$max_samples_per_allele = 6; # total number of samples (i.e. individuals) in the alignment
#######################################################################################
```

Let us now use the above settings to discover loci 

```Bash
scripts=/egitim/iksaglam/scripts
perl $scripts/IdentifyLoci3.pl U_OVT.novo > U_OVT.loci
```
and look at the results
```  Bash
tail U_OVT.loci | less 
```  

```
>R063751{U_OVTD04;6739;11|U_OVTD05;42075;5|U_OVTD08;43061;6|U_OVTE11;30125;6}
AAGAAATAATTAATCCCCTTGAACAACAAAGATAATTACCAATTCTAACACCTGCAACTTGGAACAAAGCTGAAGAATAA
>R063752{U_OVTD04;135968;2|U_OVTD05;10314;9|U_OVTD08;63872;5|U_OVTE09;131719;2}
AAATGTCAATAAATAAAAGATGTGACTCATAACGGAATCTGGGGACTACTGAAATGGGTGAACTCCAGTAATTTTTTTGT
>R063753{U_OVTD04;54209;5|U_OVTD05;3040;13|U_OVTE08;27482;6|U_OVTE11;65041;4}
CTTCTGCACAGTGACAGCGATTTAGGCTGAGCAAATAACGTGATGGCCTGGGCTGTTAAACTGAATAAAATTTGAACGGT
>R063754{U_OVTD04;56467;5|U_OVTD05;66468;4|U_OVTD08;86315;4|U_OVTE11;67825;4}
ATCTTCACTTCAACATTCTAGATGCCTGGAACAATCAATGGCAGGAACAAATTTGCCTGACACATCGCCACTCCTTGGGA
>R063755{U_OVTD08;30892;7|U_OVTE08;79860;3|U_OVTE09;69063;4|U_OVTE11;116063;2}
TATGTAATAACCTACTTCTTTTGGCATTGCTCTAAATCAAATACGAGTTAATTTCAAACCCTCTTCGCATTATACTAAGG
```
We can see that we have discovered 63,751 individual RAD loci.
  
Now let us change `$max_alignment_score` to `60` (meaning we are allowing two SNPs within each RAD loci) and redo the analaysis
  
```Bash
scripts=/egitim/iksaglam/scripts
perl $scripts/IdentifyLoci3.pl U_OVT.novo > U_OVT_v2.loci
```  
```Bash
tail U_OVT_v2.loci | less 
```    
```
>R048758{U_OVTD08;94045;3|U_OVTE08;32552;5|U_OVTE09;110666;2|U_OVTE11;17694;7}
CCATGAGCTGGGCAGCAGTGCAACGGTGTTAGTTAGTTCTGGTTTCGGTCTAAAAGAAAGCAGTTATGCCAAGAGATAGC
>R048759{U_OVTD05;90043;3|U_OVTD08;32614;7|U_OVTE08;84603;3|U_OVTE09;128650;2|U_OVTE11;46068;5}
GCCCGAAGATGTTGCAATACTACTTGTCCCAGGTAACTTGCTGCCGAAAATAGCTCCGTCCTAATCCCAGCTGACTGCTG
>R048760{U_OVTD04;6739;11|U_OVTD05;42075;5|U_OVTD08;43061;6|U_OVTE11;30125;6}
AAGAAATAATTAATCCCCTTGAACAACAAAGATAATTACCAATTCTAACACCTGCAACTTGGAACAAAGCTGAAGAATAA
>R048761{U_OVTD04;135968;2|U_OVTD05;10314;9|U_OVTD08;63872;5|U_OVTE09;131719;2}
AAATGTCAATAAATAAAAGATGTGACTCATAACGGAATCTGGGGACTACTGAAATGGGTGAACTCCAGTAATTTTTTTGT
>R048762{U_OVTD08;30892;7|U_OVTE08;79860;3|U_OVTE09;69063;4|U_OVTE11;116063;2}
TATGTAATAACCTACTTCTTTTGGCATTGCTCTAAATCAAATACGAGTTAATTTCAAACCCTCTTCGCATTATACTAAGG
```
Since we allowed more variability within each RAD loci the total number of loci dropped from `63,751` to `48,762`. Play around with some of the other parameters to see how they influence your results.
  
Finally let say we are comfortable with U_OVT.loci and these are the loci we want to move forward with. Let us do some cleaning up and convert our loci file into a clean looking fasta file using the custom perl script given [here](https://github.com/iksaglam/Zonguldak/blob/main/Scripts/SimplifyLoci2.pl).
  
```Bash
scripts=/egitim/iksaglam/scripts
SimplifyLoci2.pl U_OVT.loci | grep _1 -A1 --no-group-separator | sed 's/_1//' > U_OVT_ref.fasta
```
Don't forget to take a look at your nice list of RAD loci!
```Bash
less U_OVT_ref.fasta
```
```
>R000001
CTGTCAACTCCAACGACGTCCGGCGGGCTCACTTGTCCCAACTGAAACACTACAGGGTATCTTAGTTTTTGTATATTGCT
>R000002
CCATGAGCTGGGCAGCAGTCGTCCTTGGAAGCCTAGGCCCTTCAAAGGCTGTAGTGCAAGGGGGTTAGTGTTAGTTGGAC
>R000003
TGCCCAGAACTTTGTAAACCCAGACTACTCACAGTGATGGTAACAGGTCCTTGGTTCGGACAATGCTGTAATTTAATTTC
>R000004
ATGTAAATCCTCTCGGACAGCGTCTCCGTCTGTCTACACTGTTCTCGTAACATGGAATGCACAAGAAAATCCTACGTTAC
>R000005
CCATTAGCTGGGCAGCAGTCGTCCTTGGAAGCCTAGACCCTTCAAAGGCTGTAGTGCAACAGTGTTAGTTAGTTAGTTGA
.
.
.
>R063751
AAGAAATAATTAATCCCCTTGAACAACAAAGATAATTACCAATTCTAACACCTGCAACTTGGAACAAAGCTGAAGAATAA
>R063752
AAATGTCAATAAATAAAAGATGTGACTCATAACGGAATCTGGGGACTACTGAAATGGGTGAACTCCAGTAATTTTTTTGT
>R063753
CTTCTGCACAGTGACAGCGATTTAGGCTGAGCAAATAACGTGATGGCCTGGGCTGTTAAACTGAATAAAATTTGAACGGT
>R063754
ATCTTCACTTCAACATTCTAGATGCCTGGAACAATCAATGGCAGGAACAAATTTGCCTGACACATCGCCACTCCTTGGGA
>R063755
TATGTAATAACCTACTTCTTTTGGCATTGCTCTAAATCAAATACGAGTTAATTTCAAACCCTCTTCGCATTATACTAAGG
```
