# Genome assembly by alignment/mapping of nanopore reads to a reference

In this exercise we introduce a simple pipeline for the alignment of Nanopore sequencing data to a known reference sequence. It is by no means a comprehensive guide but should be viewed as a resource for first time users.
Important aspects such as quality control of the initial dataset and the usage of quality scores in the alignment are not covered. 

Software:
- [poretools](https://github.com/arq5x/poretools) (Toolkit for working with Nanopore data)
- [LAST](http://last.cbrc.jp) (Alignment software for long reads)
- [samtools](http://www.htslib.org) (A suite of programs for interacting with HT sequencing data)
- [Tablet](https://ics.hutton.ac.uk/tablet/) (Alignment/Assembly Viewer)

The executables to all these packages should already be in your Path. In case you are wondering, they are located in /usr/local/ or /usr/local/bin.

For a short descripton of file formats see [here](https://github.com/demharters/dtc-expTechniques-nanopore/blob/master/dataFormats.md).

To find out more about the options of a programme (poretools, samtools, lastdb etc.), in the terminal type:
```
[programme] --help
or
man [programme]
```

### The Exercise
Due to storage requirements you will be working in /tmp. Go to /tmp and create a new directory. In the terminal do:
```
cd /tmp
mkdir nanoporeData
cd nanoporeData
```

## Explore your data
Poretools can help you characterise your nanopore reads.
For each barcoded sample you should create a new directory in /tmp/nanoporeData. Try to mirror the folder structure in /datasets. I.e. for barcode BC03 in run1 you would do:
```
cd /tmp/nanoporeData
mkdir run1
cd run1
mkdir BC03
cd BC03
```

If you want to exit a directory, do:
```
cd ..
```

### Note
**In the following commands, "reads_folder" needs to be replaced with your respective read folder path, for example with "/datasets/run1/BC03".**

Plot the total yield of your run:
```
poretools yield_plot \
          --plot-type reads \
          --saveas reads_yield.pdf\
          reads_folder
```
*To plot the total basepair count set --plot-type to "basepairs"*

Use software such as libreOffice to look at the output. Provided you are in the right folder, in the terminal type:
```
libreOffice reads_yield.pdf
```

Plot the distribution of read lengths:
```
poretools hist reads_folder --saveas readLength_hist1.pdf
```

You can modify four plot by adding the following options:
```
poretools hist --min-length 1000 --max-length 10000 reads_folder --saveas readLength_hist2.pdf
poretools hist --num-bins 20 --max-length 10000 read_folder --saveas readLength_hist3.pdf
```

Look at the throughput of each of the pores:
```
poretools occupancy reads_folder --saveas occupancy.pdf
```

For more commands see the [poretools documentation](https://poretools.readthedocs.io/en/latest/content/examples.html).


## Alignment

##### Reference Sequences
An alignment requires one ore more reference sequences. We provided some relevant references for you under /datasets/refSequences.

For **plasmid alignments** use 'allRef.fasta'. This file contains reference sequences for genes gacS, GFP, penicillin amidase and a gentamycin resistance gene and the plasmid backbones pUC18 and pME6010. The alignment software will attempt to align your reads against all those sequences. At the end of this exercise you should know which of these genes and plasmids are present and whether there are any insertions, deletions or other mutations.

For the **genomic samples** use the file 'p_fluorescens.fasta'.

Before we continue with the alignments, copy the provided directory of reference files under /datasets into /tmp/nanoporeData:
```
cd /tmp/nanoporeData
cp -r /datasets/refSequences .
```

##### Step 1: Extract reads as FASTA
The nanopore basecaller (Epi2Me) returns each of the basecalled reads as individual fast5 files. These are stored in binary format. We need to use **poretools** to extract the reads from the fast5 folder and store them in a single fasta file.

Go to a directory such as /tmp/nanoporeData/BC03 and do:
```
poretools fasta /datasets/run1/BC03 > reads.fasta
```

Note, that we could also have used the fastq format which includes the quality scores and may help with the alignment step. However, for simplicity we will use the fasta format.


##### Script
The following commands can be run with [this script](https://github.com/demharters/dtc-expTechniques-nanopore/blob/master/alignmentScript.sh). Replace the values for "barcode" and "reference" accordingly.

*I suggest you go through the commands one by one first.*

```
#/bin/bash

#barcode="BC03"
#reference="allRef"
barcode="BC06"
reference="p_fluorescens"
lastdb -Q 0 ../../refSequences/$reference".lastindex" ../../refSequences/$reference".fasta"
lastal -s 2 -T 0 -Q 0 -a 1 ../../refSequences/$reference".lastindex" $barcode"_reads.fasta" > $barcode"_reads_aligned.maf"
maf-convert sam $barcode"_reads_aligned.maf" > $barcode"_reads_aligned.sam"
samtools faidx ../../refSequences/$reference".fasta"
samtools view -b -S -t ../../refSequences/$reference".fasta.fai" -o $barcode"_reads_aligned.bam" $barcode"_reads_aligned.sam"
samtools sort $barcode"_reads_aligned.bam" $barcode"_reads_aligned.sorted"
samtools index $barcode"_reads_aligned.sorted.bam"
```

Place the script under /tmp/nanoporeData/ and if located in /tmp/nanoporeData/run1/BC03 run it so:

```
../../alignmentScript.sh
```

You may have to give it execution permissions first, like so:

```
chmod +x ../../alignmentScript.sh
```

##### Step 2: Index reference file
Alignment tools like to index their files to improve computational efficiency. Provided you are working in a directory such as /tmp/nanoporeData/run1/BC03, and the reference file is named 'reference.fasta', the command would look like this:

```
lastdb -Q 0 ../../refSequences/reference.lastindex ../../refSequences/reference.fasta
```
```
*HELP*
-Q	input format: 0=fasta, 1=fastq
```

This will generate a set of files in /tmp/nanoporeData/refSequences, with different extensions (.ssp, .tis, .sds, .des, .prj, .suf, .bck). These will be used by the alignment software in subsequent steps.

##### Step 3: Alignment
Align the extracted reads to the reference sequence with the following command:

```
lastal -s 2 -T 0 -Q 0 -a 1 ../../refSequences/reference.lastindex reads.fasta > reads_aligned.maf
```
```
*HELP*
-s	0=reverse, 1=forward, 2=both
-T	type of alignment: 0=local, 1=overlap
-Q	input format: 0=fasta, 1=fastq
-a	gap penalty
```

##### Step 4: Generate Genome Viewer friendly alignment
Convert your alignment to the .sam format with **maf-convert.py**:

```
maf-convert sam reads_aligned.maf > reads_aligned.sam
```

Index your reference file (samtools has its own referencing scheme ..):

```
samtools faidx ../../refSequences/reference.fasta
```

Compress .sam to .bam:

```
samtools view -b -S -t ../../refSequences/reference.fasta.fai -o reads_aligned.bam reads_aligned.sam
```
```
*HELP*
-b	output BAM
-S	input SAM
-t	reference index file
-o	output file name
```

Sort your alignment by genome location to allow for pile-up:

```
samtools sort reads_aligned.bam reads_aligned.sorted
```

Index the sorted alignment (required to view the alignment in the Tablet alignment viewer)
```
samtools index reads_aligned.sorted.bam
```

##### Step 5: Alignment Visualisation

- Open the Tablet alignment viewer. In the terminal type:
```
tablet
```

- Go to >Open Assembly.
- Select your sorted alignment file.
- Select your reference.
- In the left column under contigs select your alignment.
- You should now be able to see your assembly.
- To see the coverage, go to > Advanced > Coverage.

![Tablet alignment](https://github.com/demharters/assemblyTutorial/blob/master/figures/tablet.png)
**Figure 1** An alignment of Nanopore reads of E.coli as viewed in Tablet.


### Classify plasmid variants
For this part use the data provided in /datasets/backup. Note, that we were using an older nanopore chemistry and basecaller. To ensure high accuracy we performed 2D sequencing. So when you extract your fasta sequences with poretools do something like:

```
poretools fasta --type 2D /datasets/backup/BC04 > BC04_reads.fasta
```

The 5 barcoded samples BC04-BC8 correspond to the following plasmids:

- P1. pGFP GFPmut3.1 in pME6010 – TetR
- P2. pME3258 GacS from P. fluor CHAO, in pME6010 - TetR
- P3. pME3281 GacS wt in miniTn7/pUC18 - GmR
- P4. pME3283 GacS H294R in miniTn7/pUC18 - GmR
- P5. pME6191 GacS Δ76 in miniTn7/pUC18 - GmR

Can you figure out which barcode matches which plasmid? Apply what you have learnt. You may want to add more sequences to your reference file.

[This paper](http://apsjournals.apsnet.org/doi/abs/10.1094/MPMI.2003.16.7.634) will help you identify the different gac variants:

* Zuber, Sophie, et al. "GacS sensor domains pertinent to the regulation of exoproduct formation and to the biocontrol potential of Pseudomonas fluorescens CHA0." Molecular plant-microbe interactions 16.7 (2003): 634-644.

Feel free to collaborate.

