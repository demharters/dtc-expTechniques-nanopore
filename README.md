# dtc-expTechniques-nanopore

## The Exercise
[In this exercise](https://github.com/demharters/dtc-expTechniques-nanopore/blob/master/exercise.md) you will learn about nanopore read characterisation, genome alignment and alignment analysis.

### Linux
This exercise will require Linux (Ubuntu). 

**Useful commands:**

|Command|Action|Example|
|---|---|---|
|mkdir|make directory|mkdir folderName|
|cd|change directory<br/>go up one directory|cd folderName<br/>cd ..|
|pwd|print working directory|pwd|
|rm|remove file<br/>remove folder|rm fileName<br/>rm -r folderName|
|mv|move or rename file or folder|mv oldFolderName newFolderName|
|cp|copy file|cp file fileCopy|
|cp -r|copy folder|cp -r folder folderCopy|
|ls|list folder contents<br/>list folder contents and show most recent file at bottom|ls folderName<br/>ls -ltr|

*When you type a command in the terminal, try pressing "TAB" while you're writing to see if software or file name exists (tab completion).*

* Popular tools for opening text files in Linux are **vim** (Terminal), **less** (Terminal), **nano** (Terminal), **gedit** (GUI) and **libreoffice** (GUI).
* Tools for opening pdfs and images in Ubuntu include among others **evince**, **okular** and **libreoffice**.


### Sequencing Data
Each pair of you was assigned one DNA sample and one of 9 barcodes. As there were 18 pairs in total, each DNA sample and barcode was shared between two pairs. After you gave us your barcoded libraries on Tuesday we quantified them. Using those concentrations we split your samples up into two sets of 9 samples, pooled them and sequenced them on two independent Oxford Nanopore flow cells. The respective datasets are located under **/datasets/run1** and **/datasets/run2**. Unfortunately, not all samples were sucessfully sequenced. This could be due to a variety of reasons.

The data should already be on your machines. You can download the same data from [here](https://figshare.com/s/1ef7923d6717ba8172ec).

For groups with no reads or a low number of reads, look at the samples of other pairs. The following table shows where the sequencing data can be found.

|Group|Barcode|Run|Type|
|---|---|---|---|
|5|BC03|1|Plasmid|
|9|BC05|1|Plasmid|
|10|BC05|2|Plasmid|
|11|BC06|2|Genomic|
|12|BC06|1|Genomic|
|13|BC07|2|Genomic|
|14|BC07|1|Genomic|
|17|BC09|1|Genomic|

Follow the steps in the [exercise](https://github.com/demharters/dtc-expTechniques-nanopore/blob/master/exercise.md) to characterise the sequencing data, perform alignments and identify the barcoded samples.
