#/bin/bash

barcode="BC03"
reference="allRef"
lastdb -Q 0 ../../refSequences/$reference".lastindex" ../../refSequences/$reference".fasta"
lastal -s 2 -T 0 -Q 0 -a 1 ../../refSequences/$reference".lastindex" $barcode"_reads.fasta" > $barcode"_reads_aligned.maf"
maf-convert sam $barcode"_reads_aligned.maf" > $barcode"_reads_aligned.sam"
samtools faidx ../../refSequences/$reference".fasta"
samtools view -b -S -t ../../refSequences/$reference".fasta.fai" -o $barcode"_reads_aligned.bam" $barcode"_reads_aligned.sam"
samtools sort $barcode"_reads_aligned.bam" $barcode"_reads_aligned.sorted"
samtools index $barcode"_reads_aligned.sorted.bam"
