#!/bin/bash
#Enter the fastq file pathways and output directory (make sure these directories exist and are empty!)
source ~/.bashrc # load $PATH to all the software/tools required in the pipeline

# Define input, output, and reference direction. Modify according to your own pathway
Input_DIR=/full/path/to/Illumina_RNA_Reads

Output_DIR=/full/path/to/RNAseq_output

ref_DIR=/full/path/to/ref_genome



cd $Output_DIR
mkdir clean_seq
mkdir bowtie2_index



# Adapter-trimming.  According to FastQC, minimal trimming should happen as almost no adapter sequences were detected

cd $Input_DIR
for fq in *_R1_001.fastq.gz
do bbduk.sh in1=$fq in2=${fq%_R1_001.fastq.gz}_R2_001.fastq.gz out1=$Output_DIR/clean_seq/trim_${fq%_R1_001.fastq.gz}_R1.fq out2=$Output_DIR/clean_seq/trim_${fq%_R1_001.fastq.gz}_R2.fq ref=/Users/xinyuewang/software/bbmap/resources/adapters.fa ktrim=N k=23 mink=11 hdist=1 tpe tbo

done

# Quality trimming
$Output_DIR/clean_seq
for fq in trim_*.fastq.gz
do bbduk.sh -Xmx1g in=$fq out=clean_${fq} qtrim=rl trimq=10
done

# remove phix
cd $Output_DIR/clean_seq
for fq in clean_*_R1_001.fastq.gz
bbduk.sh in1=$fq in2=${fq%_R1_001.fastq.gz}_R2_001.fastq.gz 
out1=$Output_DIR/clean_seq/filtered_${fq%_R1_001.fastq.gz}_R1.fq out2=$Output_DIR/clean_seq/filtered_${fq%_R1_001.fastq.gz}_R2.fq 
ref=phix.fa k=31 hdist=1 stats=stats.txt
done

# use bowtie2 to build reference file 

bowtie2-build -f $ref_DIR/WT_6206_genome.fasta $Output_DIR/bowtie2_index/mtb 

# align the QCed RNAseq reads to the mtb genome and create .bam files

cd $Output_DIR/clean_seq
mkdir bam_sam_wt6206

for i in clean_*_R1.fq
do bowtie2 -x $Output_DIR/bowtie2_index/mtb -1 $i -2 ${i%_R1.fq}_R2.fq --un $Output_DIR/clean_seq/bam_sam_wt6206/${i%_R1.fq}_unaligned.fq -p 4 -S $Output_DIR/clean_seq/bam_sam_wt6206/${i%_R1.fq}.sam

done


# creat bam and bai files
cd $Output_DIR/clean_seq/bam_sam_wt6206
for s in *.sam
do samtools sort -O bam -@4 $s -o ${s%.sam}_bowtie2.bam 
done

for b in *.bam
do samtools index $b  ${b}.bai
done

# process the downstream in R


