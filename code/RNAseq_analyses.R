#### RNAseq analyses
### load packages + install ones that aren't  already installed in the libraries

packs_required<-c("BiocManager","pacman")
if_insta_packs<-packs_required %in% rownames(installed.packages())
if (any(if_insta_packs == FALSE)) {
  install.packages(packs_required[!if_insta_packs])
}

# load CRAN and Bioconductor packages
pacman::p_load(DESeq2,GenomicAlignments,GenomicRanges,Rsamtools,GenomicFeatures,ReportingTools,Gviz,RColorBrewer,ggplot2,tidyverse,ggpubr,ggstance,ggrepel)

# load bam and bai files
PATH<-"path/to/bam_and_bai_files"

# load bam and bai files
bai_files = list.files(path =PATH ,pattern="*.bai")
bam_files<-gsub(".bai","",bai_files) # define bam files based on the bai file names

# generate count table

Bamfiles_m<-NULL
# set a progress bar
pb = txtProgressBar(min = 0, max = length(bam_files), initial = 0)
stepi=1

for (f in bam_files) {
  temp_n<-gsub("_bowtie2.bam","",f)
  tempf<-BamFile(f, index=paste0(f,".bai"), yieldSize=NA_integer_, obeyQname=FALSE,
                 asMates=FALSE, qnamePrefixEnd=NA, qnameSuffixStart=NA)
  setTxtProgressBar(pb,stepi)
  stepi<-stepi+1
  Bamfiles_m[[temp_n]]<-tempf
}
# quickBamFlagSummary(tempf)

## We specify the files using the BamFileList() function (this groups multiple bamfiles)
# In this instance we are specifying both conditions (white clover, WC) and (Broth, B)

Bamfiles_l<- BamFileList(Bamfiles_m)

Bamfiles_l

# modify pathway and load reference gff files and features
gff_file1 <- makeTxDbFromGFF("path/to/ref.gff", format = "gff3")
feature_table<-read.delim("path/to/ref_feature_table.csv", sep = ",")

# Make a list of genes, the following line produces a GRangesList of all transcripts grouped by gene

Glist <- genes(gff_file1)
Glist

# Perform overlap queries between reads and genomic features
seq_count_WT_rescreen <- summarizeOverlaps(features=Glist, reads=Bamfiles_l,  mode = "Union")

# sanity check
colSums(assay(seq_count_WT_rescreen))

#Check read counts to the first 100 genes
head(assay(seq_count_WT_rescreen), 100)
seq_count_WT_rescreen_copy<-as.data.frame(assay(seq_count_WT_rescreen))
seq_count_WT_rescreen_copy$geneid<-rownames(seq_count_WT_rescreen_copy)

# add Strain metadata to the counts table, please modify the into = c("names") to suit dataset
metaRNA<-data.frame(SeqID=colnames(seq_count_WT_rescreen)) %>% separate(col = SeqID, sep = "_", into = c("a","strain","replicate","RNAseqID"), remove = FALSE ) %>% select(-a) 

rownames(metaRNA)<-metaRNA$SeqID
metaRNA$strain<-factor(metaRNA$strain)

# set Wt as the control # change the strain name to your own reference where applicable
metaRNA$strain<- relevel(metaRNA$strain, ref = "Wt") 

colData(seq_count_WT_rescreen)<- DataFrame(metaRNA)
colData(seq_count_WT_rescreen)$strain

# Calculate differential expression using the normalised dataset, using strain as the comparisons
dds_seq_count_WT_rescreen<- DESeqDataSet(seq_count_WT_rescreen,design = ~ strain,)


# Calculate the size factor for counts normalization

dds_seq_count_WT_rescreen_normalised <- estimateSizeFactors(dds_seq_count_WT_rescreen)

sizeFactors(dds_seq_count_WT_rescreen_normalised)

# DESeq2 doesn’t actually use normalized counts, rather it uses the raw counts and models the normalization inside the Generalized Linear Model (GLM). These normalized counts will be useful for downstream visualization of results, but cannot be used as input to DESeq2 or any other tools that peform differential expression analysis which use the negative binomial model.

# DESeq2 performs an internal normalization (i.e. median of ratios method) where geometric mean is calculated for each gene across all samples. The counts for a gene in each sample is then divided by this mean. The median of these ratios in a sample is the size factor for that sample. This procedure corrects for library size and RNA composition bias, which can arise for example when only a small number of genes are very highly expressed in one experiment condition but not in the other.

# convert to normalized count table 
Norm_rnaseq_table<-as.data.frame(counts(dds_seq_count_WT_rescreen_normalised, normalize = TRUE))


# Futher modification or reshape where needed then export table as needed


