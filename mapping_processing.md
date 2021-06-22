# 1. Mapping and processing WGS data

![image](https://user-images.githubusercontent.com/64896369/122912564-5d3c2800-d350-11eb-8d10-d90331ea9adc.png)

## 1.1. Splitting input file (if needed)
#### Separating the read's information

_Extracting reads' ID_
```
awk '1==NR%4' <accession.fastq> > accession_ID
```
_Extracting reads' sequence_
```
awk '2==NR%4' <accession.fastq> > accession_SEQ
```
_Extracting reads' quality scores_
```
awk '0==NR%4' <accession.fastq> > accession_QS
```

#### Merging the extracted information into a single table file
```
paste <accession_ID> <accession_SEQ> <accession_QS> > TABLE 
```
#### Selecting IDs with pairs
 The reads ID information is used to select reads that have pairs and compose a list of their IDs using a simple bash command.
```
cat <accession_ID> | awk '{print $2}' | cut -f1 -d# | sort | uniq -c | awk '$1==2 {print $2}'  > IDs_with_pairs
```

#### Writting paired reads into their respective forward/reverse files

 The table and list are used by an in-house perl script (*select_paired_reads_frantz.pl*) that works by identifying the IDs with pairs in the table, checking the read’s orientation and writing its respective information into the respective output file.
```
perl select_paired_reads_frantz.pl <TABLE>  <IDs_with_pairs> <R1.fastq> <R2.fastq>
```

## 1.2. FastQC: basic quality control analysis

 Assessing reads' quality, length distribution and adapter presence.
```
fastqc <R1.fastq> <R2.fastq>
```

## 1.3. Adapter Trimming (if present)

  Trimming adapter sequences if their presence was detected in FastQC analysis.
 ```
 flexbar -r <R1.fastq> -p <R2.fastq> [-aa or --adapter-preset] TruSeq [-ap or --adapter-pair-overlap] ON
 ```
 
 ## 1.4. Quality Filtering 
 
  Removing low quality reads using this Prinseq-lite command.
 ```
 prinseq-lite.pl [-fastq input_fastq_file] <R1.fastq> [-min_len n] [-max_len n] [-min_qual_score n]
 prinseq-lite.pl [-fastq input_fastq_file] <R2.fastq> [-min_len n] [-max_len n] [-min_qual_score n]
 ```
 
 ## 1.5. Reads Pairing
 
  #### Selecting IDs with pairs
  
  These bash commands were used to get the reads IDs from the filtered forward and reverse files and select the ones that still have pairs.
  
  _Extracting the forward and reverse IDs_
  ```
	awk '1==NR%4' <R1_filtered.fastq>  | awk -F '#' '{print $1}' > R1_IDs
	awk '1==NR%4' <R2_filtered.fastq>  | awk -F '#' '{print $1}' > R2_IDs
  ```
  _Selecting IDs with pairs_
  ```
	cat <R1_IDs> <R2_IDs>  | sort | uniq -c | awk '$1==2 {print $2}' > IDs_with_pairs
  ```
  #### Extracting reads' information

  _Extracting reads' ID_
  ```
  awk '1==NR%4' <R1_filtered.fastq> > R1_accession_ID
  awk '1==NR%4' <R2_filtered.fastq> > R2_accession_ID
  ```
  _Extracting reads' sequence_
  ```
  awk '2==NR%4' <R1_filtered.fastq> > R1_accession_SEQ
  awk '2==NR%4' <R2_filtered.fastq> > R2_accession_SEQ
  ```
  _Extracting reads' quality scores_
  ```
  awk '0==NR%4' <R1_filtered.fastq> > R1_accession_QS
  awk '0==NR%4' <R2_filtered.fastq> > R2_accession_QS
  ```
  #### Merging the extracted information into tables
  ```
  paste <R1_accession_ID> <R1_accession_SEQ> <R1_accession_QS> > R1_filtered_TABLE
  paste <R2_accession_ID> <R2_accession_SEQ> <R2_accession_QS> > R2_filtered_TABLE
  ```
  
  #### Writting paired reads into their respective forward/reverse files
  
   The tables and lists are used by a perl in-house script (*select_paired_good_reads_frantz.pl*), that selects paired reads from the tables based on the provided list of paired IDs.
	
  ```
  perl select_paired_good_reads_frantz.pl <IDs_with_pairs> <R1_filtered_TABLE> <R1_filtered_TABLE_paired.fastq>
  perl select_paired_good_reads_frantz.pl <IDs_with_pairs> <R2_filtered_TABLE> <R2_filtered_TABLE_paired.fastq>
  ```
  
  ## 1.6. Indexing reference genome
  
  Firstly, the reference genome must be used to create the BWT indexes.
  ```
	bwa index <genome.fa>
  ```
  
  ## 1.7.  Mapping to Reference Genome
   The BWA mem command uses seeding alignments and considers the maximal number of exact matches (MEMs). The seeds are then extended using the Smith-Waterman algorithm to produce an alignment in SAM file format.
	
  ```
  bwa mem <indexed_genome.fa> <R1_filtered_TABLE_paired.fastq> <R2_filtered_TABLE_paired.fastq> > mapped.sam
  ```
  
  ## 1.8. Add Read Groups
  Read Groups are added using Picard’s “AddOrReplaceReadGroups” feature. Since this feature requires a .bam file as an input it is firstly required to convert the alignment’s output file from .sam to .bam.
	
  ```
  samtools view [-Sb .sam to .bam] [-o output] alignment.bam alignment.sam
  ```
  
  ```
	java -jar picard.jar AddOrReplaceReadGroups [-I --input] alignment.bam [-O --output] alignmentRG.bam [-RGID] N [-RGLB] libN [-RGPL] ILLUMINA [-RGPU] unitN [-RGSM] 20
  ```
  
  ## 1.9. Sorting the BAM file
  
  Sorting the alignment file is then required before it is possible to remove PCR duplicates by using Picard’s “MarkDuplicates” feature.
	
  ```
  samtools sort <alignmentRG.bam> <alignmentRG_sorted.bam>
  ```
  
  ## 1.10. Mark and Remove Duplicated 
  
  ```
	java -jar picard.jar MarkDuplicates [-I –input] alignmentRG_sorted.bam [-O --output] alignmentRG_sorted_duplicates_removed.bam [-M --metrics-file] alignmentRG_sorted_removed_dup_metrics.txt [--remove-sequencing-duplicates] true
  ```

  ## 1.11. Collecting Mapping Statistics
  
  ```
  samtools flagstats alignmentRG_sorted_duplicates_removed.bam
  ```
  
  ## 1.12. Merging Same-sample Runs
  
  Occasionaly, samples can be deposited using multiple runs. Therefore, for these samples a final merging step is required.
	
  ```
  samtools merge <merged_sample.bam> <in1.bam> … <inN.bam>
  ```
