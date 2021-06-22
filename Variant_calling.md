
# 2. Variant Calling

![image](https://user-images.githubusercontent.com/64896369/122931392-c2017d80-d364-11eb-8bd4-df5c0b053ecd.png)


## 2.1. DNAnexus' SAMtools Variant Caller app

Variant calling was performed in a cloud-based platform for data analysis, DNA nexus.
Samples can be processed using DNAnexus’ application “SAMtools Variant Caller” v1.0.6. This application allows variant calling in “batch” (several samples can be called at the same time). The application is responsible for indexing the .bam and .fa files, building the mpileup data structure, and call variants.

_Indexing mapped sample_

```
samtools index <mapped_sample.bam>
```

_Indexing Reference Genome_

```
samtools faidx <genome.fa>
```

_Building the mpileup data structure_

```
samtools mpileup -f <genome.fa> <mapped_sample.bam> -u[uncompressed format] [-C capQcoef] [-q minMapQ] [-Q minBaseQ] [-D per-sample read depth] [-S strand bias]
```

_Calling variants_

```
bcftools view [-vg gziped VCF] 
```

## 2.2. Removing INDELs

Variants regarding insertions and deletions should be removed when working with only single nucleotide variations (SNVs).

```
 bcftools view -i 'TYPE="snp"' <variants_sample.vcf.gz> > <variants_sample_indel_rmv.vcf>
```

## 2.2. SNP filtering

SNPs with low coverage and representation should be removed in order to remove the false positives frequency.

```
bcftools view <variants_sample_indel_rmv.vcf> | vcfutils.pl varFilter [-a minimum number of alternate bases] [-d minimum read depth] > <variants_sample_indel_rmv_filtered.vcf>
```

