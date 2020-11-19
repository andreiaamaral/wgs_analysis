# wgs_analysis

## Participants in the hackathon (14): 
                                
                                Andreia Amaral

                               Juliana Menezes
                               Pedro Sá
                               
                               Lucio Marcello 
                               
                               Leila Mansouri-CRG
                               
                               Suzanne Jin|BovReg|CRG
                               
                               Meenu Bhati | BovReg | ETH
                               
                               Kylie Munyard Curtin Uni/ Fr-AgENCODE
                               
                               Cervin Guyomar - FAANG - INRAE Toulouse
                               
                               Zeinab Manzari 
                               
                               Phanindra-SLU-FAANG
                               
                               Philippe BARDOU
                               
                               Alexandre Gilardet 
                               
                               Jose Espinosa
                               
                               
     
     
     
## WORKFLOW                         
                               
          
          
 ![](https://github.com/andreiaamaral/wgs_analysis/blob/main/Slide1.jpg)
                                         
                               
     
## Tools for the workflow and manual instructions

### PICARD
1- create genome dictionary

java -jar /data/software/picard/build/libs/picard.jar CreateSequenceDictionary \ 
      R=GCF_002310715.1_ASM231071v1_genomic.fna \ 
      O=rGCF_002310715.1_ASM231071v1_genomic.dict 
      
2- java -jar /data/software/picard/build/libs/picard.jar MarkDuplicates \
I=E2_alnRG_sorted.bam \
O=E2_alnRGsorted_duplicates_removed.bam \
M=E2_alnRGsorted_duplicates_removed_dup_metrics.txt \
REMOVE_SEQUENCING_DUPLICATES=true

3- After having all individual VCFs merged them in a single VCF
   
   java -jar /data/software/picard/build/libs/picard.jar MergeVcfs \
          I=E1_output_filtered.vcf.gz\
          I=E2_output_filtered.vcf.gz \
          O=E12output_variants.vcf.gz

### SAMTOOLS

1. Genome indexing

samtools faidx GCF_002310715.1_ASM231071v1_genomic.fna

2. convert SAM to BAM
samtools view -S -b E2-aln_1223.sam.sam > E2-aln_1223.bam

### GATK

1- Variant calling

gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R GCF_002310715.1_ASM231071v1_genomic.fna \
   -I E2_alnRGsorted_duplicates_removed.bam \
   -O E2_output.vcf.gz \
   -bamout E2_bamout.bam
   
2. Variant filter
   
   gatk VariantFiltration \
   -R GCF_002310715.1_ASM231071v1_genomic.fna \
   -V E2_output.vcf.gz \
   -O E2_output_filtered.vcf.gz \
   --filter-name "my_filter1" \
   --filter-expression "QUAL < 0 || DP<10.0 || MQ < 30.00 || SOR > 10.000 || QD < 2.00 || QD> 5.00|| FS > 200.000 || ReadPosRankSum < -20.000 || ReadPosRankSum > 20.000"  
   
   




     
     
