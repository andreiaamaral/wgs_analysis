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

### SAMTOOLS

1. Genome indexing

samtools faidx GCF_002310715.1_ASM231071v1_genomic.fna

2. convert SAM to BAM
samtools view -S -b E2-aln_1223.sam.sam > E2-aln_1223.bam
     
     
