
# 4. Population Fitness Analysis

![workflow_phylogenetic_tree_final_heritage-Page-2](https://user-images.githubusercontent.com/64896369/123460606-c9808b00-d5df-11eb-88d8-029fa2242af6.png)

## 4.1. Differentiation and Introgression

### Calculating Fst and θπ density 

![pop_structure_analysis-Page-5 (1)](https://user-images.githubusercontent.com/64896369/123659609-c9280000-d82a-11eb-811d-2d03c177ce88.png)


Fst (Fixation Index) can be used to measure differentiation between populations, due to genetic structure variations.

θπ (Nucleotide Diversity) can be used to measure the degree of variation within a population.


#### Renaming sample VCFs (If needed)
Multiple times, VCF files refer to their sample name as only the number "20". Further use of multiple samples in the same VCF file requires the sample name to be changed in order to preserve the identity of the sample. Furthermore, BCFTools will not be able to merge multiple samples with the same name.

You can use the following command to check the VCF's sample name:

```
bcftools view -h input_file.vcf.gz | tail -n1
```


The following *for* loop command uses BCFTools *reheader* feature to change the sample name to the corresponding file name. Furthermore, it makes use of tabix to index the renamed file.

```
for file in $(ls *filtered.vcf.gz); do bcftools reheader -s <(echo ${file/.vcf.gz}) -o ${file/.vcf.gz/.renamed.vcf.gz} $file; tabix -fp vcf ${file/.vcf.gz/.renamed.vcf.gz}; done
```


#### Creating multi-sample VCF file

To calculate nucleotide diversity within a population it is necessary to create a merged file of samples from that population.
To calculate the fixation index between two populations it is necessary to create a merged file with all the existing samples from both populations.

For calculating pop1 and pop2 θπ: 
```
bcftools merge [-Oz compressed VCF output] [-o output] pop1.vcf.gz Input1.vcf.gz ... InputN.vcf.gz

bcftools merge [-Oz compressed VCF output] [-o output] pop2.vcf.gz Input1.vcf.gz ... InputN.vcf.gz
```

For calculating pop1 vs. pop2 Fst:
```
bcftools merge [-Oz compressed VCF output] [-o output] pop1_vs_pop2.vcf.gz Input1.vcf.gz ... InputN.vcf.gz
```


#### Indexing multi-sample VCF file

```
tabix -p vcf <file>
```

#### Calculating θπ : 

```
vcftools --vcf popN_autosomes.vcf [--window-pi] --out popN_pi
```

#### Calculating Fst:

```
vcftools --vcf popN_autosomes.vcf [--fst-window-size] [--weir-fst-pop pop1] [--weir-fst-pop pop2] [--out pop1_vs_pop2_fst]
```

### R code for plotting θπ and Fst



#### Loading VCFTools output file
```
fst.file <- read.table("pop1_vs_pop2_fst.windowed.weir.fst", header = T, sep = ",")
pi.fileN <- read.table("popN_pi.windowed.weir.pi", header = T, sep = ",")
```

#### Selecting autosomal windows

For the purpose of our experiment, only autosomal SNP information is considered for analysis. Therefore, windows from sexual chromosomes, Mitochondria and SNPs present in unplaced scaffolds should be removed. 
```
fst.autosomes <- subset(fst.file, CHROM %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18"))

pi.autosomesN <- subset(pi.fileN, CHROM %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18"))
```

#### Calculating log2(θπpop1 / θπpop2)


Calculating log2(θπpop1 / θπpop2) is important in our analysis because it will allow us to identify genomic regions in  **pop2** that have low levels of genetic diversity when compared to **pop1**. These regions that display low levels of genetic diversity may be conserved regions, holding genes that can be associated with desirable features.

*Subsetting θπ values according to chromosome*

In order to generate a table with values for pop1 and pop2 θπ we need first to subset θπ values according to chromosome. This is necessary because the position at each chromosome is reset and merging all chromosomes at once would fail since R wouldn't be able to merge multiple bins with the same positions. 

Pop1:
```
pi.pop1.chr1 <- subset(pi.pop1, CHROM == "1")

pi.pop1.chrN <- subset(pi.pop1, CHROM == "N")
```

Pop2:
```
pi.pop2.chr1 <- subset(pi.pop2, CHROM == "1")

pi.pop2.chrN <- subset(pi.pop2, CHROM == "N")
```

*Merging tables according to chromosome*

```
pi.pop1.pop2.chr1 <- merge(pi.pop1.chr1, pi.pop2.chr1, by = "BIN_START", All=FALSE)

pi.pop1.pop2.chrN <- merge(pi.pop1.chrN, pi.pop2.chrN, by = "BIN_START", All=FALSE)
```

*Bind chromosome-wide tables into a single table*

```
pi.pop1.pop2 <- rbind(pi.pop1.pop2.chr1, ..., pi.pop1.pop2.chrN)
```

*Calculate log2(θπpop1 / θπpop2)*

```
pi.pop1.pop2$ratio <- log2(pi.pop1.pop2$PI.y/pi.pop1.pop2$PI.x)
```

*Subsetting θπ and Fst tables according to chromosome*

This processe needs to be repeated in order to join Fst and θπ ratio values in the same table.

θπ ratio:
```
pi.pop1.pop2.chr1 <- subset(pi.pop1, CHROM == "1")

pi.pop1.pop2.chrN <- subset(pi.pop1, CHROM == "N")
```

Fst:
```
fst.pop1.pop2.chr1 <- subset(fst.pop1.pop2, CHROM == "1")

fst.pop1.pop2.chrN <- subset(fst.pop1.pop2, CHROM == "N")
```

*Merging tables according to chromosome*

```
pi.fst.chr1 <- merge(pi.pop1.pop2.chr1, fst.pop1.pop2.chr1, by = "BIN_START", All=FALSE)

pi.fst.chrN <- merge(pi.pop1.pop2.chrN, fst.pop1.pop2.chrN, by = "BIN_START", All=FALSE)
```

*Bind chromosome-wide tables into a single table*

```
pi.fst <- rbind(pi.fst.chr1, ..., pi.fst.chrN)
```

#### Plotting θπ and Fst density plots

```
library(tidyverse)

pi.fst$Legend <- "pop1 vs. pop2"
```

```
fst.density <- ggplot(pi.fst, aes(MEAN_FST, fill=Legend) 
      + geom_histogram(alpha = 0.6, aes(y = ..density..), position = 'identity', bins = 200,color= "black", show.legend = F) 
      + labs(x = expression(("F"["st"])), y = "Density") 
      + xlim(c(-0.5, 1.0)) + theme_classic() 
      + theme(axis.line = element_line(colour = "black", size = 1.3, linetype = "solid"), text=element_text(size=20, face = "bold", color = "black"))
```

```
pi.density <- ggplot(pi.fst, aes(ratio, fill=Legend)) 
      + geom_histogram(alpha = 0.6, aes(y = ..density..), position = 'identity', bins = 200,color= "black", show.legend = F) 
      + labs(x = expression(log[2](theta[pi~"Pop1"]/theta[pi~"Pop2"])), y = "Density") 
      + xlim(c(-10, 10)) + theme_classic() 
      + theme(axis.line = element_line(colour = "black", size = 1.3, linetype = "solid"), text=element_text(size=20, face = "bold", color = "black"))
```

![pop_structure_analysis-Page-6](https://user-images.githubusercontent.com/64896369/123674218-42c6ea80-d839-11eb-9594-30b94ffb8218.png)


### Building Fst and θπ manhattan plots

Plotting Fst and θπ values according to position can be complicated since the position is reset in each chromosome. Therefore we can creat a new column in our dataframe (**bins**) according to the bin position along the genome, whithout being reset in each chromosome.
```
pi.fst$bins <- 1:nrow(pi.fst)
```

Furthermore, our dataframe can be used to group bins by each chromosome in order to facilitate chromosome identification.
```
pos <- pi.fst %>%
    + group_by(CHROM) %>% 
    + summarize(avg = round(mean(bins))) %>%
    + pull(avg)
```

Plotting θπ manhattan plot:

```
pi <- ggplot(pi.fst, aes(x=bins, y=ratio, color=factor(CHROM))) 
      + geom_point(alpha=0.75, show.legend = F) + theme_classic() 
      + scale_x_discrete(limits = pos, labels = unique(pi.fst$CHROM)) 
      + scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(pi.fst$CHROM)))) 
      + xlab("Chromosome") + ylab(expression(log[2](theta[pi~"Pop1"]/theta[pi~"Pop2"]))) 
      + theme(axis.text = element_text(size = 26, face = "bold"), axis.title=element_text(size=14,face="bold")) 
      + scale_y_continuous(position = "right")
```

Plotting Fst manhattan plot:

```
fst <- ggplot(pi.fst, aes(x=bins, y=MEAN_FST, color=factor(CHROM))) 
      + geom_point(alpha=0.75, show.legend = F) + theme_classic() 
      + scale_x_discrete(limits = pos, labels = unique(pi.fst$CHROM)) 
      + scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(pi.fst$CHROM)))) 
      + xlab("Chromosome") + ylab(expression("F"["st"])) 
      + theme(axis.text = element_text(size = 30, face = "bold"), axis.title=element_text(size=30,face="bold"))
```

### Building Fst and θπ Scatter plot

Determine the 95% quantile windows with the highest Fst and θπ ratio values:
```
fst.quantile <- quantile(pi.fst$MEAN_FST, 0.95)
pi.quantile <- quantile(pi.fst$ratio, 0.95)
```

Extracting 95% quantile windows:
```
pi.fst.95 <- subset(pi.fst, ratio > pi.quantile & MEAN_FST > fst.quantile)
```

```
scatter <- ggplot() + geom_point(data = pi.fst, aes(x=ratio, y=MEAN_FST)) 
      + geom_point(data = pi.fst, aes(x=ratio, y=MEAN_FST), color = "red") + theme(legend.position = "none") 
      + theme_bw() + ylab(expression("F"["st"])) + xlab(expression(log[2](theta[pi~"Pop1"]/theta[pi~"Pop2"]))) 
      + geom_vline(xintercept = pi.quantile, linetype = "dashed") 
      + geom_hline(yintercept = fst.quantile, linetype = "dashed")
```

### Testing Equality of means (Mann-Whitney U test):

It is important to understand if there is any significant difference between the two sets of data of Whole genome Fst and θπ values compared to the extract 95% outliet values. The Mann-Whitney U test is usefull in these situations.

Fst:
```
wilcox.test(pi.fst$MEAN_FST, pi.fst.95$MEAN_FST)
```

θπ:
```
wilcox.test(pi.fst$ratio, pi.fst.95$ratio)
```

Extract 95% quantile windows position (Chromosome, Window Start and Window End):

```
pi.fst.95.pos <- pi.fst.95[c("CHROM", "BIN_START", "BIN_END")]

write.csv(pi.fst.95.pos, "C:/path-to-folder/pi_fst_95.csv", row.names = FALSE)
```


### Identifying Genes in the Outlier Regions

```
library(biomaRt)
```

The following commands can be used to check the available reference genomes and their versions.
```
ensembl = useMart("ensembl")
listDatasets(ensembl, verbose = FALSE)

ensembl = useMart("ensembl", dataset = "ref_gen")
```

Extracting genes IDs:
```
genes_list <- cbind(
    pi.fst.95.pos,
    genes = apply(pi.fst.95.pos, 1, function(i){
    x <- getBM(attributes=c("external_gene_name"),
    filters = c("chromosome_name" , "start", "end"),
    values = list(i[1], i[2], i[3]),
    mart = ensembl)
# keeping only 3 genes, as output is too long.
x <- head(x, 3)
# return genes, comma separated
paste(x$external_gene_name, collapse = ",")
})
)
```

Renaming dataframe's columns:
```
genes_list_rename <- genes_list %>% rename(CHROM = V1, BIN_START = V2, BIN_END = V3, GENES = genes)
```

Extracting dataframe:
```
write.csv(genes_list_rename, "C:/path-to-folder/pop1_vs_pop2_95_outliers_genes.csv", row.names = FALSE)
```


### Calculating the Inbreeding Coefficient (FIS)

![pop_structure_analysis-Page-7](https://user-images.githubusercontent.com/64896369/123676094-6854f380-d83b-11eb-8d41-582deed5e576.png)

#### Renaming sample VCFs (If needed)
Multiple times, VCF files refer to their sample name as only the number "20". Further use of multiple samples in the same VCF file requires the sample name to be changed in order to preserve the identity of the sample. Furthermore, BCFTools will not be able to merge multiple samples with the same name.

You can use the following command to check the VCF's sample name:

```
bcftools view -h input_file.vcf.gz | tail -n1
```


The following *for* loop command uses BCFTools *reheader* feature to change the sample name to the corresponding file name. Furthermore, it makes use of tabix to index the renamed file.

```
for file in $(ls *filtered.vcf.gz); do bcftools reheader -s <(echo ${file/.vcf.gz}) -o ${file/.vcf.gz/.renamed.vcf.gz} $file; tabix -fp vcf ${file/.vcf.gz/.renamed.vcf.gz}; done
```


#### Creating multi-sample VCF file

```
bcftools merge [-Oz compressed VCF output] [-o output] output_file.vcf.gz Input1.vcf.gz ... InputN.vcf.gz
```

#### Indexing and BGZiping multi-sample VCF file

```
tabix -p vcf output_file.vcf.gz
```

#### Selecting autosomes for analysis

For the purpose of our experiment, only autosomal SNP information is considered for analysis. Therefore, SNP data from sexual chromosomes, Mitochondria and SNPs present in unplaced scaffolds should be removed. 

```
tabix -h output_file.vcf.gz [Chr1 ... ChrN] > output_file_autosomes.vcf
```

#### Compressing Autosomal VCF to a BGZip format

```
bgzip output_file_autosomes.vcf
```

#### Calculating Observed and Expected Heterozygosity

```
plink [--vcf output_file_autosomes.vcf.gz] [--double-id] [--hardy] [--allow-extra-chr] [--out popX]
```

PLINK command will calculate expected and observed heterozygosity for each *loci*. Furthermore, these values can be used to calculate mean values and use them to calculate the inbreeding coeficient.

```
cat popX.hwe | awk '{print $7}' | grep -v 'O(HET)' | awk '{sum+=$1} END {print "Mean O(HET) = ",sum/NR}' 
```

```
cat popX.hwe | awk '{print $8}' | grep -v 'E(HET)' | awk '{sum+=$1} END { print "Mean E(HET) = ",sum/NR}' 
```

#### Calculating FIS

Inbreeding coeficient can be calculated using the following command:

FIS=(HET(E)-HET(O))/(HET(E))
