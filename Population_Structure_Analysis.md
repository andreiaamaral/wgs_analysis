
# 3. Population Structure Analysis

![pop_structure_analysis](https://user-images.githubusercontent.com/64896369/123308599-c623cb80-d51b-11eb-8bec-9912ee11b989.png)


## 3.1. Principal Component Analysis


![pop_structure_analysis-Page-2](https://user-images.githubusercontent.com/64896369/123317329-caed7d00-d525-11eb-84fb-01526497a94d.png)


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

#### Pruning redundant SNP data

SNP pruning is an important step in several population structure analysis, where multiple SNPs are removed from the dataset, as it is very common to find extensive loci redundancy, considering many pairs of SNPs may display high levels of linkage disequilibrium (LD).
```
plink --vcf output_file_autosomes.vcf [--double-id] [--allow-extra-chr] [--set-missing-var-ids @:#] [--indep-pairwise <window size> <step size> <r^2 threshold>] [--out autosomes_pruned_data]
```

#### Principal Component Analysis

Based on the unpruned VCF file, __PLINK__ can extract the list of SNPs that remained after pruning and build two output files: a *.eigenvec* and *.eigenval*.

```
plink --vcf output_file_autosomes.vcf [--double-id] [--allow-extra-chr] [--set-missing-var-ids @:#] [--extract autosomes_pruned_data.prune.in] [--make-bed[ [--pca] [--out  autosomes_pca]
```

### R-code for detecting significant PCs and plotting

#### Identifying significant PCs

A Tracy-Widom test can be applied to the outputed .eigenval vector containing the principal components' eigenvalues. Here we determine the significant PCs considering a 99% confidence level and analysing a total of 20 generated PCs, using *tw()* function from the __AssocTests__ package.

```
library(AssocTests)

eigenval <- scan("autosomes_pca.eigenval")
tw(eigenvalues = eigenval, eigenL = 20, criticalpoint = 2.0234)
```

#### Plotting the PCA


##### Loading *.eigenvec* and *.eigenval* data structures
```
library(tidyverse)


pca <- read_table2("./autosomes_pca.eigenvec", col_names = FALSE)
eigenval <- scan("./autosomes_pca.eigenval")
```

##### Reshaping the principal components matrix 

```
pca <- pca[,-1]
names(pca)[1] <- "samples"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
```

#### Introducing population information 

It is very important to group samples by their population, and introduce information about which samples belong to which population before plotting. Here we introduced a *.txt* file with each sample population's name, one per line, matching the order of the samples in the *.eigenvec* file.

```
pops <- read.table("samples_populations", header = F)
pops <- as.vector(pops)
pca <- as.tibble(data.frame(pca, pops))
```

#### Plotting PC's explained variance

Plotting the explained variance of each PC is very important so that we may have an ideia regarding how well can the PCA represent our data.

```
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
var_plot <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity") + ylab("Percentage variance explained") + theme_light()
cumsum(pve$pve)
```

### Plotting our data
```
pca_plot <- ggplot(pca, aes(PC1, PC2, col = V1)) + geom_point(size = 3) 
            +  coord_equal() + theme_light() +  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) 
            + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

```
## 3.2. Phylogenetic Analysis

![workflow_phylogenetic_tree_final_heritage (4)](https://user-images.githubusercontent.com/64896369/123458953-8c1afe00-d5dd-11eb-9b53-ac56602e99fc.png)


### 3.2.1. Phylogenetic Analysis - Maternal Heritage

#### Indexing VCF files

```
tabix -p vcf output_file.vcf.gz
```

#### Selecting Mitochondrial SNPs

```
tabix output_file.vcf.gz [MT chr] [-h keep header] > output_file_chrMT.vcf
```

#### Compressing Mitochondrial VCF to a BGZip format

```
bgzip output_file_chrMT.vcf
```

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
bcftools merge [-Oz compressed VCF] [-o output] <multi_sample_chrMT.vcf.gz> <input1> ... <inputN>
```

#### Creating Multi-sample MUSCLE alignment


![pop_structure_analysis-Page-4](https://user-images.githubusercontent.com/64896369/123457770-2da15000-d5dc-11eb-9f42-729c5fc29368.png)


VCF-kit is a collection of utilities that allows multiple sequence alignment (using MUSCLE algorithm) of multiple samples' SNP data. This alignment can further be used to plot a phylogenetic analysis, using maximum likelyhood method in MEGA.

```
vk phylo fasta <multi_sample_chrMT.vcf.gz> > Groenen_Frantz_Ganda_alignment.fasta
```

VCF-kit is also able to plot a phylogenetic analysis using the neighbor joining method without requiring a previous alignment. The newick file format can then be loaded into a online tool (https://itol.embl.de/upload.cgi).

```
vk phylo tree nj <multi_sample_chrMT.vcf.gz> > multi_sample_chrMT.newick
```


## 3.3. Pairwise Fst

![pop_structure_analysis-Page-3](https://user-images.githubusercontent.com/64896369/123321533-fe7ed600-d52a-11eb-8d1c-88ccbeae09e1.png)


#### Renaming sample VCFs (If needed)
Multiple times, VCF files refer to their sample name as only the number "20". Further use of multiple samples in the same VCF file requires the sample name to be changed in order to preserve the identity of the sample. Furthermore, BCFTools will not be able to merge multiple samples with the same name.

You can use the following command to check the VCF's sample name:

```
bcftools view -h input_file.vcf.gz | tail -n1
```


The following command uses BCFTools *reheader* feature to change the sample name to the corresponding file name. Furthermore, it makes use of tabix to index the renamed file.

```
for file in $(ls *filtered.vcf.gz); do bcftools reheader -s <(echo ${file/.vcf.gz}) -o ${file/.vcf.gz/.renamed.vcf.gz} $file; tabix -fp vcf ${file/.vcf.gz/.renamed.vcf.gz}; done
```


#### Creating multi-sample VCF file

In this multi-sample VCF, samples from the two populations we want to compare must be present.

```
bcftools merge [-Oz compressed VCF output] [-o output] output_file.vcf.gz Input1.vcf.gz ... InputN.vcf.gz
```

#### Indexing and BGZiping multi-sample VCF file

```
tabix -p vcf output_file.vcf.gz
```

#### Calculating Fst 

VCFTools requires, not only the merged multi-sample VCF but two files, *pop1* and *pop2*, with the names of the samples belonging to each of the populations in study. 

```
vcftools --gzvcf output_file_autosomes.vcf [--fst-window-size] [--weir-fst-pop pop1] [--weir-fst-pop pop2] --out pop1_vs_pop2
```

### R code for calculating pairwise mean Fst

#### Loading VCFTools output file
```
fst.file <- read.table("pop1_vs_pop2.windowed.weir.fst", header = T)
```

#### Selecting autosomal windows

For the purpose of our experiment, only autosomal SNP information is considered for analysis. Therefore, Fst windows from sexual chromosomes, Mitochondria and SNPs present in unplaced scaffolds should be removed. 
```
fst.autosomes <- subset(fst.file, CHROM %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18"))
```

**Fst values have biological significance between 0 and 1. Eventhough, tools like VCFTools may calculate Fst values between -1 and 1. There are different ways to deal with Fst negative values.**

#### Changing Negative values to 0
```
fst.autosomes <- mutate(fst.autosomes, MEAN_FST = if_else(MEAN_FST < 0, 0, MEAN_FST))
```
or

#### Considering only Non-Negative values
```
fst.autosomes <- subset(fst.autosomes, CHROM >= 0)
```

#### Calculating Mean Fst
```
mean(fst.autosomes$MEAN_FST)
```

#### Creating heatmap with Fst mean values

Input should be loaded in this format:
```
X           Y           Fst
pop1        pop2        value1
pop1        pop3        value2
pop2        pop3        value3
```

#### R code 

```
library(tidyverse)

pairwise <- read.table("pairwise_mean_fst.csv", header = T, sep = ";")

ggplot(pairwise, aes(x=X, y=Y, fill= Fst)) + geom_tile() + geom_text(aes(label = round(Fst, 2)), size = 5) + scale_fill_gradient(low = "red", high = "#07D72A") + theme_classic() + xlab("") + ylab("") + theme(axis.text = element_text(size = 14, face = "bold"))
```
