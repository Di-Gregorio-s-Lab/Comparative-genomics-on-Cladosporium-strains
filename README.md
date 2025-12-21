# Comparative genomics on *Cladosporium* strains
This pipeline contains the codes used in the paper "Comparative genomics of the genus *Cladosporium* seemingly identified lack of selective pressure in soils historically polluted by hexachlorocyclohexane and polychlorobiphenyls".

## Data accessibility
Raw data and outputs are pubblicly available.
Raw genomic (g)DNA sequences and assembled genomes are pubblicly available under the NCBI BioProject [PRJNA1348484](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1348484).\
Raw sequences were obtained via PacBio HiFi circular consensus sequencing.\
Quality and completeness of the assemblies are reported in the supplementary materials of the paper (scrivere supplementary material come link).

## Pipeline overview
The following steps were used in this work:
- [genome assembly with HiFiasm v0.16.0 and Minimap2 v2.30](#Genome-assembly-with-HiFiasm)
- [genome quality check with Quast v5.3 and BUSCO v6.0.0](#Genome-quality-check-with-Quast-and-BUSCO)
- [taxonomic identification with local and online BLAST v2.17.0](#Taxonomic-identification-with-BLAST)
- [UPGMA tree with hclust and ggtree v3.12.0 on R v4.4](#UPGMA-tree-with-hclust-and-ggtree)
- [gene prediction with Braker3 v3.0.8](#Gene-prediction-with-Braker3)
- [functional annotation with Diamond v2.0.15 against the Swissprot database](#Functional-annotation-with-Diamond-+-GO-and-CAZyme-annotation)
- [Gene Ontology annotation with the Uniprot Retrive/ID mapping tool](#Functional-annotation-with-Diamond-+-GO-and-CAZyme-annotation)
- [CAZyme annotation with dbcan3 v5.1.2](#Functional-annotation-with-Diamond-+-GO-and-CAZyme-annotation)
- [PCoA and heatmaps with ggplot2 v3.5.2 and pheatmap v1.0.13 on R](#PCoA-and-heatmaps-with-ggplot2-and-pheatmap)

## Bioinformatic pipelines in Bash
To better navigate in the multiple analyses performed in thiese pipelines, a folders architecture was created as follows:
- "main" folder, containing all the sub-folders for this analysis
- "work" folder, containing raw and processed data for each pipeline
- "results" folder, containing the final output of each pipeline
- "scripts" folder, containing all the scripts
An overview of the folders architecture is reported here:
inserire immagine folders

### Genome assembly with HiFiasm

### Genome quality check with Quast and BUSCO

### Taxonomic identification with BLAST

### Gene prediction with Braker3

### Functional annotation with Diamond + GO and CAZyme annotation

## Bioinformatic pipelines in R
### UPGMA tree with hclust and ggtree
This pipeline was modified starting from [Brandon Güell et al.](https://fuzzyatelin.github.io/bioanth-stats/module-24/module-24.html)

Required packages and set working directory:

```r
library(adegenet)
library(ape)
library(ggplot2)
library(ggtree)
library(phangorn)
library(apex)

setwd('absolute_path_to_folder/main')
getwd()
```

Inport previously obtained [multi-FASTA sequences](#Taxonomic-identification-with-BLAST):

```r
files <- list.files(
  'work/tree/fasta',
  pattern = '\\.fasta$',
  full.names = TRUE
)

multifasta <- read.multiFASTA(files = files)
getNumSequences(multifasta)
```

Pipeline:

```r
#Visualize sequence polymorphisms
phydat_multifasta <- multidna2multiphyDat(multifasta)
concat_multifasta <- concatenate(phydat_multifasta)
image(concat_multifasta)

#Plot distances in a heatmap
dnabin_multifasta <- as.DNAbin(concat_multifasta)

dist_matrix <- dist.dna(dnabin_multifasta, model = "TN93")
length(dist_matrix) 

dist_plot <- as.data.frame(as.matrix(dist_matrix))
table.paint(dist_plot, cleg=0, clabel.row=.5, clabel.col=.5) 

#UPGMA clustering with hclust
h_cluster <- hclust(dist_matrix, method = "average", members = NULL) 

#check if UPGMA is appropriate
tre <- as.phylo(hclust(dist_matrix,method = "average"))
x <- as.vector(dist_matrix)
y <- as.vector(as.dist(cophenetic(tre)))
plot(x, y, xlab="original pairwise distances", ylab="pairwise distances on the tree",
 main="Is UPGMA appropriate?", pch=20, col=transp("black",.1), cex=3)
abline(lm(y~x), col="red")
cor(x,y)^2

#plot phylogenetic tree
phylogenetic_tree <- ggtree(h_cluster) + xlim(-0.4, .35) + theme_tree2() +
# Cladosporium strains isolated in this work are highlighted in bold (+ italic)
  geom_tiplab(fontface = 4, family = "serif",
              aes(subset = (node %in% c(9,8,15,10,12,13)))) +
# other strains are in italic
  geom_tiplab(fontface = 3, family = "serif",
              aes(subset = (node %in% c(14,11,17,16,5,4,6,7,18,2,1,3,19))))

#save plot
#ggsave(filename = "cladosporium_tree.png", plot = phylogenetic_tree,
path = 'results', width = 17, height =10, units = "cm")

```


### PCoA and heatmaps with ggplot2 and pheatmap

```r
```
