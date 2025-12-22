# Comparative genomics on *Cladosporium* strains
This pipeline contains the codes used in the paper "Comparative genomics of the genus *Cladosporium* seemingly identified lack of selective pressure in soils historically polluted by hexachlorocyclohexane and polychlorobiphenyls".

## Data accessibility
Raw genomic (g)DNA sequences and assembled genomes are pubblicly available under the NCBI BioProject [PRJNA1348484](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1348484).\
Raw sequences were obtained via PacBio HiFi circular consensus sequencing.\
Quality and completeness of the assemblies are reported in the supplementary materials of the paper (scrivere supplementary material come link).

## Pipeline overview
The following steps were used in this work:
- [genome assembly with HiFiasm v0.16.0 and Minimap2 v2.30](#Genome-assembly-with-HiFiasm-and-Minimap)
- [genome quality check with Quast v5.3 and BUSCO v6.0.0](#Genome-quality-check-with-Quast-and-BUSCO)
- [taxonomic identification with local and online BLAST v2.17.0](#Taxonomic-identification-with-BLAST)
- [UPGMA tree with hclust and ggtree v3.12.0 on R v4.4](#UPGMA-tree-with-hclust-and-ggtree)
- [gene prediction with Braker3 v3.0.8](#Gene-prediction-with-Braker3)
- [functional annotation with Diamond v2.0.15 against the Swissprot database](#Functional-annotation-with-Diamond-+-CAZyme-annotation)
- [CAZyme annotation with dbcan3 v5.1.2](#Functional-annotation-with-Diamond-+-CAZyme-annotation)
- [Gene Ontology annotation with the Uniprot Retrive/ID mapping tool](#Gene-Ontology-(GO)-and-Enzyme-Comissions-(EC)-annotation)
- [Miltiple sequence alignment with ClustalW](#Miltiple-sequence-alignment-with-ClustalW)
- [PCoA and heatmaps with ggplot2 v3.5.2 and pheatmap v1.0.13 on R](#PCoA-and-heatmaps-with-ggplot2-and-pheatmap)

To better navigate in the multiple analyses performed in thiese pipelines, a folders architecture was created as follows:
- "main" folder, containing all the sub-folders for this analysis
  - "work" (sub)folder, containing raw and processed data for each pipeline
  - "results" (sub)folder, containing the final output of each pipeline
  - "scripts" (sub)folder, containing all the scripts
An overview of the folders architecture is reported here:\
inserire immagine folders

## Bioinformatic pipelines in Bash
Bioinformatic pipelines in Bash were performed in [Grex HPC](https://um-grex.github.io/grex-docs/grex/).\
[Singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html) was adopted for containerization. \
An example on how singularity sif files were installed is reported in the [genome assembly section](#Genome-assembly-with-HiFiasm-and-Minimap)

### Required files

### Genome assembly with HiFiasm and Minimap
A genome assembly step was performed using [HiFiasm](https://github.com/chhylp123/hifiasm), a tool built for PacBio libraries capable of handling high heterozygosis and eukaryote genomes.

```bash
cd path_to/main

# singularity and samtools are already present in Grex
module load singularity
module load samtools

cd work/sif
#install the sif for hifiasm
singularity search hifiasm
#copy and paste the path to hifiasm_latest.sif
singularity pull hifiasm_latest.sif library://pastepath/hifiasm_latest.sif
cd path_to/main

#perform genome assembly with HiFiasm (25 threads)
ls_names=($(cat names_bash.txt)) &&
MAX_names=$["$(cat names_bash.txt | wc -l)" - 1] &&
\
for i in $(seq 0 $MAX_names);
do
singularity exec work/sif/hifiasm-nf_latest.sif hifiasm \
-o work/genomes/${ls_names[$i]}_assembly work/raw_reads/${ls_names[$i]}_reads.fastq \
--primary -t 25
done

#we obtain a .gfa, to be converted in .fa
ls_names=($(cat names_bash.txt)) &&
MAX_names=$["$(cat names_bash.txt | wc -l)" - 1] &&
\
for i in $(seq 0 $MAX_names);
do
awk '/^S/{print ">"$2\n"$3}' work/genomes/${ls_names[$i]}_assembly.p_ctg.gfa | fold > work/genomes/CCD_clado-${ls_names[$i]}.p_ctg.fa
done

#remove unassembled sequences (length < 25,000 bp)
#install bbmap (not present in singularity)
cd work/sif
wget https://sourceforge.net/projects/bbmap/files/latest/download -O BBTools.tar.gz
tar -xzf BBTools.tar.gz
cd path_to/main

#use reformat in bbmap
ls_names=($(cat names_bash.txt)) &&
MAX_names=$["$(cat names_bash.txt | wc -l)" - 1] &&
\
for i in $(seq 0 $MAX_names);
do
work/sif/bbmap/reformat.sh in=work/genomes/CCD_clado-${ls_names[$i]}.p_ctg.fa \
out=work/genomes/CCD_clado-${ls_names[$i]}_filtered.p_ctg.fa minlength=25000
done
```

We obtained assembled genomes in FASTA format, ready for downstream analyses.

The genome assembly of strain "O" had an abnormal size. This might derive from the presence of contaminant sequences.

To solve this issue, raw reads of strain "O" were first aligned to the assembled genome of strain "F105".
Aligned reads, free of contamination, were then used for genome alignment in HiFiasm.

```bash
cd path_to/main

#use minimap2 to align raw DNA sequences of strain O to the previously assembled genome of strain F105
singularity exec work/sif/minimap2_v2.30.sif minimap2 \
work/genomes/CCD_clado-f105_filtered.p_ctg.fa work/raw_sequences/O_SRR35901019.fastq \
--sam-hit-only -x map-hifi -t 25 > work/minimap_align/O_minimapalign.sam

#samtools does not support paths. Moving working directory.
cd work/minimap_align/
#convert .sam in .fastq
samtools fastq O_minimapalign.sam > O_reads.fastq
#return back to "main" folder
cd path_to/main
```

It is now possible to use the filtered reads "O_reads.fastq" to perform a genome assembly using HiFiasm, as previously described.

### Genome quality check with Quast and BUSCO
To assess genome quality and completeness, the tools [Quast](https://github.com/ablab/quast) and [BUSCO](https://busco.ezlab.org/) were used.

```bash
#Quast was already present in Grex
module load Quast

ls_names=($(cat names_bash.txt)) &&
MAX_names=$["$(cat names_bash.txt | wc -l)" - 1] &&
\
for i in $(seq 0 $MAX_names);
do
quast work/genomes/${ls_names[$i]}_filtered.p_ctg.fa -o work/quast
done

#running BUSCO from singularity
ls_names=($(cat names_bash.txt)) &&
MAX_names=$["$(cat names_bash.txt | wc -l)" - 1] &&
\
for i in $(seq 0 $MAX_names);
do
singularity exec work/sif/busco_v6.0.0.sif busco \
 -i work/genomes/${ls_names[$i]}_filtered.p_ctg.fa -o work/BUSCO/${ls_names[$i]} -m genome -l -c 25
done
```

### Taxonomic identification with BLAST

```bash
```

### Inport of additional Cladosporium genomes
Additional Cladosporium genomes were inported following the [NCBI tutorial](https://github.com/ncbi/datasets)

The [codes_additional.txt and names_additional.txt files](#Required-files) are required:

```bash
ls_codes=($(cat codes_additional.txt)) &&
ls_names=($(cat names_additional.txt)) &&
MAX_names=$["$(cat names_additional.txt | wc -l)" - 1] &&
\
for i in $(seq 0 $MAX_names);
do
alias datasets=work/genomes
datasets download genome accession ${ls_codes[$i]} --filename ${ls_names[$i]}.zip
done

#unzip files
ls_names=($(cat names_additional.txt)) &&
MAX_names=$["$(cat names_additional.txt | wc -l)" - 1] &&
\
for i in $(seq 0 $MAX_names);
do
echo "start ${ls_names[$i]}"
unzip ${ls_names[$i]}.zip
echo "end ${ls_names[$i]}"
done
```
Downloaded genomes were moved to the folder "genomes". \
From here on, the file names_bash.txt was updated to contain also the additional strains.

### Gene prediction with Braker3

```bash
```

### Functional annotation with Diamond + CAZyme annotation
Functional annotation was performed on the protein sequences obtained from [Braker3](#Gene-prediction-with-Braker3).
[Diamond](https://github.com/bbuchfink/diamond) is an effective alignment tool and accepts only protein sequences as input.

A Diamond database was created starting from the [Swissprot protein database](https://www.uniprot.org/uniprotkb?query=reviewed:true). The database was downloaded in the folder "work/database" under the name "uniprot_sprot.fasta.gz".

```bash
cd path_to/main

#diamond was already present in Grex
module load diamond

#make diamond database
diamond makedb -d swissprot -p 10 --in 'work/database/uniprot_sprot.fasta.gz' -d 'work/database'

#run diamond
ls_names=($(cat names_bash.txt)) &&
MAX_names=$["$(cat names_bash.txt | wc -l)" - 1] &&
\
for i in $(seq 0 $MAX_names);
do
diamond blastp --db 'database/swissprot.dmnd' -q 'work/braker_out/${ls_names[$i]}/braker.aa' \
--out 'work/annotation/${ls_names[$i]}_diamond.txt' -p 10 --sensitive --max-target-seqs 1
done
```

It is important to specify the option "--max-target-seqs 1" as each predicted gene must have only one annotation. \
We obtained a series of diamond.txt files containing both annotated protein names and their relative Uniprot identifiers.
These identifiers can be used to obtain Gene Ontology (GO) terms and Enzyme Commission (EC) numbers.

A section covernig how to [obtain GO and EC annotations] is reported in the following section.

## Bioinformatic pipelines in R
### Required files
Prepare a "metadata" file containing all the grouping factors of your samples.\
Here you can find an example of a metadata.xlsx file:\
head metadata

Prepare a "names_R" file containing a code to recognize your annotation files and the species name of each strain.\
Here you can find an example of a names.txt file:\
head names



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
  'work/tree',
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
ggsave(filename = "cladosporium_tree.png", plot = phylogenetic_tree,
path = 'results', width = 17, height =10, units = "cm")
```

### Gene Ontology (GO) and Enzyme Comissions (EC) annotation
With the use of the Swissprot protein database, it is possible to retrive both Gene Ontology (GO) and Enzyme Commission (EC) annotations, among others.

Required libraries and set working directory:

```r
library(stringr)

setwd('absolute_path_to_folder/main')
getwd()
```

Inport [names_R](#Required-files) and and previously obtained [diamond files](#Functional-annotation-with-Diamond-+-CAZyme-annotation):

```r
names <- read.table("names_R.txt", sep = ';', fill = T)

for (i in c(1:length(names[,1]))) {
  df_strain <- names[i,1]
  temp_df <- read.table(paste("work/annotation/", df_strain, "_diamond.tsv", sep = ""),
              sep = '\"', fill = T, quote="")
#subset annotations with <30% similarity
  temp_df_subset <- subset(temp_df, temp_df$V3 >= 30)
#split the ID column by "|" and retrive only the UniProt codes
  temp_list <- strsplit(temp_df_subset$V2, "[|]")
  temp_df2 <- do.call(rbind.data.frame, temp_list) 
  temp_list2 <- list(temp_df2[,3])
#save Uniprot ID list
  data.table::fwrite(temp_list2, file = paste("work/annotation/", df_strain, "_uniprot_list.txt", sep = ""))
#create and order and unique ID df with frequencies of each ID
  df_unique <- as.data.frame(table(temp_df2[,3]))
  df_unique <- df_unique[order(df_unique$Var1),]
#
  assign(paste(df_strain,"unique", sep = "_"), df_unique)
}
```

An UniProt ID list for each strain is saved in "work/annotation/..._uniprot_list.txt".
Each of these lists must be uploaded in the [UniProt batch retrival/ID/ tool](https://www.ebi.ac.uk/training/online/courses/uniprot-exploring-protein-sequence-and-functional-info/how-to-use-uniprot-tools-clone/batch-retrieval-id-mapping/). 

In this online tool, columns of interest can be selected in the "configure columns" section.
In this work, the position of certain columns is important for downstream analysis:
- Column 6 - Protein names, which will become V9
- Column 8 - Gene Ontology (GO), which will become V11

Once obtained, the tables must be saved with the "download" button by selecting the "tsv" and "uncompressed" format. 
All saved tables are named "idmapping.tsv" and can be inported on R to obtain the complete annotation table.

```r
for (i in c(1:length(names[,1]))) {
  df_strain <- names[i,1]
#inport output di Uniptot
  temp_df <- read.table('Funghi_braker/Braker_out/Braker_F32/idmapping.tsv',
            sep = '\t', header = TRUE, fill = T, , quote="")
#quote = "" is necessary because the character ' is present
#fill =T because there are some empty cells
#order dataframe
  temp_df <- temp_df[order(temp_df$From),]
  assign(paste(df_strain, "idmapping", sep = "_"), temp_df)
}

#check if both dataframes are correctly ordered
for (i in c(1:length(names[,1]))) {
  df_strain <- names[i,1]
  unique <- get(paste(df_strain, "_unique"))
  idmapping <- get(paste(df_strain, "_idmapping"))
  df_strain
  table(unique$Var1 == idmapping$From)
 }

#if all dataframes are correct, merge them
for (i in c(1:length(names[,1]))) {
  df_strain <- names[i,1]
  temp_df <- cbind(get(paste(df_strain, "_unique")),
                   get(paste(df_strain, "_idmapping")))
  write.table(temp_df, file = paste("work/annotation/". df_strain, "_annot.tsv", sep = ""), row.names = F)
 }
```

We obtained complete annotation tables that can be used for downstream analyses. \
An example of "annot.tsv" table is reported below: \
annot.tsv

### Miltiple sequence alignment with ClustalW
This analysis was performed with the online tool [ClustalW](https://www.genome.jp/tools-bin/clustalw).

Multi-FASTA files for each sequence of interest were used as input.

Multi-FASTA files were obtained in R, starting from the output of [GO functional annotation](#Gene-Ontology-(GO)-and-Enzyme-Comissions-(EC)-annotation).

Required libraries and set working directory:

```r
library(dplyr)
library(stringr)
library(tidyr)
library(tidytext)
library(readxl)

#avoid scientific numbering
options(scipen=999)

setwd('absolute_path_to_folder/main')
getwd()
```

Inport [names_R](#Required-files) and and previously obtained [annotation files](#Gene-Ontology-(GO)-and-Enzyme-Comissions-(EC)-annotation):

```r
names <- read.table("names_R.txt", sep = ';', fill = T)
sum_df <- data.frame()

#two different annotation files are required:
#1. Diamond annotation
#2. Uniprot Gene Ontology annotation

#In the second table, due to the presence of all common separators in the main text of the annotation files, " was used to separate values in the table. 
for (i in c(1:length(names[,1]))) {
  df_strain <- names[i,1]
  temp_df <-  read.table(paste("work/annotation/", df_strain, "_annot.tsv", sep = ""),
                 sep = '\"', fill = T, quote="")
#
  second_df <-  read.table(paste("work/annotation/", df_strain, "_diamond.tsv", sep = ""),
                 sep = '\"', fill = T, quote="")
#
  assign(paste(df_strain,"complete", sep = "_"), temp_df0)
  assign(paste(df_strain, "diamond", sep = "_", second_df)
}
```

Pipeline (still requires optimization):

```r
# repeat the pipeline for these four functions and for each strain
#.*GO:0018583 biphenil diol
#.*GO:0018784 haloacid
#.*GO:0018786 haloalkane
#.*GO:0018666 2,4-dichlorophenol 6-monooxygenase

subset(F32_complete, grepl(".*GO:0018583" ,F32_complete$V11))
#insert the value of V2 (i.e. the Uniprot identifier of the function) in the following function:
subset(F32_diamond, grepl(".*TFDB_CUPPJ", F32_diamond$V2))

#obtain the gene ID and search it in the fasta file (you can use the notes app)
#copy the selected fasta in the multi-fasta file and rename the header to specify function and strain
```

In ClustalW, the options "output format = clustal", "slow/accurate" and "protein" were selected.

Copy and paste the alignment scores in a separate note and format it like a tsv file.\
An example of a tsv of the alignment scores is reported here:\
alignment scores

The alignments section in the ClustalW output is also very informative as it highlights the matches and mismatches of your sequences.

### PCoA and heatmaps with ggplot2 and pheatmap
This pipeline was modified starting from a [CD Genomics tutorial](https://bioinfo.cd-genomics.com/resource-pcoa-analysis-with-R.html).

Required packages and set working directory:

```r
library(vegan)
library(ggplot2)
library(ggforce)
library(dplyr)
library(stringr)
library(tidyr)
library(tidytext)
library(readxl)
library(pheatmap)
library(ggplotify)
library(ggvegan)
library(dendextend)
library(rstatix)

#avoid scientific numbering
options(scipen=999)
setwd('absolute_path_to_folder/main')
getwd()
```

Inport [metadata](#Required-files), [names_R](#Required-files) and previously obtained [annotation files](#Gene-Ontology-(GO)-and-Enzyme-Comissions-(EC)-annotation):

```r
metadata <- as.data.frame(as.matrix(read_excel("work/annotation/metadata.xlsx")))
metadata <- metadata[order(metadata$species),]

names <- read.table("names_R.txt", sep = ';', fill = T)
sum_df <- data.frame()

#Due to the presence of all common separators in the main text of the annotation files, " was used to separate values in the table.
for (i in c(1:length(names[,1]))) {
  df_strain <- names[i,1]
  temp_df0 <-  read.table(paste("work/annotation/", df_strain, "_annot.tsv", sep = ""),
                 sep = '\"', fill = T, quote="")
  temp_df <- as.data.frame(unlist(str_split(as.vector(temp_df0$V11[-1]), ";")))
  temp_df$strain <- names[i,2]
  sum_df <- rbind(sum_df, temp_df)
  assign(paste(df_strain,"complete", sep = "_"), temp_df0)
  assign(df_strain, temp_df)
}

assign("GO_df", sum_df)
colnames(GO_df) <- c("GO", "strain")
GO_df <- subset(GO_df, GO_df$GO!="")
```

Pipeline:

```r
GO_df$Freq <- 1

#obtain the overall number of GOs in each strain
GO_sum <- with(GO_df, aggregate(Freq, by = list(strain), FUN = "sum"))
colnames(GO_sum) <- c("strain",  "freq")

#obtain the number of each GO term in each strain
Df <- with(GO_df, aggregate(Freq, by = list(strain, GO), FUN = "sum"))
colnames(Df) <- c("strain", "GO",  "freq")

#transform the data from long to wide format
df_wide <- reshape(Df, idvar = "strain", timevar = "GO", direction = "wide")
df_wide[is.na(df_wide)] <- 0
rownames(df_wide) <- df_wide$strain
df_wide <- subset(df_wide, select=-c(strain))
colnames(df_wide) <- sub("freq\\. *", "", colnames(df_wide))

#normalize GO frequencies as "GO per million"
ifelse(unique(rownames(df_wide)) == GO_sum$strain,
       df_wide <- df_wide/GO_sum$freq * 1000000,
       print("ERROR")
)

#check that all strains are present in the metadata file
table(rownames(df_wide) == metadata$species)
```

PCoA on overall GO terms:

```r
#obtain Bray-Curtis distances
dist_clado_wide <- vegdist(df_wide, method = "bray")
pcoa_clado_wide <- cmdscale(dist_clado_wide, eig = TRUE, k = 3)

#obtain position of the points and variance covered for each PCoA axis
points0c_wide <- as.data.frame(pcoa_clado_wide$points)
colnames(points0c_wide) <- c("PCoA1", "PCoA2","PCoA3")
pointsc_wide <- cbind(points0c_wide, metadata)
variancec_wide <- round(100 * pcoa_clado_wide$eig / sum(pcoa_clado_wide$eig), 2)

#plot the PCoA
pcoa_clado_view_wide <- ggplot(pointsc_wide, aes(x = PCoA1, y = PCoA2, colour = isolation)) +
  geom_point(size = 3) + scale_color_manual(values = c("lightblue4","red3")) +
#  geom_text(aes(label = ID_codes, hjust = -0.75, vjust = -0.75, size = 12)) +
  stat_ellipse(level = 0.95, linetype = 2, alpha = 0.3, linewidth = 1) +
  labs(x = paste0("PCoA1 (", variancec_wide[1], "%)"),
       y = paste0("PCoA2 (", variancec_wide[2], "%)")
  ) +
  theme_test() + 
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 14))

pcoa_clado_view_1_3_wide <- ggplot(pointsc_wide, aes(x = PCoA1, y = PCoA3, colour = isolation)) +
  geom_point(size = 3) + scale_color_manual(values = c("lightblue4","red3")) +
#  geom_text(aes(label = ID_codes, hjust = -0.75, vjust = -0.75, size = 12)) +
  stat_ellipse(level = 0.95, linetype = 2, alpha = 0.3, linewidth = 1) +
  labs(x = paste0("PCoA1 (", variancec_wide[1], "%)"),
       y = paste0("PCoA3 (", variancec_wide[3], "%)")) + 
  theme_test() +
  theme(axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 14))

#save plots
ggsave("pcoa_clado_GO1_2.jpg",plot = pcoa_clado_view_wide,
path = "results/", width = 6, height = 5.5)
ggsave("pcoa_clado_GO1_3.jpg",plot = pcoa_clado_view_1_3_wide,
path = "results/", width = 7.5, height = 5.5)

#the two graphs were combined using paint.net, which was also used to insert clear labels for each dot
```

Multiple heatmaps were created with almost identical pipelines.
Here, two examples are provided.

Heatmap on Halogenated Organic Compounds putative degrading genes:

```r
#select GOs of interest
haloalkane <- subset(Df, grepl(".*GO:0018786", Df$GO))
haloalkane$GO <- haloalkane$GO[1]
haloacid <- subset(Df, grepl(".*GO:0018784", Df$GO))
haloacid$GO <- haloacid$GO[1]
biphenyl <- subset(Df, grepl(".*GO:0018583", Df$GO))
biphenyl$GO <- biphenyl$GO[1]
dichloro <- subset(Df, grepl(".*GO:0018666", Df$GO))
dichloro$GO <- dichloro$GO[1] 
cyclohexadiene <- subset(Df, grepl(".*GO:0018502", Df$GO))
cyclohexadiene$GO <- cyclohexadiene$GO[1] 

#Laccase numbers were obtained with CAZyme annotation (group AA1)
lac <- read_csv("work/annotation/AA1.csv", show_col_types = F, )

laccase <- data.frame(strain = lac$strain,
                      GO =  "Laccase [CAZyme:AA1]",
                      freq = lac$freq)

#Cytochrome P450 numbers were obtained with protein names
cyp <- c()
cyp_strain <- c()

for (i in c(1:length(names[,1]))) {
  df_strain <- paste(names[i,1], "complete", sep = "_")
  cyp_temp <- dim(subset(get(df_strain),
                  grepl(".*P450",
                  get(df_strain)$V9)))[1]
  strain_temp <- names[i,2]
  cyp <- c(cyp, cyp_temp)
  cyp_strain <- c(cyp_strain, strain_temp)
}

cyp450 <- data.frame(strain = cyp_strain,
                      GO =  "Cytochrome P450 [EC 1.14.x.x]",
                      freq = cyp) 

#GO numbers were not normalized since we are looking at specific functions
degradation <- rbind(haloacid,haloalkane,biphenyl,dichloro,cyclohexadiene,cyp450,laccase)
long_degradation <- with(degradation, aggregate(freq, by=list(strain, GO), FUN=sum))
colnames(long_degradation) <- c("strain", "GO", "freq")

#transform the data from long to wide format
wide_degradation <- reshape(long_degradation, idvar = "strain", timevar = "GO",
direction = "wide")
rownames(wide_degradation) <- wide_degradation$strain
wide_degradation <- wide_degradation[order(wide_degradation$strain),]
wide_degradation <- subset(wide_degradation, select=-c(strain))
wide_degradation[is.na(wide_degradation)] <- 0

#check that all strains are present in the metadata file
table(rownames(wide_degradation) == metadata$species)

heat_clado_deg <- wide_degradation
colnames(heat_clado_deg) <- sub("freq\\. *", "", colnames(heat_clado_deg))
heat_clado_deg_t <- t(heat_clado_deg)

#calculate Bray-Curtis distance
dist_deg_rows <- vegdist(heat_clado_deg_t, method = "bray")
dist_deg_cols <- vegdist(heat_clado_deg, method = "bray")

#attach metadata
dex_clado <- data.frame(row.names = metadata$species, "Isolation source" = metadata$isolation)

#plot
heat_clado_deg <- as.ggplot(pheatmap(heat_clado_deg_t, scale="row", angle_col = 45, , display_numbers = heat_clado_deg_t,
                                clustering_distance_cols = dist_deg_cols, clustering_distance_rows = dist_deg_rows,
                           annotation_col = dex_clado, annotation_colors=list('Isolation.source'=c('polluted soil'="red3",
                                                                                                   'non-polluted matrices'="lightblue4")),
                           number_color = "black", fontsize_number = 10,
                                color=colorRampPalette(c("navy", "white", "red"))(50))) + theme(plot.margin=unit(c(0,0,0,1.25), 'cm'))

#save plot
ggsave("heat_clado_deg.jpg",plot = heat_clado_deg,
path = "results/", width = 14, height = 5)
```

Heatmap of [ClustalW results](#Miltiple-sequence-alignment-with-ClustalW):

```r
#showing the pipeline with haloacid dehalogenase
#other heatmaps were obtained following the same steps
haloacid <- read.table("work/clustal/Haloacid_dehalogenase_clustal_score.txt",
sep = ',', fill = T, header = T)

#obtain all comparisons between strains
strain_order <- c("C. europaeum F32","C. xylophilum F42","C. cladosporioides F51",
                  "C. pseudocladosporioides F105","C. cladosporioides LV","C. pseudocladosporioides O",
                  "C. velox C4","C. sphaerospermum R2409","C. rectoides A179",
                  "C. halotolerans T138_S3","C. allicinum IBT_42152","C. cladosporioides SYC63",
                  "C. cladosporioides ACCC_36060","C. cladosporioides ZJUC171","C. fusiforme IBT_42164",
                  "C. inversicolor IBT_42153","C. pseudocladosporioides ZJUC127","C. sp SL-16",
                  "C. xylophilum F42_1","C. xylophilum F42_2",
                  "C. velox C4_1","C. velox C4_2",
                  "C. sp SL-16_1","C. sp SL-16_2",
                  "C. pseudocladosporioides F105_1","C. pseudocladosporioides F105_2",
                  "C. pseudocladosporioides O_1","C. pseudocladosporioides O_2",
                  "C. cladosporioides SYC63_1","C. cladosporioides SYC63_2",
                  "C. cladosporioides ZJUC171_1","C. cladosporioides ZJUC171_2",
                  "C. fusiforme IBT_42164_1","C. fusiforme IBT_42164_2",
                  "C. pseudocladosporioides ZJUC127_1","C. pseudocladosporioides ZJUC127_2")
against_self <-  data.frame(strain1 = strain_order,
                            strain2 = strain_order,
                            score = 100)
haloacid2 <- data.frame(strain1 = haloacid$strain2,
                        strain2 = haloacid$strain1,
                        score = haloacid$score)

haloacid <- rbind(haloacid, haloacid2, subset(against_self,
 against_self$strain1 %in% haloacid$strain1 |
 against_self$strain1 %in% haloacid$strain2))

#remove strains with problematic reads
haloacid <- subset(haloacid, haloacid$strain1 != "C. cladosporioides ZJUC171" &
                     haloacid$strain2 != "C. cladosporioides ZJUC171" &
                     haloacid$strain1 != "C. velox C4" &
                     haloacid$strain2 != "C. velox C4")

#transform the data from long to wide format
wide_haloacid <- reshape(haloacid, idvar = "strain1", timevar = "strain2", direction = "wide")

rownames(wide_haloacid) <- wide_haloacid$strain1
wide_haloacid <- wide_haloacid[order(wide_haloacid$strain1),]
wide_haloacid <- subset(wide_haloacid, select=-c(strain1))
wide_haloacid[is.na(wide_haloacid)] <- 0

heat_clado_haloacid <- wide_haloacid
heat_clado_haloacid_t <- t(heat_clado_haloacid)

#calculate Bray-Curtis distance
dist_haloacid_rows <- vegdist(heat_clado_haloacid_t, method = "bray")
dist_haloacid_cols <- vegdist(heat_clado_haloacid, method = "bray")

#attach metadata
meta_haloacid <- subset(metadata, metadata$species %in% haloacid$strain1)
dex_clado <- data.frame(row.names = meta_haloacid$species, "Isolation source" = meta_haloacid$isolation)

#plot heatmap
heat_clado_haloacid <- as.ggplot(pheatmap(heat_clado_haloacid_t, angle_col = 45, #display_numbers = heat_clado_haloacid_t,
                                          clustering_distance_cols = dist_haloacid_cols, clustering_distance_rows = dist_haloacid_rows,
                                          annotation_col = dex_clado, annotation_colors=list('Isolation.source'=c('polluted soil'="red3",
                                                                                                                  'non-polluted matrices'="lightblue4")),
                                          #number_color = "black", fontsize_number = 10,
                                          color=colorRampPalette(c("navy", "white", "red"))(50),
                                          breaks = seq(from=60, to=100,by=((100-60)/50)))) +
  theme(plot.margin=unit(c(0,0,0,2), 'cm'))

#save plot
ggsave("heat_haloacid.jpg",plot = heat_clado_haloacid,
path = "results/", width = 10.5, height = 5.5)
```



