# Comparative genomics on *Cladosporium* strains
This pipeline contains the codes used in the paper "Comparative genomics of the genus *Cladosporium* seemingly identified lack of selective pressure in soils historically polluted by hexachlorocyclohexane and polychlorobiphenyls".

## Data accessibility
Raw data and outputs are pubblicly available.
Raw genomic (g)DNA sequences and assembled genomes are pubblicly available under the NCBI BioProject [PRJNA1348484](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1348484).\
Raw sequences were obtained via PacBio HiFi circular consensus sequencing.\
Quality and completeness of the assemblies are reported in the supplementary materials of the paper (scrivere supplementary material come link).

## Pipeline overview
The following steps were used in this work:
- genome assembly with HiFiasm v0.16.0 and Minimap2 v2.30
- genome quality check with BUSCO v6.0.0 and Quast v5.3
- taxonomic identification with local and online BLAST v2.17.0
- UPGMA tree using hclust and ggtree v3.12.0 on R v4.4
- gene prediction with Braker3 v3.0.8
- functional annotation with Diamond v2.0.15 against the Swissprot database
- Gene Ontology annotation with the Uniprot Retrive/ID mapping tool
- CAZyme annotation with dbcan3 v5.1.2
- PCoA and heatmaps with ggplot v3.5.2 and pheatmap v1.0.13 on R

