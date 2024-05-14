#!/usr/bin/bash -l


# Obtaining up-to-date genome information


# NCBI RefSeq Genome-based GO annotations Updated 12-18-2023
NCBI_URL="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/Bemisia_tabaci/latest_assembly_versions/GCF_001854935.1_ASM185493v1/GCF_001854935.1_ASM185493v1_gene_ontology.gaf.gz"


# UniProt Protein-based GO annotations Updated
UNIPROT_URL="https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/5124817.B_tabaci_1.goa"
UNIPROT_cBTQ1_URL="https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/370882.C_endosymbiont_cBtQ1_of_Bemisia_tabaci.goa"
UNIPROT_A_URL="https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/4481243.A_endosymbiont_of_Bemisia_tabaci_Q2.goa"


# Download remote files
wget -P genome $NCBI_URL