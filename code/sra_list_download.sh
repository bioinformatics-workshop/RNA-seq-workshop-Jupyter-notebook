#!/bin/bash

# NCBI tools
module load ncbi_edirect

# Get all SRA runs for a given BioProject
# https://ncbi-hackathons.github.io/EDirectCookbook/

DB="sra"
BIOPROJ="PRJNA950346"
OUT_DIR="raw"

# Search for sample metadata
esearch -db bioproject -query $BIOPROJ | \
elink -target sra | efetch -format native -mode xml | \
xtract -pattern EXPERIMENT_PACKAGE \
-block SAMPLE -element VALUE -block RUN_SET -element PRIMARY_ID >> \
$OUT_DIR/${BIOPROJ}.metadata.tmp
