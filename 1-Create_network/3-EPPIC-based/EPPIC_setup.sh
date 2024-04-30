#!/usr/bin/bash

## download the git source
git clone https://github.com/eppic-team/eppic


## setup the databases
cd ../
cd Databases
mkdir Uniprot100
cd Uniprot100
wget https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz
gunzip -c uniref100.fasta.gz | makeblastdb -dbtype prot -logfile makeblastdb.log -parse_seqids -out uniref100.fasta -title uniref100.fasta


## mmseq2 software via conda
conda install -c conda-forge -c bioconda mmseqs2
# OR
wget https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz
tar xvfz mmseqs-linux-avx2.tar.gz
export PATH=$(pwd)/mmseqs/bin/:$PATH

## clustalomega
wget http://prdownloads.sourceforge.net/argtable/argtable2-13.tar.gz
tar -xf argtable2-13.tar.gz
./configure --prefix=/home/newazkha/common_bin/argtable2-13

wget http://www.clustal.org/omega/clustal-omega-1.2.4.tar.gz
tar -xf clustal-omega-1.2.4.tar.gz
./configure CFLAGS='-I/home/newazkha/common_bin/argtable2-13/src' LDFLAGS='-L/home/newazkha/common_bin/argtable2-13/lib'
make
make install prefix=/home/newazkha/common_bin/clustal-omega-1.2.4

## EPIIC sofwtare
wget http://eppic-web.org/downloads/eppic.zip
unzip eppic.zip
vim ~/.eppic.conf
## Add
# BLAST_DB_DIR
# BLAST_DB
# MMSEQS_BIN
# BLASTP_BIN
# CLUSTALO_BIN
