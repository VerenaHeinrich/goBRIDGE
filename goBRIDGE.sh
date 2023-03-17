#!/bin/bash

# custom paths:
ROOT="/home/verena_laupert/"
DATA_DIR=$ROOT"DATA/"
PIPELINE_DIR=$ROOT"BRIDGE_PIPELINE/"
META_DIR=$PIPELINE_DIR"META/"
REFERENCE_DIR=$ROOT"REFERENCE/"
RESULTS_DIR=$ROOT"RESULTS/"
FIGURES_DIR=$ROOT"ANALYSIS/FIGURES/"

# STAR version:
STAR=$ROOT"STAR-2.7.10a/source/STAR"

# general parameter:
THREADS=24
CPU=16
ulimit -n 2048

# load functions:
source $PIPELINE_DIR/functions.sh

help ()
{
   echo "Run Pipeline for BRIDGE samples."
   echo
   echo "options:"
   echo "h     Print this Help."
   echo "f     Runs fastqc on fastq files given a list with accessions (e.g. BRM1). Each line one accession."
   echo "c     Counts reads in fastq files given a list with accessions. See option 'f'."
   echo "p     Counts pattern/Barcodes in fastq files given a list with accessions and a list with pattern."    
   echo "s     Run STARsolo given a list of accessions (see option 'f') and a species (e.g. 'human' or 'rat')."
   echo
}

while getopts "hfcps:" option; do
   case $option in
      h) # display Help
         help
         exit;;
      f) # runs fastqc on fastq files given a list with accessions:
         run_fastqc $DATA_DIR $META_DIR/accessions
         ;;
      c) # count reads in fastq files given a list with accessions
         count_reads_in_fastq $DATA_DIR $META_DIR/accessions
         ;;
      p) # count pattern in fastq files given a list with accessions and a whitelist
         count_pattern_in_reads $DATA_DIR $META_DIR'whitelist96' $META_DIR/accessions
         Rscript R/summarize_pattern.R $META_DIR/accessions $DATA_DIR/ $FIGURES_DIR/
         ;;
      s) # run STAR solo
         SPECIES="$OPTARG"
         run_STARsolo $DATA_DIR $RESULTS_DIR $REFERENCE_DIR $SPECIES $META_DIR/accessions
         Rscript R/summarize_STARlog.R $META_DIR/accessions $RESULTS_DIR $FIGURES_DIR/
         ;;
      \?) # Invalid option
         echo "Error: Invalid option"
         exit;;
   esac
done
