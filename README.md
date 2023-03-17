goBRIDGE.sh: simple pipleline to process drug-seq data

get help message with: bash goBRIDGE.sh -h
>  Run Pipeline for BRIDGE samples.

>   options:
>   h     Print this Help.
>   f     Runs fastqc on fastq files given a list with accessions
>   c     Counts reads in fastq files given a list with accessions
>   p     Counts pattern/Barcodes in fastq files given a list with accessions and a list with pattern.
>   s     Run STARsolo given a list of accessions and a species (e.g. 'human' or 'rat').

> IMPORTANT:
adapt list of accessions: META/accessions. Each line one accession (e.g. BRM1).Scripts will go through every accession listed here.
adapt settings for STAR: META/STAR.param.GENERAL.sh