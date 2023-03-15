# goBRIDGE.sh: simple pipleline to process drug-seq data

# get help message with: bash goBRIDGE.sh -h
#   Run Pipeline for BRIDGE samples.
#
#   options:
#   h     Print this Help.
#   f     Runs fastqc on fastq files given a list with accessions (e.g. BRM1). Each line one accession.
#   c     Counts reads in fastq files given a list with accessions. See option 'f'.
#   p     Counts pattern/Barcodes in fastq files given a list with accessions and a list with pattern.
#   s     Run STARsolo given a list of accessions (see option 'f') and a species (e.g. 'human' or 'rat').
