#!/usr/bin/env Rscript
args <-  commandArgs(trailingOnly=TRUE)
accession_list <- args[1]
data_dir <- args[2]
out <- args[3]

#accession_list <- "/home/verena_laupert/DATA/accessions"
#data_dir='~/DATA/'
# out <- '~/ANALYSIS/FIGURES/'

library(dplyr)
getCurrentFileLocation <-  function()
{
  this_file <- commandArgs() %>%
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
  if (length(this_file)==0)
  {
    this_file <- getSourceEditorContext()$path
  }
  return(dirname(this_file))
}
this.dir <- getCurrentFileLocation()
source(paste0(this.dir, '/functions.R'))
summarize_pattern_count(data_dir, out, accession_list)
