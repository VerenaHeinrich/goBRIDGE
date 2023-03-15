library(reshape2)
library(ggplot2)
library(readxl)

# own color palette (Bay colors):
bay_colors <- function(palette = 'help'){
  color <-  c('Blue', 'Green', 'Purple')
  hues <- c('Bright', 'Core', 'Mid', 'Dark')
  D <- data.frame(col=apply(expand.grid(color, hues), 1, function(x) paste(x[2], x[1], sep=" ")))
  D$hex <- c('#00BCFF',
             '#89D329',
             '#FF3162',
             '#0091DF',
             '#66B512',
             '#D30F4B',
             '#00617F',
             '#2B6636',
             '#624963',
             '#10384F',
             '#004422',
             '#443247')
  
  if(palette == 'help'){
    c(paste0(unique(gsub('.* ','',D$col)),'s'),
      unique(gsub(' .*','',D$col)),
      'Mixed')
  }
  else if(palette == 'Blues'){
    D$hex[grep('Blue',D$col)]
  }else if(palette == 'Greens'){
    D$hex[grep('Green',D$col)]
  }else if(palette == 'Purples'){
    D$hex[grep('Purple',D$col)]
  }else if(palette == 'Bright'){
    D$hex[grep('Bright',D$col)]
  }else if(palette == 'Core'){
    D$hex[grep('Core',D$col)]
  }else if(palette == 'Mid'){
    D$hex[grep('Mid',D$col)]
  }else if(palette == 'Dark'){
    D$hex[grep('Dark',D$col)]
  }else if(palette == 'Mixed'){
    D$hex[c(1,3,5,7,9,11,2,4,6,8,10,12)]
  }
}

# read list (each entry, one row); return data.frame
read_lists <- function(list, col.label=NULL){
  d <- read.table(list)
  if(! is.null(col.label)) colnames(d) <- col.label
  return(d)
}

# read star log file:
get_star_qc_table <- function(file){
  this <- scan(file=file, what='character',allowEscapes=T, sep='\n')
  this <- this[grepl('\\|', this)]
  this <- gsub(' {2,}','',this)
  this <- lapply(this, function(x) unlist(strsplit(x, '\\|\t')))
  this <- data.frame(t(matrix(unlist(this), nrow=length(this), byrow=TRUE)))
  colnames(this) <- gsub(' $','',this[1,])
  this[2,] <- gsub('%','',this[2,])
  this <- this[-1,]
  return(this)
}

# get a glimpse of STAR ouput:
summarize_STAR <- function(dir, out, accession_list){
  accessions <- read_lists(accession_list, 'accession')
  
  df <- NULL
  for(acc in accessions$accession){
    file <- list.files(paste0(dir, acc, '/STARsolo'), pattern=paste0('Log.final.out'), full.names = T, recursive = T)
    this <- get_star_qc_table(file)
    this <- this[,-c(1:4)]
    this <- this[-grep('Average input read length|number|Number of splices|Number of reads mapped|Number of reads unmapped|splices|chimeric|Deletion|Insertion', colnames(this))]
    this$acc <- acc
    if(is.null(df)) df <- this
    else df <- rbind(df, this)
  }
  melt <- reshape2::melt(df, id.vars = 'acc')
  melt$value <- as.numeric(melt$value)
  
  gg <- ggplot(melt, aes(x=acc, y=value, fill=acc))+
    geom_bar(stat='identity') +
    facet_wrap(.~variable, scales = 'free', )+
    theme_bw()+
    scale_fill_manual(values=bay_colors('Mixed'))+
    theme(axis.text.x = element_text(size=15, angle=90, hjust=0.5, vjust=0),
          axis.text.y = element_text(size=15),
          axis.title = element_text(size=18),
          legend.position = 'none')+
    xlab('')
  
    ggsave(gg, device='png', width = 300, height = 220, units='mm',
    file=paste0(out, paste(unlist(accessions), collapse='_'), '.STAR.summary.png'))
}

# visualize pattern/BC count in reads:
summarize_pattern_count <- function(data_dir, out, accession_list){
  library(reshape2)
  accessions <- read_lists(accession_list, 'accession')
  
  df <- NULL
  for(acc in accessions$accession){
    this <- read.table(paste0(data_dir, acc, '/pattern_counts.txt'), header=T)
    if(is.null(df)) df <- this
    else df <- rbind(df, this)
  }
  
  melt <- reshape2::melt(df)
  melt$run <- gsub('_.*','',basename(melt$sample))
  melt$pair <- ifelse(grepl('_R1_',melt$sample), 'R1', 'R2')
  
  gg <- ggplot(melt, aes(x=pattern, y=value, col=run, shape=pair)) +
    geom_point() +
    theme_classic()+
    theme(axis.text.x = element_text(angle=90, hjust = 1),
          axis.ticks.x = element_blank(),
          legend.position = 'right')+
    ylab('# reads') + xlab('') +
    scale_colour_manual(values=bay_colors('Mixed'))+
    scale_y_continuous(
      trans = "log10"
    )
  ggsave(gg, device='png', width = 300, height = 100, units='mm',
         file=paste0(out, paste(unlist(accessions), collapse='_'), '.BC.counts.png'))
}

# read excel file
read_plate_design(file, BC){
  d <- read_excel(file, col_names=T)
  d <- d[,-1]
  bc <-  read.table(BC)
  
  df <- NULL
  counter=0
  # create table:
  for(r in 1:8){
    for(c in 1:12){
      counter <- counter+1
      this <- data.frame(BAY_ID=as.character(d[r,c]),
                         BC=bc$V1[counter])
      if(is.null(df)) df <- this
      else df <- rbind(df, this)
    }
  }
  
  # add replicate information:
  df$rep <- 1
  tab <- table(df$BAY_ID)
  for(this in names(tab)){
    df$rep[df$BAY_ID==this] <- seq(tab[this])
  }
    
}
