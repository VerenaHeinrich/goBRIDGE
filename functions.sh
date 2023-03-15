
read_accession ()
{
  local accession_list=$1
  
  i=0
  while read this
    do 
    array[ $i ]="$this"        
    (( i=$i+1 ))
    done < $accession_list
  echo "${array[@]}" 
}

count_reads_in_fastq () 
{
  local dir_to_fastq=$1
  local accession_list=$2
  local accession_list=( $(read_accession $accession_list) )
   
  for acc in "${accession_list[@]}"
  do
    echo $acc
    
    OUT=$dir_to_fastq$acc"/read_counts.txt"
    if [ -f "OUT" ]; then
      rm $OUT
    fi
    echo -e 'sample\tnumber' | tee $OUT
  
    for i in $dir_to_fastq$acc/*fastq.gz; 
      do 
        n=$(($(zless -S "$i" |wc -l) / 4 ))
        echo -e $i'\t' $n | tee -a $OUT
      done
    echo "results are here: "$OUT
  done
}

count_pattern_in_reads () 
{
  local dir_to_fastq=$1
  local pattern_file=$2
  local accession_list=$3
  local accession_list=( $(read_accession $accession_list) )
  
  for acc in "${accession_list[@]}"
  do
    echo $acc
    
    OUT=$dir_to_fastq$acc"/pattern_counts.txt"
    if [ -f "OUT" ]; then
      rm $OUT
    fi
    echo -e 'pattern\tsample\tnumber' | tee $OUT
  
    for i in $dir_to_fastq$acc/*fastq.gz; 
      do 
        while read this
        do 
          pattern=$(echo $this | sed "s/ .*//")
          n=$(zless -S $i | grep $pattern | wc -l) 
          echo -e $pattern'\t'$i'\t'$n | tee -a $OUT
        done < $pattern_file
      done
      echo "results are here: "$OUT
  done
}

run_fastqc ()
{
  local dir_to_fastq=$1
  local dir_to_results=$2
  local accession_list=$5
  local accession_list=( $(read_accession $accession_list) )
   
  for acc in "${accession_list[@]}"
  do
    OUT_PREFIX=$dir_to_results$acc"/"
    OUT=$OUT_PREFIX/fastqc/
    if [ ! -d "$OUT_PREFIX" ]; then
      mkdir $OUT_PREFIX
      mkdir $OUT
    fi
    
    echo $acc
    fastq_files=`find $dir_to_fastq -type f -name "${acc}*"|sort`
    fastqc $fastq_files 
    
    echo "results are here: "$OUT
  done
}

run_STARsolo () 
{
  local dir_to_fastq=$1
  local dir_to_results=$2
  local dir_to_reference=$3
  local species=$4
  local whitelist=$5
  local accession_list=$6
  local accession_list=( $(read_accession $accession_list) )
  
  if [ ! -d "$dir_to_reference$species" ]; then
    echo "$dir_to_reference$species does not exist."
    exit
  fi
  
  if [[ $species == "human" ]]; then
  source $META_DIR/STAR.param.human.sh
  fi

  for acc in "${accession_list[@]}"
  do
    echo $acc
    
    OUT_PREFIX=$dir_to_results$acc"/"
    OUT=$OUT_PREFIX/STARsolo/
    if [ ! -d "$OUT_PREFIX" ]; then
      mkdir $OUT_PREFIX
      mkdir $OUT
    fi
    
    fastq_files=`find $dir_to_fastq -type f -name "${acc}*"|sort`
    fastq_array=(${fastq_files//:/ })  
  echo $fastq_files
    $STAR \
    --outSAMattributes NH HI AS nM MD \
    --runThreadN $THREADS \
    --genomeDir $REFERENCE \
    --sjdbGTFfile $GTF \
    --outSAMunmapped Within \
    KeepPairs \
    --outSAMtype BAM \
    SortedByCoordinate \
    --outSAMorder Paired \
    --limitBAMsortRAM 60413900847 \
    --readFilesCommand gunzip -c \
    --outFileNamePrefix $OUT \
    --readFilesIn ${fastq_array[1]} ${fastq_array[0]} \
    --soloType CB_UMI_Simple \
    --soloStrand Forward \
    --soloCBstart 25 \
    --soloCBlen 10 \
    --soloUMIstart 17 \
    --soloUMIlen 8 \
    --soloBarcodeReadLength 0 \
    --soloFeatures Gene GeneFull SJ Velocyto \
    --soloCBwhitelist $whitelist
    
    echo "results are here: "$OUT
  done
}