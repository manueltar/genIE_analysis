#!/usr/bin/env bash
 
 
MASTER_ROUTE=$1
output_dir=$MASTER_ROUTE

input_bwa_bam=$2
REFERENCE=$3
max_samclip=$4
sample=$5 


mem=$6
pc=$7
queue=$8

TOF_file=$9
replicates_file=${10}
tag=${11}
name=${12}



#### CONDA

source /software/hgi/installs/anaconda3/etc/profile.d/conda.sh


conda_samclip=$(echo "/nfs/team151/software/samclip/")

conda deactivate

conda activate $conda_samclip

output_dir=$(echo "$MASTER_ROUTE""/")


#### SAMCLIP STEP



type=$(echo "samclip_step")

outfile_samclip=$(echo "$output_dir""outfile""_""$type""_""$sample"".out")
name_samclip=$(echo "$type""_""$sample""_job")

output_bwa_bam_filtered=$(echo "$output_dir""$sample""_filtered.bam")


bsub -G team151 -o $outfile_samclip -q normal -n$pc -J $name_samclip -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
"samtools view -h $input_bwa_bam|samclip --ref $REFERENCE --max $max_samclip|samtools sort > $output_bwa_bam_filtered"

#### INDEX


type=$(echo "samclip_index_step")

outfile_index=$(echo "$output_dir""outfile""_""$type""_""$sample"".out")
name_index=$(echo "$type""_""$sample""_job")


bsub -G team151 -o $outfile_index -w"done($name_samclip)" -q normal -n$pc -J $name_index -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
"samtools index $output_bwa_bam_filtered"


#### COUNT


type=$(echo "samclip_count_step")

outfile_count=$(echo "$output_dir""outfile""_""$type""_""$sample"".out")
name_count=$(echo "$type""_""$sample""_job")

OUTPUT_counts_filtered=$(echo "$output_dir""$sample""_count"".txt")


bsub -G team151 -o $outfile_count -w"done($name_index)" -q normal -n$pc -J $name_count -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
"samtools view $output_bwa_bam_filtered|cut -f3|sort|uniq -c|sort -rg > $OUTPUT_counts_filtered"


#### ToF and replicates

echo -e "$name\t$sample\t$tag\t$output_bwa_bam_filtered" >> $replicates_file
echo -e "$OUTPUT_counts_filtered" >> $TOF_file
