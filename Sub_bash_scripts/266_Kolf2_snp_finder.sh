#!/bin/bash
 
################################################################################################################################## Dependencies:
 

#bcftools=/nfs/users/nfs_m/mt19/sOFTWARE/bcftools-1.6/bcftools

#################################################################################################### Script
 

global=$1
input_variants=$2
master_tbi=$3
type=$4


output_file=$(echo "$global""/""$type""_SNPS"".txt")
stderr_file=$(echo "$global""/""$type""_outfile"".out")
touch $stderrfile
echo -n "" > $stderrfile

name_file=$(echo "$type")

 
# bsub -G team151 -o $stderr_file -M 4000 -J $name__file  -R"select[mem>4000] rusage[mem=4000] span[hosts=1]" -n1 -q normal -- \


bsub -G team151 -o $stderr_file -M 4000  -J $name_file -R"select[mem>=4000] rusage[mem=4000] span[hosts=1]" -n1 -q normal -- \
"bcftools view -R $input_variants $master_tbi > $output_file"
