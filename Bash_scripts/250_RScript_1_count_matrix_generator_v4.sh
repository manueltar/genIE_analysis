#!/bin/bash
 
# /software/R-3.6.1/bin/Rscript
Rscript=/software/R-4.1.0/bin/Rscript
script_113=/nfs/users/nfs_m/mt19/Scripts/R/113_my_own_functions_partII_v3_ATAC.R
script_113=/nfs/users/nfs_m/mt19/Scripts/R/113_my_own_functions_partII_v3.R
script_248=/nfs/users/nfs_m/mt19/Scripts/Wraper_scripts/248_CIGAR_extractor.sh
#script_114=/nfs/users/nfs_m/mt19/Scripts/R/114_my_own_functions_genIE_index_demultiplexation_v2_ATAC.
script_114=/nfs/users/nfs_m/mt19/Scripts/R/114_my_own_functions_genIE_index_demultiplexation_v2.R

 

################################################################################################################################## Dependencies: 

MASTER_ROUTE=$1
REFERENCE=$2
PREFIXES_Table_EXPANDED=$3
output=$4
route_fastq=$5
parameter=$6
input_regions=$7
batch=$8

index_check_summary=$(echo "$MASTER_ROUTE""index_check_summary.txt")
# /lustre/scratch115/teams/soranzo/projects/genIE_analysis/ReSeq/BWA/

#/nfs/users/nfs_m/mt19/Scripts/Wraper_scripts/200_CREATION.sh


#-w\"done(Rscript_1_NO_FLASH)\"

touch $output
echo -n "" > $output

echo -e "#!/bin/bash"  >> $output
echo -e "\n\n"  >> $output


output_dir=$(echo "$MASTER_ROUTE""NO_FLASH")
ToF_NO_FLASH=$(echo "$output_dir""/""ToF"".txt")
replicates_NO_FLASH=$(echo "$output_dir""/""replicates"".tsv")
counter_NO_FLASH=$(echo "$parameter""_""$batch""_NO_FLASH")

stderrfile_NO_FLASH_113=$(echo "$output_dir""/""Rscript_113.out")

touch $stderrfile_NO_FLASH_113
echo -n "" > $stderrfile_NO_FLASH_113

stderrfile_NO_FLASH_114=$(echo "$output_dir""/""Rscript_114.out")

touch $stderrfile_NO_FLASH_114
echo -n "" > $stderrfile_NO_FLASH_114

outfile_Script_248=$(echo "$output_dir""/""Script_248.out")

touch $outfile_Script_248
echo -n "" > $outfile_Script_248



echo -e "cd $output_dir/" >> $output
echo -e "bsub -G team151 -o $output_dir/Rscript_113.out -M 4000  -J Rscript_113_$counter_NO_FLASH -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
echo -e "\"$Rscript $script_113 --file $ToF_NO_FLASH --fasta $REFERENCE --prefixes $PREFIXES_Table_EXPANDED --type NO_FLASH\"" >> $output

echo -e "bsub -G team151 -o $output_dir/Script_248.out -M 16000 -J Script_248_$counter_NO_FLASH -R\"select[mem>=16000] rusage[mem=16000] span[hosts=1]\" -n4 -q normal -- \\" >> $output
echo -e "\"bash $script_248 $replicates_NO_FLASH $MASTER_ROUTE NO_FLASH\"" >> $output

echo -e "bsub -G team151 -o $output_dir/Rscript_114.out -M 4000  -J Rscript_114_$counter_NO_FLASH -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
echo -e "\"$Rscript $script_114 --prefixes_EXPERIMENT $index_check_summary --prefixes_Eve $PREFIXES_Table_EXPANDED --route_fastq $route_fastq --out INDEXES_TABLE.txt\"" >> $output


echo -e "\n\n"  >> $output




#######################################################################################################################################################################################################

output_dir_filtered=$(echo "$MASTER_ROUTE""NO_FLASH""_FILTERED")
ToF_NO_FLASH_filtered=$(echo "$output_dir_filtered""/""ToF_filtered"".txt")
replicates_NO_FLASH_filtered=$(echo "$output_dir_filtered""/""replicates_filtered"".tsv")
counter_NO_FLASH_FILTERED=$(echo "$parameter""_""$batch""_NO_FLASH_FILTERED")


stderrfile_NO_FLASH_filtered_113=$(echo "$output_dir_filtered""/""Rscript_113.out")

touch $stderrfile_NO_FLASH_filtered_113
echo -n "" > $stderrfile_NO_FLASH_filtered_113

stderrfile_NO_FLASH_filtered_114=$(echo "$output_dir_filtered""/""Rscript_114.out")

touch $stderrfile_NO_FLASH_filtered_114
echo -n "" > $stderrfile_NO_FLASH_filtered_114

outfile_Script_248=$(echo "$output_dir_filtered""/""Script_248.out")

touch $outfile_Script_248
echo -n "" > $outfile_Script_248




echo -e "cd $output_dir_filtered/" >> $output
echo -e "bsub -G team151 -o $output_dir_filtered/Rscript_113.out -M 4000  -J Rscript_113_$counter_NO_FLASH_FILTERED -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
echo -e "\"$Rscript $script_113 --file $ToF_NO_FLASH_filtered --fasta $REFERENCE --prefixes $PREFIXES_Table_EXPANDED --type NO_FLASH_FILTERED\"" >> $output

echo -e "bsub -G team151 -o $output_dir_filtered/Script_248.out -M 16000 -J Script_248_$counter_NO_FLASH_FILTERED -R\"select[mem>=16000] rusage[mem=16000] span[hosts=1]\" -n4 -q normal -- \\" >> $output
echo -e "\"bash $script_248 $replicates_NO_FLASH_filtered $MASTER_ROUTE NO_FLASH_FILTERED\"" >> $output

echo -e "bsub -G team151 -o $output_dir_filtered/Rscript_114.out -M 4000  -J Rscript_114_$counter_NO_FLASH_FILTERED -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
echo -e "\"$Rscript $script_114 --prefixes_EXPERIMENT $index_check_summary --prefixes_Eve $PREFIXES_Table_EXPANDED --route_fastq $route_fastq --out INDEXES_TABLE.txt\"" >> $output



echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output


#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################


output_dir_FLASH=$(echo "$MASTER_ROUTE""FLASH")
ToF_FLASH=$(echo "$output_dir_FLASH""/""ToF"".txt")
replicates_FLASH=$(echo "$output_dir_FLASH""/""replicates"".tsv")
counter_FLASH=$(echo "$parameter""_""$counter""_FLASH")


stderrfile_FLASH_filtered_113=$(echo "$output_dir_FLASH""/""Rscript_113.out")

touch $stderrfile_FLASH_filtered_113
echo -n "" > $stderrfile_FLASH_filtered_113

stderrfile_FLASH_filtered_114=$(echo "$output_dir_FLASH""/""Rscript_114.out")

touch $stderrfile_FLASH_filtered_114
echo -n "" > $stderrfile_FLASH_filtered_114

outfile_Script_248=$(echo "$output_dir_FLASH""/""Script_248.out")

touch $outfile_Script_248
echo -n "" > $outfile_Script_248




echo -e "cd $output_dir_FLASH/" >> $output
echo -e "bsub -G team151 -o $output_dir_FLASH/Rscript_113.out -M 4000  -J Rscript_113_$counter_FLASH -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
echo -e "\"$Rscript $script_113 --file $ToF_FLASH --fasta $REFERENCE --prefixes $PREFIXES_Table_EXPANDED --type FLASH\"" >> $output

echo -e "bsub -G team151 -o $output_dir_FLASH/Script_248.out -M 16000 -J Script_248_$counter_FLASH -R\"select[mem>=16000] rusage[mem=16000] span[hosts=1]\" -n4 -q normal -- \\" >> $output
echo -e "\"bash $script_248 $replicates_FLASH $MASTER_ROUTE FLASH\"" >> $output

echo -e "bsub -G team151 -o $output_dir_FLASH/Rscript_114.out -M 4000  -J Rscript_114_$counter_FLASH -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
echo -e "\"$Rscript $script_114 --prefixes_EXPERIMENT $index_check_summary --prefixes_Eve $PREFIXES_Table_EXPANDED --route_fastq $route_fastq --out INDEXES_TABLE.txt\"" >> $output


echo -e "\n\n"  >> $output

echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output


######################################################################################################################################################################################################################

output_dir_FLASH_filtered=$(echo "$MASTER_ROUTE""FLASH""_FILTERED")
ToF_FLASH_filtered=$(echo "$output_dir_FLASH_filtered""/""ToF_filtered"".txt")
replicates_FLASH_filtered=$(echo "$output_dir_FLASH_filtered""/""replicates_filtered"".tsv")
counter_FLASH_FILTERED=$(echo "$parameter""_""$batch""_FLASH_FILTERED")


stderrfile_FLASH_filtered_113=$(echo "$output_dir_FLASH_filtered""/""Rscript_113.out")

touch $stderrfile_FLASH_filtered_113
echo -n "" > $stderrfile_FLASH_filtered_113

stderrfile_FLASH_filtered_114=$(echo "$output_dir_FLASH_filtered""/""Rscript_114.out")

touch $stderrfile_FLASH_filtered_114
echo -n "" > $stderrfile_FLASH_filtered_114

outfile_Script_248=$(echo "$output_dir_FLASH_filtered""/""Script_248.out")

touch $outfile_Script_248
echo -n "" > $outfile_Script_248



echo -e "cd $output_dir_FLASH_filtered/" >> $output
echo -e "bsub -G team151 -o $output_dir_FLASH_filtered/Rscript_113.out -M 4000  -J Rscript_113_$counter_FLASH_FILTERED -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
echo -e "\"$Rscript $script_113 --file $ToF_FLASH_filtered --fasta $REFERENCE --prefixes $PREFIXES_Table_EXPANDED --type FLASH_FILTERED\"" >> $output

echo -e "bsub -G team151 -o $output_dir_FLASH_filtered/Script_248.out -M 16000 -J Script_248_$counter_FLASH_FILTERED -R\"select[mem>=16000] rusage[mem=16000] span[hosts=1]\" -n4 -q normal -- \\" >> $output
echo -e "\"bash $script_248 $replicates_FLASH_filtered $MASTER_ROUTE FLASH_FILTERED\"" >> $output

echo -e "bsub -G team151 -o $output_dir_FLASH_filtered/Rscript_114.out -M 4000  -J Rscript_114_$counter_FLASH_FILTERED -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
echo -e "\"$Rscript $script_114 --prefixes_EXPERIMENT $index_check_summary --prefixes_Eve $PREFIXES_Table_EXPANDED --route_fastq $route_fastq --out INDEXES_TABLE.txt\"" >> $output


echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output

echo "bash $output"






