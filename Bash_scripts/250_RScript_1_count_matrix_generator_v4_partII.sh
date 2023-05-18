#!/bin/bash
   
# /software/R-3.6.1/bin/Rscript
Rscript=/software/R-4.1.0/bin/Rscript

script_116=/nfs/users/nfs_m/mt19/Scripts/R/116_QC_proceser_v3.R
 
 
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


stderrfile_NO_FLASH_116=$(echo "$output_dir""/""Rscript_116.out")

touch $stderrfile_NO_FLASH_116
echo -n "" > $stderrfile_NO_FLASH_116

Count_NO_FLASH_matrix=$(echo "$MASTER_ROUTE""NO_FLASH""/""Count_NO_FLASH_matrix.txt")
NO_FLASH_CIGAR_mismatch_CONDENSED=$(echo "$MASTER_ROUTE""NO_FLASH_CIGAR_mismatch_CONDENSED.txt")
NO_FLASH_TOTAL_READS=$(echo "$MASTER_ROUTE""NO_FLASH_TOTAL_READS.txt")
DEMULTIPLEX_RESULT=$(echo "$MASTER_ROUTE""NO_FLASH""/""INDEXES_TABLE.txt")
out=$(echo  "$MASTER_ROUTE""NO_FLASH""/")

echo -e "bsub -G team151 -o $output_dir/Rscript_116.out -M 4000 -J Rscript_116_$counter_NO_FLASH -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
echo -e "\"$Rscript $script_116 --COUNT_MATRIX $Count_NO_FLASH_matrix --mismatch_CONDENSED $NO_FLASH_CIGAR_mismatch_CONDENSED --TOTAL_READS $NO_FLASH_TOTAL_READS --DEMULTIPLEX_RESULT $DEMULTIPLEX_RESULT --type NO_FLASH --out $out\"" >> $output

echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output

type=$(echo "Maximum_match""_""$counter_NO_FLASH")
outfile_Maximum_match=$(echo "$output_dir""/""outfile""_""$type"".out")
touch $outfile_Maximum_match
echo -n "" > $outfile_Maximum_match
name_Maximum_match=$(echo "$type""_job")


R_script_Maximum_match=/nfs/users/nfs_m/mt19/Scripts/R/116_QC_proceser_v3_Maximum_match.R

echo -e "bsub -G team151 -o $outfile_Maximum_match -M 4000 -w\"done(Rscript_116_$counter_NO_FLASH)\" -J $name_Maximum_match -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
echo -e "\"$Rscript $R_script_Maximum_match --COUNT_MATRIX $Count_NO_FLASH_matrix --mismatch_CONDENSED $NO_FLASH_CIGAR_mismatch_CONDENSED --TOTAL_READS $NO_FLASH_TOTAL_READS --DEMULTIPLEX_RESULT $DEMULTIPLEX_RESULT --type NO_FLASH --out $out\"" >> $output


echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output


type=$(echo "SOFT""_""$counter_NO_FLASH")
outfile_SOFT=$(echo "$output_dir""/""outfile""_""$type"".out")
touch $outfile_SOFT
echo -n "" > $outfile_SOFT
name_SOFT=$(echo "$type""_job")


R_script_SOFT=/nfs/users/nfs_m/mt19/Scripts/R/116_QC_proceser_v3_SOFT.R

echo -e "bsub -G team151 -o $outfile_SOFT -M 4000 -w\"done(Rscript_116_$counter_NO_FLASH)\" -J $name_SOFT -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
echo -e "\"$Rscript $R_script_SOFT --COUNT_MATRIX $Count_NO_FLASH_matrix --mismatch_CONDENSED $NO_FLASH_CIGAR_mismatch_CONDENSED --TOTAL_READS $NO_FLASH_TOTAL_READS --DEMULTIPLEX_RESULT $DEMULTIPLEX_RESULT --type NO_FLASH --out $out\"" >> $output

echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output


type=$(echo "HARD""_""$counter_NO_FLASH")
outfile_HARD=$(echo "$output_dir""/""outfile""_""$type"".out")
touch $outfile_HARD
echo -n "" > $outfile_HARD
name_HARD=$(echo "$type""_job")


R_script_HARD=/nfs/users/nfs_m/mt19/Scripts/R/116_QC_proceser_v3_HARD.R

echo -e "bsub -G team151 -o $outfile_HARD -M 4000 -w\"done(Rscript_116_$counter_NO_FLASH)\" -J $name_HARD -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
echo -e "\"$Rscript $R_script_HARD --COUNT_MATRIX $Count_NO_FLASH_matrix --mismatch_CONDENSED $NO_FLASH_CIGAR_mismatch_CONDENSED --TOTAL_READS $NO_FLASH_TOTAL_READS --DEMULTIPLEX_RESULT $DEMULTIPLEX_RESULT --type NO_FLASH --out $out\"" >> $output

echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output


type=$(echo "DEL""_""$counter_NO_FLASH")
outfile_DEL=$(echo "$output_dir""/""outfile""_""$type"".out")
touch $outfile_DEL
echo -n "" > $outfile_DEL
name_DEL=$(echo "$type""_job")


R_script_DEL=/nfs/users/nfs_m/mt19/Scripts/R/116_QC_proceser_v3_DEL.R

echo -e "bsub -G team151 -o $outfile_DEL -M 4000 -w\"done(Rscript_116_$counter_NO_FLASH)\" -J $name_DEL -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
echo -e "\"$Rscript $R_script_DEL --COUNT_MATRIX $Count_NO_FLASH_matrix --mismatch_CONDENSED $NO_FLASH_CIGAR_mismatch_CONDENSED --TOTAL_READS $NO_FLASH_TOTAL_READS --DEMULTIPLEX_RESULT $DEMULTIPLEX_RESULT --type NO_FLASH --out $out\"" >> $output

echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output


type=$(echo "INS""_""$counter_NO_FLASH")
outfile_INS=$(echo "$output_dir""/""outfile""_""$type"".out")
touch $outfile_INS
echo -n "" > $outfile_INS
name_INS=$(echo "$type""_job")


R_script_INS=/nfs/users/nfs_m/mt19/Scripts/R/116_QC_proceser_v3_INS.R

echo -e "bsub -G team151 -o $outfile_INS -M 4000 -w\"done(Rscript_116_$counter_NO_FLASH)\" -J $name_INS -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
echo -e "\"$Rscript $R_script_INS --COUNT_MATRIX $Count_NO_FLASH_matrix --mismatch_CONDENSED $NO_FLASH_CIGAR_mismatch_CONDENSED --TOTAL_READS $NO_FLASH_TOTAL_READS --DEMULTIPLEX_RESULT $DEMULTIPLEX_RESULT --type NO_FLASH --out $out\"" >> $output

echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output


type=$(echo "PRINTER""_""$counter_NO_FLASH")
outfile_PRINTER=$(echo "$output_dir""/""outfile""_""$type"".out")
touch $outfile_PRINTER
echo -n "" > $outfile_PRINTER
name_PRINTER=$(echo "$type""_job")


R_script_PRINTER=/nfs/users/nfs_m/mt19/Scripts/R/116_QC_proceser_v3_PRINTER.R

echo -e "bsub -G team151 -o $outfile_PRINTER -M 4000 -w\"done($name_INS) && done($name_DEL) && done($name_HARD) && done($name_SOFT) && done($name_Maximum_match)\" -J $name_PRINTER -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
echo -e "\"$Rscript $R_script_PRINTER --type NO_FLASH --out $out\"" >> $output


echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output




#######################################################################################################################################################################################################

output_dir_filtered=$(echo "$MASTER_ROUTE""NO_FLASH""_FILTERED")
ToF_NO_FLASH_filtered=$(echo "$output_dir_filtered""/""ToF_filtered"".txt")
replicates_NO_FLASH_filtered=$(echo "$output_dir_filtered""/""replicates_filtered"".tsv")
counter_NO_FLASH_FILTERED=$(echo "$parameter""_""$batch""_NO_FLASH_FILTERED")


stderrfile_NO_FLASH_filtered_116=$(echo "$output_dir_filtered""/""Rscript_116.out")

touch $stderrfile_NO_FLASH_filtered_116
echo -n "" > $stderrfile_NO_FLASH_filtered_116

stderrfile_NO_FLASH_filtered_117=$(echo "$output_dir_filtered""/""Rscript_117.out")



Count_NO_FLASH_FILTERED_matrix=$(echo "$MASTER_ROUTE""NO_FLASH_FILTERED""/""Count_NO_FLASH_FILTERED_matrix.txt")
NO_FLASH_FILTERED_CIGAR_mismatch_CONDENSED=$(echo "$MASTER_ROUTE""NO_FLASH_FILTERED_CIGAR_mismatch_CONDENSED.txt")
NO_FLASH_FILTERED_TOTAL_READS=$(echo "$MASTER_ROUTE""NO_FLASH_FILTERED_TOTAL_READS.txt")
DEMULTIPLEX_RESULT=$(echo "$MASTER_ROUTE""NO_FLASH_FILTERED""/""INDEXES_TABLE.txt")
out=$(echo  "$MASTER_ROUTE""NO_FLASH_FILTERED""/")

echo -e "bsub -G team151 -o $output_dir_filtered/Rscript_116.out -M 4000 -J Rscript_116_$counter_NO_FLASH_FILTERED -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
echo -e "\"$Rscript $script_116 --COUNT_MATRIX $Count_NO_FLASH_FILTERED_matrix --mismatch_CONDENSED $NO_FLASH_FILTERED_CIGAR_mismatch_CONDENSED --TOTAL_READS $NO_FLASH_FILTERED_TOTAL_READS --DEMULTIPLEX_RESULT $DEMULTIPLEX_RESULT --type NO_FLASH_FILTERED --out $out\"" >> $output

echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output

type=$(echo "Maximum_match""_""$counter_NO_FLASH_FILTERED")
outfile_Maximum_match=$(echo "$output_dir_filtered""/""outfile""_""$type"".out")
touch $outfile_Maximum_match
echo -n "" > $outfile_Maximum_match
name_Maximum_match=$(echo "$type""_job")

R_script_Maximum_match=/nfs/users/nfs_m/mt19/Scripts/R/116_QC_proceser_v3_Maximum_match.R

echo -e "bsub -G team151 -o $outfile_Maximum_match -M 4000 -w\"done(Rscript_116_$counter_NO_FLASH_FILTERED)\" -J $name_Maximum_match -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
echo -e "\"$Rscript $R_script_Maximum_match --COUNT_MATRIX $Count_NO_FLASH_FILTERED_matrix --mismatch_CONDENSED $NO_FLASH_FILTERED_CIGAR_mismatch_CONDENSED --TOTAL_READS $NO_FLASH_FILTERED_TOTAL_READS --DEMULTIPLEX_RESULT $DEMULTIPLEX_RESULT --type NO_FLASH_FILTERED --out $out\"" >> $output

echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output


type=$(echo "SOFT""_""$counter_NO_FLASH_FILTERED")
outfile_SOFT=$(echo "$output_dir_filtered""/""outfile""_""$type"".out")
touch $outfile_SOFT
echo -n "" > $outfile_SOFT
name_SOFT=$(echo "$type""_job")


R_script_SOFT=/nfs/users/nfs_m/mt19/Scripts/R/116_QC_proceser_v3_SOFT.R

echo -e "bsub -G team151 -o $outfile_SOFT -M 4000 -w\"done(Rscript_116_$counter_NO_FLASH_FILTERED)\" -J $name_SOFT -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
echo -e "\"$Rscript $R_script_SOFT --COUNT_MATRIX $Count_NO_FLASH_FILTERED_matrix --mismatch_CONDENSED $NO_FLASH_FILTERED_CIGAR_mismatch_CONDENSED --TOTAL_READS $NO_FLASH_FILTERED_TOTAL_READS --DEMULTIPLEX_RESULT $DEMULTIPLEX_RESULT --type NO_FLASH_FILTERED --out $out\"" >> $output

echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output


type=$(echo "HARD""_""$counter_NO_FLASH_FILTERED")
outfile_HARD=$(echo "$output_dir_filtered""/""outfile""_""$type"".out")
touch $outfile_HARD
echo -n "" > $outfile_HARD
name_HARD=$(echo "$type""_job")


R_script_HARD=/nfs/users/nfs_m/mt19/Scripts/R/116_QC_proceser_v3_HARD.R

echo -e "bsub -G team151 -o $outfile_HARD -M 4000 -w\"done(Rscript_116_$counter_NO_FLASH_FILTERED)\" -J $name_HARD -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
echo -e "\"$Rscript $R_script_HARD --COUNT_MATRIX $Count_NO_FLASH_FILTERED_matrix --mismatch_CONDENSED $NO_FLASH_FILTERED_CIGAR_mismatch_CONDENSED --TOTAL_READS $NO_FLASH_FILTERED_TOTAL_READS --DEMULTIPLEX_RESULT $DEMULTIPLEX_RESULT --type NO_FLASH_FILTERED --out $out\"" >> $output

echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output


type=$(echo "DEL""_""$counter_NO_FLASH_FILTERED")
outfile_DEL=$(echo "$output_dir_filtered""/""outfile""_""$type"".out")
touch $outfile_DEL
echo -n "" > $outfile_DEL
name_DEL=$(echo "$type""_job")


R_script_DEL=/nfs/users/nfs_m/mt19/Scripts/R/116_QC_proceser_v3_DEL.R

echo -e "bsub -G team151 -o $outfile_DEL -M 4000 -w\"done(Rscript_116_$counter_NO_FLASH_FILTERED)\" -J $name_DEL -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
echo -e "\"$Rscript $R_script_DEL --COUNT_MATRIX $Count_NO_FLASH_FILTERED_matrix --mismatch_CONDENSED $NO_FLASH_FILTERED_CIGAR_mismatch_CONDENSED --TOTAL_READS $NO_FLASH_FILTERED_TOTAL_READS --DEMULTIPLEX_RESULT $DEMULTIPLEX_RESULT --type NO_FLASH_FILTERED --out $out\"" >> $output

echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output


type=$(echo "INS""_""$counter_NO_FLASH_FILTERED")
outfile_INS=$(echo "$output_dir_filtered""/""outfile""_""$type"".out")
touch $outfile_INS
echo -n "" > $outfile_INS
name_INS=$(echo "$type""_job")


R_script_INS=/nfs/users/nfs_m/mt19/Scripts/R/116_QC_proceser_v3_INS.R

echo -e "bsub -G team151 -o $outfile_INS -M 4000 -w\"done(Rscript_116_$counter_NO_FLASH_FILTERED)\" -J $name_INS -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
echo -e "\"$Rscript $R_script_INS --COUNT_MATRIX $Count_NO_FLASH_FILTERED_matrix --mismatch_CONDENSED $NO_FLASH_FILTERED_CIGAR_mismatch_CONDENSED --TOTAL_READS $NO_FLASH_FILTERED_TOTAL_READS --DEMULTIPLEX_RESULT $DEMULTIPLEX_RESULT --type NO_FLASH_FILTERED --out $out\"" >> $output

echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output


type=$(echo "PRINTER""_""$counter_NO_FLASH_FILTERED")
outfile_PRINTER=$(echo "$output_dir_filtered""/""outfile""_""$type"".out")
touch $outfile_PRINTER
echo -n "" > $outfile_PRINTER
name_PRINTER=$(echo "$type""_job")


R_script_PRINTER=/nfs/users/nfs_m/mt19/Scripts/R/116_QC_proceser_v3_PRINTER.R

echo -e "bsub -G team151 -o $outfile_PRINTER -M 4000 -w\"done($name_INS) && done($name_DEL) && done($name_HARD) && done($name_SOFT) && done($name_Maximum_match)\" -J $name_PRINTER -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
echo -e "\"$Rscript $R_script_PRINTER --type NO_FLASH_FILTERED --out $out\"" >> $output


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
counter_FLASH=$(echo "$parameter""_""$batch""_FLASH")


stderrfile_FLASH_filtered_116=$(echo "$output_dir_FLASH""/""Rscript_116.out")

touch $stderrfile_FLASH_filtered_116
echo -n "" > $stderrfile_FLASH_filtered_116



Count_FLASH_matrix=$(echo "$MASTER_ROUTE""FLASH""/""Count_FLASH_matrix.txt")
FLASH_CIGAR_mismatch_CONDENSED=$(echo "$MASTER_ROUTE""FLASH_CIGAR_mismatch_CONDENSED.txt")
FLASH_TOTAL_READS=$(echo "$MASTER_ROUTE""FLASH_TOTAL_READS.txt")
DEMULTIPLEX_RESULT=$(echo "$MASTER_ROUTE""FLASH""/""INDEXES_TABLE.txt")
out=$(echo "$MASTER_ROUTE""FLASH""/")

echo -e "bsub -G team151 -o $output_dir_FLASH/Rscript_116.out -M 4000 -J Rscript_116_$counter_FLASH -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
echo -e "\"$Rscript $script_116 --COUNT_MATRIX $Count_FLASH_matrix --mismatch_CONDENSED $FLASH_CIGAR_mismatch_CONDENSED --TOTAL_READS $FLASH_TOTAL_READS --DEMULTIPLEX_RESULT $DEMULTIPLEX_RESULT --type FLASH --out $out\"" >> $output

echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output

type=$(echo "Maximum_match""_""$counter_FLASH")
outfile_Maximum_match=$(echo "$output_dir_FLASH""/""outfile""_""$type"".out")
touch $outfile_Maximum_match
echo -n "" > $outfile_Maximum_match
name_Maximum_match=$(echo "$type""_job")

R_script_Maximum_match=/nfs/users/nfs_m/mt19/Scripts/R/116_QC_proceser_v3_Maximum_match.R

echo -e "bsub -G team151 -o $outfile_Maximum_match -M 4000 -w\"done(Rscript_116_$counter_FLASH)\" -J $name_Maximum_match -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
echo -e "\"$Rscript $R_script_Maximum_match --COUNT_MATRIX $Count_FLASH_matrix --mismatch_CONDENSED $FLASH_CIGAR_mismatch_CONDENSED --TOTAL_READS $FLASH_TOTAL_READS --DEMULTIPLEX_RESULT $DEMULTIPLEX_RESULT --type FLASH --out $out\"" >> $output

echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output


type=$(echo "SOFT""_""$counter_FLASH")
outfile_SOFT=$(echo "$output_dir_FLASH""/""outfile""_""$type"".out")
touch $outfile_SOFT
echo -n "" > $outfile_SOFT
name_SOFT=$(echo "$type""_job")


R_script_SOFT=/nfs/users/nfs_m/mt19/Scripts/R/116_QC_proceser_v3_SOFT.R

echo -e "bsub -G team151 -o $outfile_SOFT -M 4000 -w\"done(Rscript_116_$counter_FLASH)\" -J $name_SOFT -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
echo -e "\"$Rscript $R_script_SOFT --COUNT_MATRIX $Count_FLASH_matrix --mismatch_CONDENSED $FLASH_CIGAR_mismatch_CONDENSED --TOTAL_READS $FLASH_TOTAL_READS --DEMULTIPLEX_RESULT $DEMULTIPLEX_RESULT --type FLASH --out $out\"" >> $output

echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output


type=$(echo "HARD""_""$counter_FLASH")
outfile_HARD=$(echo "$output_dir_FLASH""/""outfile""_""$type"".out")
touch $outfile_HARD
echo -n "" > $outfile_HARD
name_HARD=$(echo "$type""_job")


R_script_HARD=/nfs/users/nfs_m/mt19/Scripts/R/116_QC_proceser_v3_HARD.R

echo -e "bsub -G team151 -o $outfile_HARD -M 4000 -w\"done(Rscript_116_$counter_FLASH)\" -J $name_HARD -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
echo -e "\"$Rscript $R_script_HARD --COUNT_MATRIX $Count_FLASH_matrix --mismatch_CONDENSED $FLASH_CIGAR_mismatch_CONDENSED --TOTAL_READS $FLASH_TOTAL_READS --DEMULTIPLEX_RESULT $DEMULTIPLEX_RESULT --type FLASH --out $out\"" >> $output

echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output


type=$(echo "DEL""_""$counter_FLASH")
outfile_DEL=$(echo "$output_dir_FLASH""/""outfile""_""$type"".out")
touch $outfile_DEL
echo -n "" > $outfile_DEL
name_DEL=$(echo "$type""_job")


R_script_DEL=/nfs/users/nfs_m/mt19/Scripts/R/116_QC_proceser_v3_DEL.R

echo -e "bsub -G team151 -o $outfile_DEL -M 4000 -w\"done(Rscript_116_$counter_FLASH)\" -J $name_DEL -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
echo -e "\"$Rscript $R_script_DEL --COUNT_MATRIX $Count_FLASH_matrix --mismatch_CONDENSED $FLASH_CIGAR_mismatch_CONDENSED --TOTAL_READS $FLASH_TOTAL_READS --DEMULTIPLEX_RESULT $DEMULTIPLEX_RESULT --type FLASH --out $out\"" >> $output

echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output


type=$(echo "INS""_""$counter_FLASH")
outfile_INS=$(echo "$output_dir_FLASH""/""outfile""_""$type"".out")
touch $outfile_INS
echo -n "" > $outfile_INS
name_INS=$(echo "$type""_job")


R_script_INS=/nfs/users/nfs_m/mt19/Scripts/R/116_QC_proceser_v3_INS.R

echo -e "bsub -G team151 -o $outfile_INS -M 4000 -w\"done(Rscript_116_$counter_FLASH)\" -J $name_INS -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
echo -e "\"$Rscript $R_script_INS --COUNT_MATRIX $Count_FLASH_matrix --mismatch_CONDENSED $FLASH_CIGAR_mismatch_CONDENSED --TOTAL_READS $FLASH_TOTAL_READS --DEMULTIPLEX_RESULT $DEMULTIPLEX_RESULT --type FLASH --out $out\"" >> $output

echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output


type=$(echo "PRINTER""_""$counter_FLASH")
outfile_PRINTER=$(echo "$output_dir_FLASH""/""outfile""_""$type"".out")
touch $outfile_PRINTER
echo -n "" > $outfile_PRINTER
name_PRINTER=$(echo "$type""_job")


R_script_PRINTER=/nfs/users/nfs_m/mt19/Scripts/R/116_QC_proceser_v3_PRINTER.R

echo -e "bsub -G team151 -o $outfile_PRINTER -M 4000 -w\"done($name_INS) && done($name_DEL) && done($name_HARD) && done($name_SOFT) && done($name_Maximum_match)\" -J $name_PRINTER -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
echo -e "\"$Rscript $R_script_PRINTER --type FLASH --out $out\"" >> $output


echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output


######################################################################################################################################################################################################################

output_dir_FLASH_filtered=$(echo "$MASTER_ROUTE""FLASH""_FILTERED")
ToF_FLASH_filtered=$(echo "$output_dir_FLASH_filtered""/""ToF_filtered"".txt")
replicates_FLASH_filtered=$(echo "$output_dir_FLASH_filtered""/""replicates_filtered"".tsv")
counter_FLASH_FILTERED=$(echo "$parameter""_""$batch""_FLASH_FILTERED")


stderrfile_FLASH_filtered_116=$(echo "$output_dir_FLASH_filtered""/""Rscript_116.out")

touch $stderrfile_FLASH_filtered_116
echo -n "" > $stderrfile_FLASH_filtered_116

 

Count_FLASH_FILTERED_matrix=$(echo "$MASTER_ROUTE""FLASH_FILTERED""/""Count_FLASH_FILTERED_matrix.txt")
FLASH_FILTERED_CIGAR_mismatch_CONDENSED=$(echo "$MASTER_ROUTE""FLASH_FILTERED_CIGAR_mismatch_CONDENSED.txt")
FLASH_FILTERED_TOTAL_READS=$(echo "$MASTER_ROUTE""FLASH_FILTERED_TOTAL_READS.txt")
DEMULTIPLEX_RESULT=$(echo "$MASTER_ROUTE""FLASH_FILTERED""/""INDEXES_TABLE.txt")

out=$(echo "$MASTER_ROUTE""FLASH_FILTERED""/")

echo -e "bsub -G team151 -o $output_dir_FLASH_filtered/Rscript_116.out -M 4000  -J Rscript_116_$counter_FLASH_FILTERED -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
echo -e "\"$Rscript $script_116 --COUNT_MATRIX $Count_FLASH_FILTERED_matrix --mismatch_CONDENSED $FLASH_FILTERED_CIGAR_mismatch_CONDENSED --TOTAL_READS $FLASH_FILTERED_TOTAL_READS --DEMULTIPLEX_RESULT $DEMULTIPLEX_RESULT --type FLASH_FILTERED --out $out\"" >> $output

echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output

type=$(echo "Maximum_match""_""$counter_FLASH_FILTERED")
outfile_Maximum_match=$(echo "$output_dir_FLASH_filtered""/""outfile""_""$type"".out")
touch $outfile_Maximum_match
echo -n "" > $outfile_Maximum_match
name_Maximum_match=$(echo "$type""_job")

R_script_Maximum_match=/nfs/users/nfs_m/mt19/Scripts/R/116_QC_proceser_v3_Maximum_match.R

echo -e "bsub -G team151 -o $outfile_Maximum_match -M 4000 -w\"done(Rscript_116_$counter_FLASH_FILTERED)\" -J $name_Maximum_match -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
echo -e "\"$Rscript $R_script_Maximum_match --COUNT_MATRIX $Count_FLASH_FILTERED_matrix --mismatch_CONDENSED $FLASH_FILTERED_CIGAR_mismatch_CONDENSED --TOTAL_READS $FLASH_FILTERED_TOTAL_READS --DEMULTIPLEX_RESULT $DEMULTIPLEX_RESULT --type FLASH_FILTERED --out $out\"" >> $output

echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output

type=$(echo "SOFT""_""$counter_FLASH_FILTERED")
outfile_SOFT=$(echo "$output_dir_FLASH_filtered""/""outfile""_""$type"".out")
touch $outfile_SOFT
echo -n "" > $outfile_SOFT
name_SOFT=$(echo "$type""_job")


R_script_SOFT=/nfs/users/nfs_m/mt19/Scripts/R/116_QC_proceser_v3_SOFT.R

echo -e "bsub -G team151 -o $outfile_SOFT -M 4000 -w\"done(Rscript_116_$counter_FLASH_FILTERED)\" -J $name_SOFT -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
echo -e "\"$Rscript $R_script_SOFT --COUNT_MATRIX $Count_FLASH_FILTERED_matrix --mismatch_CONDENSED $FLASH_FILTERED_CIGAR_mismatch_CONDENSED --TOTAL_READS $FLASH_FILTERED_TOTAL_READS --DEMULTIPLEX_RESULT $DEMULTIPLEX_RESULT --type FLASH_FILTERED --out $out\"" >> $output

echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output


type=$(echo "HARD""_""$counter_FLASH_FILTERED")
outfile_HARD=$(echo "$output_dir_FLASH_filtered""/""outfile""_""$type"".out")
touch $outfile_HARD
echo -n "" > $outfile_HARD
name_HARD=$(echo "$type""_job")


R_script_HARD=/nfs/users/nfs_m/mt19/Scripts/R/116_QC_proceser_v3_HARD.R

echo -e "bsub -G team151 -o $outfile_HARD -M 4000 -w\"done(Rscript_116_$counter_FLASH_FILTERED)\" -J $name_HARD -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
echo -e "\"$Rscript $R_script_HARD --COUNT_MATRIX $Count_FLASH_FILTERED_matrix --mismatch_CONDENSED $FLASH_FILTERED_CIGAR_mismatch_CONDENSED --TOTAL_READS $FLASH_FILTERED_TOTAL_READS --DEMULTIPLEX_RESULT $DEMULTIPLEX_RESULT --type FLASH_FILTERED --out $out\"" >> $output

echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output


type=$(echo "DEL""_""$counter_FLASH_FILTERED")
outfile_DEL=$(echo "$output_dir_FLASH_filtered""/""outfile""_""$type"".out")
touch $outfile_DEL
echo -n "" > $outfile_DEL
name_DEL=$(echo "$type""_job")


R_script_DEL=/nfs/users/nfs_m/mt19/Scripts/R/116_QC_proceser_v3_DEL.R

echo -e "bsub -G team151 -o $outfile_DEL -M 4000 -w\"done(Rscript_116_$counter_FLASH_FILTERED)\" -J $name_DEL -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
echo -e "\"$Rscript $R_script_DEL --COUNT_MATRIX $Count_FLASH_FILTERED_matrix --mismatch_CONDENSED $FLASH_FILTERED_CIGAR_mismatch_CONDENSED --TOTAL_READS $FLASH_FILTERED_TOTAL_READS --DEMULTIPLEX_RESULT $DEMULTIPLEX_RESULT --type FLASH_FILTERED --out $out\"" >> $output

echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output


type=$(echo "INS""_""$counter_FLASH_FILTERED")
outfile_INS=$(echo "$output_dir_FLASH_filtered""/""outfile""_""$type"".out")
touch $outfile_INS
echo -n "" > $outfile_INS
name_INS=$(echo "$type""_job")


R_script_INS=/nfs/users/nfs_m/mt19/Scripts/R/116_QC_proceser_v3_INS.R

echo -e "bsub -G team151 -o $outfile_INS -M 4000 -w\"done(Rscript_116_$counter_FLASH_FILTERED)\" -J $name_INS -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
echo -e "\"$Rscript $R_script_INS --COUNT_MATRIX $Count_FLASH_FILTERED_matrix --mismatch_CONDENSED $FLASH_FILTERED_CIGAR_mismatch_CONDENSED --TOTAL_READS $FLASH_FILTERED_TOTAL_READS --DEMULTIPLEX_RESULT $DEMULTIPLEX_RESULT --type FLASH_FILTERED --out $out\"" >> $output

echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output


type=$(echo "PRINTER""_""$counter_FLASH_FILTERED")
outfile_PRINTER=$(echo "$output_dir_FLASH_filtered""/""outfile""_""$type"".out")
touch $outfile_PRINTER
echo -n "" > $outfile_PRINTER
name_PRINTER=$(echo "$type""_job")


R_script_PRINTER=/nfs/users/nfs_m/mt19/Scripts/R/116_QC_proceser_v3_PRINTER.R

echo -e "bsub -G team151 -o $outfile_PRINTER -M 4000 -w\"done($name_INS) && done($name_DEL) && done($name_HARD) && done($name_SOFT) && done($name_Maximum_match)\" -J $name_PRINTER -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
echo -e "\"$Rscript $R_script_PRINTER --type FLASH_FILTERED --out $out\"" >> $output

echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output
echo -e "\n\n"  >> $output

echo "bash $output"






