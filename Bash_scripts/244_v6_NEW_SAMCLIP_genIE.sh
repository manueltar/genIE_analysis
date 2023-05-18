#!/bin/bash
  

################################################################################################################################## Dependencies: 

#samtools=/nfs/users/nfs_m/mt19/sOFTWARE/samtools-1.6/bin/samtools
python=/usr/bin/python2.7
parser_flash_output=/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/src/getFlashOutputDetails.py
barcode_checker=/nfs/users/nfs_m/mt19/Scripts/Wraper_scripts/247_indexes_fastq_check_v2.sh



#REFERENCE=/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/ReSeq/genIE.fasta
#REFERENCE=/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/ReSeq/genIE_NO_haplotype2.fasta
REFERENCE=$2

MASTER_ROUTE=$3
# /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/ReSeq/BWA/

#rm -rf $MASTER_ROUTE/
#mkdir -p $MASTER_ROUTE/

 
output=$4
#/nfs/users/nfs_m/mt19/Scripts/Wraper_scripts/200_CREATION.sh

parameter=$6

FASTQ_ROUTE=$7

 
touch $output
echo -n "" > $output

echo -e "#!/bin/bash"  >> $output
echo -e "\n\n"  >> $output





##################### Parser of barcodes

path=$(pwd)

input_file=$(echo "$path""/""$1")

stderr=$(echo "$MASTER_ROUTE""barcode_checker.out")
#echo "--->$input_file"

echo -e "bsub -G team151 -o $stderr -M 4000 -J barcode_checker -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
echo -e "\"bash $barcode_checker $input_file $MASTER_ROUTE\"" >> $output


# #############################################################################################
 




output_dir=$(echo "$MASTER_ROUTE""NO_FLASH")
output_dir_filtered=$(echo "$MASTER_ROUTE""NO_FLASH""_FILTERED")

 # rm -rf $output_dir/
  #mkdir -p $output_dir/

#  rm -rf $output_dir_filtered/
#  mkdir -p $output_dir_filtered/


ToF_NO_FLASH=$(echo "$output_dir""/""ToF"".txt")

echo $ToF_NO_FLASH


touch $ToF_NO_FLASH
echo -n "" > $ToF_NO_FLASH


ToF_NO_FLASH_filtered=$(echo "$output_dir_filtered""/""ToF_filtered"".txt")

echo $ToF_NO_FLASH_filtered


touch $ToF_NO_FLASH_filtered
echo -n "" > $ToF_NO_FLASH_filtered


replicates_NO_FLASH=$(echo "$output_dir""/""replicates"".tsv")

echo $replicates_NO_FLASH


touch $replicates_NO_FLASH
echo -n "" > $replicates_NO_FLASH

echo -e "name\treplicate\ttype\tbam"  >> $replicates_NO_FLASH

replicates_NO_FLASH_filtered=$(echo "$output_dir_filtered""/""replicates_filtered"".tsv")

echo $replicates_NO_FLASH_filtered


touch $replicates_NO_FLASH_filtered
echo -n "" > $replicates_NO_FLASH_filtered

echo -e "name\treplicate\ttype\tbam"  >> $replicates_NO_FLASH_filtered

end_of_file=0

echo "FILE------------------->:$1"

while [[ $end_of_file == 0 ]]
do
  read -r line
  end_of_file=$?

echo "LINE:$line"

if [[ "$line" =~ [^[:space:]] ]]; then

    
    a=($(echo "$line" | tr "\t" '\n'))
    sample=${a[2]}
    Read_1=${a[9]}
    Read_2=${a[10]}
    replicate=${a[2]}
    type=${a[6]}


    echo "DEF:Sample: $sample\t$Read_1\t$Read_2\n"

#    exit
    
    if [ $sample != "replicate_full" ]; then # Don't use the header line NOT very elegant


    b=($(echo "$sample" | tr "_" '\n'))
    name=${b[0]}
    # replicate=${b[1]}
    # letter=$(echo "$replicate"|sed -r 's/[0-9]+\-[0-9]+\-//g')
    # letter=$(echo "$letter"|sed -r 's/[0-9]+/DNA/g')
    # replicate=sed -e 's/[0-9]+//g'


#	exit
     # if [ $name == "EROS" ]; then

     #     Read_1=${a[8]}
     # 	 Read_2=${a[9]}

     # fi

	echo "DEF:Sample: $sample"
	echo "name: $name"
	echo "replicate: $replicate"
	echo "type: $type"


	echo "Read_1: $Read_1"
	echo "Read_2: $Read_2"

#        exit

     if [ $replicate != "NA" ]; then # Don't inexisting alignments

     # if [ $name == "EROS" ]; then
 
#exit

     # fi


     OUTPUT_BWA=$output_dir/$sample.sam
     read_group_id=$(echo "'@RG"'\\t'"ID:1"'\\t'"SM:$sample"'\\t'"PL:Illumina"'\\t'"LB:1"'\\t'"PU:1'")

     echo $read_group_id
     echo "------------------------------------------------------------------------------------->$OUTPUT_BWA"
 
     counter=$(echo "$sample""_""$parameter")

 
     echo -e "bsub -G team151 -o $output_dir/BWA_$counter.out -M 8000  -J BWA_$counter -R\"select[mem>=8000] rusage[mem=8000] span[hosts=1]\" -n2 -q normal -- \\" >> $output
     echo -e "\"bwa mem -M -R \\" >> $output
     echo -e "$read_group_id \\" >> $output
     echo -e "-t 2 \\" >> $output
     echo -e "'-O 24,48 -E 1 -A 4 -B 16 -T 70 -k 19 -w 200 -d 600 -L 20 -U 40' \\" >> $output
     echo -e "$REFERENCE \\" >> $output
     echo -e "$Read_1 \\" >> $output
     echo -e "$Read_2 > $OUTPUT_BWA\""  >> $output

 
     # step 2 convert sam to bam

     OUTPUT_BWA_bam=$output_dir/$sample.bam

     echo -e "bsub -G team151 -o $output_dir/BWA_$counter.out -M 4000 -w\"done(BWA_$counter)\" -J SAM_TO_BAM_$counter -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
     echo -e "\"samtools view -b $OUTPUT_BWA >  $OUTPUT_BWA_bam\"" >> $output

     # step 2.5 clear SAM

     echo -e "bsub -G team151 -o $output_dir/BWA_$counter.out -M 4000 -w\"done(SAM_TO_BAM_$counter)\" -J clear_1_$counter -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
     echo -e "\"rm $OUTPUT_BWA\"" >> $output

     # step 3 sort bam by coord

      OUTPUT_BWA_bam_sorted=$(echo "$output_dir""/""$sample""_sorted"".bam")

     echo -e "bsub -G team151 -o $output_dir/BWA_$counter.out -M 4000 -w\"done(clear_1_$counter)\" -J SORT_BAM_$counter -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
     echo -e "\"samtools sort -o $OUTPUT_BWA_bam_sorted $OUTPUT_BWA_bam\"" >> $output

     # step 4 INDEX

    echo -e "bsub -G team151 -o $output_dir/BWA_$counter.out -M 4000 -w\"done(SORT_BAM_$counter)\" -J INDEX_$counter -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
    echo -e "\"samtools index $OUTPUT_BWA_bam_sorted\"" >> $output

     echo -e "$name\t$sample\t$type\t$OUTPUT_BWA_bam_sorted" >> $replicates_NO_FLASH
 
     # step 5 clear 2

     echo -e "bsub -G team151 -o $output_dir/BWA_$counter.out -M 4000 -w\"done(INDEX_$counter)\" -J clear_2_$counter -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
     echo -e "\"rm $OUTPUT_BWA_bam\"" >> $output

     # step 6 filter soft clipped bigger than XXX

     bash_samclip=/nfs/users/nfs_m/mt19/Scripts/Wraper_scripts/244_v6_sub_bash.sh


     echo -e "bsub -G team151 -o $output_dir_filtered/BWA_$counter.out -M 4000 -w\"done(INDEX_$counter)\" -J samclip_$counter -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
#     echo -e "bsub -G team151 -o $output_dir_filtered/BWA_$counter.out -M 4000  -J samclip_$counter -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
     echo -e "\"bash $bash_samclip $output_dir_filtered $OUTPUT_BWA_bam_sorted $REFERENCE $parameter $sample 4000 1 normal $ToF_NO_FLASH_filtered $replicates_NO_FLASH_filtered $type $name\"" >> $output



     
     # step 11 samtools cgr counts

     OUTPUT_counts=$(echo "$output_dir""/""$sample""_count"".txt")

     echo  $OUTPUT_counts >> $ToF_NO_FLASH


     echo -e "bsub -G team151 -o $output_dir/count_$counter.out -M 4000 -w\"done(clear_2_$counter)\" -J count_$counter -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
#     echo -e "bsub -G team151 -o $output_dir/count_$counter.out -M 4000 -J count_$counter -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
     echo -e "\"samtools view $OUTPUT_BWA_bam_sorted |cut -f3|sort|uniq -c|sort -rg > $OUTPUT_counts\"" >> $output

     echo -e "\n\n"  >> $output
     echo -e "\n\n"  >> $output

#     exit
     
     
    fi
    fi


fi

done < "$1"


#exit

#### Here is the place for $4 the flash input file


output_dir_FLASH=$(echo "$MASTER_ROUTE""FLASH")
output_dir_FLASH_filtered=$(echo "$MASTER_ROUTE""FLASH""_FILTERED")

 rm -rf $output_dir_FLASH/
 mkdir -p $output_dir_FLASH/

echo -e "cd $output_dir_FLASH/" >> $output
echo -e "echo -e \"File name\tOutput name\tTotal pairs\tCombined pairs\tUncombined pairs\tPercent combined\tMin overlap\tMax overlap\tMax mismatch dens\tWarning\" > flash_output.summary.tsv" >> $output


 rm -rf $output_dir_FLASH_filtered/
 mkdir -p $output_dir_FLASH_filtered/

echo -e "cd $output_dir_FLASH_filtered/" >> $output
echo -e "echo -e \"File name\tOutput name\tTotal pairs\tCombined pairs\tUncombined pairs\tPercent combined\tMin overlap\tMax overlap\tMax mismatch dens\tWarning\" > flash_output.summary.tsv" >> $output

echo -e "cd $FASTQ_ROUTE" >> $output

ToF_FLASH=$(echo "$output_dir_FLASH""/""ToF"".txt")

echo $ToF_FLASH


touch $ToF_FLASH
echo -n "" > $ToF_FLASH


ToF_FLASH_filtered=$(echo "$output_dir_FLASH_filtered""/""ToF_filtered"".txt")

echo $ToF_FLASH_filtered


touch $ToF_FLASH_filtered
echo -n "" > $ToF_FLASH_filtered


replicates_FLASH=$(echo "$output_dir_FLASH""/""replicates"".tsv")

echo $replicates_FLASH


touch $replicates_FLASH
echo -n "" > $replicates_FLASH

echo -e "name\treplicate\ttype\tbam"  >> $replicates_FLASH

replicates_FLASH_filtered=$(echo "$output_dir_FLASH_filtered""/""replicates_filtered"".tsv")

echo $replicates_FLASH_filtered


touch $replicates_FLASH_filtered
echo -n "" > $replicates_FLASH_filtered

echo -e "name\treplicate\ttype\tbam"  >> $replicates_FLASH_filtered


end_of_file2=0

while [[ $end_of_file2 == 0 ]]
do
  read -r line
  end_of_file2=$?

#echo "LINE:$line"

if [[ "$line" =~ [^[:space:]] ]]; then

    a=($(echo "$line" | tr "\t" '\n')) # check
    sample=${a[0]} # check
    Read_1=${a[1]} # check
    Read_2=${a[2]} #check
    read_length=${a[3]}
    read_length_string=$(echo "-r ""$read_length")
    fragment_size=${a[4]}
    fragment_size_string=$(echo "-f ""$fragment_size")
    fragment_sd=${a[5]}
    fragment_sd_string=$(echo "-s ""$fragment_sd")
    min_overlap=${a[6]}
    min_overlap_string=$(echo "-m ""$min_overlap")
    
    max_overlap_string=$(echo "-M ""65")


    max_mismatch_dens=${a[7]}
    max_mismatch_dens_string=$(echo "-x ""$max_mismatch_dens")
    flash_output=${a[8]}

 
    check=$(echo "$sample"|sed -r 's/_.+//g')


    # STEP 0 ELIMINATE CARRIAGE RETURN FROM R!!!!

    flash_output=($(echo "$flash_output" | tr "\r\t" '\n'))


    if [ $sample != "replicate_full" ]; then # Don't use the header line NOT very elegant


	name=$(echo "$sample"|sed -r 's/_.+$//g')
	echo "name:$name"
	replicate=$(echo "$sample"|sed -r 's/^[^_]+_[^_]+_//g')
	echo "replicate:$replicate"
#	type=$(echo "$replicate"|sed -r 's/\[0-9\]+//g')
	type=$(echo "$replicate"|sed -r 's/C.+/cDNA/g')

	echo "type:$type"

#	exit
	
#	type=$(echo "$replicate"|sed -r 's/G/gDNA/g')


	echo "------------>$line"

	if [ $type != "cDNA" ];
	then

	    type=$(echo "$replicate"|sed -r 's/G.+/gDNA/g')
    	    echo "type:$type"
#           exit

	fi


	
     
    if [ $Read_1 != "NA" ]; then 

    counter=$(echo "$sample""_""$parameter")

    # if [ $name == "EROS" ]; then 

#    	exit
    # fi

   echo "-----------$check----------------->$sample\t$type\t$Read_1\t$Read_2"


 #    step 1 FLASH convert the fq files

    
    if [ $check == "SAV1" ];
    then

	min_overlap_string=$(echo "-m ""1")
	fragment_size_string=$(echo "-f ""274")
	max_overlap_string=$(echo "-M ""100")

	

	
	     echo -e "bsub -G team151 -o $output_dir_FLASH/FLASH_$counter.out -M 4000 -J FLASH_$counter -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
	     echo -e "\"flash2 -z --allow-outies --output-prefix=$sample $read_length_string $fragment_size_string $min_overlap_string $max_overlap_string $Read_1 $Read_2\"" >> $output

#	     exit

     else

     echo -e "bsub -G team151 -o $output_dir_FLASH/FLASH_$counter.out -M 4000 -J FLASH_$counter -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
     echo -e "\"flash2 -z --allow-outies --output-prefix=$sample $read_length_string $fragment_size_string $fragment_sd_string $min_overlap_string $max_mismatch_dens_string $Read_1 $Read_2\"" >> $output	 

    fi

 

     echo -e "bsub -G team151 -o $output_dir_FLASH/FLASH_parse_$counter.out -M 4000 -w\"done(FLASH_$counter)\" -J FLASH_parse_$counter -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
     echo -e "\"$python $parser_flash_output --file $output_dir_FLASH/FLASH_$counter.out|sed '1d' >> flash_output.summary.tsv\"" >> $output

	

     ### STEP 2  BWA 

     OUTPUT_BWA=$output_dir_FLASH/$sample.sam
     read_group_id=$(echo "'@RG"'\\t'"ID:1"'\\t'"SM:$sample"'\\t'"PL:Illumina"'\\t'"LB:1"'\\t'"PU:1'")

#     echo $read_group_id


     echo -e "bsub -G team151 -o $output_dir_FLASH/BWA_$counter.out -M 8000 -w\"done(FLASH_$counter)\" -J BWA_FLASH_$counter -R\"select[mem>=8000] rusage[mem=8000] span[hosts=1]\" -n2 -q normal -- \\" >> $output
     echo -e "\"bwa mem -M -R \\" >> $output
     echo -e "$read_group_id \\" >> $output
     echo -e "-t 2 \\" >> $output
     echo -e "'-O 24,48 -E 1 -A 4 -B 16 -T 70 -k 19 -w 200 -d 600 -L 20 -U 40' \\" >> $output
     echo -e "$REFERENCE \\" >> $output
     echo -e "$flash_output > $OUTPUT_BWA\""  >> $output

         # step 2 convert sam to bam

     OUTPUT_BWA_bam=$output_dir_FLASH/$sample.bam

     echo -e "bsub -G team151 -o $output_dir_FLASH/BWA_$counter.out -M 4000 -w\"done(BWA_FLASH_$counter)\" -J FLASH_SAM_TO_BAM_$counter -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
     echo -e "\"samtools view -b $OUTPUT_BWA >  $OUTPUT_BWA_bam\"" >> $output

     # clear SAM

     echo -e "bsub -G team151 -o $output_dir_FLASH/BWA_$counter.out -M 4000 -w\"done(FLASH_SAM_TO_BAM_$counter)\" -J c_1F_$counter -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
     echo -e "\"rm $OUTPUT_BWA\"" >> $output

     # step 3 sort bam by coord

      OUTPUT_BWA_bam_sorted=$(echo "$output_dir_FLASH""/""$sample""_sorted"".bam")

     echo -e "bsub -G team151 -o $output_dir_FLASH/BWA_$counter.out -M 4000 -w\"done(c_1F_$counter)\" -J FLASH_SORT_BAM_$counter -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
     echo -e "\"samtools sort -o $OUTPUT_BWA_bam_sorted $OUTPUT_BWA_bam\"" >> $output

      echo -e "$name\t$sample\t$type\t$OUTPUT_BWA_bam_sorted" >> $replicates_FLASH


     # step 4 INDEX

     echo -e "bsub -G team151 -o $output_dir_FLASH/BWA_$counter.out -M 4000 -w\"done(FLASH_SORT_BAM_$counter)\" -J FLASH_INDEX_$counter -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
     echo -e "\"samtools index $OUTPUT_BWA_bam_sorted\"" >> $output

     # step 5 clear

     echo -e "bsub -G team151 -o $output_dir_FLASH/BWA_$counter.out -M 4000 -w\"done(FLASH_INDEX_$counter)\" -J c_2F_$counter -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
     echo -e "\"rm $OUTPUT_BWA_bam\"" >> $output

     # step 6 filter soft clipped bigger than XXX

     bash_samclip=/nfs/users/nfs_m/mt19/Scripts/Wraper_scripts/244_v6_sub_bash.sh
#     output_dir_FLASH_filtered=$(echo "$output_dir_FLASH_filtered""/")     

     echo -e "bsub -G team151 -o $output_dir_FLASH_filtered/BWA_$counter.out -M 4000 -w\"done(FLASH_INDEX_$counter)\" -J FLASH_samclip_$counter -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
#     echo -e "bsub -G team151 -o $output_dir_FLASH_filtered/BWA_$counter.out -M 4000  -J samclip_$counter -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
     echo -e "\"bash $bash_samclip $output_dir_FLASH_filtered $OUTPUT_BWA_bam_sorted $REFERENCE $parameter $sample 4000 1 normal $ToF_FLASH_filtered $replicates_FLASH_filtered $type $name\"" >> $output



     # step 11 count

     OUTPUT_counts_FLASH=$(echo "$output_dir_FLASH""/""$sample""_count"".txt")

     echo  $OUTPUT_counts_FLASH >> $ToF_FLASH


     echo -e "bsub -G team151 -o $output_dir_FLASH/count_$counter.out -M 4000 -w\"done(c_2F_$counter)\" -J FLASH_count_$counter -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
#     echo -e "bsub -G team151 -o $output_dir_FLASH/count_$counter.out -M 4000  -J FLASH_count_$counter -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q normal -- \\" >> $output
     echo -e "\"samtools view $OUTPUT_BWA_bam_sorted |cut -f3|sort|uniq -c|sort -rg > $OUTPUT_counts_FLASH\"" >> $output

     echo -e "\n\n"  >> $output
     echo -e "\n\n"  >> $output

#     exit
    fi
    fi

fi
 
done < "$5"




echo "bash $output"
