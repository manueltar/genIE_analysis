#!/bin/bash

master_file=$1
MASTER_ROUTE=$2
reference=$3
out_master_file=$4
out_master_script=$5

mem=$6
pc=$7
queue=$8


now="$(date +'%d_%m_%Y')"
now2="$(date +"%T")"
now4="_"
now3=$now$now4$now2

java=/software/jre1.8.0_131/bin/java
GATK=/nfs/team151/software/GATK_3.8/GenomeAnalysisTK.jar

 
input_vcf_dir=$(echo "$MASTER_ROUTE""vcf_files""/")
input_intervals_dir=$(echo "$MASTER_ROUTE""/""intervals_files""/")
output_dir=$(echo "$MASTER_ROUTE""/""ALT_FA_sequences""/")

rm -rf $output_dir
mkdir -p $output_dir

# out_master_file=$(echo "$MASTER_ROUTE""/""$type""/""$type""_MASTER_TILES_PLUS_REF_AND_NCGR_SUBS_PLUS_ALT_fa.tsv")
# out_master_script=$(echo "$MASTER_ROUTE""/""$type""/""$type""_MASTER_TILES_PLUS_REF_AND_NCGR_SUBS_PLUS_ALT_fa.sh")

touch $out_master_script
echo -n "" > $out_master_script

touch $out_master_file
echo -n "" > $out_master_file
 
echo "1-START_printing_ALT_sequences:\t$1\t$2\t$3\n"
echo "$now3"

end_of_file=0

while [[ $end_of_file == 0 ]]
do
  read -r line
  end_of_file=$?

#echo "LINE:$line"

if [[ "$line" =~ [^[:space:]] ]]; then


    a=($(echo "$line" | tr "\t" '\n'))
    VAR=${a[0]}
    FLAG=${a[1]}
    vcf_file=${a[2]}
    interval_file=${a[3]}


    if [ $VAR == "VAR" ]; then

	NEW_line=$(echo -e $line"\t""ALT_fa_file") 
#	echo "HEADER:$NEW_line:HEADER"
	echo $NEW_line  >> $out_master_file
	
    fi
    
    if [ $VAR != "VAR" ]; then # Don't use the header line NOT very elegant

        echo "DEF:VAR: $VAR"
	echo "DEF:FLAG: $FLAG"
        echo "DEF:vcf_file: $vcf_file"
        echo "DEF:input_vcf_dir: $input_vcf_dir"
        echo "DEF:interval_file: $interval_file"

#	exit
#	vcf=$(echo "$input_vcf_dir$vcf_file")
#	intervals=$(echo "$input_intervals_dir""$interval_file")

#	output_file=($(echo "$vcf_file" |sed -r 's/$input_vcf_dir//g'))

 #       echo "1:output_file: $output_file"

	#exit
	output_file=($(echo "$vcf_file" | tr "\.vcf" "_ALT"))

	echo "2:output_file: $output_file"

#         exit
	NEW_line=$(echo -e $line"\t"$output_file".fa")
#	echo "DEF:NEW_line: $NEW_line"
	echo $NEW_line  >> $out_master_file
	name=$(echo "ALT_fa_""$vcf_file")
	
	output_file=$(echo "$output_dir""$output_file"".fa")

	echo "3:output_file: $output_file"

#        exit



	



#	stdout_file=($(echo "$vcf_file" | tr "\.vcf" "_ALT"))
        stdout_file=$(echo "$output_dir""FASTA_ALTERNATE_""output.out")

	touch $stdout_file
	echo -n "" > $stdout_file

	echo "bsub -G team151 -o $stdout_file -M $mem  -J $name -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $out_master_script
	echo "\"$java -jar $GATK -T FastaAlternateReferenceMaker -R $reference -o $output_file -L $(echo "$input_intervals_dir$interval_file") -V $(echo "$input_vcf_dir$vcf_file")\"" >> $out_master_script

#	exit

    fi


fi

done < "$1"

echo "2-FINISH master FILE:$1\t$2\t$3\n"
echo "$now3"

bash $out_master_script

echo "3-FINISH RUN FILE:$1\t$2\t$3\n"
echo "$now3"
