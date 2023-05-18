#!/bin/bash
  
  
################################################################################################################################## Dependencies: 

samtools=/nfs/users/nfs_m/mt19/sOFTWARE/samtools-1.6/bin/samtools
python=/usr/bin/python2.7
parser_flash_output=/lustre/scratch115/teams/soranzo/projects/genIE_analysis/src/getFlashOutputDetails.py

MASTER_ROUTE=$2

#file=$3

output=$(echo "$MASTER_ROUTE""index_check_summary.txt")





touch $output
echo -n "" > $output


output_dir=$(echo "$MASTER_ROUTE")

cd $output_dir
echo -e "reads\tindexes_string\tfile" >> $output


end_of_file=0

while [[ $end_of_file == 0 ]]
do
  read -r line
  end_of_file=$?

#echo "LINE:$line"

if [[ "$line" =~ [^[:space:]] ]]; then

    
    a=($(echo "$line" | tr "\t" '\n'))
    sample=${a[2]}
    name=${a[0]}
    Read_1=${a[9]}
    Read_2=${a[10]}

    if [ $name == "EROS" ]; then 

	  Read_1=${a[8]}
	  Read_2=${a[9]}

    fi
	echo "DEF:Sample: $sample"
	echo "Read_1: $Read_1"
	echo "Read_2: $Read_2"    

    # if [ $name != "EROS" ]; then 

    # 	exit
    # fi

    if [ $sample != "replicate_full" ]; then # Don't use the header line NOT very elegant
    if [ $Read_1 != "NA" ]; then # Don't use the header line NOT very elegant
#	echo "Hello world"

	
#	exit

	zcat $Read_1|grep "^@"|awk -F"N:0:" '{print $2}'|sort|uniq -c|sort -rg|head -5|awk -v pat=$Read_1 -F" " '{print $1"\t"$2"\t"pat}' >> index_check_summary.txt
	zcat $Read_2|grep "^@"|awk -F"N:0:" '{print $2}'|sort|uniq -c|sort -rg|head -5|awk -v pat=$Read_2 -F" " '{print $1"\t"$2"\t"pat}' >> index_check_summary.txt

#	exit

    fi
    fi


fi

done < "$1"
