#!/bin/bash
  
 
################################################################################################################################## Dependencies: 

samtools=/nfs/users/nfs_m/mt19/sOFTWARE/samtools-1.6/bin/samtools
python=/usr/bin/python2.7
parser_flash_output=/lustre/scratch115/teams/soranzo/projects/genIE_analysis/src/getFlashOutputDetails.py


MASTER_ROUTE=$2

type=$3
# /lustre/scratch115/teams/soranzo/projects/genIE_analysis/ReSeq/BWA/
output0=$(echo "$MASTER_ROUTE""/""$type""_""TOTAL_READS.txt")
output=$(echo "$MASTER_ROUTE""/""$type""_""CIGAR_summary.txt")
output2=$(echo "$MASTER_ROUTE""/""$type""_""CIGAR_SORT.txt")
output3=$(echo "$MASTER_ROUTE""/""$type""_""CIGAR_mismatch.txt")
output4=$(echo "$MASTER_ROUTE""/""$type""_""CIGAR_mismatch_CONDENSED.txt")
#/nfs/users/nfs_m/mt19/Scripts/Wraper_scripts/200_CREATION.sh


touch $output0
echo -n "" > $output0


touch $output
echo -n "" > $output

touch $output2
echo -n "" > $output2

touch $output3
echo -n "" > $output3

touch $output4
echo -n "" > $output4



output_dir=$(echo "$MASTER_ROUTE")

cd $output_dir

echo -e "name\tSample_Name\tTOTAL_reads" >> $output0
echo -e "name\tCIGAR" >> $output
echo -e "name\tCIGAR" >> $output2
echo -e "name\tSample_Name\tCIGAR\tNM\treads" >> $output4



end_of_file=0

while [[ $end_of_file == 0 ]]
do
  read -r line
  end_of_file=$?

echo "LINE:$line"

if [[ "$line" =~ [^[:space:]] ]]; then

    
    a=($(echo "$line" | tr "\t" '\n'))
    name=${a[0]}
    sample=${a[1]}
    bam=${a[3]}
    
    if [ $name != "name" ]; then # Don't use the header line NOT very elegant

	echo "DEF:name: $name"
	echo "DEF:Sample: $sample"
	echo "DEF:bam: $bam"

#	echo "ME"
#	exit
	samtools view $bam|wc -l|awk -v pat=$sample '{print pat"\t"$0}'|awk -v pat=$name -F"\t" '{print pat"\t"$0}' >> $output0
#	exit
	samtools view $bam|cut -f6|sort|uniq -c|sort -rg|head -10|awk -v pat=$name -F" " '{print pat"\t"$2}' >> $output
#exit


    fi


fi

done < "$1"

#exit

####################################################################################


awk 'NR < 2{ next } {FS="\t"; print $1"\t"$2}' $output|sort|uniq -c|sort -rg|awk -F" " '{print $2"\t"$3}'|sort -k1 > $output2


#exit

###################################################################################


end_of_file=0

while [[ $end_of_file == 0 ]]
do
  read -r line
  end_of_file=$?

#echo "LINE:$line"

if [[ "$line" =~ [^[:space:]] ]]; then


    a=($(echo "$line" | tr "\t" '\n'))
    name=${a[0]}
    sample=${a[1]}
    bam=${a[3]}

    if [ $name != "name" ]; then # Don't use the header line NOT very elegant

        echo "DEF_MASTER:name: $name"
        echo "DEF_MASTER:Sample: $sample"
        echo "DEF_MASTER:bam: $bam"

	$label=$(echo $name"\t"$sample)

        echo "ME"

        #samtools view $bam|cut -f6|sort|uniq -c|sort -rg|head -10|awk -v pat=$name -F" " '{print pat"\t"$2}' >> $output

	end_of_file2=0

	while [[ $end_of_file2 == 0 ]]
	do
	    read -r line2
	    end_of_file2=$?

            #echo "LINE2:$line2"

	    if [[ "$line2" =~ [^[:space:]] ]]; then


		a=($(echo "$line2" | tr "\t" '\n'))
		name2=${a[0]}
		CIGAR=${a[1]}

		if [ $name2 != "name" ]; then # Don't use the header line2 NOT very elegant

		    echo "DEF:name2: $name2"
		    echo "DEF:CIGAR: $CIGAR"
		    echo "ME"

		    if [ $name2 == $name ]; then

#			samtools view $bam|awk -F"\t" '{print $6"\t"$12}' 
			samtools view $bam|awk -v pat=$CIGAR -F"\t" '$6 == pat {print $6"\t"$12}'|awk -v pat=$sample -F"\t" '{print pat"\t"$1"\t"$2}'|awk -v pat=$name -F"\t" '{print pat"\t"$1"\t"$2"\t"$3}' >> $output3



		    fi
		fi

	   fi

	done < "$output2"



    fi


fi

done < "$1"

####################################################################################

awk -F"\t" '{print $1" "$2" "$3" "$4}' $output3|sort|uniq -c|sort -rg|sort -k3|awk -F" " '{print $2"\t"$3"\t"$4"\t"$5"\t"$1}' >> $output4

rm $output3


 
