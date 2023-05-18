#!/bin/bash
 
Rscript=/software/R-4.1.0/bin/Rscript

MASTER_ROUTE=$1
mem=$2
pc=$3
queue=$4
output=$5

touch $output
echo -n "" > $output

echo -e "#!/bin/bash"  >> $output
echo -e "\n\n"  >> $output




output_dir=$(echo "$MASTER_ROUTE")




echo "#############################################################################################################################################################">> $output
echo "#################################### PRINT BED and FILE FOR SNP search  ######################################################################################################"  >> $output

type=$(echo "bed_and_SNP_search")
outfile_bed_and_SNP_search=$(echo "$output_dir""outfile""_""$type"".out")
touch $outfile_bed_and_SNP_search
echo -n "" > $outfile_bed_and_SNP_search
name_bed_and_SNP_search=$(echo "$type""_job")

Rscript_bed_and_SNP_search=$(echo "/nfs/users/nfs_m/mt19/Scripts/R/239_genIE_fasta_builder_INITIAL_STEPS.R")
span=$(echo '250')
input_Eve=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/updated_input_corrected.csv")
tag=$(echo "round3_HL60_THP1")

echo "bsub -G team151 -o $outfile_bed_and_SNP_search -M $mem  -J $name_bed_and_SNP_search -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_bed_and_SNP_search \\" >> $output
echo "--input_Eve $input_Eve \\" >> $output
echo "--span $span \\" >> $output
echo "--type $tag --out $output_dir\"" >> $output




echo "#################################### bash search on HPSI0114i-kolf_2.wgs.gatk.haplotype_caller.20170425.genotypes.vcf.gz  ######################################################################################################"  >> $output

type=$(echo "bash_search_on_HPSI")
outfile_bash_search_on_HPSI=$(echo "$output_dir""outfile""_""$type"".out")
touch $outfile_bash_search_on_HPSI
echo -n "" > $outfile_bash_search_on_HPSI
name_bash_search_on_HPSI=$(echo "$type""_job")


bash_script_bash_search_on_HPSI=$(echo "/nfs/users/nfs_m/mt19/Scripts/Wraper_scripts/266_Kolf2_snp_finder.sh") 
input_bash_script_bash_search_on_HPSI=$(echo "$output_dir""genERA_VAR_POST_HF1_2_3_region_file_FROM_TO.txt")
HPSI_file=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/reference_files/HPSI0114i-kolf_2.wgs.gatk.haplotype_caller.20170425.genotypes.vcf.gz")


echo "bsub -G team151 -o $outfile_bash_search_on_HPSI -M $mem -w\"done($name_bed_and_SNP_search)\" -J $name_bash_search_on_HPSI -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"bash $bash_script_bash_search_on_HPSI $output_dir $input_bash_script_bash_search_on_HPSI $HPSI_file $tag\"" >> $output


echo "#################################### bedtools get REF fasta  ######################################################################################################"  >> $output

REFERENCE=$(echo "/lustre/scratch115/resources/ref/Homo_sapiens/GRCh37_53/Homo_sapiens.GRCh37.dna.all.fa")

outfile_get_REF_fasta=$(echo "$output_dir""outfile_get_REF_fasta""_""$type"".out")
touch $outfile_get_REF_fasta
echo -n "" > $outfile_get_REF_fasta
name_get_REF_fasta=$(echo "$type""_name_get_REF_fasta")


bedfile=$(echo "$output_dir""$tag""_TILES"".bed")

REF_fasta=$(echo "$output_dir""$tag""_REF"".fasta")


echo "bsub -G team151 -o $outfile_get_REF_fasta -M $mem  -w\"done($name_bed_and_SNP_search)\" -J $name_get_REF_fasta -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"bedtools getfasta -fi $REFERENCE -fo $REF_fasta -bed $bedfile\"" >> $output

echo "#############################################################################################################################################################">> $output
echo "#################################### vcf step  ######################################################################################################"  >> $output

type=$(echo "vcf_step")
outfile_vcf_step=$(echo "$output_dir""outfile""_""$type"".out")
touch $outfile_vcf_step
echo -n "" > $outfile_vcf_step
name_vcf_step=$(echo "$type""_job")

Rscript_vcf_step=$(echo "/nfs/users/nfs_m/mt19/Scripts/R/240_genIE_fasta_builder_vcf_step.R")
span=$(echo '250')


snps=$(echo "$output_dir""$tag""_SNPS.txt")
SNPS_parsed=$(echo "$output_dir""$tag""_SNPS_parsed.txt")
input_corrected_unblanked=$(echo "$output_dir""input_corrected_unblanked.tsv")
tag=$(echo "round3_HL60_THP1")
REF_FASTA=$(echo "$output_dir""$tag""_REF.fasta")
bedfile=$(echo "$output_dir""$tag""_TILES.bed")


echo "bsub -G team151 -o $outfile_vcf_step -M $mem -w\"done($name_bash_search_on_HPSI) &&  done($name_get_REF_fasta)\" -J $name_vcf_step -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_vcf_step \\" >> $output
echo "--snps $snps \\" >> $output
echo "--SNPS_parsed $SNPS_parsed \\" >> $output
echo "--input_corrected_unblanked $input_corrected_unblanked \\" >> $output
echo "--span $span \\" >> $output
echo "--REF_FASTA $REF_FASTA \\" >> $output
echo "--bedfile $bedfile \\" >> $output
echo "--type $tag --out $output_dir\"" >> $output

echo "#############################################################################################################################################################">> $output
echo "#################################### FASTA ALTERNATE  ######################################################################################################"  >> $output

##### step IIc ALT fasta

type=$(echo "LAUNCH_FASTA_ALTERNATE")
outfile_LAUNCH_FASTA_ALTERNATE=$(echo "$output_dir""outfile""_""$type"".out")
touch $outfile_LAUNCH_FASTA_ALTERNATE
echo -n "" > $outfile_LAUNCH_FASTA_ALTERNATE
name_LAUNCH_FASTA_ALTERNATE=$(echo "$type""_job")

master_file=$(echo "$output_dir""$tag""_vcf_interval_corresp.tsv")
FASTA_ALTERNATE=$(echo "/nfs/users/nfs_m/mt19/Scripts/Wraper_scripts/91_GATK_ALTERNATE_genIE.sh")
out_master_file=$(echo "$output_dir""$tag""_vcf_interval_corresp_PLUS_ALT_FA.tsv")
out_master_script=$(echo "/nfs/users/nfs_m/mt19/Scripts/Wraper_scripts/291_CREATION_FASTA_ALTERNATE.sh")



echo "bsub -G team151 -o $outfile_LAUNCH_FASTA_ALTERNATE -w\"done($name_vcf_step)\" -J $name_LAUNCH_FASTA_ALTERNATE -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -M $mem -n$pc -q $queue -- \\" >> $output
echo "\"bash $FASTA_ALTERNATE $master_file $output_dir $REFERENCE $out_master_file $out_master_script $mem $pc $queue\"" >> $output


echo "bash $output"
