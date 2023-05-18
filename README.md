
# Original fastq files are in /nfs/team151/WetLab/GenIE/Sequencing/

# find SNPs from Kolf2 and create the fasta alignment reference

$ bash ~/Scripts/Wraper_scripts/291_genIE_fasta_meta_flash_wrangler.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/ 4\
000 1 normal ~/Scripts/Wraper_scripts/291_CREATION.sh

/software/R-4.1.0/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/R/241_genIE_fasta_builder_REAL_TILE_to_seq_name.R --corresp_file round3_HL60_THP1_vcf_interval_c\
orresp_PLUS_ALT_FA.tsv --intermediate_table round3_HL60_THP1_intermediate_table.tsv --master_file /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_a\
nalysis/round3_HL60_THP1/round3_HL60_THP1_MASTER_FILE.tsv --input_Eve /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/updated_input_correc\
ted.csv --type round3_HL60_THP1 --out /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/

$ bwa index -a bwtsw round3_wt_reference.fasta
$ samtools faidx round3_wt_reference.fasta

$ bwa index -a bwtsw NEW_reference.fasta
$ samtools faidx NEW_reference.fasta


# Build the meta and flash files

/software/R-4.1.0/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/R/242_meta_and_flash_files_builder.R --master_file /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/round3_HL60_THP1_MASTER_FILE.tsv --manifest /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/HL60/manifest_HL60.csv --equivalence /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/HL60/equivalence_HL60.csv --FILES /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/HL60/fastq/File_list.txt --master_path /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/HL60/fastq/ --PREFIXES_Table_EXPANDED /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/HL60/HL60_PREFIXES_Table_EXPANDED.txt --input_regions /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/round3_HL60_THP1_input_regions.tsv  --type HL60 --out /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/HL60/


/software/R-4.1.0/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/R/242_meta_and_flash_files_builder.R --master_file /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/round3_HL60_THP1_MASTER_FILE.tsv --manifest /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/THP1/manifest_THP1.csv --equivalence /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/THP1/equivalence_THP1.csv --FILES /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/THP1/fastq/File_list.txt --master_path /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/THP1/fastq/ --PREFIXES_Table_EXPANDED /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/THP1/THP1_PREFIXES_Table_EXPANDED.txt --input_regions /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/round3_HL60_THP1_input_regions.tsv  --type THP1 --out /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/THP1/

/software/R-4.1.0/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/R/242_meta_and_flash_files_builder.R --master_file /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/round3_HL60_THP1_MASTER_FILE.tsv --manifest /nfs/team151/WetLab/GenIE/MiSeq\ Sample\ sheets/14_07_2020_miseqwalkup275_Genie2_2.csv --equivalence ~/analysis_genIE/round2/Genie_20_21_MANUEL_LIST.csv --FILES /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/genIE_kolf_round275_post_covid/fastq/File_list.txt  --master_path /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/genIE_kolf_round275_post_covid/fastq/ --PREFIXES_Table_EXPANDED /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/genIE_kolf_round275_post_covid/round275_PREFIXES_Table_EXPANDED.txt --input_regions /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/round3_HL60_THP1_input_regions.tsv --type round275 --out /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/genIE_kolf_round275_post_covid/


##### K562 under construction


/nfs/users/nfs_m/mt19/analysis_genIE/K562/Initial.csv

/software/R-4.1.0/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/R/242_meta_and_flash_files_builder.R --master_file /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/round3_HL60_THP1_MASTER_FILE.tsv --manifest /nfs/users/nfs_m/mt19/analysis_genIE/K562/Initial.csv --equivalence /nfs/users/nfs_m/mt19/analysis_genIE/'Genera round 1 gene info simplified_MTS_v3_corrected_SAV1.csv' --FILES /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/K562/fastq/File_list.txt --master_path /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/K562/fastq/ --PREFIXES_Table_EXPANDED /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/K562/K562_PREFIXES_Table_EXPANDED.txt --input_regions /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/round3_HL60_THP1_input_regions.tsv --type K562 --out /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/K562/


#### round2 under construction



/software/R-4.1.0/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/R/242_meta_and_flash_files_builder.R --master_file /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/round3_HL60_THP1_MASTER_FILE.tsv --manifest /nfs/team151/WetLab/GenIE/MiSeq\ Sample\ sheets/11_02_2020_miseqwalkup257_Geni e2_0_and_2_1_final.csv --equivalence /nfs/users/nfs_m/mt19/analysis_genIE/round2/Genie_20_21_MANUEL_LIST.csv --FILES /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round2_exp0/fastq/File_list.txt --master_path /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round2_exp0/fastq/ - \PREFIXES_Table_EXPANDED /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round2_exp0/round2_exp0_PREFIXES_Table_EXPANDED.txt --input_regions /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/round3_HL60_THP1_input_regions.tsv --type round2_exp0 --out /lustre/scratch 123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round2_exp0/


/software/R-4.1.0/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/R/242_meta_and_flash_files_builder.R --master_file /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/round3_HL60_THP1_MASTER_FILE.tsv --manifest /nfs/team151/WetLab/GenIE/MiSeq\ Sample\ sheets/11_02_2020_miseqwalkup257_Geni e2_0_and_2_1_final.csv --equivalence /nfs/users/nfs_m/mt19/analysis_genIE/round2/Genie_20_21_MANUEL_LIST.csv --FILES /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round2_exp1/fastq/File_list.txt --master_path /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round2_exp1/fastq/ --PREFIXES_Table_EXPANDED /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round2_exp1/round2_exp1_PREFIXES_Table_EXPANDED.txt --input_regions /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/round3_HL60_THP1_input_regions.tsv --type round2_exp1 --out /lustre/scratch 123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round2_exp1/


#### round1 under construction


/software/R-4.1.0/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/R/242_meta_and_flash_files_builder.R --master_file /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/round3_HL60_THP1_MASTER_FILE.tsv --manifest /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round1_ReSeq /24_07_19_changed.csv --equivalence /nfs/users/nfs_m/mt19/analysis_genIE/'Genera round 1 gene info simplified_MTS_v3_corrected_SAV1.csv' --FILES /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round1_ReSeq/fastq/File_list.txt --master_path /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_a nalysis/round1_ReSeq/fastq/ --PREFIXES_Table_EXPANDED /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round1_ReSeq/round1_ReSeq_PREFIXES_Table_EXPANDED.txt --input_regions /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/round3_HL60_THP1_input_regions.tsv --type round1_ReSeq --out /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round1_ReSeq/

# Alignment to reference

# Alignment


bash ~/Scripts/Wraper_scripts/244_v6_NEW_SAMCLIP_genIE.sh  HL60_genie.meta.tsv /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/NEW_reference.fasta /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/HL60/BWA_SAMCLIP_10/ ~/Scripts/Wraper_scripts/244_CREATION_HL60.\
sh HL60_flash_input.tsv 10  /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/HL60/fastq/

nohup bash /nfs/users/nfs_m/mt19/Scripts/Wraper_scripts/244_CREATION_HL60.sh > my.log 2>&1 &


bash ~/Scripts/Wraper_scripts/244_v6_NEW_SAMCLIP_genIE.sh THP1_genie.meta.tsv /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/NEW_reference.fasta /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/THP1/BWA_SAMCLIP_10/ ~/Scripts/Wraper_scripts/244_CREATION_THP1.s\
h THP1_flash_input.tsv 10  /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/THP1/fastq/

nohup bash /nfs/users/nfs_m/mt19/Scripts/Wraper_scripts/244_CREATION_THP1.sh > my.log 2>&1 &

bash ~/Scripts/Wraper_scripts/244_v6_NEW_SAMCLIP_genIE.sh  round275_genie.meta.tsv /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/NEW_reference.fasta /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/genIE_kolf_round275_post_covid/BWA_SAMCLIP_10/ ~/Scripts/Wra\
per_scripts/244_CREATION_275.sh  round275_flash_input.tsv 10 /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/genIE_kolf_round275_post_covid/fastq/

nohup bash /nfs/users/nfs_m/mt19/Scripts/Wraper_scripts/244_CREATION_275.sh > my.log 2>&1 &

#### K562 under construction FALSE FRIEND mkdir BWA_SAMCLIP_10

bash ~/Scripts/Wraper_scripts/244_v6_NEW_SAMCLIP_genIE.sh  K562_genie.meta.tsv /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/NEW_reference.fasta /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/K562/BWA_SAMCLIP_10/ ~/Scripts/Wraper_scripts/244_CREATION_K562.\
sh  K562_flash_input.tsv 10 /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/K562/fastq/

nohup bash /nfs/users/nfs_m/mt19/Scripts/Wraper_scripts/244_CREATION_K562.sh > my.log 2>&1 &

#### round2_exp0 round2_exp1

bash ~/Scripts/Wraper_scripts/244_v6_NEW_SAMCLIP_genIE.sh  round2_exp1_genie.meta.tsv /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/NEW_reference.fasta /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round2_exp1/BWA_SAMCLIP_10/ ~/Scripts/Wraper_scripts/244_\
CREATION_round2_exp1.sh round2_exp1_flash_input.tsv 10  /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round2_exp1/fastq/

nohup bash /nfs/users/nfs_m/mt19/Scripts/Wraper_scripts/244_CREATION_round2_exp1.sh > my.log 2>&1 &

bash ~/Scripts/Wraper_scripts/244_v6_NEW_SAMCLIP_genIE.sh  round2_exp0_genie.meta.tsv /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/NEW_reference.fasta /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round2_exp0/BWA_SAMCLIP_10/ ~/Scripts/Wraper_scripts/244_\
CREATION_round2_exp0.sh round2_exp0_flash_input.tsv 10  /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round2_exp0/fastq/

nohup bash /nfs/users/nfs_m/mt19/Scripts/Wraper_scripts/244_CREATION_round2_exp0.sh > my.log 2>&1 &

#### round1_ReSeq

bash ~/Scripts/Wraper_scripts/244_v6_NEW_SAMCLIP_genIE.sh  round1_ReSeq_genie.meta.tsv /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/NEW_reference.fasta /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round1_ReSeq/BWA_SAMCLIP_10/ ~/Scripts/Wraper_scripts/24\
4_CREATION_round1_ReSeq.sh round1_ReSeq_flash_input.tsv 10  /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round1_ReSeq/fastq/

nohup bash /nfs/users/nfs_m/mt19/Scripts/Wraper_scripts/244_CREATION_round1_ReSeq.sh > my.log 2>&1 &

# QC part I

bash ~/Scripts/Wraper_scripts/250_RScript_1_count_matrix_generator_v4.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/HL60/BWA_SAMCLIP_10/ /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/NEW_reference.fasta  /lustre/scratch123/hgi/mdt1/teams/soranzo/projec\
ts/genIE_analysis/HL60/HL60_PREFIXES_Table_EXPANDED.txt ~/Scripts/Wraper_scripts/250_15_CREATION_HL60.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/HL60/fastq/ 15 /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/NEW_input_regions.tsv HL60

bash ~/Scripts/Wraper_scripts/250_RScript_1_count_matrix_generator_v4.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/THP1/BWA_SAMCLIP_10/ /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/NEW_reference.fasta  /lustre/scratch123/hgi/mdt1/teams/soranzo/projec\
ts/genIE_analysis/THP1/THP1_PREFIXES_Table_EXPANDED.txt ~/Scripts/Wraper_scripts/250_15_CREATION_THP1.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/THP1/fastq/ 15 /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/NEW_input_regions.tsv THP1

bash ~/Scripts/Wraper_scripts/250_RScript_1_count_matrix_generator_v4.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/genIE_kolf_round275_post_covid/BWA_SAMCLIP_10/ /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/NEW_reference.fasta  /lustre/scratch123/hgi\
/mdt1/teams/soranzo/projects/genIE_analysis/genIE_kolf_round275_post_covid/round275_PREFIXES_Table_EXPANDED.txt ~/Scripts/Wraper_scripts/250_15_CREATION_round275.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/genIE_kolf_round275_post_covid/fastq/ 15 /lustre/scratch123/hgi/mdt1/teams/soranz\
o/projects/genIE_analysis/round3_HL60_THP1/NEW_input_regions.tsv round275

#### K562 under construction FALSE FRIEND genIE_kolf_round275_post_covid

bash ~/Scripts/Wraper_scripts/250_RScript_1_count_matrix_generator_v4.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/K562/BWA_SAMCLIP_10/ /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/NEW_reference.fasta  /lustre/scratch123/hgi/mdt1/teams/soranzo/projec\
ts/genIE_analysis/K562/K562_PREFIXES_Table_EXPANDED.txt ~/Scripts/Wraper_scripts/250_15_CREATION_K562.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/K562/fastq/ 15 /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/NEW_input_regions.tsv K562

#### round2_exp0 & round2_exp1


bash ~/Scripts/Wraper_scripts/250_RScript_1_count_matrix_generator_v4.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round2_exp0/BWA_SAMCLIP_10/ /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/NEW_reference.fasta  /lustre/scratch123/hgi/mdt1/teams/soranzo\
/projects/genIE_analysis/round2_exp0/round2_exp0_PREFIXES_Table_EXPANDED.txt ~/Scripts/Wraper_scripts/250_15_CREATION_round2_exp0.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round2_exp0/fastq/ 15 /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/NEW_inpu\
t_regions.tsv round2_exp0

bash ~/Scripts/Wraper_scripts/250_RScript_1_count_matrix_generator_v4.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round2_exp1/BWA_SAMCLIP_10/ /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/NEW_reference.fasta  /lustre/scratch123/hgi/mdt1/teams/soranzo\
/projects/genIE_analysis/round2_exp1/round2_exp1_PREFIXES_Table_EXPANDED.txt ~/Scripts/Wraper_scripts/250_15_CREATION_round2_exp1.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round2_exp1/fastq/ 15 /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/NEW_inpu\
t_regions.tsv round2_exp1

#### round1_ReSeq


bash ~/Scripts/Wraper_scripts/250_RScript_1_count_matrix_generator_v4.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round1_ReSeq/BWA_SAMCLIP_10/ /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/NEW_reference.fasta  /lustre/scratch123/hgi/mdt1/teams/soranz\
o/projects/genIE_analysis/round1_ReSeq/round1_ReSeq_PREFIXES_Table_EXPANDED.txt ~/Scripts/Wraper_scripts/250_15_CREATION_round1_ReSeq.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round1_ReSeq/fastq/ 15 /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/NEW\
_input_regions.tsv round1_ReSeq

# QC part II

bash ~/Scripts/Wraper_scripts/250_RScript_1_count_matrix_generator_v4_partII.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/HL60/BWA_SAMCLIP_10/ /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/NEW_reference.fasta  /lustre/scratch123/hgi/mdt1/teams/soranzo\
/projects/genIE_analysis/HL60/HL60_PREFIXES_Table_EXPANDED.txt ~/Scripts/Wraper_scripts/250_15_CREATION_HL60_partII.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/HL60/fastq/ 15 /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/NEW_input_regions.tsv HL60

bash ~/Scripts/Wraper_scripts/250_RScript_1_count_matrix_generator_v4_partII.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/THP1/BWA_SAMCLIP_10/ /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/NEW_reference.fasta  /lustre/scratch123/hgi/mdt1/teams/soranzo\
/projects/genIE_analysis/THP1/THP1_PREFIXES_Table_EXPANDED.txt ~/Scripts/Wraper_scripts/250_15_CREATION_THP1_partII.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/THP1/fastq/ 15 /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/NEW_input_regions.tsv THP1

bash ~/Scripts/Wraper_scripts/250_RScript_1_count_matrix_generator_v4_partII.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/genIE_kolf_round275_post_covid/BWA_SAMCLIP_10/ /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/NEW_reference.fasta  /lustre/scratch\
123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/genIE_kolf_round275_post_covid/round275_PREFIXES_Table_EXPANDED.txt ~/Scripts/Wraper_scripts/250_10_CREATION_round275_partII.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/genIE_kolf_round275_post_covid/fastq/ 10 /lustre/scratch123/hgi/mdt\
1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/NEW_input_regions.tsv round275

bash ~/Scripts/Wraper_scripts/250_RScript_1_count_matrix_generator_v4_partII.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/K562/BWA_SAMCLIP_10/ /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/NEW_reference.fasta  /lustre/scratch123/hgi/mdt1/teams/soranzo\
/projects/genIE_analysis/K562/K562_PREFIXES_Table_EXPANDED.txt ~/Scripts/Wraper_scripts/250_10_CREATION_K562_partII.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/K562/fastq/ 10 /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/NEW_input_regions.tsv K562

bash ~/Scripts/Wraper_scripts/250_RScript_1_count_matrix_generator_v4_partII.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round2_exp0/BWA_SAMCLIP_10/ /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/NEW_reference.fasta  /lustre/scratch123/hgi/mdt1/teams/\
soranzo/projects/genIE_analysis/round2_exp0/round2_exp0_PREFIXES_Table_EXPANDED.txt ~/Scripts/Wraper_scripts/250_15_CREATION_round2_exp0_partII.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round2_exp0/fastq/ 15 /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60\
_THP1/NEW_input_regions.tsv round2_exp0

bash ~/Scripts/Wraper_scripts/250_RScript_1_count_matrix_generator_v4_partII.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round2_exp1/BWA_SAMCLIP_10/ /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/NEW_reference.fasta  /lustre/scratch123/hgi/mdt1/teams/\
soranzo/projects/genIE_analysis/round2_exp1/round2_exp1_PREFIXES_Table_EXPANDED.txt ~/Scripts/Wraper_scripts/250_15_CREATION_round2_exp1_partII.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round2_exp1/fastq/ 15 /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60\
_THP1/NEW_input_regions.tsv round2_exp1

bash ~/Scripts/Wraper_scripts/250_RScript_1_count_matrix_generator_v4_partII.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round1_ReSeq/BWA_SAMCLIP_10/ /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/NEW_reference.fasta  /lustre/scratch123/hgi/mdt1/teams\
/soranzo/projects/genIE_analysis/round1_ReSeq/round1_ReSeq_PREFIXES_Table_EXPANDED.txt ~/Scripts/Wraper_scripts/250_15_CREATION_round1_ReSeq_partII.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round1_ReSeq/fastq/ 15 /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3\
_HL60_THP1/NEW_input_regions.tsv round1_ReSeq


# rgenie package


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round2_exp0/BWA_SAMCLIP_10/FLASH_FILTERED/rgenie_round2_exp0_10_10.out -M 16000 -R"select[mem >16000] rusage[mem=16000] span[hosts=1]" -n4 -q normal "/software/R-4.1.0/bin/Rscript  /nfs/users/nfs_m/mt19/Scripts/R/New_genIE_script_for_R_v4.R --regions /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/NEW_input_regions.tsv --replicates /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round2_exp0/BWA_SAMCLIP_10/FLASH_FILTERED/replicates_AFTER_QC.tsv --required_match_left_THRESHOLD 10 --required_match_right_THRESHOLD 10 --dropouts_ALL_PLOTS NA --QC_PASS C2CD5,KLF16,SLC9A3R1,STAT6 --type FLASH_FILTERED_AFTER_QC --out /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round2_exp0/BWA_SAMCLIP_10/FLASH_FILTERED/ --input /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/input_corrected_unblanked.tsv"


bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/K562/BWA_SAMCLIP_10/FLASH_FILTERED/rgenie_K562_10_10.out -M 16000 -R"select[mem >16000] rusage[mem=16000] span[hosts=1]" -n4 -q normal "/software/R-4.1.0/bin/Rscript  /nfs/users/nfs_m/mt19/Scripts/R/New_genIE_script_for_R_v4.R --regions /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/NEW_input_regions.tsv --replicates /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/K562/BWA_SAMCLIP_10/FLASH_FILTERED/replicates_AFTER_QC.tsv --required_match_left_THRESHOLD 10 --required_match_right_THRESHOLD 10 --dropouts_ALL_PLOTS UGCG --QC_PASS BID,BRAP,CUX1,C2CD5,DNMT1,EPB41,EROS,FOXP1,NBN,SH2B3,TNRC6A,DLG4,UBAC1,ZNRF3 --type FLASH_FILTERED_AFTER_QC --out /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/K562/BWA_SAMCLIP_10/FLASH_FILTERED/ --input /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/input_corrected_unblanked.tsv"

bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/HL60/BWA_SAMCLIP_10/FLASH_FILTERED/rgenie_HL60_10_10.out -M 16000 -R"select[mem >16000] rusage[mem=16000] span[hosts=1]" -n4 -q normal "/software/R-4.1.0/bin/Rscript  /nfs/users/nfs_m/mt19/Scripts/R/New_genIE_script_for_R_v4.R --regions /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/NEW_input_regions.tsv --replicates /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/HL60/BWA_SAMCLIP_10/FLASH_FILTERED/replicates_AFTER_QC.tsv --required_match_left_THRESHOLD 10 --required_match_right_THRESHOLD 10 --dropouts_ALL_PLOTS NA --QC_PASS BID,CEP104,C2CD5,DNMT1,EPB41,FOXP1,FUT8,LUC7L,PLEKHG3,SH2B3,UGCG --type FLASH_FILTERED_AFTER_QC --out /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/HL60/BWA_SAMCLIP_10/FLASH_FILTERED/ --input /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/input_corrected_unblanked.tsv"

bsub -G team151 -o /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/genIE_kolf_round275_post_covid/BWA_SAMCLIP_10/FLASH_FILTERED/rgenie_round275_10_10.out -M 16000 -R"select[mem >16000] rusage[mem=16000] span[hosts=1]" -n4 -q normal "/software/R-4.1.0/bin/Rscript  /nfs/users/nfs_m/mt19/Scripts/R/New_genIE_script_for_R_v4.R --regions /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/NEW_input_regions.tsv --replicates /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/genIE_kolf_round275_post_covid/BWA_SAMCLIP_10/FLASH_FILTERED/replicates_AFTER_QC.tsv --required_match_left_THRESHOLD 10 --required_match_right_THRESHOLD 10 --dropouts_ALL_PLOTS NA --QC_PASS C2CD5,DNMT1,CUX1,EPB41,FOXP1,FUT8,SH2B3,STAT6 --type FLASH_FILTERED_AFTER_QC --out /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/genIE_kolf_round275_post_covid/BWA_SAMCLIP_10/FLASH_FILTERED/ --input /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/round3_HL60_THP1/input_corrected_unblanked.tsv"

# GLOBAL_ANALYSIS

bash ~/Scripts/Wraper_scripts/294_genIE_GLOBAL_ANALYSIS.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/GLOBAL_ANALYSIS/ 4000 1 normal ~/Scripts/Wraper_scripts/294_CREATION.sh

# EXPRESSION PART


/software/R-4.1.0/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/R/233_EXP_kolf2_v2.R --STAR_kolf2_1 /nfs/team151/WetLab/GenIE/AlignmentSTRINGENT_NOgDNAIndex_output_ERR1203463.isoforms.results --STAR_kolf2_2 /nfs/team151/WetLab/GenIE/AlignmentSTRINGENT_NOgDNAIndex_output_ERR1243466.isoforms.results --Equiv.table /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/BP_Iso_Reanalysis/OLD/Homo_sapiens.GRCh37.87_Transcripts_table.txt --SELECTED_GENES DNMT1,S1PR2 --VEP_route /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/GLOBAL_ANALYSIS/AS_genIE/VEP_CSQ_csv_tables/  --type Kolf2_EXP --out /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/GLOBAL_ANALYSIS/
