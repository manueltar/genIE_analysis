
suppressMessages(library("plyr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("data.table", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("crayon", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("withr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("ggplot2", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))
library("farver", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("labeling", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
suppressMessages(library("optparse", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("dplyr", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("withr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("backports", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("broom", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("rstudioapi", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("cli", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("tzdb", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("tidyverse", lib.loc="/nfs/team151/software/manuel_R_libs_4_1//"))
suppressMessages(library("BiocGenerics", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("S4Vectors", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("IRanges", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("GenomeInfoDb", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("GenomicRanges", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("Biobase", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("AnnotationDbi", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("GenomicFeatures", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("OrganismDbi", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
library("XVector", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
suppressMessages(library("Biostrings", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))

suppressMessages(library("vcfR", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))


opt = NULL

parse_266 = function(option_list)
{
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  
  #### Read SNPS ----
  
  SNPS<-fread(opt$snps, header=T, stringsAsFactors = F, skip = "CHROM")
  
  cat("SNPS_\n")
  str(SNPS)
  cat("\n")
  
  # quit(status=1)
  
  # if VAR is present in Kolf2 has a , in ALT
  
  Kolf2_proficient2<-SNPS[grep("\\,", SNPS$ALT),]
  
  cat("Kolf2_proficient2_\n")
  str(Kolf2_proficient2)
  cat("\n")
  
  Kolf2_proficient2$FlagHOM<-"NA"
  Kolf2_proficient2$FlagHOM[grep("1\\/1",Kolf2_proficient2$'HPSI0114i-kolf_2')]<-"HOM_KOLF2"
  Kolf2_proficient2$FlagHOM[grep("0\\/1",Kolf2_proficient2$'HPSI0114i-kolf_2')]<-"HET_KOLF2"
  Kolf2_proficient2$FlagHOM[grep("1\\/0",Kolf2_proficient2$'HPSI0114i-kolf_2')]<-"HET_KOLF2"
  Kolf2_proficient2$FlagHOM[grep("0\\/0",Kolf2_proficient2$'HPSI0114i-kolf_2')]<-"REF_KOLF2"
  
  
  
  
  Kolf2_proficient2$FlagVAR_TYPE<-"Kolf2"
  
  Kolf2_proficient2_nonREF<-Kolf2_proficient2[(Kolf2_proficient2$FlagHOM != "REF_KOLF2"),]
  
  cat("Kolf2_proficient2_nonREF_\n")
  str(Kolf2_proficient2_nonREF)
  cat("\n")
  
  colnames(Kolf2_proficient2_nonREF)[1]<-"chrom"
  Kolf2_proficient2_nonREF$chrom<-paste('chr',Kolf2_proficient2_nonREF$chrom,sep='')
  colnames(Kolf2_proficient2_nonREF)[2]<-"pos"
  colnames(Kolf2_proficient2_nonREF)[3]<-"id"
  colnames(Kolf2_proficient2_nonREF)[4]<-"ref"
  colnames(Kolf2_proficient2_nonREF)[5]<-"alt"
  colnames(Kolf2_proficient2_nonREF)[6]<-"QUAL"
  colnames(Kolf2_proficient2_nonREF)[7]<-"FILTER"
  colnames(Kolf2_proficient2_nonREF)[8]<-"INFO"
  colnames(Kolf2_proficient2_nonREF)[9]<-"DUMMY"
  colnames(Kolf2_proficient2_nonREF)[10]<-"PHASE"
  Kolf2_proficient2_nonREF$alt<-gsub(",<NON_REF>","",Kolf2_proficient2_nonREF$alt)
  
  
  cat("Kolf2_proficient2_nonREF_\n")
  str(Kolf2_proficient2_nonREF)
  cat("\n")
  
  #### SAVE 3 ----
  
  setwd(out)
  
  filename<-paste(type,"_SNPS_parsed",".txt",sep='')
  
  write.table(Kolf2_proficient2_nonREF, 
              file=filename, 
              sep="\t", quote=F, row.names = F)
  
  # quit(status=1)
  
}

Print_vcf = function(option_list)
{
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### READ and transform span ----
  
  span = opt$span
  
  cat("OUT_\n")
  cat(sprintf(as.character(span)))
  cat("\n")
  
  #### Read snps_parsed ----
  
  SNPS_parsed<-read.table(file=opt$SNPS_parsed, header=T, stringsAsFactors = F)
  
  indexes_interest<-c(which(colnames(SNPS_parsed) == "chrom"),
                      which(colnames(SNPS_parsed) == "pos"),
                      which(colnames(SNPS_parsed) == "id"),
                      which(colnames(SNPS_parsed) == "ref"),
                      which(colnames(SNPS_parsed) == "alt"),
                      which(colnames(SNPS_parsed) == "FlagHOM"),
                      which(colnames(SNPS_parsed) == "FlagVAR_TYPE"))
  
  SNPS_parsed<-unique(SNPS_parsed[,indexes_interest])
  
  cat("SNPS_parsed_\n")
  str(SNPS_parsed)
  cat("\n")
  
  #### Read input_corrected ----
  
  input_corrected_unblanked<-read.table(file=opt$input_corrected_unblanked, header=T, stringsAsFactors = F)
  
  cat("input_corrected_unblanked\n")
  str(input_corrected_unblanked)
  cat("\n")
  
  
  
 
  #### Read bedfile ----
  
  BED<-fread(opt$bedfile, header=F, stringsAsFactors = F,sep="\t")
  
  BED<-BED[,-c(5:6)]
  colnames(BED)<-c("chr","START","END","NAME")
  
  cat("BED_\n")
  str(BED)
  cat("\n")
  
  #quit(status=1)
  ##### READ REF FILE & extract EROS REF and ALT----
  
 
  
  fastaFile<-readDNAStringSet(file=opt$REF_FASTA)
  
  seq_name = names(fastaFile)
  sequence = paste(fastaFile)
  REF_fasta <- data.frame(seq_name, sequence, stringsAsFactors = F)
  
  cat("REF_fasta_\n")
  str(REF_fasta)
  cat("\n")
  
  REF_fasta$chr<-gsub(":.+$","",REF_fasta$seq_name)
  
  REF_fasta$START<-gsub("^[^:]+:","",REF_fasta$seq_name)
  REF_fasta$START<-gsub("-.+$","",REF_fasta$START)
  
  REF_fasta$END<-gsub("^[^-]+-","",REF_fasta$seq_name)
  
  cat("REF_fasta_\n")
  str(REF_fasta)
  cat("\n")
  
  REF.df<-merge(REF_fasta,
                BED,
                by=c("chr","START","END"))
  
  cat("REF.df_\n")
  str(REF.df)
  cat("\n")
  
  REF.df$VAR<-gsub("__.+$","",REF.df$NAME)
  REF.df$HGNC<-gsub(".+__","",REF.df$NAME)
  
  REF.df$chr<-gsub("Chr","chr",REF.df$chr)
  REF.df$START<-as.integer(REF.df$START)
  REF.df$END<-as.integer(REF.df$END)
  
  
  cat("REF.df_\n")
  str(REF.df)
  cat("\n")
  
  REF.df_EROS<-REF.df[which(REF.df$HGNC == "EROS"),]
  
  cat("REF.df_EROS_\n")
  str(REF.df_EROS)
  cat("\n")
  
  total_length<-nchar(REF.df_EROS$sequence)
  
  REL_POS_VAR<-total_length-span
  
  cat("REL_POS_VAR:")
  cat(sprintf(as.character(REL_POS_VAR)))
  cat("\n")
  
  ref_CTRL<-substr(REF.df_EROS$sequence, REL_POS_VAR,REL_POS_VAR)
  
  cat("ref_CTRL:")
  cat(sprintf(as.character(ref_CTRL)))
  cat("\n")
  
  alt_CTRL<-"NA"
  
  if(ref_CTRL == "C")
  {
    alt_CTRL<-"T"
  }
  if(ref_CTRL == "T")
  {
    alt_CTRL<-"C"
  }
  if(ref_CTRL == "G")
  {
    alt_CTRL<-"A"
  }
  if(ref_CTRL == "A")
  {
    alt_CTRL<-"G"
  }
  
  ref<-ref_CTRL
  alt<-alt_CTRL
  pos<-REF.df_EROS$END-span
  
  variant_CTRL<-paste(REF.df_EROS$chr,pos,ref,alt, sep="_")
  
  
  REF.df$VAR[(REF.df$HGNC == "EROS")]<-variant_CTRL
  
  REF.df_EROS<-REF.df[which(REF.df$HGNC == "EROS"),]
  
  cat("REF.df_EROS_\n")
  str(REF.df_EROS)
  cat("\n")
  
  #### VARS LOOP ----
  
  VARS<-unique(input_corrected_unblanked$VAR)
  
  cat("VARS_\n")
  str(VARS)
  cat("\n")
  
  list_VARS<-list()
  
  for(i in 1:length(VARS))
  {
    VAR_sel<-VARS[i]
    
    cat("------------------------------------->\t")
    cat(sprintf(as.character(VAR_sel)))
    cat("\n")
    
  
    
    input_corrected_unblanked_sel<-input_corrected_unblanked[which(input_corrected_unblanked$VAR == VAR_sel),]
    
    # cat("input_corrected_unblanked_sel_\n")
    # str(input_corrected_unblanked_sel)
    # cat("\n")
    
    chr_sel<-input_corrected_unblanked_sel$chr
    chr_sel<-as.character(paste('chr',chr_sel,sep=''))
    pos_sel<-input_corrected_unblanked_sel$pos
    ref_sel<-input_corrected_unblanked_sel$ref
    alt_sel<-input_corrected_unblanked_sel$alt
    rsid_sel<-input_corrected_unblanked_sel$rsid
    HGNC_sel<-input_corrected_unblanked_sel$HGNC
    
    ### span
    
    limit_upstream<-pos_sel-span
    
    # cat("limit_upstream>\t")
    # cat(sprintf(as.character(limit_upstream)))
    # cat("\n")
    
    limit_downstream<-pos_sel+span
    
    # cat("limit_downstream>\t")
    # cat(sprintf(as.character(limit_downstream)))
    # cat("\n")
    
    #### select table
    
    VAR_sel_table<-as.data.frame(cbind(chr_sel,pos_sel,rsid_sel),
                                 stringsAsFactors = FALSE)
    colnames(VAR_sel_table)<-c("chr","pos","id")
    VAR_sel_table$pos<-as.integer(VAR_sel_table$pos)
    
   
    
    VAR_sel_table$VAR<-VAR_sel
    VAR_sel_table$HGNC<-HGNC_sel
    
    VAR_sel_table$vcf_VAR<-VAR_sel
    
    
    # cat("VAR_sel_table_\n")
    # str(VAR_sel_table)
    # cat("\n")
    
    if(VAR_sel == "Control")
    {
      
      VAR_sel_table<-as.data.frame(cbind(variant_CTRL),
                                   stringsAsFactors = FALSE)
      
      colnames(VAR_sel_table)<-"VAR"
      
      VAR_sel_table$chr<-gsub("_.+$","",VAR_sel_table$VAR)
      VAR_sel_table$pos<-gsub("^chr[^_]+_","",VAR_sel_table$VAR)
      VAR_sel_table$pos<-gsub("_.+$","",VAR_sel_table$pos)
      
      # VAR_sel_table$ref<-gsub("^chr[^_]+_[^_]+_","",VAR_sel_table$VAR)
      # VAR_sel_table$ref<-gsub("_.+$","",VAR_sel_table$ref)
      # 
      # VAR_sel_table$alt<-gsub("^chr[^_]+_[^_]+_[^_]+_","",VAR_sel_table$VAR)
      
      VAR_sel_table$id<-"Control"
      
      VAR_sel_table$HGNC<-HGNC_sel
      
      # cat("VAR_sel_table_\n")
      # str(VAR_sel_table)
      # cat("\n")
      
      VAR_sel_table$vcf_VAR<-variant_CTRL
      
      
      
      
      # cat("VAR_sel_table_\n")
      # str(VAR_sel_table)
      # cat("\n")
      # 
      # 
      # 
      # quit(status=1)
      
      #VAR_sel<-variant_CTRL
    }
    
    # 
    # accomp.variants<-as.data.frame(unique(SNPS_parsed[(SNPS_parsed$chr == chr_sel),]),
    #                                stringsAsFactors = FALSE)
    # 
    # cat("accomp.variants_\n")
    # str(accomp.variants)
    # cat("\n")
    # 
    
    
    accomp.variants<-as.data.frame(unique(SNPS_parsed[(SNPS_parsed$chr == chr_sel &
                                                         SNPS_parsed$pos >= limit_upstream &
                                                         SNPS_parsed$pos <= limit_downstream),]),
                                   stringsAsFactors = FALSE)
    
    # cat("accomp.variants_\n")
    # str(accomp.variants)
    # cat("\n")
    # 
    # quit(status = 1)
    
    if(dim(accomp.variants)[1] >0)
    {
      accomp.variants$RSid<-as.character(paste(accomp.variants$id,
                                               accomp.variants$FlagHOM,
                                               accomp.variants$FlagVAR_TYPE,
                                               accomp.variants$alt,
                                               sep="_"))
      
      
      # accomp.variants$FlagVAR_TYPE,
      
      
      accomp.variants_table<-as.data.frame(cbind(accomp.variants$chrom,accomp.variants$pos,accomp.variants$RSid),
                                           stringsAsFactors = FALSE)
      
      colnames(accomp.variants_table)<-c("chr","pos","id")
      
      accomp.variants_table$pos<-accomp.variants_table$pos
      
      accomp.variants_table$VAR<-VAR_sel
      accomp.variants_table$HGNC<-HGNC_sel
      accomp.variants_table$vcf_VAR<-paste(accomp.variants_table$chr,accomp.variants_table$pos,accomp.variants$ref, accomp.variants$alt,sep="_")
      
      # cat("accomp.variants_table_\n")
      # str(accomp.variants_table)
      # cat("\n")
      
      DEF<-rbind(VAR_sel_table,accomp.variants_table)
      
      # cat("DEF_HELLO_WORLD\n")
      # str(DEF)
      # cat("\n")
      # 
      # quit(status=1)
      
    }else{
      
      
      DEF<-VAR_sel_table
      
      # cat("DEF_\n")
      # str(DEF)
      # cat("\n")
    }
    
    list_VARS[[i]]<-DEF
    
    
  }#i vars
  
  
  TABLE_vcf = unique(as.data.frame(data.table::rbindlist(list_VARS, fill = T)))
  
  cat("TABLE_vcf\n")
  str(TABLE_vcf)
  cat("\n")
  # quit(status = 1)
  
  TABLE_vcf_EROS = TABLE_vcf[which(TABLE_vcf$HGNC == "EROS"),]
  
  cat("TABLE_vcf_EROS\n")
  str(TABLE_vcf_EROS)
  cat("\n")
  
  
  

  #### merge with REF file ----
  
  TABLE_vcf<-merge(TABLE_vcf,
                        REF.df,
                        by=c("chr","HGNC","VAR"))
  
  cat("TABLE_vcf\n")
  str(TABLE_vcf)
  cat("\n")
  
  colnames(TABLE_vcf)[which(colnames(TABLE_vcf) == "seq_name")]<-"REF_interval"
  colnames(TABLE_vcf)[which(colnames(TABLE_vcf) == "sequence")]<-"REF_sequence"
  
  
  
  # quit(status = 1)
  
  #### LOOP HGNC PRINT ----
  
  vcf_header<-c('##fileformat=VCFv4.2',
                '##FILTER=<ID=PASS,Description="All filters passed">',
                '##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">',
                '##FILTER=<ID=LowQual,Description="Low quality">',
                '##GVCFBlock0-1=minGQ=0(inclusive),maxGQ=1(exclusive)',
                '##GVCFBlock1-2=minGQ=1(inclusive),maxGQ=2(exclusive)',
                '##GVCFBlock10-11=minGQ=10(inclusive),maxGQ=11(exclusive)',
                '##GVCFBlock11-12=minGQ=11(inclusive),maxGQ=12(exclusive)',
                '##GVCFBlock12-13=minGQ=12(inclusive),maxGQ=13(exclusive)',
                '##GVCFBlock13-14=minGQ=13(inclusive),maxGQ=14(exclusive)',
                '##GVCFBlock14-15=minGQ=14(inclusive),maxGQ=15(exclusive)',
                '##GVCFBlock15-16=minGQ=15(inclusive),maxGQ=16(exclusive)',
                '##GVCFBlock16-17=minGQ=16(inclusive),maxGQ=17(exclusive)',
                '##GVCFBlock17-18=minGQ=17(inclusive),maxGQ=18(exclusive)',
                '##GVCFBlock18-19=minGQ=18(inclusive),maxGQ=19(exclusive)',
                '##GVCFBlock19-20=minGQ=19(inclusive),maxGQ=20(exclusive)',
                '##GVCFBlock2-3=minGQ=2(inclusive),maxGQ=3(exclusive)',
                '##GVCFBlock20-21=minGQ=20(inclusive),maxGQ=21(exclusive)',
                '##GVCFBlock21-22=minGQ=21(inclusive),maxGQ=22(exclusive)',
                '##GVCFBlock22-23=minGQ=22(inclusive),maxGQ=23(exclusive)',
                '##GVCFBlock23-24=minGQ=23(inclusive),maxGQ=24(exclusive)',
                '##GVCFBlock24-25=minGQ=24(inclusive),maxGQ=25(exclusive)',
                '##GVCFBlock25-26=minGQ=25(inclusive),maxGQ=26(exclusive)',
                '##GVCFBlock26-27=minGQ=26(inclusive),maxGQ=27(exclusive)',
                '##GVCFBlock27-28=minGQ=27(inclusive),maxGQ=28(exclusive)',
                '##GVCFBlock28-29=minGQ=28(inclusive),maxGQ=29(exclusive)',
                '##GVCFBlock29-30=minGQ=29(inclusive),maxGQ=30(exclusive)',
                '##GVCFBlock3-4=minGQ=3(inclusive),maxGQ=4(exclusive)',
                '##GVCFBlock30-31=minGQ=30(inclusive),maxGQ=31(exclusive)',
                '##GVCFBlock31-32=minGQ=31(inclusive),maxGQ=32(exclusive)',
                '##GVCFBlock32-33=minGQ=32(inclusive),maxGQ=33(exclusive)',
                '##GVCFBlock33-34=minGQ=33(inclusive),maxGQ=34(exclusive)',
                '##GVCFBlock34-35=minGQ=34(inclusive),maxGQ=35(exclusive)',
                '##GVCFBlock35-36=minGQ=35(inclusive),maxGQ=36(exclusive)',
                '##GVCFBlock36-37=minGQ=36(inclusive),maxGQ=37(exclusive)',
                '##GVCFBlock37-38=minGQ=37(inclusive),maxGQ=38(exclusive)',
                '##GVCFBlock38-39=minGQ=38(inclusive),maxGQ=39(exclusive)',
                '##GVCFBlock39-40=minGQ=39(inclusive),maxGQ=40(exclusive)',
                '##GVCFBlock4-5=minGQ=4(inclusive),maxGQ=5(exclusive)',
                '##GVCFBlock40-41=minGQ=40(inclusive),maxGQ=41(exclusive)',
                '##GVCFBlock41-42=minGQ=41(inclusive),maxGQ=42(exclusive)',
                '##reference=file:///lustre/scratch116/vr/projects/hipsci/vrpipe/refs/human/ncbi37/hs37d5.fa',
                '##bcftools_concatVersion=1.3+htslib-1.3',
                '##bcftools_concatCommand=concat -a -D -f /lustre/scratch116/vr/projects/hipsci/vrpipe/a/8/4/3e0e3e374e1265c870395af411fa3/3053477/1_bcftools_concat/merge_list.txt',
                '##bcftools_annotateVersion=1.3+htslib-1.3',
                '##bcftools_annotateCommand=annotate -c CHROM,POS,ID,REF,ALT -a /lustre/scratch116/vr/projects/hipsci/vrpipe/refs/human/ncbi37/resources/annots-rsIDs-dbSNPv138.2014-04-02.tab.gz',
                '##bcftools_viewVersion=1.3.1+htslib-1.3.2',
                '##bcftools_viewCommand=view -Oz -o /lustre/scratch116/vr/projects/hipsci/releases/data/wgs/gatk_calls/HPSI0114i-kolf_2/HPSI0114i-kolf_2.wgs.gatk.haplotype_caller.20170425.genotypes.vcf.gz',
                '##bcftools_viewVersion=1.6+htslib-1.6',
                '##bcftools_viewCommand=view -R genERA_VAR_POST_HF1_2_3_region_file_FROM_TO.txt /lustre/scratch115/teams/soranzo/projects/genIE_analysis/reference_files/HPSI0114i-kolf_2.wgs.gatk.haplotype_caller.20170425.genotypes.vcf.gz; Date=Tue Oct 1 17:20:56 2019',
                paste('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO', sep="\t"))
    
    
  VAR_vector<-unique(TABLE_vcf$VAR)
  
  cat("VAR_vector_\n")
  str(VAR_vector)
  cat("\n")
  
  list_VAR_vector<-list()
  
  path_vcf<-paste(out,'vcf_files', sep='')
  
  if (file.exists(path_vcf)){
    
  } else {
    dir.create(file.path(path_vcf))
    
    
  }
  
  path_intervals<-paste(out,'intervals_files', sep='')
  
  if (file.exists(path_intervals)){
    
  } else {
    dir.create(file.path(path_intervals))
    
    
  }
  
  list_files<-list()
  list_final_table<-list()
  
  for(i in 1:length(VAR_vector))
  {
    VAR_sel<-VAR_vector[i]
    
    cat("------------------------------------->\t")
    cat(sprintf(as.character(VAR_sel)))
    cat("\n")
    
    TABLE_vcf_VAR_sel<-TABLE_vcf[which(TABLE_vcf$VAR == VAR_sel),]
    
    cat("TABLE_vcf_VAR_sel\n")
    str(TABLE_vcf_VAR_sel)
    cat("\n")
    
    HGNC_sel<-unique(TABLE_vcf_VAR_sel$HGNC)
    
    cat("------------------------------------->\t")
    cat(sprintf(as.character(HGNC_sel)))
    cat("\n")
    
    chr<-gsub("_.+$","",TABLE_vcf_VAR_sel$vcf_VAR)
    
    pos<-gsub("^[^_]+_","",TABLE_vcf_VAR_sel$vcf_VAR)
    
    pos<-as.integer(gsub("_.+$","",pos))
    ref<-gsub("^[^_]+_[^_]+_","",TABLE_vcf_VAR_sel$vcf_VAR)
    ref<-gsub("_.+$","",ref)
    alt<-gsub("^[^_]+_[^_]+_[^_]+_","",TABLE_vcf_VAR_sel$vcf_VAR)
    
    # cat("--VCF fields:--->\t")
    # cat(sprintf(as.character(chr)))
    # cat("\t")
    # cat(sprintf(as.character(pos)))
    # cat("\t")
    # cat(sprintf(as.character(ref)))
    # cat("\t")
    # cat(sprintf(as.character(alt)))
    # cat("\n")
    
    vcf_df<-as.data.frame(cbind(chr,pos,paste(chr,pos,ref,alt,sep="_"),ref,alt), stringsAsFactors=F)
    vcf_df$QUAL<-'100'
    vcf_df$FILTER<-'PASS'
    vcf_df$INFO<-'.'
    vcf_df$chr<-gsub("^chr","Chr",vcf_df$chr)
    
    cat("vcf_df\n")
    str(vcf_df)
    cat("\n")
    
    vcf_df_ordered<-vcf_df[order(vcf_df$pos, decreasing = F),]
    
    cat("vcf_df_ordered\n")
    str(vcf_df_ordered)
    cat("\n")
    
    if(dim(TABLE_vcf_VAR_sel)[1] == 1)
    {
      TABLE_vcf_VAR_sel$FLAG<-"NA"
      
      setwd(path_vcf)
      
      
      filename_vcf<-paste(paste(HGNC_sel,VAR_sel,sep="__"),".vcf",sep='')
      unlink(filename_vcf)
      
      # vector_files_vcf[i]<-filename_vcf
      
      
      
      #### vcf print
      
      cat(vcf_header,
          file=filename_vcf,sep="\n")
      write.table(vcf_df_ordered,
                  file=filename_vcf, sep="\t",append = T,
                  quote=F,col.names = F, row.names = F, eol="\n")
      
      ### interval file
      
      setwd(path_intervals)
      filename_interval<-paste(paste(HGNC_sel,VAR_sel,sep="__"),".intervals",sep='')
      unlink(filename_interval)
   
      interval_sel<-unique(TABLE_vcf_VAR_sel$REF_interval)
      cat(interval_sel,file=filename_interval,sep="\n")
      
      
      #master_df<-as.data.frame(cbind(VAR_sel,TABLE_vcf_VAR_sel$FLAG,paste(path_vcf,filename_vcf,sep="/"),paste(path_intervals,filename_interval,sep="/")), stringsAsFactors=F)
      master_df<-as.data.frame(cbind(VAR_sel,TABLE_vcf_VAR_sel$FLAG,filename_vcf,filename_interval), stringsAsFactors=F)
      
      
      colnames(master_df)<-c("VAR","FLAG","vcf_file","interval_file")
      cat("master_df\n")
      str(master_df)
      cat("\n")
      
      list_files[[i]]<-master_df
      
      # quit(status=1)
      
      #file_df<-as
      
      list_final_table[[i]]<-TABLE_vcf_VAR_sel
        
    }else{
      
      if(dim(TABLE_vcf_VAR_sel)[1] == 2)
      {
        TABLE_vcf_VAR_sel_accompanying_variants<-TABLE_vcf_VAR_sel[which(TABLE_vcf_VAR_sel$vcf_VAR != VAR_sel),]
        
        cat("TABLE_vcf_VAR_sel_accompanying_variants\n")
        str(TABLE_vcf_VAR_sel_accompanying_variants)
        cat("\n")
        
        FLAG<-TABLE_vcf_VAR_sel_accompanying_variants$id
        
        cat("------------------------------------->\t")
        cat(sprintf(as.character(FLAG)))
        cat("\n")
        
        FLAG<-gsub("^[^_]+_","",FLAG)
        FLAG<-gsub("_.+$","",FLAG)
        
        cat("------------------------------------->\t")
        cat(sprintf(as.character(FLAG)))
        cat("\n")
        
        if(FLAG == "HET")
        {
          TABLE_vcf_VAR_sel$FLAG<-"1_HET_polymorphism"
          
          #### vcf print
          
          setwd(path_vcf)
          
          filename_vcf_double_ALT<-paste(paste(HGNC_sel,VAR_sel,"double_ALT",sep="__"),".vcf",sep='')
          unlink(filename_vcf_double_ALT)
          
          
          cat(vcf_header,
              file=filename_vcf_double_ALT,sep="\n")
          write.table(vcf_df_ordered,
                      file=filename_vcf_double_ALT, sep="\t",append = T,
                      quote=F,col.names = F, row.names = F, eol="\n")
          
        
          
          ###
          
          vcf_df_ordered_HET_polymorphism_ALT<-vcf_df_ordered[which(vcf_df_ordered$V3%in%TABLE_vcf_VAR_sel_accompanying_variants$vcf_VAR),]
          unlink(vcf_df_ordered_HET_polymorphism_ALT)
          
          
          cat("vcf_df_ordered_HET_polymorphism_ALT\n")
          str(vcf_df_ordered_HET_polymorphism_ALT)
          cat("\n")
          
          filename_vcf_HET_polymorphism_ALT<-paste(paste(HGNC_sel,VAR_sel,"HET_polymorphism_ALT",sep="__"),".vcf",sep='')
          
          cat(vcf_header,
              file=filename_vcf_HET_polymorphism_ALT,sep="\n")
          write.table(vcf_df_ordered_HET_polymorphism_ALT,
                      file=filename_vcf_HET_polymorphism_ALT, sep="\t",append = T,
                      quote=F,col.names = F, row.names = F, eol="\n")
          
          ###
          
          vcf_df_ordered_VAR_ALT<-vcf_df_ordered[which(vcf_df_ordered$V3 == VAR_sel),]
          unlink(vcf_df_ordered_VAR_ALT)
          
          
          cat("vcf_df_ordered_VAR_ALT\n")
          str(vcf_df_ordered_VAR_ALT)
          cat("\n")
          
          filename_vcf_VAR_ALT<-paste(paste(HGNC_sel,VAR_sel,"VAR_ALT",sep="__"),".vcf",sep='')
          
          cat(vcf_header,
              file=filename_vcf_VAR_ALT,sep="\n")
          write.table(vcf_df_ordered_VAR_ALT,
                      file=filename_vcf_VAR_ALT, sep="\t",append = T,
                      quote=F,col.names = F, row.names = F, eol="\n")
          
          ### interval file
          
          filename_interval<-paste(paste(HGNC_sel,VAR_sel,sep="__"),".intervals",sep='')
          setwd(path_intervals)
          unlink(filename_interval)
          
          
          interval_sel<-unique(TABLE_vcf_VAR_sel$REF_interval)
          
          cat(interval_sel,file=filename_interval,sep="\n")
          
          
          
          master_df<-as.data.frame(cbind(VAR_sel,TABLE_vcf_VAR_sel$FLAG,
                                         rbind(filename_vcf_double_ALT,filename_vcf_HET_polymorphism_ALT,filename_vcf_VAR_ALT),
                                         filename_interval), stringsAsFactors=F)
          
          colnames(master_df)<-c("VAR","FLAG","vcf_file","interval_file")
          cat("master_df\n")
          str(master_df)
          cat("\n")
          
          list_files[[i]]<-master_df
          list_final_table[[i]]<-TABLE_vcf_VAR_sel
          
          # quit(status=1)
          
          
        }else{
          
          if(FLAG == "HOM")
          {
            TABLE_vcf_VAR_sel$FLAG<-"1_HOM_polymorphism"
            
            #### vcf print
            
            setwd(path_vcf)
            
            filename_vcf_double_ALT<-paste(paste(HGNC_sel,VAR_sel,"double_ALT",sep="__"),".vcf",sep='')
            unlink(filename_vcf_double_ALT)
            
            
            cat(vcf_header,
                file=filename_vcf_double_ALT,sep="\n")
            write.table(vcf_df_ordered,
                        file=filename_vcf_double_ALT, sep="\t",append = T,
                        quote=F,col.names = F, row.names = F, eol="\n")
            
            
            
            ###
            
            vcf_df_ordered_HOM_polymorphism_ALT<-vcf_df_ordered[which(vcf_df_ordered$V3%in%TABLE_vcf_VAR_sel_accompanying_variants$vcf_VAR),]
            unlink(vcf_df_ordered_HOM_polymorphism_ALT)
            
            
            cat("vcf_df_ordered_HOM_polymorphism_ALT\n")
            str(vcf_df_ordered_HOM_polymorphism_ALT)
            cat("\n")
            
            filename_vcf_HOM_polymorphism_ALT<-paste(paste(HGNC_sel,VAR_sel,"HOM_polymorphism_ALT",sep="__"),".vcf",sep='')
            
            cat(vcf_header,
                file=filename_vcf_HOM_polymorphism_ALT,sep="\n")
            write.table(vcf_df_ordered_HOM_polymorphism_ALT,
                        file=filename_vcf_HOM_polymorphism_ALT, sep="\t",append = T,
                        quote=F,col.names = F, row.names = F, eol="\n")
            
            ### interval file
            
            filename_interval<-paste(paste(HGNC_sel,VAR_sel,sep="__"),".intervals",sep='')
            setwd(path_intervals)
            unlink(filename_interval)
            
            
            interval_sel<-unique(TABLE_vcf_VAR_sel$REF_interval)
            
            cat(interval_sel,file=filename_interval,sep="\n")
            
            
            
            master_df<-as.data.frame(cbind(VAR_sel,TABLE_vcf_VAR_sel$FLAG,
                                           rbind(filename_vcf_double_ALT,filename_vcf_HOM_polymorphism_ALT),
                                           filename_interval), stringsAsFactors=F)
            
            colnames(master_df)<-c("VAR","FLAG","vcf_file","interval_file")
            cat("master_df\n")
            str(master_df)
            cat("\n")
            
            list_files[[i]]<-master_df
            list_final_table[[i]]<-TABLE_vcf_VAR_sel
            
            #quit(status=1)
            
            
          }else{
            
            cat("WARNING: unknown flag\t")
            quit(status=1)
          }
        }
      }else{
       
        
        cat("WARNING: >1 variants haplotype\t")
        cat("\n")
        
        TABLE_vcf_VAR_sel_accompanying_variants<-TABLE_vcf_VAR_sel[which(TABLE_vcf_VAR_sel$vcf_VAR != VAR_sel),]
        
        cat("TABLE_vcf_VAR_sel_accompanying_variants\n")
        str(TABLE_vcf_VAR_sel_accompanying_variants)
        cat("\n")
        
        FLAG<-TABLE_vcf_VAR_sel_accompanying_variants$id
        
        cat("------------------------------------->\t")
        cat(sprintf(as.character(FLAG)))
        cat("\n")
        
        FLAG<-gsub("^[^_]+_","",FLAG)
        FLAG<-gsub("_.+$","",FLAG)
        
        cat("------------------------------------->\t")
        cat(sprintf(as.character(FLAG)))
        cat("\n")
        
        unique_FLAG<-unique(FLAG)
        
        if(unique_FLAG == "HOM")
        {
          TABLE_vcf_VAR_sel$FLAG<-">1_HOM_polymorphism"
          
          
          #### vcf print
          
          setwd(path_vcf)
          
          filename_vcf_double_ALT<-paste(paste(HGNC_sel,VAR_sel,"ALL_ALT",sep="__"),".vcf",sep='')
          unlink(filename_vcf_double_ALT)
          
          
          cat(vcf_header,
              file=filename_vcf_double_ALT,sep="\n")
          write.table(vcf_df_ordered,
                      file=filename_vcf_double_ALT, sep="\t",append = T,
                      quote=F,col.names = F, row.names = F, eol="\n")
          
          
          
          ###
          
          vcf_df_ordered_HOM_polymorphism_ALT<-vcf_df_ordered[which(vcf_df_ordered$V3%in%TABLE_vcf_VAR_sel_accompanying_variants$vcf_VAR),]
          unlink(vcf_df_ordered_HOM_polymorphism_ALT)
          
          
          cat("vcf_df_ordered_HOM_polymorphism_ALT\n")
          str(vcf_df_ordered_HOM_polymorphism_ALT)
          cat("\n")
          
          filename_vcf_HOM_polymorphism_ALT<-paste(paste(HGNC_sel,VAR_sel,"HOM_polymorphism_ALT",sep="__"),".vcf",sep='')
          
          cat(vcf_header,
              file=filename_vcf_HOM_polymorphism_ALT,sep="\n")
          write.table(vcf_df_ordered_HOM_polymorphism_ALT,
                      file=filename_vcf_HOM_polymorphism_ALT, sep="\t",append = T,
                      quote=F,col.names = F, row.names = F, eol="\n")
          
          ### interval file
          
          filename_interval<-paste(paste(HGNC_sel,VAR_sel,sep="__"),".intervals",sep='')
          setwd(path_intervals)
          unlink(filename_interval)
          
          
          interval_sel<-unique(TABLE_vcf_VAR_sel$REF_interval)
          
          cat(interval_sel,file=filename_interval,sep="\n")
          
          
          
          master_df<-as.data.frame(cbind(VAR_sel,TABLE_vcf_VAR_sel$FLAG,
                                         rbind(filename_vcf_double_ALT,filename_vcf_HOM_polymorphism_ALT),
                                         filename_interval), stringsAsFactors=F)
          
          colnames(master_df)<-c("VAR","FLAG","vcf_file","interval_file")
          cat("master_df\n")
          str(master_df)
          cat("\n")
          
          list_files[[i]]<-master_df
          list_final_table[[i]]<-TABLE_vcf_VAR_sel
          
          # quit(status=1)
          
          
        }else{
          
          cat("WARNING: >1 mix of HOMS and HETS \t")
          quit(status=1)
        }
        
        
        
       
        
      
      }
    }
  }#i
  
  
  master_df_DEF = unique(as.data.frame(data.table::rbindlist(list_files, fill = T)))
  
  cat("master_df_DEF\n")
  str(master_df_DEF)
  cat("\n")
  
  
  intermediate_table = unique(as.data.frame(data.table::rbindlist(list_final_table, fill = T)))
  
  cat("intermediate_table\n")
  str(intermediate_table)
  cat("\n")
  
  
  # quit(status = 1)
  
  #### SAVE ----
  
  setwd(out)
  
  write.table(intermediate_table, file=paste(type,"_intermediate_table",".tsv", sep=''), sep="\t", row.names = F, quote=F)
  
  write.table(master_df_DEF, file=paste(type,"_vcf_interval_corresp",".tsv", sep=''), sep="\t", row.names = F, quote=F)
  
  
  
  
}

printList = function(l, prefix = "    ") {
  list.df = data.frame(val_name = names(l), value = as.character(l))
  list_strs = apply(list.df, MARGIN = 1, FUN = function(x) { paste(x, collapse = " = ")})
  cat(paste(paste(paste0(prefix, list_strs), collapse = "\n"), "\n"))
}

#### main script ----

main = function() {
  cmd_line = commandArgs()
  cat("Command line:\n")
  cat(paste(gsub("--file=", "", cmd_line[4], fixed=T),
            paste(cmd_line[6:length(cmd_line)], collapse = " "),
            "\n\n"))
  option_list <- list(
    make_option(c("--snps"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--SNPS_parsed"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--input_corrected_unblanked"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--REF_FASTA"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--bedfile"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--span"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="filename", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
    
  )
  parser = OptionParser(usage = "140__Rscript_v106.R
                        --subset type
                        --TranscriptEXP FILE.txt
                        --cadd FILE.txt
                        --ncboost FILE.txt
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
 parse_266(opt)
 Print_vcf(opt)
  
}


###########################################################################

system.time( main() )
