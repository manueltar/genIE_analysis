



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
suppressMessages(library("seqinr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("seqRFLP", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("splitstackshape", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))



opt = NULL

generate_fasta_reference_FOR_DESIGN = function(option_list)
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

  #### Read corresp_file ----
  
  corresp_file<-read.table(file=opt$corresp_file, sep=" ", header=T, stringsAsFactors = F)
  
  cat("corresp_file\n")
  cat(str(corresp_file))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(corresp_file$FLAG))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(corresp_file$FLAG)))))
  cat("\n")
  
  
  
  #### Read intermediate_table ----
  
  intermediate_table<-read.table(file=opt$intermediate_table, sep="\t", header=T, stringsAsFactors = F)
  
  cat("intermediate_table\n")
  cat(str(intermediate_table))
  cat("\n")
  
  check<-intermediate_table[-which(intermediate_table$VAR%in%corresp_file$VAR),]
  
  cat("check\n")
  cat(str(check))
  cat("\n")
  
 # quit(status=1)  
  
 target_table<-intermediate_table[which(intermediate_table$VAR == intermediate_table$vcf_VAR),]
 
 cat("target_table\n")
 cat(str(target_table))
 cat("\n")
 
 # quit(status=1)
  
  #### LOOP to read the ALT fasta files ----
  
 VARS<-unique(target_table$VAR)
 
 cat("VARS_\n")
 cat(str(VARS))
 cat("\n")
 
 list_files<-list()
 
  
 for(i in 1:length(VARS))
 {
   
   VAR_sel<-VARS[i]
   
   cat("--------->\t")
   cat(sprintf(as.character(VAR_sel)))
   cat("\n")
   
   target_table_sel<-target_table[which(target_table$VAR == VAR_sel),]
  
   cat("target_table_sel\n")
   cat(str(target_table_sel))
   cat("\n")
   
   indx.int<-c(which(colnames(target_table_sel) == "chr"),which(colnames(target_table_sel) == "HGNC"),which(colnames(target_table_sel) == "VAR"),
               which(colnames(target_table_sel) == "id"),which(colnames(target_table_sel) == "START"),which(colnames(target_table_sel) == "END"),
               which(colnames(target_table_sel) == "FLAG"),which(colnames(target_table_sel) == "REF_sequence"))
   
   
   Gather_table<-unique(target_table_sel[,indx.int])
   
   cat("Gather_table\n")
   cat(str(Gather_table))
   cat("\n")
   
     
   corresp_file_sel<-corresp_file[which(corresp_file$VAR == VAR_sel),]
   
   cat("corresp_file_sel\n")
   cat(str(corresp_file_sel))
   cat("\n")
   
   FLAG<-unique(target_table_sel$FLAG)
   
   cat("---->\t")
   cat(sprintf(as.character(FLAG)))
   cat("\n")
   
   pivotal_point<-sum(is.na(FLAG))
   
   cat("---->\t")
   cat(sprintf(as.character(pivotal_point)))
   cat("\n")
   
   path_intervals<-paste(out,'ALT_FA_sequences', sep='')
   
   if (file.exists(path_intervals)){
     setwd(path_intervals)
   } else {
     dir.create(file.path(path_intervals))
     setwd(path_intervals)
     
   }
   
   
   if(pivotal_point == 1)
   {
     ALT_fa_file<-corresp_file_sel$ALT_fa_file
     
     fastaFile<-readDNAStringSet(file=ALT_fa_file)

     # cat("fastaFile\n")
     # cat(str(fastaFile))
     # cat("\n")

     seq_name = names(fastaFile)
     sequence = paste(fastaFile)
     ALT_fasta <- data.frame(seq_name, sequence, stringsAsFactors = F)
     
     cat("ALT_fasta\n")
     cat(str(ALT_fasta))
     cat("\n")
     
     #### FASTA ALTERNATe 0-based???? 
     
     ALT_fasta$sequence<-substr(ALT_fasta$sequence,2,nchar(ALT_fasta$sequence))
     
     cat("ALT_fasta\n")
     cat(str(ALT_fasta))
     cat("\n")
     
     Gather_table$ALT_sequence<-ALT_fasta$sequence
     
     cat("Gather_table\n")
     cat(str(Gather_table))
     cat("\n")
     
     list_files[[i]]<-Gather_table
     
     cat("---nchar_REF->\t")
     cat(sprintf(as.character(nchar(Gather_table$REF_sequence))))
     cat("\n")
     cat("---nchar_ALT->\t")
     cat(sprintf(as.character(nchar(Gather_table$ALT_sequence))))
     cat("\n")
     
      # quit(status=1)
   }else{
     
     if(FLAG == "1_HET_polymorphism")
     {
       
       indx.double_alt<-grep("double_ALT_ALT",corresp_file_sel$ALT_fa_file)
       
       cat("--------->\t")
       cat(sprintf(as.character(indx.double_alt)))
       cat("\n")
       
       # quit(status=1)
       
       double_alt<-corresp_file_sel$ALT_fa_file[indx.double_alt]
       
       fastaFile<-readDNAStringSet(file=double_alt)
       
       # cat("fastaFile\n")
       # cat(str(fastaFile))
       # cat("\n")
       
       seq_name = names(fastaFile)
       sequence = paste(fastaFile)
       ALT_fasta <- data.frame(seq_name, sequence, stringsAsFactors = F)
       
       cat("ALT_fasta\n")
       cat(str(ALT_fasta))
       cat("\n")
       
       ALT_fasta$sequence<-substr(ALT_fasta$sequence,2,nchar(ALT_fasta$sequence))
       
       cat("ALT_fasta\n")
       cat(str(ALT_fasta))
       cat("\n")
       
       Gather_table$ALT_sequence<-ALT_fasta$sequence
       
       cat("Gather_table\n")
       cat(str(Gather_table))
       cat("\n")
       
       HET_polymorphism_ALT<-corresp_file_sel$ALT_fa_file[grep("HET_polymorphism_ALT",corresp_file_sel$ALT_fa_file)]
       
       fastaFile<-readDNAStringSet(file=HET_polymorphism_ALT)
       
       # cat("fastaFile\n")
       # cat(str(fastaFile))
       # cat("\n")
       
       seq_name = names(fastaFile)
       sequence = paste(fastaFile)
       REF_fasta <- data.frame(seq_name, sequence, stringsAsFactors = F)
       
       cat("REF_fasta\n")
       cat(str(REF_fasta))
       cat("\n")
       
       REF_fasta$sequence<-substr(REF_fasta$sequence,2,nchar(REF_fasta$sequence))
       
       cat("REF_fasta\n")
       cat(str(REF_fasta))
       cat("\n")
       
       Gather_table$REF_sequence<-REF_fasta$sequence
       
       cat("Gather_table\n")
       cat(str(Gather_table))
       cat("\n")
       
       list_files[[i]]<-Gather_table
       
       # quit(status=1)
       
     }else{
       
       if(FLAG == "1_HOM_polymorphism")
       {
         indx.double_alt<-grep("double_ALT_ALT",corresp_file_sel$ALT_fa_file)
         
         cat("--------->\t")
         cat(sprintf(as.character(indx.double_alt)))
         cat("\n")
         
         # quit(status=1)
         
         double_alt<-corresp_file_sel$ALT_fa_file[indx.double_alt]
         
         fastaFile<-readDNAStringSet(file=double_alt)
         
         # cat("fastaFile\n")
         # cat(str(fastaFile))
         # cat("\n")
         
         seq_name = names(fastaFile)
         sequence = paste(fastaFile)
         ALT_fasta <- data.frame(seq_name, sequence, stringsAsFactors = F)
         
         cat("ALT_fasta\n")
         cat(str(ALT_fasta))
         cat("\n")
         
         ALT_fasta$sequence<-substr(ALT_fasta$sequence,2,nchar(ALT_fasta$sequence))
         
         cat("ALT_fasta\n")
         cat(str(ALT_fasta))
         cat("\n")
         
         Gather_table$ALT_sequence<-ALT_fasta$sequence
         
         cat("Gather_table\n")
         cat(str(Gather_table))
         cat("\n")
         
         HOM_polymorphism_ALT<-corresp_file_sel$ALT_fa_file[grep("HOM_polymorphism_ALT",corresp_file_sel$ALT_fa_file)]
         
         fastaFile<-readDNAStringSet(file=HOM_polymorphism_ALT)
         
         # cat("fastaFile\n")
         # cat(str(fastaFile))
         # cat("\n")
         
         seq_name = names(fastaFile)
         sequence = paste(fastaFile)
         REF_fasta <- data.frame(seq_name, sequence, stringsAsFactors = F)
         
         cat("REF_fasta\n")
         cat(str(REF_fasta))
         cat("\n")
         
         REF_fasta$sequence<-substr(REF_fasta$sequence,2,nchar(REF_fasta$sequence))
         
         cat("REF_fasta\n")
         cat(str(REF_fasta))
         cat("\n")
         
         Gather_table$REF_sequence<-REF_fasta$sequence
         
         cat("Gather_table\n")
         cat(str(Gather_table))
         cat("\n")
         
         list_files[[i]]<-Gather_table
         
         # quit(status=1)
       }else{
         
         if(FLAG == ">1_HOM_polymorphism")
         {
           indx.ALL_ALT<-grep("ALL_ALT",corresp_file_sel$ALT_fa_file)
           
           cat("--------->\t")
           cat(sprintf(as.character(indx.ALL_ALT)))
           cat("\n")
           
           # quit(status=1)
           
           ALL_ALT<-corresp_file_sel$ALT_fa_file[indx.ALL_ALT]
           
           fastaFile<-readDNAStringSet(file=ALL_ALT)
           
           # cat("fastaFile\n")
           # cat(str(fastaFile))
           # cat("\n")
           
           seq_name = names(fastaFile)
           sequence = paste(fastaFile)
           ALT_fasta <- data.frame(seq_name, sequence, stringsAsFactors = F)
           
           cat("ALT_fasta\n")
           cat(str(ALT_fasta))
           cat("\n")
           
           ALT_fasta$sequence<-substr(ALT_fasta$sequence,2,nchar(ALT_fasta$sequence))
           
           cat("ALT_fasta\n")
           cat(str(ALT_fasta))
           cat("\n")
           
           Gather_table$ALT_sequence<-ALT_fasta$sequence
           
           cat("Gather_table\n")
           cat(str(Gather_table))
           cat("\n")
           
           HOM_polymorphism_ALT<-corresp_file_sel$ALT_fa_file[grep("HOM_polymorphism_ALT",corresp_file_sel$ALT_fa_file)]
           
           fastaFile<-readDNAStringSet(file=HOM_polymorphism_ALT)
           
           # cat("fastaFile\n")
           # cat(str(fastaFile))
           # cat("\n")
           
           seq_name = names(fastaFile)
           sequence = paste(fastaFile)
           REF_fasta <- data.frame(seq_name, sequence, stringsAsFactors = F)
           
           cat("REF_fasta\n")
           cat(str(REF_fasta))
           cat("\n")
           
           REF_fasta$sequence<-substr(REF_fasta$sequence,2,nchar(REF_fasta$sequence))
           
           cat("REF_fasta\n")
           cat(str(REF_fasta))
           cat("\n")
           
           Gather_table$REF_sequence<-REF_fasta$sequence
           
           cat("Gather_table\n")
           cat(str(Gather_table))
           cat("\n")
           
           list_files[[i]]<-Gather_table
           
           # quit(status=1)
           
         }else{
           
           cat("WARNING:UNKNOWn FLAG\n")
           quit(status=1)
         }
       }
     }
     
     
   }
 }#i
 
 master_df_DEF = unique(as.data.frame(data.table::rbindlist(list_files, fill = T)))
 
 master_df_DEF$seq_name<-paste(paste(master_df_DEF$VAR,master_df_DEF$HGNC, sep="__"),"haplotype_REF",sep="_")
 
 cat("master_df_DEF\n")
 str(master_df_DEF)
 cat("\n")
 
 #### FASTA Printing & copy for alignment----
 
 indx.fasta<-c(which(colnames(master_df_DEF) == "seq_name"),which(colnames(master_df_DEF) == "REF_sequence"))
 
 
 Merge2_TO_FASTA<-unique(master_df_DEF[,indx.fasta])
 
 cat("Merge2_TO_FASTA\n")
 str(Merge2_TO_FASTA)
 cat("\n")
 
 df.fasta = dataframe2fas(Merge2_TO_FASTA)
 
 # Here we print the fasta file
 
 # writeLines(df.fasta, sep ="\n", paste(type,'_alignment_check','.fasta',sep=''))
 # 
 # Merge2_TO_FASTA<-unique(master_df_DEF[,indx.fasta])
 # 
 # df.fasta = dataframe2fas(Merge2_TO_FASTA, file="df.fasta")
 
 # Here we print the fasta file
 setwd(out)
 writeLines(df.fasta, sep ="\n", paste(type,'_reference_FOR_DESIGN','.fasta',sep=''))
 
  
   #### SAVE MASTER TABLE ----
  
  
  
  setwd(out)
  
  filename<-paste(type,"_MASTER_FILE",".tsv", sep='')
  write.table(master_df_DEF, file=filename, sep="\t", row.names = F, quote=F)

 
  
  
  
}


generate_input_regions_file_and_reference_FASTA = function(option_list)
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
  
  #### Read master file ----
  
  master_file<-as.data.frame(fread(opt$master_file, header=T, sep="\t"), stringsAsFactors = F)
  
  master_file$pos.VAR<-gsub("^[^_]+_","",master_file$VAR)
  master_file$pos.VAR<-gsub("_.+$","",master_file$pos.VAR)
  master_file$relPOS<-as.integer(master_file$pos.VAR)-as.integer(master_file$START +1) +1 # necessary correction
  master_file$ref.VAR<-gsub("^[^_]+_[^_]+_","",master_file$VAR)
  master_file$ref.VAR<-gsub("_.+$","",master_file$ref.VAR)
  master_file$alt.VAR<-gsub("^[^_]+_[^_]+_[^_]+_","",master_file$VAR)
  
  
  
  
  
  master_file$RevComp_sequence<-revComp(master_file$REF_sequence)
  
  cat("master_file_\n")
  str(master_file)
  cat("\n")
  
  #### Read PRIMER_TABLE & TRANSFORM ----
  
  PRIMER_TABLE<-read.csv(opt$input_Eve, header=T, stringsAsFactors = F)
  
  # cat("PRIMER_TABLE_\n")
  # str(PRIMER_TABLE)
  # cat("\n")
  
  # quit(status = 1)
  
  PRIMER_TABLE_NO_NA<-PRIMER_TABLE[(PRIMER_TABLE$VAR != ""),]
  
  # cat("PRIMER_TABLE_NO_NA_0\n")
  # str(PRIMER_TABLE_NO_NA)
  # cat("\n")
  
  
  ##### NO "" HGNC -----
  
  PRIMER_TABLE_NO_NA<-PRIMER_TABLE_NO_NA[(PRIMER_TABLE_NO_NA$HGNC != ""),]
  
  cat("PRIMER_TABLE_NO_NA_1\n")
  str(PRIMER_TABLE_NO_NA)
  cat("\n")
  
  #### correct eros pos mistake
  
  PRIMER_TABLE_NO_NA$pos[PRIMER_TABLE_NO_NA$HGNC == "EROS"]<-"80407100"
  
  #### LOOP Gather ----
  
  VARS<-unique(master_file$VAR)
  
  cat("VARS_1\n")
  str(VARS)
  cat("\n")
  
  # vector.samples<-unique(as.numeric(as.character(PRIMER_TABLE_NO_NA$sample)))
  # SAMPLES<-sort(vector.samples, decreasing = T)
  # 
  # cat("--->SAMPLES\n")
  # cat(sprintf(as.character(SAMPLES)))
  # cat("\n")
  # 
  # #quit(status =1)
  # 
  Gather<- data.frame(matrix(vector(), 0, 14,
                             dimnames=list(c(),
                                           c("HGNC","rsid","strand",
                                             "FWD","REV","protospacer",
                                             "amplicon","amplicon_length",
                                             "EDIT_site","CUT_site",
                                             "obs","VAR","ref_string","alt_string"
                                           ))),
                      stringsAsFactors=F)
  
  
  for(i in 1:length(VARS))
  {
    
    VAR_sel<-VARS[i]
    
    cat("--------->\t")
    cat(sprintf(as.character(VAR_sel)))
    cat("\n")
    
    master_file_sel<-master_file[which(master_file$VAR == VAR_sel),]
    
    cat("master_file_sel\n")
    cat(str(master_file_sel))
    cat("\n")
    
    
    PRIMER_TABLE_NO_NA_sel<-PRIMER_TABLE_NO_NA[which(PRIMER_TABLE_NO_NA$VAR%in%master_file_sel$VAR &
                                                       PRIMER_TABLE_NO_NA$HGNC%in%master_file_sel$HGNC),]
    
    cat("--->PRIMER_TABLE_NO_NA_sel_\n")
    str(PRIMER_TABLE_NO_NA_sel)
    cat("\n")
    
    if(dim(PRIMER_TABLE_NO_NA_sel)[1] >0)
    {
      
    }else{
      
      if(VAR_sel == "chr17_80407100_A_G") # EROS MOCK VAR
      {
        PRIMER_TABLE_NO_NA_sel<-PRIMER_TABLE_NO_NA[which(PRIMER_TABLE_NO_NA$HGNC%in%master_file_sel$HGNC),]
        
        cat("--->PRIMER_TABLE_NO_NA_sel_1\n")
        str(PRIMER_TABLE_NO_NA_sel)
        cat("\n")
        
        FWD_vector<-as.character(PRIMER_TABLE_NO_NA_sel$Fwd.primer)
        
        cat("FWD_vector1\n")
        cat(sprintf(as.character(FWD_vector)))
        cat("\n")
        
        REV_vector<-as.character(PRIMER_TABLE_NO_NA_sel$Rev.Primer)
        
        cat("REV_vector1\n")
        cat(sprintf(as.character(REV_vector)))
        cat("\n")
        
        PRIMER_TABLE_NO_NA_sel$Rev.Primer<-FWD_vector
        PRIMER_TABLE_NO_NA_sel$Fwd.primer<-REV_vector
        
        # PRIMER_TABLE_NO_NA_sel<-PRIMER_TABLE_NO_NA_sel[1,]
        # 
        cat("--->PRIMER_TABLE_NO_NA_sel_2\n")
        str(PRIMER_TABLE_NO_NA_sel)
        cat("\n")
        
        
        
      }else{
        
        cat("--->UNKNOWN ERROR\n")
        break
        
      }
    }
    
   
    cat("--->PRIMER_TABLE_NO_NA_sel_\n")
    str(PRIMER_TABLE_NO_NA_sel)
    cat("\n")
    
    HGNC<-as.character(unique(PRIMER_TABLE_NO_NA_sel$HGNC))
    
    cat("HGNC\n")
    cat(sprintf(as.character(HGNC)))
    cat("\n")
    
    FWD<-unique(PRIMER_TABLE_NO_NA_sel$Fwd.primer)
    
    cat("FWD1\n")
    cat(sprintf(as.character(FWD)))
    cat("\n")
    
    
    FWD<-gsub(" ","",FWD)
    
    cat("FWD2\n")
    cat(sprintf(as.character(FWD)))
    cat("\n")
    
    REV<-unique(PRIMER_TABLE_NO_NA_sel$Rev.Primer)
    REV<-gsub(" ","",REV)
    REV.rev.comp<-revComp(REV)
    VAR<-unique(as.character(PRIMER_TABLE_NO_NA_sel$VAR))
    rsid<-unique(as.character(PRIMER_TABLE_NO_NA_sel$RSid))
    
    protospacer<-unique(PRIMER_TABLE_NO_NA_sel$Protospacer)
    protospacer<-gsub(" ","",protospacer)
    protospacer.rev.comp<-revComp(protospacer)
    
    check1<-data.frame(cbind(HGNC,FWD,REV,REV.rev.comp,VAR,protospacer,protospacer.rev.comp),
                       stringsAsFactors = F)
    
    cat("--->check1\n")
    str(check1)
    cat("\n")
    
    # if(VAR_sel == "chr17_80407100_A_G") # EROS MOCK VAR
    # {
    #   quit(status=1)
    # }
    # 
    
    if(!is.na(FWD) & FWD != "")
    {
      FASTA.index<-grep(FWD,master_file_sel$REF_sequence)
      cat("------------------------------------------------------->FASTA.index\n")
      cat(sprintf(as.character(FASTA.index)))
      cat("\n")
      cat(sprintf(as.character(length(FASTA.index))))
      cat("\n")
      
      
      if(length(FASTA.index)<1)
      {
        strand.FASTA<-"-"
        FASTA.index<-grep(FWD,master_file_sel$RevComp_sequence)
        cat("----------NEGATIVE STRAND----------------->FASTA.index\n")
        cat(sprintf(as.character(FASTA.index)))
        cat("\n")
        cat(sprintf(as.character(length(FASTA.index))))
        cat("\n")
        
        if(length(FASTA.index)<1)
        {
          
          rev_FWD<-revComp(FWD)
          FASTA.index<-grep(rev_FWD,master_file_sel$REF_sequence)
          cat("----------NEGATIVE STRAND2----------------->FASTA.index\n")
          cat(sprintf(as.character(FASTA.index)))
          cat("\n")
          cat(sprintf(as.character(length(FASTA.index))))
          cat("\n")
          
          cat("FWD\n")
          cat(sprintf(as.character(FWD)))
          cat("\n")
          
          cat("rev_FWD\n")
          cat(sprintf(as.character(rev_FWD)))
          cat("\n")
          
          cat("REV\n")
          cat(sprintf(as.character(REV)))
          cat("\n")
          
          cat("REV.rev.comp\n")
          cat(sprintf(as.character(REV.rev.comp)))
          cat("\n")
          
          cat("master_file_sel$REF_sequence\n")
          cat(sprintf(as.character(master_file_sel$REF_sequence)))
          cat("\n")
          
          cat("master_file_sel$RevComp_sequence\n")
          cat(sprintf(as.character(master_file_sel$RevComp_sequence)))
          cat("\n")
          
          # quit(status=1)
          
          
        }else{
          master_file_sel_sel<-master_file_sel[(FASTA.index),]
          rsid.FASTA<-unique(as.character(master_file_sel_sel$id))
          VAR.FASTA<-unique(as.character(master_file_sel_sel$VAR))
          relPOS.FASTA<-unique(as.character(master_file_sel_sel$relPOS))
          ref.FASTA<-revComp(unique(as.character(master_file_sel_sel$ref.VAR)))
          alt.FASTA<-revComp(unique(as.character(master_file_sel_sel$alt.VAR)))
          sequence<-master_file_sel_sel$RevComp_sequence
          
        }
        
      }else{
        
        strand.FASTA<-"+"
        master_file_sel_sel<-master_file_sel[(FASTA.index),]
        rsid.FASTA<-unique(as.character(master_file_sel_sel$id))
        VAR.FASTA<-unique(as.character(master_file_sel_sel$VAR))
        relPOS.FASTA<-unique(as.character(master_file_sel_sel$relPOS))
        ref.FASTA<-unique(as.character(master_file_sel_sel$ref.VAR))
        alt.FASTA<-unique(as.character(master_file_sel_sel$alt.VAR))
        sequence<-master_file_sel_sel$REF_sequence
      }
      
      
      
      
      check2<-data.frame(cbind(FASTA.index,rsid.FASTA,VAR.FASTA,
                               relPOS.FASTA,ref.FASTA,alt.FASTA,strand.FASTA,
                               sequence),
                         stringsAsFactors = F)
      
      cat("--->check2\n")
      str(check2)
      cat("\n")
      
      if(rsid.FASTA == "Control")
      {
        rsid <-"POSITIVECTRL"       
      }
      
      # quit(status=1)
      
      if(VAR_sel == VAR.FASTA) # check rsid
      {
        cat("--->HELLO_WORLD\n")
        
        
        #### Get the amplicon
        
        ## Break FWD
        Break.FWD<-cSplit(check2,"sequence", sep = FWD, fixed =T, drop=T)
        
        cat("--->Break.FWD\n")
        str(Break.FWD)
        cat("\n")
        
        Break.FWD.REV<-cSplit(Break.FWD,"sequence_2", sep = REV.rev.comp, 
                              fixed =T, drop=T)
        
        cat("--->Break.FWD.REV\n")
        str(Break.FWD.REV)
        cat("\n")
        
        offset_FWD<-nchar(as.character(Break.FWD.REV$sequence_1))
        
        cat("offset_FWD--->\t")
        cat(sprintf(as.character(offset_FWD)))
        cat("\n")
        
        
        relPOS.amplicon<-as.numeric(relPOS.FASTA)-offset_FWD
        
        cat("relPOS.amplicon--->\t")
        cat(sprintf(as.character(relPOS.amplicon)))
        cat("\n")
        
        #nchar(as.character(Break.FWD.REV$sequence_2_2))
        
        Break.FWD.REV$amplicon<-paste(FWD,Break.FWD.REV$sequence_2_1,REV.rev.comp, sep='')
        
        cat("--->Break.FWD.REV\n")
        str(Break.FWD.REV)
        cat("\n")
        
        
        
        amplicon_length<-nchar(as.character(Break.FWD.REV$amplicon))
        
        cat("amplicon_length--->\t")
        cat(sprintf(as.character(amplicon_length)))
        cat("\n")
        
        
        #### string ref and string alt
        
        length.prior<-relPOS.amplicon - 1
        
        
        cat("length.prior--->\t")
        cat(sprintf(as.character(length.prior)))
        cat("\n")
        
        dash_prior_SNP<-paste(rep("-",length.prior),collapse='')
        
        
        length.post<-unique(amplicon_length) - as.numeric(relPOS.amplicon)
        
        cat("length.post--->\t")
        cat(sprintf(as.character(length.post)))
        cat("\n")
        
        
        dash_post_SNP<-paste(rep("-",length.post),collapse='')
        
        ref_string<-paste(dash_prior_SNP,ref.FASTA,dash_post_SNP, sep='')
        
        cat("ref_string--->\t")
        cat(sprintf(as.character(ref_string)))
        cat("\n")
        
        #nchar(ref_string) == amplicon_length
        
        alt_string<-paste(dash_prior_SNP,alt.FASTA,dash_post_SNP, sep='')
        
        cat("alt_string--->\t")
        cat(sprintf(as.character(alt_string)))
        cat("\n")
        
        #nchar(alt_string) == amplicon_length
        
      # quit(status=1)
        
        
        
        ####  LOCATE CUT SITE ----
        
        ### Case 1 protospacer binds like REV primer
        
        case1<-length(grep(protospacer.rev.comp,Break.FWD.REV$amplicon))
        
        cat("case1--->\t")
        cat(sprintf(as.character(case1)))
        cat("\n")
        
        ### Case 2 protospacer binds like FWD primer
        
        case2<-sum(grep(protospacer,Break.FWD.REV$amplicon))
        
        cat("case2--->\t")
        cat(sprintf(as.character(case2)))
        cat("\n")
        
        if(case1 >= 1 & case2 ==0)
        {
          Break.FWD.REV<-cSplit(Break.FWD.REV,"amplicon", sep = protospacer.rev.comp, fixed =T, drop=F)
          
          # nchar(Break.FWD.REV$amplicon)
          # nchar(as.character(Break.FWD.REV$amplicon_1))
          # nchar(as.character(Break.FWD.REV$amplicon_2))
          
          CUT_relPOS<-unique(nchar(as.character(Break.FWD.REV$amplicon_1))+3)
          
          A<-as.data.frame(cbind(HGNC,rsid,strand.FASTA,FWD,REV,protospacer,
                                 Break.FWD.REV$amplicon,amplicon_length,
                                 relPOS.amplicon,CUT_relPOS,
                                 "CUT_SITE_protospacer_binds_like_REV_primer",
                                 VAR.FASTA,ref_string,alt_string))
          
          
         
          colnames(A)<-colnames(Gather)
          
          cat("A1\t")
          cat(str(A))
          cat("\n")
          
          
          Gather<-rbind(Gather,A)
          
         
          
        }else{
          
          if(case1 == 0 & case2 >= 1)
          {
            
            Break.FWD.REV<-cSplit(Break.FWD.REV,"amplicon", sep = protospacer, fixed =T, drop=F)
            
            # nchar(Break.FWD.REV$amplicon)
            # nchar(as.character(Break.FWD.REV$amplicon_1))
            # nchar(as.character(Break.FWD.REV$amplicon_2))
            
            CUT_relPOS<-nchar(as.character(Break.FWD.REV$amplicon_1))+nchar(protospacer)-2
            A<-as.data.frame(cbind(HGNC,rsid,strand.FASTA,FWD,REV,protospacer,
                                   Break.FWD.REV$amplicon,amplicon_length,
                                   relPOS.amplicon,CUT_relPOS,
                                   "CUT_SITE_protospacer_binds_like_FWD_primer",
                                   VAR.FASTA,ref_string,alt_string))
            
            
            colnames(A)<-colnames(Gather)
            
            cat("A1\t")
            cat(str(A))
            cat("\n")
            
            Gather<-rbind(Gather,A)
            
          }else{
            
            CUT_relPOS<-"NA"
            A<-as.data.frame(cbind(HGNC,rsid,strand.FASTA,FWD,REV,protospacer,
                                   Break.FWD.REV$amplicon,amplicon_length,
                                   relPOS.amplicon,CUT_relPOS,
                                   "not_found_cut_site",
                                   VAR.FASTA,ref_string,alt_string))
            
            
            colnames(A)<-colnames(Gather)
            
            cat("A3\t")
            cat(str(A))
            cat("\n")
            
            Gather<-rbind(Gather,A)
            
            
            
          }
        }
        
        # quit(status=1)
      } #VAR_sel == VAR.FASTA
      else{
        
        A<-as.data.frame(cbind(HGNC,rsid,strand.FASTA,FWD,REV,protospacer,"NA","NA","NA",
                               "NA","not_matching rsids",VAR.FASTA,"NA","NA"))
        
       
        
        colnames(A)  <-colnames(Gather)
        
        cat("A4\t")
        cat(str(A))
        cat("\n")
        
        Gather<-rbind(Gather,A)
      }
      
    }#!is.na(FWD)
    else{
      
      A<-as.data.frame(cbind(HGNC,rsid,"NA",FWD,REV,protospacer,"NA","NA","NA",
                             "NA","not_found","NA","NA","NA"))
      
     
      
      colnames(A)  <-colnames(Gather)
      
      cat("A5\t")
      cat(str(A))
      cat("\n")
      
      Gather<-rbind(Gather,A)
    }
    
    
    cat("--->Gather\n")
    str(Gather)
    cat("\n")

    # quit(status =1)
    
    #quit(status =1)
  }#i
  
  #### temporary solution for TNRC6A ----
  
  cat("Gather_1\n")
  str(Gather)
  cat("\n")
  
  Gather$start<-1
  Gather$CUT_site[(Gather$HGNC == "TNRC6A")]<-"144"
  
  cat("Gather_2\n")
  str(Gather)
  cat("\n")
  
  cat("master_file--->\t")
  cat(str(master_file))
  cat("\n")
  
  
  #### SAVE middle step ----
  
  common_colnames<-colnames(Gather)[which(colnames(Gather)%in%colnames(master_file))]
  
  cat("common_colnames--->\t")
  cat(sprintf(as.character(common_colnames)))
  cat("\n")
  
  Gather<-merge(Gather,
                master_file,
                by=common_colnames,
                all.x=T)
  
  cat("Gather--->\t")
  cat(str(Gather))
  cat("\n")
  
  
  
  # quit(status = 1)
  
  setwd(out)
  
  filename<-paste(type,"_input_regions_EXPANDED",".tsv",sep="")
  
  write.table(Gather, 
              file=filename, 
              sep="\t", quote=F, row.names = F)
  
  #### SAVE input regions
  
  cat("colnames(Gather)--->\t")
  cat(sprintf(as.character(colnames(Gather))))
  cat("\n")
  
  cat("colnames(master_file)--->\t")
  cat(sprintf(as.character(colnames(master_file))))
  cat("\n")
  
  
  setwd(out)
  
 
  
  # quit(status = 1)
  
  indexes.int<-c(which(colnames(Gather) == "HGNC"),
                 which(colnames(Gather) == "seq_name"),
                 which(colnames(Gather) == "start"),
                 which(colnames(Gather) == "amplicon_length"),
                 which(colnames(Gather) == "CUT_site"),
                 which(colnames(Gather) == "EDIT_site"),
                 which(colnames(Gather) == "ref_string"),
                 which(colnames(Gather) == "alt_string"),
                 which(colnames(Gather) == "amplicon"))
  
  region_file<-unique(Gather[,indexes.int])
  
  region_file2<-region_file
  
  colnames(region_file2)<-c("name","sequence_name","start",
                            "end","cut_site","highlight_site",
                            "wt_allele_profile","hdr_allele_profile",
                            "ref_sequence")
  
  filename<-paste(type,"_input_regions.tsv",sep="")
  
  write.table(region_file2, file=filename,
              sep="\t", quote=F,
              row.names = F)
  
  
  #### FASTA Printing & copy for alignment----
  
  # quit(status=1)
  
  indx.fasta<-c(which(colnames(Gather) == "seq_name"),which(colnames(Gather) == "amplicon"))
  
  
  Merge2_TO_FASTA<-unique(Gather[,indx.fasta])
  
  Merge2_TO_FASTA_NO_NA<-Merge2_TO_FASTA[!is.na(Merge2_TO_FASTA$seq_name),]
  
  cat("Merge2_TO_FASTA_NO_NA\n")
  str(Merge2_TO_FASTA_NO_NA)
  cat("\n")
  
  df.fasta = dataframe2fas(Merge2_TO_FASTA_NO_NA)
  
  # Here we print the fasta file
  
  # writeLines(df.fasta, sep ="\n", paste(type,'_alignment_check','.fasta',sep=''))
  # 
  # Merge2_TO_FASTA_NO_NA<-unique(Gather[,indx.fasta])
  # 
  # df.fasta = dataframe2fas(Merge2_TO_FASTA_NO_NA, file="df.fasta")
  
  # Here we print the fasta file
  setwd(out)
  writeLines(df.fasta, sep ="\n", paste(type,'_reference','.fasta',sep=''))
  
  cat("THE END\n")
  
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
    make_option(c("--corresp_file"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--intermediate_table"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--master_file"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--input_Eve"), type="character", default=NULL,
                metavar="FILE.txt",
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
  
  
 generate_fasta_reference_FOR_DESIGN(opt)
 generate_input_regions_file_and_reference_FASTA(opt)
  
}


###########################################################################

system.time( main() )
