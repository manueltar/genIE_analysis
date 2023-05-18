#!/usr/bin/env Rscript

#### libraries ----
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
suppressMessages(library("vroom", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))



#### functions ----

opt = NULL

options(warn=1)


runCounts = function(option_list)
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
  
  
  files.df = as.data.frame(fread(opt$files, sep="\t", header=F), stringsAsFactors=F)
  colnames(files.df)<-"route"
  
  cat("files.df_\n")
  str(files.df)
  cat("\n")
  
  
  
  ### substitute the route and get the file name
  
  A<-unlist(strsplit(as.character(files.df$route), "/"))
  indx<-grep("count", A)
  B<-A[indx]
  
  cat("Manuel_2_\n")
  cat(sprintf(as.character(B)))
  cat("\n")
  
  # quit(status=1)
  
  ### Read the fasta file
  
  fastaFile <- read.fasta(file=opt$fasta)

  cat("fastaFile_\n")
  str(fastaFile)
  cat("\n")


  seq_name = fastaFile[seq(1, length(fastaFile), by=2)]
  sequence = fastaFile[seq(2, length(fastaFile), by=2)]
  df <- data.frame(seq_name, sequence, stringsAsFactors = F)

  colnames(df)<-c("seq_name","sequence")
  
  df$seq_name<-gsub(">","",df$seq_name)

  cat("files.df_\n")
  str(files.df)
  cat("\n")
  cat("df$seq_name\n")
  cat(sprintf(as.character(df$seq_name)))
  cat("\n")
  
  #### LOOP TO READ THE FILES ----
  
  seq_names_from_fasta<-c(as.character(df$seq_name),"UNDETERMINED")
  
  cat("seq_names_from_fasta\n")
  cat(sprintf(as.character(seq_names_from_fasta)))
  cat("\n")
  cat(sprintf(as.character(length(seq_names_from_fasta))))
  cat("\n")
  
  # Gather<- data.frame(matrix(vector(), 0, sum(1+length(seq_names_from_fasta)),
  #                            dimnames=list(c(), 
  #                                          c("sample",seq_names_from_fasta
  #                                          ))),
  #                     stringsAsFactors=F)
  # 
  # factor_order<-factor(colnames(Gather),
  #                      levels=colnames(Gather),
  #                      ordered = T)
  # 
  # cat("Gather\n")
  # cat(str(Gather))
  # cat("\n")
  
  list_gather<-list()
  
  TOTAL_files<-unique(B)
  
  cat("TOTAL_files\n")
  cat(sprintf(TOTAL_files))
  cat("\n")
  
  for(i in 1:length(TOTAL_files))
  {
    file_sel<-TOTAL_files[i]
    
    FLAG_exist<-file.exists(file_sel)

    cat("FLAG_exist/file_sel\n")
    cat(sprintf(as.character(FLAG_exist)))
    cat("\n")
    cat(sprintf(as.character(file_sel)))
    cat("\n")
    
    #quit(status = 1)
    
    if(FLAG_exist == TRUE)
    {
      info = file.info(file_sel)
      FLAG_empty <-(info$size == 0)
      
      cat("FLAG_empty?\n")
      cat(sprintf(as.character(FLAG_empty)))
      cat("\n")
      cat(sprintf(as.character(info)))
      cat("\n")
      
      #quit(status = 1)
      
      if(FLAG_empty != TRUE)
      {
        
        cat("------------------------> FILE exists & not empty\n")
        cat("\n")
        
        
        file_table<-as.data.frame(fread(file=file_sel, sep=" ", header=F, stringsAsFactors = F))
        
        colnames(file_table)<-c("reads","seq_name")
        
        # cat("file_table_PRE\n")
        # str(file_table)
        # cat("\n")
        
        ### Substitute *  
        
        file_table$seq_name[which(file_table$seq_name == '*')]<-"UNDETERMINED"
        
        # cat("file_table_POST\n")
        # str(file_table)
        # cat("\n")
        
        file_table$VAR<-gsub("\\|.+$","",file_table$seq_name)
        file_table$rsid<-gsub("^[^\\|]+\\|","",file_table$seq_name)
        file_table$rsid<-gsub("\\|.+$","",file_table$rsid)
        file_table$HGNC<-gsub("^[^\\|]+\\|[^\\|]+\\|","",file_table$seq_name)
        file_table$HGNC<-gsub("\\|.+$","",file_table$HGNC)
        file_table$START<-gsub("^[^\\|]+\\|[^\\|]+\\|[^\\|]+\\|","",file_table$seq_name)
        file_table$START<-gsub("\\|.+$","",file_table$START)
        file_table$STOP<-gsub("^[^\\|]+\\|[^\\|]+\\|[^\\|]+\\|","",file_table$seq_name)
        file_table$STOP<-gsub("\\|.+$","",file_table$STOP)
        
        file_table$amplicon<-paste(file_table$HGNC,file_table$VAR, sep="__")
        file_table$amplicon[(file_table$HGNC == "UNDETERMINED")]<-"UNDETERMINED"
          
        file_table$file<-file_sel
        file_table$sample<-gsub("_count.txt","",file_table$file)
       
        cat("file_table_2\n")
        str(file_table)
        cat("\n")
        
        indx.int<-c(which(colnames(file_table) == "amplicon"),
                    which(colnames(file_table) == "file"),
                    which(colnames(file_table) == "sample"),
                    which(colnames(file_table) == "reads"))
        
        
        file_table_subset<-unique(file_table[,indx.int])
        
        # cat("file_table_subset\n")
        # str(file_table_subset)
        # cat("\n")
        # 
        # cat("file_table\n")
        # str(file_table)
        # cat("\n")
        
        # if(file_table$HGNC == "NA")
        # {
        #   
        #   quit(status = 1)
        # }
        
        # 
        
                
        list_gather[[i]]<-file_table_subset
       
      }
      
      
    }# FLAG_exist == TRUE
    
  } #i
  
  Gather<-unique(as.data.frame(data.table::rbindlist(list_gather, fill = T)))
  
  cat("Gather_FINAL\n")
  cat(str(Gather))
  cat("\n")
  
  #### check GA2
  
  # quit(status = 1)
  ####
  
  amplicon_order<-sort(unique(Gather$amplicon))
  
  cat("amplicon_order\n")
  cat(str(amplicon_order))
  cat("\n")
  
  
  Gather$amplicon<-factor(Gather$amplicon,
                          levels=amplicon_order,
                          ordered=T)
  
  cat("Gather_2\n")
  cat(str(Gather))
  cat("\n")
  cat(sprintf(as.character(names(summary(amplicon_order)))))
  cat("\n")
  cat(sprintf(as.character(summary(amplicon_order))))
  cat("\n")
  
  
  Gather_wide<-as.data.frame(pivot_wider(Gather[order(Gather$amplicon),],
                                         id_cols=c("file","sample"),
                                         names_from=amplicon,
                                         values_from=reads), stringsAsFactors=F)
  
  Gather_wide[is.na(Gather_wide)]<-0
  
  cat("Gather_wide\n")
  cat(str(Gather_wide))
  cat("\n")
  
  TOTAL_files_absent<-TOTAL_files[-which(TOTAL_files%in%Gather_wide$file)]
  
  cat("TOTAL_files_absent\n")
  cat(str(TOTAL_files_absent))
  cat("\n")
  
  if(length(TOTAL_files_absent) >0)
  {
    
    cat("ABSENT FILES\n")
    
    file_vector<-TOTAL_files_absent
    sample_vector<-gsub("_count.txt","",file_vector)
    
    
    ABSENT<- data.frame(matrix(vector(), length(file_vector), dim(Gather_wide)[2],
                               dimnames=list(c(),
                                             colnames(Gather_wide))),
                        stringsAsFactors=F)
    ABSENT$file<-file_vector
    ABSENT$sample<-sample_vector
    
    colnames(ABSENT)<-gsub("\\.","-",colnames(ABSENT))
    
    ABSENT[is.na(ABSENT)]<-0
    
    
    cat("ABSENT\n")
    cat(str(ABSENT))
    cat("\n")
    
    
    
    check.ABSENT<-colnames(ABSENT)[-which(colnames(ABSENT)%in%colnames(Gather_wide))]
    
    cat("check.ABSENT\n")
    cat(sprintf(as.character(check.ABSENT)))
    cat("\n")
    
    check.Gather_wide<-colnames(Gather_wide)[-which(colnames(Gather_wide)%in%colnames(ABSENT))]
    
    cat("check.Gather_wide\n")
    cat(sprintf(as.character(check.Gather_wide)))
    cat("\n")
    
    
    Gather_wide<-rbind(Gather_wide,ABSENT)
    
    cat("Gather_wide\n")
    cat(str(Gather_wide))
    cat("\n")
    
    # quit(status = 1)
    
    
  }
  
  HGNC_vector<-colnames(Gather_wide[-c(which(colnames(Gather_wide)  == "file"),which(colnames(Gather_wide)  == "sample"))])
  
  cat("HGNC_vector\n")
  cat(sprintf(as.character(HGNC_vector)))
  cat("\n")
  
  HGNC_vector<-gsub("__.+","",HGNC_vector)
  
  cat("HGNC_vector_2\n")
  cat(sprintf(as.character(HGNC_vector)))
  cat("\n")
  
  ### sample and name
  
  
  Gather_wide$sample<-gsub(paste(HGNC_vector,collapse="|"),"",Gather_wide$sample)
  Gather_wide$sample<-gsub("^_","",Gather_wide$sample)
  
  
  cat("Gather_wide$sample\n")
  cat(sprintf(as.character(Gather_wide$sample)))
  cat("\n")
  
  # quit(status=1)
  
  
  Gather_wide$name<-gsub("_.+$","",Gather_wide$sample)
  
  cat("Gather_wide$name\n")
  cat(sprintf(as.character(Gather_wide$name)))
  cat("\n")
  
  
  Gather_wide$Experiment<-gsub("^[^_]+_","",Gather_wide$sample)
  Gather_wide$Experiment<-paste("Experiment_",gsub("_.+$","",Gather_wide$Experiment),sep='')
  
  cat("Gather_wide$Experiment\n")
  cat(sprintf(as.character(Gather_wide$Experiment)))
  cat("\n")
  
  
 
  
  Gather_wide$type<-gsub("^[^_]+_[^_]+_","",Gather_wide$sample)
  Gather_wide$type<-gsub("^[^_]+_","",Gather_wide$type)
  
  Gather_wide$type<-gsub("[0-9]+$","",Gather_wide$type)
  
  cat("Gather_wide$type\n")
  cat(sprintf(as.character(Gather_wide$type)))
  cat("\n")
  

  Gather_wide$Replicate<-gsub("^[^_]+_[^_]+_","",Gather_wide$sample)
  Gather_wide$Replicate<-as.integer(gsub("^[^0-9]+","",Gather_wide$Replicate))
  
  cat("Gather_wide$Replicate\n")
  cat(sprintf(as.character(Gather_wide$Replicate)))
  cat("\n")
 
  
  cat("Gather_wide_3\n")
  cat(str(Gather_wide))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Gather_wide$type))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Gather_wide$type)))))
  cat("\n")
  
  
  Gather_wide$type[which(Gather_wide$type == "Apcr")]<-"ATAC"
  Gather_wide$type[which(Gather_wide$type == "Gpcr")]<-"gDNA"
  Gather_wide$type[which(Gather_wide$type == "gpcr")]<-"gDNA"
  
  Gather_wide$type<-factor(Gather_wide$type,
                           levels=c("gDNA","ATAC"),
                           ordered=T)
  
  cat("Gather_wide_4\n")
  cat(str(Gather_wide))
  cat("\n")
  
  
  Gather_wide$HGNC<-gsub(paste(paste("_",Gather_wide$sample,"_count.txt",sep=''),collapse="|"),"",Gather_wide$file)
  
  cat("Gather_wide$HGNC\n")
  cat(sprintf(as.character(Gather_wide$HGNC)))
  cat("\n")
  
  # setwd(out)
  # write.table(Gather_wide,file="test.tsv", sep="\t", quote=F, row.names = F)
  
  
  
  
  Gather_wide_ordered <- Gather_wide[order(Gather_wide$name,Gather_wide$type, Gather_wide$Replicate, decreasing=T),]
  
  
  
  Gather_wide_ordered$Replicate<-paste(Gather_wide_ordered$type,Gather_wide_ordered$Replicate,sep="_")
  Gather_wide_ordered$Sample_Name<-Gather_wide_ordered$sample
  
  Gather_wide_ordered$sample<-paste(Gather_wide_ordered$HGNC,Gather_wide_ordered$sample,sep="_")
  
  
  
  
  
  sample_order<-unique(gsub("_.+","",Gather_wide_ordered$sample))
  
  
  
  # cat("sample_order_COLUMN_ORDER\n")
  # cat(sprintf(as.character(sample_order)))
  # cat("\n")
  
  cat("Gather_wide_ordered\n")
  cat(str(Gather_wide_ordered))
  cat("\n")
  
  check.Gather_wide_ordered<-Gather_wide_ordered[which(Gather_wide_ordered$HGNC == "NA"),]
  
  cat("check.Gather_wide_ordered\n")
  cat(str(check.Gather_wide_ordered))
  cat("\n")
  
  # quit(status=1)
 
  
  
  #### BIAS!! Add OLD PREFIX FOR Eve ----
  
  PREFIXES_Table_EXPANDED = fread(opt$prefixes)
  
  PREFIXES_Table_EXPANDED$sample<-paste(PREFIXES_Table_EXPANDED$HGNC, 
                                       PREFIXES_Table_EXPANDED$Sample_Name,
                                       sep="_")
  
  
  cat("PREFIXES_Table_EXPANDED_FINAL\n")
  cat(str(PREFIXES_Table_EXPANDED))
  cat("\n")
  cat("HGNC\n")
  cat(sprintf(as.character(levels(as.factor(PREFIXES_Table_EXPANDED$HGNC)))))
  cat("\n")
  
  check.PREFIXES_Table_EXPANDED<-PREFIXES_Table_EXPANDED[which(PREFIXES_Table_EXPANDED$VAR == "chr13_28674797_G_A"),]
  
  cat("check.PREFIXES_Table_EXPANDED\n")
  cat(str(check.PREFIXES_Table_EXPANDED))
  cat("\n")
  
  
  cat("Gather_wide_ordered$sample\n")
  cat(sprintf(as.character(Gather_wide_ordered$sample)))
  cat("\n")
  
  check.Gather_wide_ordered<-Gather_wide_ordered[which(Gather_wide_ordered$HGNC == "NA"),]
  
  cat("check.Gather_wide_ordered\n")
  cat(str(check.Gather_wide_ordered))
  cat("\n")
  
  
  # quit(status =1)
  
  DEF<-merge(PREFIXES_Table_EXPANDED, 
             Gather_wide_ordered[,-which(colnames(Gather_wide_ordered) == "HGNC")],
             by=c("sample","Sample_Name","Replicate","type","Experiment"),
             all.y=T)
  
  DEF$HGNC[is.na(DEF$HGNC)]<-"NA"
  
  cat("DEF_FINAL\n")
  cat(str(DEF))
  cat("\n")
  cat("HGNC\n")
  cat(sprintf(as.character(levels(as.factor(DEF$HGNC)))))
  cat("\n")
  
  check.DEF<-DEF[which(DEF$VAR == "chr13_28674797_G_A"),]
  
  cat("check.DEF\n")
  cat(str(check.DEF))
  cat("\n")
  
  # quit(status =1)
  
 #### SAVE ----
  
  setwd(out)
  
  file_name=paste("Count",opt$type,"matrix.txt",sep = "_")
  
  write.table(DEF, file=file_name, sep="\t", quote=F, row.names = F)
  
  file_name=paste("Count",opt$type,"matrix.rds",sep = "_")
  
  saveRDS(DEF, file=file_name)
  
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
    make_option(c("--files"), type="character", default=NULL, 
                metavar="FILE.tsv", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--fasta"), type="character", default=NULL, 
                metavar="FILE.fasta", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--prefixes"), type="character", default=NULL, 
                metavar="FILE.txt", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="option", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="option", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
    
      )
  
  parser = OptionParser(usage = "genie.R --files file.tsv --fasta file.fasta --prefixes file.txt -- type option",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  runCounts(opt)
  
}


###########################################################################

system.time( main() )

