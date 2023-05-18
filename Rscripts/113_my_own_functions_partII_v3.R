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

runCounts = function(option_list)
{
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
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
  
  Gather<- data.frame(matrix(vector(), 0, sum(1+length(seq_names_from_fasta)),
                             dimnames=list(c(), 
                                           c("sample",seq_names_from_fasta
                                           ))),
                      stringsAsFactors=F)
  
  factor_order<-factor(colnames(Gather),
                       levels=colnames(Gather),
                       ordered = T)
  
  TOTAL_files<-unique(B)
  
  cat("TOTAL_files\n")
  cat(sprintf(TOTAL_files))
  cat("\n")
  
  for(i in 1:length(TOTAL_files))
  {
    file_sel<-TOTAL_files[i]
    
    FLAG_exist<-file.exists(file_sel)
    
    # cat("FLAG_exist/file_sel\n")
    # cat(sprintf(as.character(FLAG_exist)))
    # cat("\n")
    # cat(sprintf(as.character(file_sel)))
    # cat("\n")
    
    #quit(status = 1)
    
    if(FLAG_exist == TRUE)
    {
      info = file.info(file_sel)
      FLAG_empty <-(info$size == 0)
      
      # cat("FLAG_empty?\n")
      # cat(sprintf(as.character(FLAG_empty)))
      # cat("\n")
      # cat(sprintf(as.character(info)))
      # cat("\n")
      
      #quit(status = 1)
      
      if(FLAG_empty != TRUE)
      {
        
        cat("------------------------> FILE exists & not empty\n")
        cat("\n")
        
        
        file_table<-as.data.frame(fread(file=file_sel, sep=" ", header=F, stringsAsFactors = F))
        
        colnames(file_table)<-c("reads","seq_name")
        
        cat("file_table_PRE\n")
        str(file_table)
        cat("\n")
        
        ### Substitute *  
        
        file_table$seq_name[which(file_table$seq_name == '*')]<-"UNDETERMINED"
        
        # cat("file_table_POST\n")
        # str(file_table)
        # cat("\n")
        
        #quit(status = 1)
        
        ## Merge tables to create a regular matrix
        
        
        PRESENT<-data.frame(matrix(as.numeric(file_table$reads), 
                                   1, length(as.character(file_table$seq_name)),
                                   dimnames=list(c(), 
                                                 c(as.character(file_table$seq_name)
                                                 ))),
                            stringsAsFactors=F)
        
        # cat("PRESENT\n")
        # str(PRESENT)
        # cat("\n")
        
        ## add 0 in the seq-names ABSENT
        
        ABSENT<-seq_names_from_fasta[-which(seq_names_from_fasta%in%colnames(PRESENT))]
        
        ABSENT.table<-data.frame(matrix(rep(0,length(ABSENT)), 1, length(ABSENT),
                                        dimnames=list(c(), 
                                                      c(ABSENT
                                                      ))),
                                 stringsAsFactors=F)
        
        
        colnames(ABSENT.table)<-gsub("X\\.","", colnames(ABSENT.table))
        
        # cat("ABSENT.table\n")
        # str(ABSENT.table)
        # cat("\n")
        
        #quit(status = 1)
        
        
        
        MERGE<-cbind(PRESENT,ABSENT.table)

        # cat("MERGE_PRE\n")
        # str(MERGE)
        # cat("\n")
        
        ## add sample
        
        MERGE$sample<-file_sel
        MERGE$sample<-gsub("_count.txt","",MERGE$sample)

        # cat("MERGE_POST\n")
        # str(MERGE)
        # cat("\n")
        
        MERGE_ordered<-MERGE[,order(factor_order)]
                
        
        
        ## Final merge with Gather
        
        Gather<-rbind(Gather,MERGE_ordered)
        
        # cat("Gather_POST\n")
        # str(Gather)
        # cat("\n")

        #quit(status = 1)

      }else{
        
        ABSENT<-seq_names_from_fasta
        
        ABSENT.table<-data.frame(matrix(rep(0,length(ABSENT)), 1, length(ABSENT),
                                        dimnames=list(c(), 
                                                      c(ABSENT
                                                      ))),
                                 stringsAsFactors=F)
        
        ABSENT.table$sample<-file_sel
        ABSENT.table$sample<-gsub("_count.txt","",ABSENT.table$sample)
        
        ABSENT.table_ordered<-ABSENT.table[,colnames(Gather)]
        
        ## Final ABSENT.table with Gather
        
        Gather<-rbind(Gather,ABSENT.table_ordered)
      }#file empty
    }else{
      
      
      ABSENT<-seq_names_from_fasta
      
      ABSENT.table<-data.frame(matrix(rep(0,length(ABSENT)), 1, length(ABSENT),
                                      dimnames=list(c(), 
                                                    c(ABSENT
                                                    ))),
                               stringsAsFactors=F)
      
      ABSENT.table$sample<-file_sel
      ABSENT.table$sample<-gsub("_count.txt","",ABSENT.table$sample)
      
      ABSENT.table_ordered<-ABSENT.table[,colnames(Gather)]
      
      ## Final ABSENT.table with Gather
      
      Gather<-rbind(Gather,ABSENT.table_ordered)
    } # file exists
    
  } #i
  
  cat("Gather_FINAL\n")
  cat(str(Gather))
  cat("\n")
  
  #### BIAS!!! Arrange and order ---
  
  Gather_arranged <- arrange(Gather,sample)
  
  cat("------------------------------------------------------------------------------------------------------->Gather_arranged_FINAL\n")
  cat(str(Gather_arranged))
  cat("\n")
  
  sample_order<-unique(gsub("_.+","",Gather_arranged$sample))
  
  column_order<-colnames(Gather_arranged)[-c(which(colnames(Gather_arranged) == "sample"),
                                             which(colnames(Gather_arranged) == "UNDETERMINED"))]
  
  cat("sample_order_COLUMN_ORDER\n")
  cat(sprintf(as.character(sample_order)))
  cat("\n")
  
  cat(sprintf(as.character(column_order)))
  cat("\n")
  
  # quit(status=1)
 
  
  column_order_2<-sort(column_order)
  cat("column_order_2_\n")
  cat(sprintf(as.character(column_order_2)))
  cat("\n")
  cat(sprintf(as.character(length(column_order_2))))
  cat("\n")
  
  column_order_3<-column_order_2
  
  column_order_3.5<-colnames(Gather_arranged)[which(colnames(Gather_arranged) == "sample")]
  
  column_order_3.75<-colnames(Gather_arranged)[which(colnames(Gather_arranged) == "UNDETERMINED")]
  
  column_order_4<-c(column_order_3.5,column_order_3,column_order_3.75)
  
  cat("column_order_4\n")
  cat(sprintf(as.character(column_order_4)))
  cat("\n")
  
  #
  
  Gather_arranged$sample<-gsub("-","_",Gather_arranged$sample)
    
  colnames(Gather_arranged)
  Gather_arranged_ordered<-Gather_arranged[,column_order_4]
  
  cat("Gather_arranged_ordered_FINAL\n")
  cat(str(Gather_arranged_ordered))
  cat("\n")
  
  # quit(status=1)
 
  
  
  #### BIAS!! Add OLD PREFIX FOR Eve ----
  
  PREFIXES_Table_EXPANDED = fread(opt$prefixes)
  
  PREFIXES_Table_EXPANDED<-PREFIXES_Table_EXPANDED[!is.na(PREFIXES_Table_EXPANDED$Sample_Name),]
  colnames(PREFIXES_Table_EXPANDED)[which(colnames(PREFIXES_Table_EXPANDED) == "hgnc")]<-"HGNC"
  
  PREFIXES_Table_EXPANDED$Sample_Name<-gsub("-","_",PREFIXES_Table_EXPANDED$Sample_Name)
  
  
  
  
  PREFIXES_Table_EXPANDED$sample<-paste(PREFIXES_Table_EXPANDED$HGNC, 
                                       PREFIXES_Table_EXPANDED$Sample_Name,
                                       sep="_")
  
  
  cat("PREFIXES_Table_EXPANDED_FINAL\n")
  cat(str(PREFIXES_Table_EXPANDED))
  cat("\n")
  cat("HGNC\n")
  cat(sprintf(as.character(levels(as.factor(PREFIXES_Table_EXPANDED$HGNC)))))
  cat("\n")
  
  
  cat("Gather_arranged_ordered$sample\n")
  cat(sprintf(as.character(Gather_arranged_ordered$sample)))
  cat("\n")
  
  
  
  DEF<-merge(PREFIXES_Table_EXPANDED, 
             Gather_arranged_ordered,
             by="sample",
             all.y=T)
  
  cat("DEF_FINAL\n")
  cat(str(DEF))
  cat("\n")
  cat("HGNC\n")
  cat(sprintf(as.character(levels(as.factor(DEF$HGNC)))))
  cat("\n")
  
  # quit(status =1)
  
 #### SAVE ----
  file_name=paste("Count",opt$type,"matrix.txt",sep = "_")
  
  write.table(DEF, file=file_name, sep="\t", quote=F, row.names = F)
  
  
  
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

