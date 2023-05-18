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

#### functions ----

opt = NULL

options(warn=1)

create_merged_table = function(option_list)
{
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  #### BIAS!! Add OLD PREFIX FOR Eve ----
  
  PREFIXES_Table_EXPANDED = read.table(opt$prefixes_Eve, sep="\t", stringsAsFactors = F, header = T)
  
  
  PREFIXES_Table_EXPANDED<-PREFIXES_Table_EXPANDED[!is.na(PREFIXES_Table_EXPANDED$Sample_Name),]
  
  # PREFIXES_Table_EXPANDED$Sample_Name<-gsub("^Genie2_","",PREFIXES_Table_EXPANDED$Sample_Name)
  
  colnames(PREFIXES_Table_EXPANDED)[which(colnames(PREFIXES_Table_EXPANDED) == "hgnc")]<-"HGNC"
  
  PREFIXES_Table_EXPANDED$Sample_Name<-gsub("-","_",PREFIXES_Table_EXPANDED$Sample_Name)
  
#  PREFIXES_Table_EXPANDED$hgnc[PREFIXES_Table_EXPANDED$hgnc == "EROS"]<-tolower(PREFIXES_Table_EXPANDED$hgnc[PREFIXES_Table_EXPANDED$HGNC == "EROS"])
  
  PREFIXES_Table_EXPANDED<-PREFIXES_Table_EXPANDED[!is.na(PREFIXES_Table_EXPANDED$Replicate),]

  cat("PREFIXES_Table_EXPANDED\n")
  cat(str(PREFIXES_Table_EXPANDED))
  cat("\n")
  cat(sprintf(as.character(PREFIXES_Table_EXPANDED$HGNC)))
  cat("\n")
  cat(sprintf(as.character(PREFIXES_Table_EXPANDED$index)))
  cat("\n")
  cat(sprintf(as.character(PREFIXES_Table_EXPANDED$index2)))
  cat("\n")
  
  
  
  PREFIXES_RUN = read.table(opt$prefixes_EXPERIMENT, sep="\t", stringsAsFactors = F, header = T)
    
  cat("PREFIXES_RUN_\n")
  cat(str(PREFIXES_RUN))
  cat("\n")
  # cat(sprintf(as.character(PREFIXES_RUN$file)))
  # cat("\n")
  # cat(sprintf(as.character(PREFIXES_RUN$indexes_string)))
  # cat("\n")
  # 
  # quit(status = 1)
  
  #### route_fastq ----
  
  route_fastq = opt$route_fastq
  
  cat("Manuel_3_\n")
  cat(sprintf(as.character(route_fastq)))
  cat("\n")
  
  #### Get sample name from file ----
  
  
  PREFIXES_RUN$Sample_Name<-gsub(route_fastq,"",PREFIXES_RUN$file)
  
  
  PREFIXES_RUN$Sample_Name<-sub("_","\\|",PREFIXES_RUN$Sample_Name)
  PREFIXES_RUN$Sample_Name<-gsub("_.+"," ",PREFIXES_RUN$Sample_Name)
  PREFIXES_RUN$Sample_Name<-sub("\\|","_",PREFIXES_RUN$Sample_Name)
  PREFIXES_RUN$Sample_Name<-sub(" ","",PREFIXES_RUN$Sample_Name)
  PREFIXES_RUN$Sample_Name<-gsub("Genie2-","",PREFIXES_RUN$Sample_Name)
  PREFIXES_RUN$Sample_Name<-gsub("_S[0-9]+","",PREFIXES_RUN$Sample_Name, perl=T)
  
  
  cat("PREFIXES_RUN$Sample_Name\n")
  cat(sprintf(as.character(PREFIXES_RUN$Sample_Name)))
  cat("\n")

  #quit(status=1)
  #### Get the MAX index per sample from top 5 ----
  
  PREFIXES_RUN_MAX<-setDT(PREFIXES_RUN)[, .SD[which.max(reads)], by=c("Sample_Name")]
  PREFIXES_RUN_MAX.df<-as.data.frame(PREFIXES_RUN_MAX)

  cat("PREFIXES_RUN_MAX_PRE\n")
  cat(str(PREFIXES_RUN_MAX.df))
  cat("\n")
  
  
  #### drop file column ----
  
  indx_del<-which(colnames(PREFIXES_RUN_MAX.df) == "file")
  PREFIXES_RUN_MAX.df<-PREFIXES_RUN_MAX.df[,-indx_del]
  
  PREFIXES_RUN_MAX.df$Sample_Name<-gsub("-","_",PREFIXES_RUN_MAX.df$Sample_Name)
  
  cat("PREFIXES_RUN_MAX_POST\n")
  cat(str(PREFIXES_RUN_MAX.df))
  cat("\n")
  
  #quit(status = 1)
  
  #### Break the index string of EXP ----
  
  PREFIXES_RUN_MAX.df<-cSplit(PREFIXES_RUN_MAX.df,
                                      "indexes_string", sep = "+", fixed =T, drop=F)
  
  colnames(PREFIXES_RUN_MAX.df)[which(colnames(PREFIXES_RUN_MAX.df) == "indexes_string_1")]<-"index2"
  colnames(PREFIXES_RUN_MAX.df)[which(colnames(PREFIXES_RUN_MAX.df) == "indexes_string_2")]<-"index"
  
  #### Do the reverse complement of index2 ----
  
  REV_COMP<-NULL
  
  for(i in 1:length(PREFIXES_RUN_MAX.df$index2))
  {
    REV_COMP[i]<-as.character(reverseComplement(DNAString(as.character(PREFIXES_RUN_MAX.df$index2)[i])))
  }
  
  PREFIXES_RUN_MAX.df$index2_rev_comp<-REV_COMP
  
  #### Do the reverse complement of index ----
  
  REV_COMP<-NULL
  
  for(i in 1:length(PREFIXES_RUN_MAX.df$index))
  {
    REV_COMP[i]<-as.character(reverseComplement(DNAString(as.character(PREFIXES_RUN_MAX.df$index)[i])))
  }
  
  PREFIXES_RUN_MAX.df$index_rev_comp<-REV_COMP
  
  cat("PREFIXES_RUN_MAX_fields\n")
  cat(sprintf(as.character(PREFIXES_RUN_MAX$Sample_Name)))
  cat("\n")
  cat(sprintf(as.character(PREFIXES_RUN_MAX.df$reads)))
  cat("\n")
  cat(sprintf(as.character(PREFIXES_RUN_MAX.df$index)))
  cat("\n")
  cat(sprintf(as.character(PREFIXES_RUN_MAX.df$index_rev_comp)))
  cat("\n")
  cat(sprintf(as.character(PREFIXES_RUN_MAX.df$index2)))
  cat("\n")
  cat(sprintf(as.character(PREFIXES_RUN_MAX.df$index2_rev_comp)))
  cat("\n")
  
  cat("PREFIXES_RUN_MAX.df_POST_POST\n")
  cat(str(PREFIXES_RUN_MAX.df))
  cat("\n")
  
 # quit(status = 1)
  
  #### Merge and print the file----
  
  cat("PREFIXES_Table_EXPANDED\n")
  cat(str(PREFIXES_Table_EXPANDED))
  cat("\n")
  
  cat("PREFIXES_RUN_MAX.df\n")
  cat(str(PREFIXES_RUN_MAX.df))
  cat("\n")
  
  file_name<-opt$out
  
  indx.check<-which(PREFIXES_Table_EXPANDED$Sample_Name%in%PREFIXES_RUN_MAX.df$Sample_Name)
  
  cat("indx.check\n")
  cat(str(indx.check))
  cat("\n")
  
  if(length(indx.check) == 0)
  {
    PREFIXES_RUN_MAX.df$Sample_Name<-paste('1_',PREFIXES_RUN_MAX.df$Sample_Name,sep='')
    
    indx.check2<-which(PREFIXES_Table_EXPANDED$Sample_Name%in%PREFIXES_RUN_MAX.df$Sample_Name)
    
    cat("indx.check2\n")
    cat(str(indx.check2))
    cat("\n")
    
    if(length(indx.check2) > 0)
    {
      check<-PREFIXES_Table_EXPANDED$Sample_Name[indx.check2]
      
      pivotal_point_0<-sum(length(check) == length(PREFIXES_Table_EXPANDED$Sample_Name))
      
      cat("pivotal_point_0\n")
      cat(str(pivotal_point_0))
      cat("\n")
      
      if(pivotal_point_0 == 0)
      {
        cat("Partial recognition\n")
        quit(status=1)
        
      }
      
    }else{
      
      cat("Error in barcodes\n")
      quit(status=1)
    }
    
    
  }
  
  # quit(status = 1)
  
  
  
  
  
  
  
  DEF<-merge(PREFIXES_Table_EXPANDED,
             PREFIXES_RUN_MAX.df,
             by="Sample_Name",
             all=T)
  
  cat("DEF\n")
  cat(str(DEF))
  cat("\n")

  # quit(status = 1)
  
  
  colnames(DEF)[which(colnames(DEF) == "index.P7")]<-"BC3.Eve"
  colnames(DEF)[which(colnames(DEF) == "index")]<-"BC5.fastq"
  
  colnames(DEF)[which(colnames(DEF) == "index.P5")]<-"BC5.Eve"
  colnames(DEF)[which(colnames(DEF) == "index2")]<-"BC3.fastq"
  colnames(DEF)[which(colnames(DEF) == "index2_rev_comp")]<-"BC3.fastq_rev_comp"
  colnames(DEF)[which(colnames(DEF) == "index_rev_comp")]<-"BC5.fastq_rev_comp"
  
  cat("DEF\n")
  cat(str(DEF))
  cat("\n")
  
  # quit(status = 1)
  
  
  indexes_reorder<-c(which(colnames(DEF)== "Sample_Name"),
                     which(colnames(DEF)== "sample"),
                     which(colnames(DEF)== "HGNC"),
                     which(colnames(DEF)== "rsid"),
                     which(colnames(DEF)== "Sample_ID"),
                     which(colnames(DEF)== "BC5.Eve"),
                     which(colnames(DEF)== "BC3.Eve"),
                     which(colnames(DEF)== "indexes_string"),
                     which(colnames(DEF)== "BC5.fastq"),
                     which(colnames(DEF)== "BC5.fastq_rev_comp"),
                     which(colnames(DEF)== "BC3.fastq"),
                     which(colnames(DEF)== "BC3.fastq_rev_comp"),
                     which(colnames(DEF)== "reads"))
  
  write.table(DEF[,indexes_reorder],file=file_name,
              sep="\t",
              quote=F,
              row.names = F)
  
}  

analysis_merged_table = function(option_list)
{
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  #### BIAS!! Add OLD PREFIX FOR Eve ----
  
  merged_table = read.table(opt$out, sep="\t", stringsAsFactors = F, header = T)
  
  # str(merged_table)
  # 
  cat("merged_table\n")
  cat(str(merged_table))
  cat("\n")
  
  #### CLASIF LOOP  initialize variables ----
  
  CLASS<-NULL
  
 
  
  for(i in 1:dim(merged_table)[1])
  {
    merged_table_sel<-merged_table[i,]
    
    cat("merged_table_sel\n")
    cat(str(merged_table_sel))
    cat("\n")
    
    Sample_Name_sel<-merged_table_sel$Sample_Name
    BC5_Eve<-merged_table_sel$BC5.Eve
    BC3_Eve<-merged_table_sel$BC3.Eve
    BC5_fastq<-merged_table_sel$BC5.fastq
    BC3_fastq<-merged_table_sel$BC3.fastq
    BC5_fastq_rev_comp<-merged_table_sel$BC5.fastq_rev_comp
    BC3_fastq_rev_comp<-merged_table_sel$BC3.fastq_rev_comp
    
    FLAG<-sum(is.na(BC5_Eve))
    
    #FLAG<-sum(Sample_Name_sel == "NA")
    
    cat("Sample_Name_sel\n")
    cat(sprintf(as.character(Sample_Name_sel)))
    cat("\n")
    cat(sprintf(as.character(BC5_Eve)))
    cat("\n")
    cat(sprintf(as.character(BC3_Eve)))
    cat("\n")
    cat(sprintf(as.character(BC5_fastq)))
    cat("\n")
    cat(sprintf(as.character(BC3_fastq)))
    cat("\n")
    cat(sprintf(as.character(BC5_fastq_rev_comp)))
    cat("\n")
    cat(sprintf(as.character(BC3_fastq_rev_comp)))
    cat("\n")
    cat(sprintf(as.character(BC5_Eve == BC5_fastq)))
    cat("\n")
    cat(sprintf(as.character(BC3_Eve == BC3_fastq_rev_comp)))
    cat("\n")
    # 
    # 
    cat(sprintf(as.character("----------------------->FLAG:\t")))
    cat(sprintf(as.character(FLAG)))
    cat("\n")
   
    if(FLAG == 0) # Do not continue in cases where everything is NA
    {
      #### algorithm of classification ----
      
      pivotal_point<-sum(BC5_Eve == BC5_fastq)
      
      cat("pivotal_point\n")
      cat(sprintf(as.character(pivotal_point)))
      cat("\n")
      
      if(!is.na(pivotal_point))
      {
        
          if(pivotal_point == 1)#(BC5_Eve == BC5_fastq)
          {
            pivotal_point2<-sum(BC3_Eve == BC3_fastq)
            
            cat("pivotal_point2\n")
            cat(sprintf(as.character(pivotal_point2)))
            cat("\n")
            
            if(pivotal_point2 ==1)
            {
              CLASS[i]<-"concordant"
              
            }else{
              
              pivotal_point3<-sum(BC3_Eve == BC3_fastq_rev_comp)
              
              cat("pivotal_point3\n")
              cat(sprintf(as.character(pivotal_point3)))
              cat("\n")
              
              if(pivotal_point3 == 1)
              {
                CLASS[i]<-"concordant_rev_comp_BC3"
                
              }else{
                
                CLASS[i]<-"discordant_BC3"
              }
            }
            
          }else{
            
            pivotal_point4<-sum(BC3_Eve == BC3_fastq)
            
            cat("pivotal_point4\n")
            cat(sprintf(as.character(pivotal_point4)))
            cat("\n")
            
            if(pivotal_point4 == 1)
            {
              pivotal_point5<-sum(BC5_Eve == BC5_fastq_rev_comp)
              
              cat("pivotal_point5\n")
              cat(sprintf(as.character(pivotal_point5)))
              cat("\n")
              
              if(pivotal_point5 == 1)
              {
                CLASS[i]<-"concordant_rev_comp_BC5"
                
              }else{
                
                CLASS[i]<-"discordant_BC5"
              }
              
            }else{
              
              cat("HELLO_WORLD\n")
              
              pivotal_point_6<-sum(BC3_Eve == BC5_fastq)
              
              cat("pivotal_point_6\n")
              cat(sprintf(as.character(pivotal_point_6)))
              cat("\n")
              
              # quit(status = 1)
              
              if(pivotal_point_6 == 1)
              {
                pivotal_point_7<-sum(BC5_Eve == BC3_fastq)
                
                cat("pivotal_point_7\n")
                cat(sprintf(as.character(pivotal_point_7)))
                cat("\n")
                
                if(pivotal_point_7 == 1)
                {
                  CLASS[i]<-"concordant_Eve_swapped_manifest"
                  
                }else{
                  
                  pivotal_point_8<-sum(BC5_Eve == BC3_fastq_rev_comp)
                  
                  cat("pivotal_point_8\n")
                  cat(sprintf(as.character(pivotal_point_8)))
                  cat("\n")
                  
                  if(pivotal_point_8 == 1)
                  {
                    CLASS[i]<-"concordant_5_concordant_rev_comp_3_Eve_swapped_manifest"
                    
                  }else{
                    
                    CLASS[i]<-"concordant_BC5_discordant_BC3_Eve_swapped_manifest"
                  }
                }
                
              }else{
                
                CLASS[i]<-"discordant_TOTAL"
              }
            }
          }
        }else{
          
          CLASS[i]<-"NA"
          
        }
    }else{
      
      cat("Error in barcodes\n")
      quit(status=1)
      
    }#!is.na(pivotal_point)
      
      
      
    
  }
  
  merged_table$CLASS<-CLASS
  
  #### print the file ----
  
  file_name<-opt$out
  
  write.table(merged_table,file=file_name,
              sep="\t",
              quote=F,
              row.names = F)
    
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
    make_option(c("--prefixes_EXPERIMENT"), type="character", default=NULL, 
                metavar="FILE.txt", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
        make_option(c("--prefixes_Eve"), type="character", default=NULL, 
                metavar="FILE.txt", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
        make_option(c("--route_fastq"), type="character", default=NULL, 
                metavar="route", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="filename", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
        
      )
  
  parser = OptionParser(usage = "114_my_own_functions_genIE_index_demultiplexation.R 
                        --prefixes_EXPERIMENT file.tsv 
                        --prefixes_Eve file.txt 
                        --route_fastq route 
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  create_merged_table(opt)
  analysis_merged_table(opt)
  
}


###########################################################################

system.time( main() )

