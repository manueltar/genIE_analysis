

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
suppressMessages(library("cowplot", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))




opt = NULL


Report_printer = function(option_list)
{
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("Manuel_4_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("Manuel_5_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### read TOTAL ----
  
  setwd(out)
  
  filename=paste(type, "_TOTAL",".rds", sep='')
  
  if (file.exists(filename)){
    
    List_TOTAL<-readRDS(file = filename)
    
  } else {
    List_TOTAL<-NULL
    
    
  }
  
  
  
  #### read HEATMAP ----
  
  filename=paste(type, "_HEATMAP",".rds", sep='')
  
  if (file.exists(filename)){
    
  List_HEATMAP<-readRDS(file = filename)
  }else{
    
    List_HEATMAP<-NULL
  }
  
  #### read Maximum_match ----
  
  filename=paste(type, "_CIGAR_Maximum_match",".rds", sep='')
  
  
  if (file.exists(filename)){
    List_Maximum_match<-readRDS(file = filename)
  }else{
    
    List_Maximum_match<-NULL
  }
  
  #### read DEL ----
  filename=paste(type, "_CIGAR_DEL",".rds", sep='')
  
  if (file.exists(filename)){
  List_DEL<-readRDS(file = filename)
  }else{
    List_DEL<-NULL
  }
  
  #### read INS ----
  
  filename=paste(type, "_CIGAR_INS",".rds", sep='')
  
  if (file.exists(filename)){
    List_INS<-readRDS(file = filename)
  }else{
    List_INS<-NULL
  }
  
  #### read SOFT ----
  
  filename=paste(type, "_CIGAR_SOFT",".rds", sep='')
  
  if (file.exists(filename)){
    List_SOFT<-readRDS(file = filename)
  }else{
    List_SOFT<-NULL
  }
  
  #### read HARD ----
  
  filename=paste(type, "_CIGAR_HARD",".rds", sep='')
  
  if (file.exists(filename)){
    List_HARD<-readRDS(file = filename)
  }else{
    List_HARD<-NULL
  }
  
 
  
  #### READ and transform COUNT_MATRIX ----
  
  setwd(out)
  
  
  filename_3<-paste(type,"_COUNT_matrix",".rds", sep='')
  
  
  COUNT_matrix = readRDS(file=filename_3)
  
  COUNT_matrix$name<-gsub("-.+","",COUNT_matrix$sample)
  
  levels_replica_order_factor<-levels(COUNT_matrix$replica_order_factor)
  
  
  cat("COUNT_matrix\n")
  cat(str(COUNT_matrix))
  cat("\n")
  
  #### LOOP GRAPH DEL ----
  
  VARS<-unique(COUNT_matrix$VAR)
  
  cat("VARS\n")
  cat(str(VARS))
  cat("\n")
  
  path_QC<-paste(out,'QC', sep='')
  
  if (file.exists(path_QC)){
    setwd(path_QC)
  } else {
    dir.create(file.path(path_QC))
    setwd(path_QC)
    
  }
  
  for(i in 1:length(VARS))
  {
    VAR_sel<-VARS[i]
    
    cat("--------------->\t")
    cat(sprintf(as.character(VAR_sel)))
    cat("\n")
    
    COUNT_matrix_sel<-COUNT_matrix[which(COUNT_matrix$VAR ==VAR_sel),]
    
    cat("COUNT_matrix_sel\n")
    cat(str(COUNT_matrix_sel))
    cat("\n")
    
    if(dim(COUNT_matrix_sel)[1]>0)
    {
      HGNC.vector<-unique(COUNT_matrix_sel$HGNC)
      
      cat("HGNC.vector\n")
      cat(str(HGNC.vector))
      cat("\n")
      
      for(k in 1:length(HGNC.vector))
      {
        HGNC_sel<-HGNC.vector[k]
        
        cat("--------------->\t")
        cat(sprintf(as.character(HGNC_sel)))
        cat("\n")
        
        
        
        COUNT_matrix_sel_HGNC_sel<-COUNT_matrix_sel[which(COUNT_matrix_sel$HGNC ==HGNC_sel),]
        
        cat("COUNT_matrix_sel_HGNC_sel\n")
        cat(str(COUNT_matrix_sel_HGNC_sel))
        cat("\n")
        
        setwd(path_QC)
        
        pdfname<-paste(type,"_Graphs_",VAR_sel,"_",HGNC_sel,".pdf", sep='')
        makepdf = TRUE
        
        if (makepdf == TRUE)
        {
          pdf ( pdfname , height=10, width=12)
        }
        
        #### plot TOTAL ----
        
        graph_TOTAL<-List_TOTAL[[VAR_sel]][[HGNC_sel]]
        
        #### plot HEATMAP ----
        
        graph_HEATMAP<-List_HEATMAP[[VAR_sel]][[HGNC_sel]]
        
        
        
        
        graph_TOTAL_and_HEATMAP<-plot_grid(graph_TOTAL, graph_HEATMAP,
                                           nrow = 2,
                                           labels = c("A","B"),
                                           label_size = 12,
                                           align = "v",
                                           rel_heights = c(1,1))
        
        print(graph_TOTAL_and_HEATMAP)
        
        #### plot Maximum_match ----
        
        if(!is.null(List_Maximum_match))
        {
          graph_Maximum_match<-List_Maximum_match[[VAR_sel]][[HGNC_sel]]
          print(graph_Maximum_match)
        
        }
        
        
        
        #### plot DEL ----
        
        if(!is.null(List_DEL))
        {
           graph_DEL<-List_DEL[[VAR_sel]][[HGNC_sel]]
           print(graph_DEL)
          
        }
        
        #### plot INS ----
        
        if(!is.null(List_INS))
        {
          graph_INS<-List_INS[[VAR_sel]][[HGNC_sel]]
          print(graph_INS)
        }
        
        #### plot SOFT ----
        
        if(!is.null(List_SOFT))
        {
          
          graph_SOFT<-List_SOFT[[VAR_sel]][[HGNC_sel]]
          print(graph_SOFT)
          
        }
        
        #### plot HARD ----
        
        if(!is.null(List_HARD))
        {
          
          graph_HARD<-List_HARD[[VAR_sel]][[HGNC_sel]]
          print(graph_HARD)
          
          
        }
        
        
        
        #### pdf printing ----
        
       
        
        
        
        
        if (makepdf == TRUE)
        {
          dev.off()
        }
        
      }
    }
  }# i VARS
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
    make_option(c("--COUNT_MATRIX"), type="character", default=NULL,
                metavar="FILE.tsv",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--mismatch_CONDENSED"), type="character", default=NULL, 
                metavar="FILE.tsv", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--TOTAL_READS"), type="character", default=NULL, 
                metavar="FILE.tsv", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--DEMULTIPLEX_RESULT"), type="character", default=NULL, 
                metavar="FILE.tsv", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="filename", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
    
  )
  
  parser = OptionParser(usage = "116_QC_proceser.R
                        --COUNT_MATRIX FILE.tsv
                        --mismatch_CONDENSED FILE.tsv 
                        --TOTAL_READS FILE.tsv 
                        --DEMULTIPLEX_RESULT FILE.tsv
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)


  
  Report_printer(opt)
 
  
}


###########################################################################

system.time( main() )
