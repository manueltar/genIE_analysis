





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
# suppressMessages(library("svglite", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
# suppressMessages(library("ggeasy", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
# suppressMessages(library("sandwich", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
# library("digest", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
# 
# library("ggforce", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
suppressMessages(library("stringr", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("seqRFLP", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))


suppressMessages(library("BiocGenerics", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("S4Vectors", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("IRanges", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("GenomeInfoDb", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("GenomicRanges", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("Biobase", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("AnnotationDbi", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("GenomicFeatures", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("OrganismDbi", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("GO.db", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("org.Hs.eg.db", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("Homo.sapiens", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("gwascat", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("rtracklayer", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("liftOver",lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))



opt = NULL

input_for_266 = function(option_list)
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
  
  cat("span_\n")
  cat(sprintf(as.character(span)))
  cat("\n")
  
    
  #### Read PRIMER_TABLE & TRANSFORM ----
  
  PRIMER_TABLE<-read.csv(opt$input_Eve, header=T, stringsAsFactors = F)
  
  cat("PRIMER_TABLE_\n")
  str(PRIMER_TABLE)
  cat("\n")
  
  # quit(status = 1)
  
  PRIMER_TABLE_NO_NA<-PRIMER_TABLE[(PRIMER_TABLE$VAR != ""),]
  
  cat("PRIMER_TABLE_NO_NA_0\n")
  str(PRIMER_TABLE_NO_NA)
  cat("\n")
  
  
  ##### NO "" HGNC -----
  
  PRIMER_TABLE_NO_NA<-PRIMER_TABLE_NO_NA[(PRIMER_TABLE_NO_NA$HGNC != ""),]
  
  cat("PRIMER_TABLE_NO_NA_1\n")
  str(PRIMER_TABLE_NO_NA)
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(PRIMER_TABLE_NO_NA$HGNC))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(PRIMER_TABLE_NO_NA$HGNC)))))
  cat("\n")
  
  # quit(status=1)

  indexes_sel<-c(which(colnames(PRIMER_TABLE_NO_NA) == "RSid"),
                 which(colnames(PRIMER_TABLE_NO_NA) == "HGNC"),
                 which(colnames(PRIMER_TABLE_NO_NA) == "VAR"),
                 which(colnames(PRIMER_TABLE_NO_NA) == "chr"),
                 which(colnames(PRIMER_TABLE_NO_NA) == "pos"),
                 which(colnames(PRIMER_TABLE_NO_NA) == "ref"),
                 which(colnames(PRIMER_TABLE_NO_NA) == "alt"))
  
  
  
  PRIMER_TABLE_NO_NA_subset<-unique(PRIMER_TABLE_NO_NA[,indexes_sel])
  colnames(PRIMER_TABLE_NO_NA_subset)[which(colnames(PRIMER_TABLE_NO_NA_subset) == "RSid")]<-"rsid"
  
  cat("PRIMER_TABLE_NO_NA_subset_1\n")
  str(PRIMER_TABLE_NO_NA_subset)
  cat("\n")
  
 
    #### file to retrieve the snps in Kolf2 ----
  
  
  indexes.int2<-c(which(colnames(PRIMER_TABLE_NO_NA_subset) == "chr"),
                  which(colnames(PRIMER_TABLE_NO_NA_subset) == "pos"))
  
  file_for_snps<-unique(PRIMER_TABLE_NO_NA_subset[,indexes.int2])
  
  file_for_snps_ordered<-file_for_snps[order(as.numeric(file_for_snps$chr), 
                                                               as.numeric(file_for_snps$pos)),]
  
  file_for_snps_ordered$START<-file_for_snps_ordered$pos-span
  
  file_for_snps_ordered$STOP<-file_for_snps_ordered$pos+span
  
  cat("file_for_snps_ordered_2\n")
  str(file_for_snps_ordered)
  cat("\n")
  
  #quit(status=1)
  
  file_name2<-"genERA_VAR_POST_HF1_2_3_region_file_FROM_TO.txt"
  
  indexes.int3<-c(which(colnames(file_for_snps_ordered) == "chr"),
                  which(colnames(file_for_snps_ordered) == "START"),
                  which(colnames(file_for_snps_ordered) == "STOP"))
  
  file_for_snps_DEF<-unique(file_for_snps_ordered[,indexes.int3])
  
  cat("file_for_snps_DEF_2\n")
  str(file_for_snps_DEF)
  cat("\n")
  
  #cat("CHROM  POS",
  #   file=file_name,sep="\n")
  
  
  ### export bed of tiles ----
  
  gr_PRIMER_TABLE_NO_NA_subset <- GRanges(
    seqnames = paste("Chr",PRIMER_TABLE_NO_NA_subset$chr,sep=''),
    ranges=IRanges(
      start=PRIMER_TABLE_NO_NA_subset$pos - span, # note the difference
      end=PRIMER_TABLE_NO_NA_subset$pos + span,
      name=paste(PRIMER_TABLE_NO_NA_subset$VAR,PRIMER_TABLE_NO_NA_subset$HGNC, sep="__")))
  
  gr_PRIMER_TABLE_NO_NA_subset_unique<-unique(gr_PRIMER_TABLE_NO_NA_subset)
  
  setwd(out)
  
  filename_20<-paste(type,"_TILES",".bed", sep='')
  
  
  export.bed(gr_PRIMER_TABLE_NO_NA_subset_unique,con=filename_20)
  
  
  
  #### SAVING ----
  
  setwd(out)
  
  write.table(file_for_snps_DEF, 
              file=file_name2, append = FALSE, sep="\t", 
              quote=F,col.names = F, row.names = F, eol="\n")
  
  write.table(PRIMER_TABLE_NO_NA_subset, 
              file="input_corrected_unblanked.tsv", sep="\t", quote=F, row.names = F)
  
#  quit(status=1)
  
 
  
}


printList = function(l, prefix = "    ") {
  list.df = data.frame(val_name = names(l), value = as.character(l))
  list_strs = apply(list.df, MARGIN = 1, FUN = function(x) { paste(x, collapse = " = ")})
  cat(paste(paste(paste0(prefix, list_strs), collapse = "\n"), "\n"))
}

getMinOverlap = function(amplicon_size, read_len) {
  minOverlap = 10
  if (2*read_len - amplicon_size < 10) {
    minOverlap = 2*read_len - amplicon_size - 1
    if (minOverlap < 4) {
      minOverlap = 4
    }
  }
  minOverlap
}


#### main script ----

main = function() {
  cmd_line = commandArgs()
  cat("Command line:\n")
  cat(paste(gsub("--file=", "", cmd_line[4], fixed=T),
            paste(cmd_line[6:length(cmd_line)], collapse = " "),
            "\n\n"))
  option_list <- list(
    make_option(c("--input_Eve"), type="character", default=NULL, 
                metavar="FILE.txt", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--span"), type="numeric", default=NULL, 
                metavar="FILE.txt", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="filename", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
    
  )
  parser = OptionParser(usage = "136_Rscript_edition_Add_Expression_List.R
                        --input_Eve FILE.txt
                        --equivalence_HL60 FILE.txt
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  input_for_266(opt)
  
}


###########################################################################

system.time( main() )
