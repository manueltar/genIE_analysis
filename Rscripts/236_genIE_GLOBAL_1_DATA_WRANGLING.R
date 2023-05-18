

suppressMessages(library("zoo", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("biomaRt", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("Sushi", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("plyr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("data.table", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("crayon", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("withr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("ggplot2", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("optparse", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("dplyr", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("withr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("backports", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("broom", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("rstudioapi", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("tzdb", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("cli", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("tidyverse", lib.loc="/nfs/team151/software/manuel_R_libs_4_1//"))
library("svglite", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("cowplot",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("digest",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("farver",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("labeling",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("ggeasy",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
suppressMessages(library("ggrepel", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
library("reshape2", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("ggrepel", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
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

options(warn=1)

data_wrangling = function(option_list)
{
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
  
  #### Read PRIMER_TABLE & TRANSFORM ----
  
  PRIMER_TABLE<-as.data.frame(fread(opt$input_Eve, header=T), stringsAsFactors = F)
  
  colnames(PRIMER_TABLE)[which(colnames(PRIMER_TABLE) == "id")]<-"rsid"
  
  
  cat("PRIMER_TABLE_\n")
  str(PRIMER_TABLE)
  cat("\n")
  
  
  
   
  indexes_sel<-c(which(colnames(PRIMER_TABLE) == "rsid"),
                 which(colnames(PRIMER_TABLE) == "HGNC"),
                 which(colnames(PRIMER_TABLE) == "VAR"))
  
  
  
  PRIMER_TABLE_subset<-unique(PRIMER_TABLE[,indexes_sel])
  
  
  
 
  cat("PRIMER_TABLE_subset\n")
  str(PRIMER_TABLE_subset)
  cat("\n")
  
  
  #### SELECTED_DELETION_FOR_CLONING ----
  
  SELECTED_DELETION_FOR_CLONING = opt$SELECTED_DELETION_FOR_CLONING

  cat("SELECTED_DELETION_FOR_CLONING\n")
  cat(sprintf(as.character(SELECTED_DELETION_FOR_CLONING)))
  cat("\n")
  
  SELECTED_DELETION_FOR_CLONING_array<-unlist(strsplit(opt$SELECTED_DELETION_FOR_CLONING, split=","))
  
  cat("SELECTED_DELETION_FOR_CLONING_array\n")
  cat(sprintf(as.character(SELECTED_DELETION_FOR_CLONING_array)))
  cat("\n")
  
  
  
  SELECTED_TARGET<-SELECTED_DELETION_FOR_CLONING_array[1]
  
  cat("SELECTED_TARGET\n")
  cat(sprintf(as.character(SELECTED_TARGET)))
  cat("\n")
  
  
  
  SELECTED_VAR<-SELECTED_DELETION_FOR_CLONING_array[2]
  
  
  cat("SELECTED_VAR\n")
  cat(sprintf(as.character(SELECTED_VAR)))
  cat("\n")
  
  
  
  SELECTION_df<-as.data.frame(cbind(SELECTED_TARGET,SELECTED_VAR))
  
  colnames(SELECTION_df)<-c("TARGET","VAR")
  
  cat("SELECTION_df\n")
  cat(str(SELECTION_df))
  cat("\n")
  
  
  SELECTION_df$chr<-gsub("_.+$","",SELECTION_df$VAR)
  SELECTION_df$pos<-gsub("^[^_]+_","",SELECTION_df$VAR)
  SELECTION_df$pos<-gsub("_.+$","",SELECTION_df$pos)
  SELECTION_df$ref<-gsub("^[^_]+_[^_]+_","",SELECTION_df$VAR)
  SELECTION_df$ref<-gsub("_.+$","",SELECTION_df$ref)
  SELECTION_df$alt<-gsub("^[^_]+_[^_]+_[^_]+_","",SELECTION_df$VAR)
  
  
  
  
  #### LiftOver 37 -> 38 ----
  
  gr_VARS <- GRanges(
    seqnames = as.character(gsub("chr","",SELECTION_df$chr)),
    ranges=IRanges(
      start=as.numeric(SELECTION_df$pos),
      end=as.numeric(SELECTION_df$pos),
      name=SELECTION_df$VAR))
  
  # cat("gr_VARS\n")
  # str(gr_VARS)
  # cat("\n")
  
  VAR_df<-data.frame(chr=as.character(paste('chr',seqnames(gr_VARS), sep='')),
                     pos=start(gr_VARS),
                     ref=SELECTION_df$ref,
                     alt=SELECTION_df$alt,
                     VAR=SELECTION_df$VAR,
                     stringsAsFactors = F)
  
  cat("VAR_df_\n")
  str(VAR_df)
  cat("\n")
  
  #path = system.file(package="liftOver", "extdata", "Hg19Tohg38.over.chain")
  ch = import.chain("/nfs/team151/software/manuel_R_ext_data_4_1/hg19ToHg38.over.chain")
  
  seqlevelsStyle(gr_VARS) = "UCSC"  # necessary
  gr_VARS38 = liftOver(gr_VARS, ch)
  gr_VARS38 = unlist(gr_VARS38)
  genome(gr_VARS38) = "hg38"
  
  if(length(gr_VARS38) >0)
  {
    
    chr_38<-as.character(seqnames(gr_VARS38))
    names_38<-as.character(names(gr_VARS38))
    
    ref_VAR38<-gsub("^chr[^_]+_[0-9]+_","",names_38)
    ref_VAR38<-gsub("_.+$","",ref_VAR38)
    
    
    # cat("ref_VAR38\n")
    # cat(sprintf(as.character(ref_VAR38)))
    # cat("\n")
    
    alt_VAR38<-gsub("^chr[^_]+_[0-9]+_[^_]+_","",names_38)
    # alt_VAR38<-gsub("_.+$","",alt_VAR38)
    
    
    # cat("alt_VAR38\n")
    # cat(sprintf(as.character(alt_VAR38)))
    # cat("\n")
    
    
    
    
    VAR_38_df<-data.frame(chr=as.character(seqnames(gr_VARS38)),
                          pos_38=start(gr_VARS38),
                          ref=ref_VAR38,
                          alt=alt_VAR38,
                          VAR=names(gr_VARS38),
                          stringsAsFactors = F)
    
    VAR_38_df$VAR_38<-paste(VAR_38_df$chr,VAR_38_df$pos_38,VAR_38_df$ref,VAR_38_df$alt,sep='_')
    
    
    cat("VAR_38_df_1\n")
    str(VAR_38_df)
    cat("\n")
    
    
    VAR_DEF_df<-unique(merge(VAR_df,
                             VAR_38_df,
                             by=c("chr","ref","alt","VAR"),
                             all=T))
    
    VAR_DEF_df$VAR[is.na(VAR_DEF_df$VAR)]<-"ABSENT"
    
    cat("VAR_DEF_df_2\n")
    str(VAR_DEF_df)
    cat("\n")
    # 
    
    
    
    df<-unique(merge(SELECTION_df,
                     VAR_DEF_df,
                     by=c("VAR","chr","pos","ref","alt"),
                     all=T))
    
    cat("df_2\n")
    str(df)
    cat("\n")
    
  }else{
    
    
    stop("NO_LIFT_OVER\n")
    
  }# length(gr_VARS38) >0
  
  
  # quit(status = 1)
  
  #### READ and transform K562_replicates ----
  
  K562_replicates<-unlist(strsplit(opt$K562_replicates, ","))
  
  cat("K562_replicates\n")
  cat(sprintf(as.character(K562_replicates)))
  cat("\n")
  
  K562_df<-as.data.frame(cbind(rep("K562",length(K562_replicates)),K562_replicates), stringsAsFactors=F)
  colnames(K562_df)<-c("Cell_Type","file")
  
  cat("K562_df\n")
  str(K562_df)
  cat("\n")
  
  K562_df$batch<-gsub('/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/',"",K562_df$file)
  
 K562_df$batch<-gsub("/.+$","",K562_df$batch)
  
  K562_df$sample<-gsub(".+/FLASH_FILTERED_AFTER_QC_","",K562_df$file)
  
  K562_df$sample<-paste(K562_df$batch,K562_df$sample, sep="_")
  
  K562_df$sample<-gsub("Rscript_|_region_stats.tsv","",K562_df$sample)
  
  
  cat("K562_df\n")
  str(K562_df)
  cat("\n")
  
  
  #### READ and transform Kolf2_replicates ----
  
  Kolf2_replicates<-unlist(strsplit(opt$Kolf2_replicates, ","))
  
  cat("Kolf2_replicates\n")
  cat(sprintf(as.character(Kolf2_replicates)))
  cat("\n")
  
  Kolf2_df<-as.data.frame(cbind(rep("Kolf2",length(Kolf2_replicates)),Kolf2_replicates), stringsAsFactors=F)
  colnames(Kolf2_df)<-c("Cell_Type","file")
  
  Kolf2_df$batch<-gsub('/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/',"",Kolf2_df$file)
  
  Kolf2_df$batch<-gsub("/.+$","",Kolf2_df$batch)
  
  Kolf2_df$sample<-gsub(".+/FLASH_AFTER_QC_|.+/FLASH_FILTERED_AFTER_QC_","",Kolf2_df$file)
  
  Kolf2_df$sample<-paste(Kolf2_df$batch,Kolf2_df$sample, sep="_")
  
  Kolf2_df$sample<-gsub("Rscript_|_region_stats.tsv","",Kolf2_df$sample)
  Kolf2_df$sample<-gsub("_region_stats.tsv","",Kolf2_df$sample)
  
  
  
  cat("Kolf2_df\n")
  str(Kolf2_df)
  cat("\n")
  
  #### READ and transform HL60_replicates ----
  
  HL60_replicates<-unlist(strsplit(opt$HL60_replicates, ","))
  
  cat("HL60_replicates\n")
  cat(sprintf(as.character(HL60_replicates)))
  cat("\n")
  
  HL60_df<-as.data.frame(cbind(rep("HL60",length(HL60_replicates)),HL60_replicates), stringsAsFactors=F)
  colnames(HL60_df)<-c("Cell_Type","file")
  
  cat("HL60_df\n")
  str(HL60_df)
  cat("\n")
  
  HL60_df$batch<-gsub('/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/',"",HL60_df$file)
  
  HL60_df$batch<-gsub("/.+$","",HL60_df$batch)
  
  HL60_df$sample<-gsub(".+/FLASH_FILTERED_AFTER_QC_","",HL60_df$file)
  
  HL60_df$sample<-paste(HL60_df$batch,HL60_df$sample, sep="_")
  
  HL60_df$sample<-gsub("Rscript_|_region_stats.tsv","",HL60_df$sample)
  
  
  cat("HL60_df\n")
  str(HL60_df)
  cat("\n")
  
  #### READ and transform THP1_replicates ----
  
  THP1_replicates<-unlist(strsplit(opt$THP1_replicates, ","))
  
  cat("THP1_replicates\n")
  cat(sprintf(as.character(THP1_replicates)))
  cat("\n")
  
  THP1_df<-as.data.frame(cbind(rep("THP1",length(THP1_replicates)),THP1_replicates), stringsAsFactors=F)
  colnames(THP1_df)<-c("Cell_Type","file")
  
  cat("THP1_df\n")
  str(THP1_df)
  cat("\n")
  
  THP1_df$batch<-gsub('/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/genIE_analysis/',"",THP1_df$file)
  
  THP1_df$batch<-gsub("/.+$","",THP1_df$batch)
  
  THP1_df$sample<-gsub(".+/FLASH_FILTERED_AFTER_QC_","",THP1_df$file)
  
  THP1_df$sample<-paste(THP1_df$batch,THP1_df$sample, sep="_")
  
  THP1_df$sample<-gsub("Rscript_|_region_stats.tsv","",THP1_df$sample)
  
  
  cat("THP1_df\n")
  str(THP1_df)
  cat("\n")
  
  #### merge & open & mine ----
  
  
  DEF<-rbind(K562_df,Kolf2_df,HL60_df,THP1_df)
  
  cat("DEF_0\n")
  str(DEF)
  cat("\n")
  
  # quit(status=1)
  
  DEF$Cell_Type<-factor(DEF$Cell_Type,
                        levels=c("Kolf2","K562","HL60","THP1"),
                        ordered=T)
  
  DEF$sample<-factor(DEF$sample,
                        levels=c("genIE_kolf_round1_reanalysis","round2_exp0_10_10","genIE_kolf_SAV1_new_design_experiment_1","genIE_kolf_round275_post_covid_10_10","K562_10_10","HL60_10_10","THP1_10_10"),
                        ordered=T)
  
  
  DEF<-DEF[order(DEF$Cell_Type,
                 DEF$sample),]
  
  row.names(DEF)<-NULL
  
  cat("DEF_1\n")
  str(DEF)
  cat("\n")
  
  
  DEF$Rep<-paste("Rep",row.names(DEF),sep="_")
  
  
  cat("DEF_2\n")
  str(DEF)
  cat("\n")
  
  
  
  
  list_gather<-list()
  
  for(i in 1:dim(DEF)[1])
  {
    
    file_sel<-DEF$file[i]
    
    
    cat("--->\t")
    cat(sprintf(as.character(file_sel)))
    cat("\n")
    
    DEF_sel<-DEF[i,]
    
    cat("DEF_sel\n")
    str(DEF_sel)
    cat("\n")
    
    Rep_sel<-DEF$Rep[i]
    
    temp<-as.data.frame(fread(file=file_sel, sep="\t", header=T), stringsAsFactors=F)
    
    cat("temp\n")
    str(temp)
    cat("\n")
    
    # quit(status=1)
    
    indx.int<-c(which(colnames(temp) == "name"),
                which(colnames(temp) == "hdr_rate_gDNA"),which(colnames(temp) == "del_rate_gDNA"),
                which(colnames(temp) == "hdr_effect"),which(colnames(temp) == "hdr_pval"),
                which(colnames(temp) == "del_effect"),which(colnames(temp) == "del_pval"))
    
    
    temp_subset<-unique(temp[,indx.int])
    
    colnames(temp_subset)[which(colnames(temp_subset) == "name")]<-"HGNC"
    
    
    temp_subset$Percentage_hdr_gDNA<-round(temp_subset$hdr_rate_gDNA*100,1)
    temp_subset$hdr_logpval<-round(-1*log10(temp_subset$hdr_pval),2)
    temp_subset$hdr_effect<-round(temp_subset$hdr_effect,3)
    
    
    temp_subset$Percentage_del_gDNA<-round(temp_subset$del_rate_gDNA*100,1)
    temp_subset$del_logpval<-round(-1*log10(temp_subset$del_pval),2)
    temp_subset$del_effect<-round(temp_subset$del_effect,3)
    
    
      
      
    cat("temp_subset_0\n")
    str(temp_subset)
    cat("\n")
    
    indx.int<-c(which(colnames(temp_subset) == "HGNC"),
                which(colnames(temp_subset) == "Percentage_hdr_gDNA"),which(colnames(temp_subset) == "Percentage_del_gDNA"),
                which(colnames(temp_subset) == "hdr_effect"),which(colnames(temp_subset) == "hdr_logpval"),
                which(colnames(temp_subset) == "del_effect"),which(colnames(temp_subset) == "del_logpval"))
    
    
    temp_subset_2<-unique(temp_subset[,indx.int])
    
    cat("temp_subset_2_0\n")
    str(temp_subset_2)
    cat("\n")
    
    
    temp_subset_2<-merge(temp_subset_2,
                       PRIMER_TABLE_subset,
                       by="HGNC",
                       all.x=T)
    
    cat("temp_subset_2_1\n")
    str(temp_subset_2)
    cat("\n")
    
    temp_subset_2$Rep<-Rep_sel
    
    cat("temp_subset_2_2\n")
    str(temp_subset_2)
    cat("\n")
    
    list_gather[[i]]<-temp_subset_2
    
    # quit(status=1)
    
    
    
    
  }#i
  
  Gather = unique(as.data.frame(data.table::rbindlist(list_gather, fill = T)))
  
  
  cat("Gather\n")
  str(Gather)
  cat("\n")
  
  
  #### melt & pivot wider ----
  
  Gather.m<-melt(Gather, id.vars=c("VAR","rsid","HGNC","Rep"))
  
  cat("Gather.m\n")
  str(Gather.m)
  cat("\n")
  
  # quit(status=1)
  
  
  # Gather_wide<-pivot_wider(Gather.m,
  #                          id_cols=c("VAR","rsid","HGNC"),
  #                          names_from=variable,
  #                          values_from=value)
  #                                       
  # cat("Gather_wide\n")
  # str(Gather_wide)
  # cat("\n")
  
  
  
  #####  combine Rep_2 and Rep_3 into 1 unique label ----
  
  ### I trust more the Rep_2 results
  
  REST.df<-Gather.m[-c(which(Gather.m$Rep == "Rep_2"),
                       which(Gather.m$Rep == "Rep_3")),]
  
  cat("REST.df\n")
  str(REST.df)
  cat("\n")
  
  
  Rep_2.df<-Gather.m[which(Gather.m$Rep == "Rep_2"),]
  
  cat("Rep_2.df\n")
  str(Rep_2.df)
  cat("\n")
  
  Rep_3.df<-Gather.m[which(Gather.m$Rep == "Rep_3"),]
  
  cat("Rep_3.df\n")
  str(Rep_3.df)
  cat("\n")
  
  
  Rep_3.df_edited<-Rep_3.df[-which(Rep_3.df$HGNC%in%Rep_2.df$HGNC),]
  
  
  cat("Rep_3.df_edited\n")
  str(Rep_3.df_edited)
  cat("\n")
  
  Rep_3.df_DEF<-rbind(Rep_3.df_edited,
                      Rep_2.df)
  
  Rep_3.df_DEF$Rep<-"Rep_3"
  
  cat("Rep_3.df_DEF\n")
  str(Rep_3.df_DEF)
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Rep_3.df_DEF$Rep))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Rep_3.df_DEF$Rep)))))
  cat("\n")
  
  
  
  
  # Rep_4.df<-Gather.m[which(Gather.m$Rep == "Rep_4"),]
  # 
  # cat("Rep_4.df\n")
  # str(Rep_4.df)
  # cat("\n")
  # 
  # Rep_5.df<-Gather.m[which(Gather.m$Rep == "Rep_5"),]
  # 
  # cat("Rep_5.df\n")
  # str(Rep_5.df)
  # cat("\n")
  # 
  # Rep_6.df<-Gather.m[which(Gather.m$Rep == "Rep_6"),]
  # 
  # cat("Rep_6.df\n")
  # str(Rep_6.df)
  # cat("\n")
  # 
  # Rep_7.df<-Gather.m[which(Gather.m$Rep == "Rep_7"),]
  # 
  # cat("Rep_7.df\n")
  # str(Rep_7.df)
  # cat("\n")
  
  Gather.m<-rbind(REST.df,Rep_3.df_DEF) 
  
  cat("Gather.m\n")
  str(Gather.m)
  cat("\n")
  
  
  # quit(status=1)
  #### SAVE ----
  
  setwd(out)
  
  filename_1<-paste("genIE_results",".rds", sep='')
  
  cat("filename_1\n")
  cat(sprintf(as.character(filename_1)))
  cat("\n")
  
  
  saveRDS(Gather.m, file= filename_1)
  

  filename_1<-paste("genIE_results",".tsv", sep='')
  
  cat("filename_1\n")
  cat(sprintf(as.character(filename_1)))
  cat("\n")
  
  
  write.table(Gather.m, sep="\t", quote=F, row.names = F, file= filename_1)
  
  filename_1<-paste("genIE_Rosetta",".rds", sep='')
  
  cat("filename_1\n")
  cat(sprintf(as.character(filename_1)))
  cat("\n")
  
  
  saveRDS(DEF, file= filename_1)
  
  
  
  

}



get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

printList = function(l, prefix = "    ") {
  list.df = data.frame(val_name = names(l), value = as.character(l))
  list_strs = apply(list.df, MARGIN = 1, FUN = function(x) { paste(x, collapse = " = ")})
  cat(paste(paste(paste0(prefix, list_strs), collapse = "\n"), "\n"))
}

#### main script ----

# HL60_replicates

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
    make_option(c("--SELECTED_DELETION_FOR_CLONING"), type="character", default=NULL, 
                metavar="FILE.txt", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--HL60_replicates"), type="character", default=NULL, 
                metavar="NEG.CTRLS", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--THP1_replicates"), type="character", default=NULL, 
                metavar="NEG.CTRLS", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--K562_replicates"), type="character", default=NULL, 
                metavar="NEG.CTRLS", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Kolf2_replicates"), type="character", default=NULL, 
                metavar="NEG.CTRLS", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Signal_Threshold"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--fdr_Threshold"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--enhancer_empirical_log_pval_Threshold"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--enhancer_pval_Threshold"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ASE_log_pval_Threshold"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--DESIGN_DROPOUTS"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ACCEPTED_REPLICATES"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--EDITION_DROPOUTS"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--EXPRESSION_DROPOUTS"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--K562_enhancer_result"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--K562_ASE_result"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Kolf2_enhancer_result"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Kolf2_ASE_result"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--HL60_enhancer_result"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--HL60_ASE_result"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--THP1_enhancer_result"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--THP1_ASE_result"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--QC_dnaREF_ALT_threshold"), type="character", default=NULL, 
                metavar="type1", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="out", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
  )
  parser = OptionParser(usage = "137_MPRA_normalization_and_filtering_Rscript_v2.R
                        --CHARAC_TABLE FILE.txt
                        --replicas charac
                        --type1 type1
                        --type2 type2
                        --barcodes_per_tile integer
                        --pvalThreshold integer
                        --FDRThreshold integer
                        --EquivalenceTable FILE.txt
                        --sharpr2Threshold charac",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
 
data_wrangling(opt)

  
  
  
}


###########################################################################

system.time( main() )
