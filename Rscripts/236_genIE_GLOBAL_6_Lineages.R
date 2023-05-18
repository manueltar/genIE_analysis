#### WARNING HARDCODED INDEXES -----


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
library("reshape2",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")


opt = NULL

options(warn = 1)

graph_function = function(option_list)
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
  

  #### READ and transform finemap_prob_Threshold ----
  
  finemap_prob_Threshold = opt$finemap_prob_Threshold
  
  cat("finemap_prob_Threshold_\n")
  cat(sprintf(as.character(finemap_prob_Threshold)))
  cat("\n")
  
  #### Read TOME_correspondence ----
  
  TOME_correspondence = read.table(opt$TOME_correspondence, sep="\t", stringsAsFactors = F, header = T)
  
  cat("TOME_correspondence_\n")
  str(TOME_correspondence)
  cat("\n")
  
  indx.int<-c(which(colnames(TOME_correspondence) == "phenotype"),
              which(colnames(TOME_correspondence) == "Lineage"))
  
  
  TOME_correspondence_subset<-unique(TOME_correspondence[,indx.int])
  
  cat("TOME_correspondence_subset_\n")
  str(TOME_correspondence_subset)
  cat("\n")
  
  #### Read genIE results ----
  
  
 
  genIE_RESULTS<-readRDS(opt$genIE_RESULTS)
  
  cat("genIE_RESULTS\n")
  str(genIE_RESULTS)
  cat("\n")
  
  #### Cell_Type colors ----
  
  Cell_Type_levels<-levels(genIE_RESULTS$Cell_Type)
  colors_Cell_Type_levels<-c('firebrick2','#32A852','#553B68','#D45E85','#6DB2EE','#62D07F','#C9244B','#87447B','#D6B8E6')
  
  df.color_Cell_Type<-as.data.frame(cbind(Cell_Type_levels,colors_Cell_Type_levels[1:length(Cell_Type_levels)]), stringsAsFactors=F)
  colnames(df.color_Cell_Type)<-c("Cell_Type","colors")
  
  cat("df.color_Cell_Type_0\n")
  cat(str(df.color_Cell_Type))
  cat("\n")
  
 
  # cat(sprintf(as.character(levels(genIE_RESULTS$variable))))
  # cat("\n")
  # cat(sprintf(as.character(levels(as.factor(genIE_RESULTS$Rep)))))
  # cat("\n")
  
  
  
  #### Read genIE_Tier_1 file ----
  
  genIE_Tier_1 = readRDS(opt$genIE_Tier_1)
  
  
  cat("genIE_Tier_1_\n")
  cat(str(genIE_Tier_1))
  cat("\n")
  
  Tier_1<-unlist(strsplit(as.character(opt$genIE_Tier_1), split='/'))
  
  # cat("Tier_1_0\n")
  # cat(sprintf(as.character(Tier_1)))
  # cat("\n")
  
  Tier_1<-gsub(".rds$","",Tier_1[length(Tier_1)])
  Tier_1<-gsub("^genIE_Tier_","",Tier_1[length(Tier_1)])
  
  
  cat("Tier_1_1\n")
  cat(sprintf(as.character(Tier_1)))
  cat("\n")
  
  # quit(status = 1)
  
  ##### genIE results ACCEPTED -----
  
 ACCEPTED_LABELS<-c("NOT_SIGNIFICANT","SIGNIFICANT")
  
  if(Tier_1 == "del")
  {
    genIE_RESULTS_ACCEPTED<-genIE_RESULTS[which(genIE_RESULTS$STATUS_del%in%ACCEPTED_LABELS),]
    
    cat("genIE_RESULTS_ACCEPTED\n")
    str(genIE_RESULTS_ACCEPTED)
    cat("\n")
    
    genIE_Tier_1_wide<-as.data.frame(pivot_wider(genIE_Tier_1,
                                                        id_cols=c("HGNC","VAR","rsid"),
                                                        names_from=Cell_Type,
                                                        values_from=STATUS_del))
    
    genIE_RESULTS_ACCEPTED_wide<-as.data.frame(pivot_wider(genIE_RESULTS_ACCEPTED,
                                                           id_cols=c("HGNC","VAR","rsid"),
                                                           names_from=Cell_Type,
                                                           values_from=STATUS_del))
    
    
    
    
    # quit(status = 1)
    
    
    
  }
  if(Tier_1 == "hdr")
  {
    genIE_RESULTS_ACCEPTED<-genIE_RESULTS[which(genIE_RESULTS$STATUS_hdr%in%ACCEPTED_LABELS),]
    
    cat("genIE_RESULTS_ACCEPTED\n")
    str(genIE_RESULTS_ACCEPTED)
    cat("\n")
    
    genIE_Tier_1_wide<-as.data.frame(pivot_wider(genIE_Tier_1,
                                                 id_cols=c("HGNC","VAR","rsid"),
                                                 names_from=Cell_Type,
                                                 values_from=STATUS_hdr))
    
    genIE_RESULTS_ACCEPTED_wide<-as.data.frame(pivot_wider(genIE_RESULTS_ACCEPTED,
                                                           id_cols=c("HGNC","VAR","rsid"),
                                                           names_from=Cell_Type,
                                                           values_from=STATUS_hdr))
    
  }
  if(Tier_1 == "del_and_hdr")
  {
    genIE_RESULTS_ACCEPTED<-genIE_RESULTS[which(genIE_RESULTS$STATUS_del%in%ACCEPTED_LABELS &
                                                genIE_RESULTS$STATUS_hdr%in%ACCEPTED_LABELS),]
    
    cat("genIE_RESULTS_ACCEPTED\n")
    str(genIE_RESULTS_ACCEPTED)
    cat("\n")
    
    genIE_Tier_1_wide<-as.data.frame(pivot_wider(genIE_Tier_1,
                                                 id_cols=c("HGNC","VAR","rsid"),
                                                 names_from=Cell_Type,
                                                 values_from=STATUS_del))
    
    genIE_RESULTS_ACCEPTED_wide<-as.data.frame(pivot_wider(genIE_RESULTS_ACCEPTED,
                                                           id_cols=c("HGNC","VAR","rsid"),
                                                           names_from=Cell_Type,
                                                           values_from=STATUS_del))
    
  }
  
  
  
  cat("genIE_Tier_1_wide\n")
  str(genIE_Tier_1_wide)
  cat("\n")
  
  cat("genIE_RESULTS_ACCEPTED_wide\n")
  str(genIE_RESULTS_ACCEPTED_wide)
  cat("\n")
  
  
  ### Read dB file ----
  
  dB = as.data.frame(fread(file=opt$dB, sep="\t", stringsAsFactors = F, header = T), stringsAsFactors =F)
  
  cat("dB_\n")
  cat(str(dB))
  cat("\n")
  
  
  dB_subset<-dB[which(dB$VAR%in%genIE_RESULTS_ACCEPTED$VAR),]
  
  cat("dB_subset_\n")
  cat(str(dB_subset))
  cat("\n")
  
  dB_subset_thresholded<-dB_subset[which(dB_subset$finemap_prob >= finemap_prob_Threshold),]
  
  cat("dB_subset_thresholded_\n")
  cat(str(dB_subset_thresholded))
  cat("\n")
  
  
 
  
  ####  MERGES  ####
  
 
  
  traits_merge<-merge(dB_subset_thresholded,
                    TOME_correspondence_subset,
                    by="phenotype",
                    all.x=T)
  
  
  
  cat("traits_merge_\n")
  cat(str(traits_merge))
  cat("\n")
  
  indx.int<-c(which(colnames(traits_merge) == "VAR"),
              which(colnames(traits_merge) == "Lineage"))
  
  PRE.Lineage<-unique(traits_merge[,indx.int])
  
  PRE.Lineage<-as.data.frame(setDT(PRE.Lineage)[,.(paste(sort(Lineage), collapse = "|")),
                                                by = VAR])
  
  colnames(PRE.Lineage)[which(colnames(PRE.Lineage) =="V1")]<-"Lineage_CLASS"
  
  
  cat("PRE.Lineage_\n")
  cat(str(PRE.Lineage))
  cat("\n")
  
  
 
  traits_merge<-merge(traits_merge,
                    PRE.Lineage,
                    by="VAR",
                    all.x=T)
  
  cat("traits_merge_\n")
  cat(str(traits_merge))
  cat("\n")
  
  #quit(status = 1)
  
  #### DEF CLASSIFICATION OF LINEAGES ----
  
  traits_merge$LINEAGE_CLASSIF_DEF<-"Pleiotropic_variants"
  
  traits_merge$LINEAGE_CLASSIF_DEF[traits_merge$Lineage_CLASS == "erythroid_lineage"]<-"Erythroid Traits"
  traits_merge$LINEAGE_CLASSIF_DEF[traits_merge$Lineage_CLASS == "mega_lineage"]<-"Megakaryocitic Traits"
  traits_merge$LINEAGE_CLASSIF_DEF[traits_merge$Lineage_CLASS == "gran_mono_lineage"]<-"GRAN-MONO Traits"
  traits_merge$LINEAGE_CLASSIF_DEF[traits_merge$Lineage_CLASS == "lymph_lineage"]<-"Lymphocitic Traits"
  
  
  
 
  # quit(status=1)
  #### list_RESULTS ----
  
  list_RESULTS<-list()
  
 
  ##### LINEAGE specificty 0 vs 1 per Cell Type ----
  
  indx.int<-c(which(colnames(traits_merge) == "VAR"),which(colnames(traits_merge) == "LINEAGE_CLASSIF_DEF"),which(colnames(traits_merge) == "Lineage_CLASS"))
  
  traits_merge_subset<-unique(traits_merge[,indx.int])
  
  cat("traits_merge_subset_1\n")
  cat(str(traits_merge_subset))
  cat("\n")
  
  
  # quit(status = 1)
  
  traits_merge_subset<-unique(traits_merge_subset[which(traits_merge_subset$LINEAGE_CLASSIF_DEF != ""),])
  
  traits_merge_subset$LINEAGE_CLASSIF_DEF<-factor(traits_merge_subset$LINEAGE_CLASSIF_DEF,
                                                         levels=c("Erythroid Traits","Megakaryocitic Traits","GRAN-MONO Traits","Lymphocitic Traits","Pleiotropic_variants"),
                                                         ordered=T)
  traits_merge_subset$Kolf2<-"NA"
  traits_merge_subset$K562<-"NA"
  traits_merge_subset$HL60<-"NA"
  traits_merge_subset$THP1<-"NA"
  
  cat("traits_merge_subset_1\n")
  cat(str(traits_merge_subset))
  cat("\n")
  cat(sprintf(as.character(names(summary(traits_merge_subset$LINEAGE_CLASSIF_DEF)))))
  cat("\n")
  cat(sprintf(as.character(summary(traits_merge_subset$LINEAGE_CLASSIF_DEF))))
  cat("\n")
  
  REP.LIN<-as.data.frame(melt(setDT(traits_merge_subset), id.vars=c("VAR","LINEAGE_CLASSIF_DEF","Lineage_CLASS")), stringsAsFactors=F)
  
  colnames(REP.LIN)[which(colnames(REP.LIN) == "variable")]<-"Cell_Type"
  colnames(REP.LIN)[which(colnames(REP.LIN) == "value")]<-"Status"
  
  
  cat("genIE_Tier_1_wide\n")
  str(genIE_Tier_1_wide)
  cat("\n")
  
  cat("REP.LIN\n")
  str(REP.LIN)
  cat("\n")
  
  genIE_Tier_1_wide_Kolf2<-genIE_Tier_1_wide[!is.na(genIE_Tier_1_wide$Kolf2),]
  
  cat("genIE_Tier_1_wide_Kolf2\n")
  str(genIE_Tier_1_wide_Kolf2)
  cat("\n")
  
  genIE_RESULTS_ACCEPTED_wide_Kolf2<-genIE_RESULTS_ACCEPTED_wide[!is.na(genIE_RESULTS_ACCEPTED_wide$Kolf2),]
  
  cat("genIE_RESULTS_ACCEPTED_wide_Kolf2\n")
  str(genIE_RESULTS_ACCEPTED_wide_Kolf2)
  cat("\n")
  
  REP.LIN$Status[which(REP.LIN$Cell_Type == "Kolf2" & REP.LIN$VAR%in%genIE_RESULTS_ACCEPTED_wide_Kolf2$VAR)]<-0
  REP.LIN$Status[which(REP.LIN$Cell_Type == "Kolf2" & REP.LIN$VAR%in%genIE_Tier_1_wide_Kolf2$VAR)]<-1
  
  genIE_Tier_1_wide_K562<-genIE_Tier_1_wide[!is.na(genIE_Tier_1_wide$K562),]
  
  cat("genIE_Tier_1_wide_K562\n")
  str(genIE_Tier_1_wide_K562)
  cat("\n")
  
  genIE_RESULTS_ACCEPTED_wide_K562<-genIE_RESULTS_ACCEPTED_wide[!is.na(genIE_RESULTS_ACCEPTED_wide$K562),]
  
  cat("genIE_RESULTS_ACCEPTED_wide_K562\n")
  str(genIE_RESULTS_ACCEPTED_wide_K562)
  cat("\n")
  
  REP.LIN$Status[which(REP.LIN$Cell_Type == "K562" & REP.LIN$VAR%in%genIE_RESULTS_ACCEPTED_wide_K562$VAR)]<-0
  REP.LIN$Status[which(REP.LIN$Cell_Type == "K562" & REP.LIN$VAR%in%genIE_Tier_1_wide_K562$VAR)]<-1
  
  genIE_Tier_1_wide_HL60<-genIE_Tier_1_wide[!is.na(genIE_Tier_1_wide$HL60),]
  
  cat("genIE_Tier_1_wide_HL60\n")
  str(genIE_Tier_1_wide_HL60)
  cat("\n")
  
  genIE_RESULTS_ACCEPTED_wide_HL60<-genIE_RESULTS_ACCEPTED_wide[!is.na(genIE_RESULTS_ACCEPTED_wide$HL60),]
  
  cat("genIE_RESULTS_ACCEPTED_wide_HL60\n")
  str(genIE_RESULTS_ACCEPTED_wide_HL60)
  cat("\n")
  
  REP.LIN$Status[which(REP.LIN$Cell_Type == "HL60" & REP.LIN$VAR%in%genIE_RESULTS_ACCEPTED_wide_HL60$VAR)]<-0
  REP.LIN$Status[which(REP.LIN$Cell_Type == "HL60" & REP.LIN$VAR%in%genIE_Tier_1_wide_HL60$VAR)]<-1
  
  
  genIE_Tier_1_wide_THP1<-genIE_Tier_1_wide[!is.na(genIE_Tier_1_wide$THP1),]
  
  cat("genIE_Tier_1_wide_THP1\n")
  str(genIE_Tier_1_wide_THP1)
  cat("\n")
  
  genIE_RESULTS_ACCEPTED_wide_THP1<-genIE_RESULTS_ACCEPTED_wide[!is.na(genIE_RESULTS_ACCEPTED_wide$THP1),]
  
  cat("genIE_RESULTS_ACCEPTED_wide_THP1\n")
  str(genIE_RESULTS_ACCEPTED_wide_THP1)
  cat("\n")
  
  REP.LIN$Status[which(REP.LIN$Cell_Type == "THP1" & REP.LIN$VAR%in%genIE_RESULTS_ACCEPTED_wide_THP1$VAR)]<-0
  REP.LIN$Status[which(REP.LIN$Cell_Type == "THP1" & REP.LIN$VAR%in%genIE_Tier_1_wide_THP1$VAR)]<-1
  

  
  REP.LIN$Status<-as.character(REP.LIN$Status)
  REP.LIN$Status<-factor(REP.LIN$Status,
                         levels=c("0","1"),
                         ordered=T)
  
  
  cat("REP.LIN_1\n")
  cat(str(REP.LIN))
  cat("\n") 
  
  
  
  cat("REP.LIN_1\n")
  cat(str(REP.LIN))
  cat("\n")
  cat(sprintf(as.character(names(summary(interaction(REP.LIN$LINEAGE_CLASSIF_DEF,REP.LIN$Status))))))
  cat("\n")
  cat(sprintf(as.character(summary(interaction(REP.LIN$LINEAGE_CLASSIF_DEF,REP.LIN$Status)))))
  cat("\n")
  
  
  
  
 
 #### First plot 0,1 breakdown -----------------------------------------------------------------------------------------------------------------------------------------------------
  
  
  MASTER_TABLE_STATUS_LINEAGE<-as.data.frame(setDT(REP.LIN)[, .N, .(Cell_Type,LINEAGE_CLASSIF_DEF,Status)], stringsAsFactors=F)
  
  colnames(MASTER_TABLE_STATUS_LINEAGE)[which(colnames(MASTER_TABLE_STATUS_LINEAGE) == "N")]<-"Associations"
  
  cat("MASTER_TABLE_STATUS_LINEAGE_\n")
  cat(str(MASTER_TABLE_STATUS_LINEAGE))
  cat("\n")
  
  
  
  
  # REP.LIN$Status[which(REP.LIN$Cell_Type == "Kolf2" & REP.LIN$VAR%in%genIE_Tier_1_wide_Kolf2$)]<-1
  
  
  MASTER_TABLE_STATUS<-as.data.frame(setDT(MASTER_TABLE_STATUS_LINEAGE)[, .(Total_associations=sum(Associations)), 
                                                                        .(Cell_Type,Status)], stringsAsFactors=F)
  
 
  
  cat("MASTER_TABLE_STATUS_\n")
  cat(str(MASTER_TABLE_STATUS))
  cat("\n")
  
  MASTER_TABLE_Lineage<-as.data.frame(setDT(MASTER_TABLE_STATUS_LINEAGE)[, .(Total_associations_per_Lineage=sum(Associations)), 
                                                                        .(Cell_Type,LINEAGE_CLASSIF_DEF)], stringsAsFactors=F)
  
  MASTER_TABLE_Lineage<-unique(MASTER_TABLE_Lineage[,-which(colnames(MASTER_TABLE_Lineage) == "Cell_Type")])
  
  MASTER_TABLE_Lineage$LINEAGE_CLASSIF_DEF<-factor(MASTER_TABLE_Lineage$LINEAGE_CLASSIF_DEF,
                                                          levels=c("Erythroid Traits","Megakaryocitic Traits","GRAN-MONO Traits","Lymphocitic Traits","Pleiotropic_variants"),
                                                          ordered=T)
  
  cat("MASTER_TABLE_Lineage_\n")
  cat(str(MASTER_TABLE_Lineage))
  cat("\n")
  
  
  # quit(status = 1)
  
  MASTER_TABLE_STATUS_LINEAGE<-as.data.frame(merge(MASTER_TABLE_STATUS_LINEAGE,
                                     MASTER_TABLE_STATUS,
                                     by=c("Cell_Type","Status")), stringsAsFactors=F)
  
 MASTER_TABLE_STATUS_LINEAGE$Perc<-round(100*(MASTER_TABLE_STATUS_LINEAGE$Associations/MASTER_TABLE_STATUS_LINEAGE$Total_associations),1)
  
  MASTER_TABLE_STATUS_LINEAGE$Cell_Type<-factor(MASTER_TABLE_STATUS_LINEAGE$Cell_Type,
                               levels=c("Kolf2","K562","HL60","THP1"),
                               ordered=T)
  MASTER_TABLE_STATUS_LINEAGE$Status<-factor(MASTER_TABLE_STATUS_LINEAGE$Status,
                            levels=c("0","1"),
                            ordered=T)
  
  MASTER_TABLE_STATUS_LINEAGE$LINEAGE_CLASSIF_DEF<-factor(MASTER_TABLE_STATUS_LINEAGE$LINEAGE_CLASSIF_DEF,
                                         levels=c("Erythroid Traits","Megakaryocitic Traits","GRAN-MONO Traits","Lymphocitic Traits","Pleiotropic_variants"),
                                         ordered=T)
  
  
  
  
  cat("MASTER_TABLE_STATUS_LINEAGE_\n")
  cat(str(MASTER_TABLE_STATUS_LINEAGE))
  cat("\n")
  
  cat(sprintf(as.character(levels(MASTER_TABLE_STATUS_LINEAGE$LINEAGE_CLASSIF_DEF))))
  cat("\n")
  
  setwd(out)
  write.table(file="test.tsv", MASTER_TABLE_STATUS_LINEAGE, sep="\t", quote=F, row.names=F)
  write.table(file="test_2.tsv", MASTER_TABLE_Lineage, sep="\t", quote=F, row.names=F)
  
  MASTER_TABLE_STATUS_LINEAGE_NO_NA<-unique(MASTER_TABLE_STATUS_LINEAGE[!is.na(MASTER_TABLE_STATUS_LINEAGE$Status),])
  
  
  cat("MASTER_TABLE_STATUS_LINEAGE_NO_NA_\n")
  cat(str(MASTER_TABLE_STATUS_LINEAGE_NO_NA))
  cat("\n")
  
  step<-round((100 - 0)/10,0)
  breaks.Rank<-seq(0,100+step,by=step)
  breaks.Rank<-seq(0,100,by=step)
  labels.Rank<-as.character(breaks.Rank)
  
    
  # scale_y_continuous(name=NULL,labels = scales::percent_format())+
    

    
    
  
  graph<-MASTER_TABLE_STATUS_LINEAGE_NO_NA %>%
    mutate(myaxis = paste0(Status, "\n", "n=", Total_associations)) %>%
    mutate(myaxis=fct_reorder(myaxis,as.numeric(Status))) %>%
    ggplot(aes(x=myaxis, y=Perc, fill=LINEAGE_CLASSIF_DEF)) +
    geom_bar(stat="identity",colour='black')+
    theme_bw()+
    theme(axis.title.y=element_text(size=24, family="sans"),
          axis.title.x=element_text(size=24, family="sans"),
          axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
          axis.text.x=element_text(angle=45,size=10, color="black", family="sans"),
          legend.title=element_text(size=16,color="black", family="sans"),
          legend.text=element_text(size=12,color="black", family="sans"))+
    scale_x_discrete(name=NULL, drop=F)+
    scale_fill_manual(values=c('#553B68','#32A852','#6DB2EE','#1877C9','#62D07F',"gray","gray","gray"),drop=F,
                      name="Totals",breaks=MASTER_TABLE_Lineage$LINEAGE_CLASSIF_DEF,
                      labels=paste(MASTER_TABLE_Lineage$LINEAGE_CLASSIF_DEF,
                                   MASTER_TABLE_Lineage$Total_associations_per_Lineage, sep =' n= '))+
    scale_y_continuous(name="Relative abundance of each class(%)",breaks=breaks.Rank,labels=labels.Rank, limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
    theme(legend.position="bottom")+
    guides(fill=guide_legend(nrow=3,byrow=TRUE)) +
    ggeasy::easy_center_title()
  
  graph <- graph + facet_grid(. ~ Cell_Type, drop=F)
  
  setwd(out)
  
  svgname<-paste("Lineage_plot_1_","INTRACELL_ALL_STATUS",".svg",sep='')
  makesvg = TRUE
  
  if (makesvg == TRUE)
  {
    ggsave(svgname, plot= graph,
           device="svg",
           height=10, width=12)
  }
  
  
  
  
  
  
  
  
  #### Second plot only active and referred to total associations-------------------------------------------------------------------------------------------------------------------------------------
  
  MASTER_TABLE_STATUS_LINEAGE<-as.data.frame(setDT(REP.LIN)[, .N, .(Cell_Type,LINEAGE_CLASSIF_DEF,Status)], stringsAsFactors=F)
  
  colnames(MASTER_TABLE_STATUS_LINEAGE)[which(colnames(MASTER_TABLE_STATUS_LINEAGE) == "N")]<-"Associations"
  
  cat("MASTER_TABLE_STATUS_LINEAGE_0\n")
  cat(str(MASTER_TABLE_STATUS_LINEAGE))
  cat("\n")
  
  
  MASTER_TABLE_STATUS<-as.data.frame(setDT(MASTER_TABLE_STATUS_LINEAGE)[, .(Total_associations=sum(Associations)), 
                                                                        .(Cell_Type)], stringsAsFactors=F)
  
  
  
  cat("MASTER_TABLE_STATUS_\n")
  cat(str(MASTER_TABLE_STATUS))
  cat("\n")
  
  
  MASTER_TABLE_STATUS_LINEAGE<-as.data.frame(merge(MASTER_TABLE_STATUS_LINEAGE,
                                                   MASTER_TABLE_STATUS,
                                                   by=c("Cell_Type")), stringsAsFactors=F)
  
  MASTER_TABLE_STATUS_LINEAGE$Perc<-round(100*(MASTER_TABLE_STATUS_LINEAGE$Associations/MASTER_TABLE_STATUS_LINEAGE$Total_associations),1) 
  cat("MASTER_TABLE_STATUS_LINEAGE_1\n")
  cat(str(MASTER_TABLE_STATUS_LINEAGE))
  cat("\n")
  
  
  ACTIVE_TOTAL_PERCENTAGES<-as.data.frame(setDT(MASTER_TABLE_STATUS_LINEAGE)[, .(Total_percentage_by_status=sum(Perc)), .(Cell_Type,Status)], stringsAsFactors=F)
  
  
  cat("ACTIVE_TOTAL_PERCENTAGES_0\n")
  cat(str(ACTIVE_TOTAL_PERCENTAGES))
  cat("\n")
  
  MASTER_TABLE_STATUS_LINEAGE<-as.data.frame(merge(MASTER_TABLE_STATUS_LINEAGE,
                                                   ACTIVE_TOTAL_PERCENTAGES,
                                                   by=c("Cell_Type","Status")), stringsAsFactors=F)
  
  
  cat("MASTER_TABLE_STATUS_LINEAGE_2\n")
  cat(str(MASTER_TABLE_STATUS_LINEAGE))
  cat("\n")
  
  MASTER_TABLE_STATUS_LINEAGE$Cell_Type<-factor(MASTER_TABLE_STATUS_LINEAGE$Cell_Type,
                                                levels=rev(c("Kolf2","K562","HL60","THP1")),
                                                ordered=T)
  MASTER_TABLE_STATUS_LINEAGE$Status<-factor(MASTER_TABLE_STATUS_LINEAGE$Status,
                                             levels=c("0","1"),
                                             ordered=T)
  
  MASTER_TABLE_STATUS_LINEAGE$LINEAGE_CLASSIF_DEF<-factor(MASTER_TABLE_STATUS_LINEAGE$LINEAGE_CLASSIF_DEF,
                                                          levels=c("Erythroid Traits","Megakaryocitic Traits","GRAN-MONO Traits","Lymphocitic Traits","Pleiotropic_variants"),
                                                          ordered=T)
  
  cat("MASTER_TABLE_STATUS_LINEAGE_3\n")
  cat(str(MASTER_TABLE_STATUS_LINEAGE))
  cat("\n")
  
  
  cat(sprintf(as.character(levels(MASTER_TABLE_STATUS_LINEAGE$LINEAGE_CLASSIF_DEF))))
  cat("\n")
  
  setwd(out)
  write.table(file="test.tsv", MASTER_TABLE_STATUS_LINEAGE, sep="\t", quote=F, row.names=F)
  
  
  MASTER_TABLE_STATUS_LINEAGE_ACTIVE<-MASTER_TABLE_STATUS_LINEAGE[which(MASTER_TABLE_STATUS_LINEAGE$Status == "1"),]
  
  
  cat("MASTER_TABLE_STATUS_LINEAGE_ACTIVE_\n")
  cat(str(MASTER_TABLE_STATUS_LINEAGE_ACTIVE))
  cat("\n")
  
  NEW_TOTAL<-as.data.frame(setDT(MASTER_TABLE_STATUS_LINEAGE_ACTIVE)[, .(Total_associations_ACTIVE=sum(Associations)), 
                                                                        .(Cell_Type)], stringsAsFactors=F)
  
  cat("NEW_TOTAL_\n")
  cat(str(NEW_TOTAL))
  cat("\n")
  
  
  MASTER_TABLE_STATUS_LINEAGE_ACTIVE<-as.data.frame(merge(MASTER_TABLE_STATUS_LINEAGE_ACTIVE,
                                                          NEW_TOTAL,
                                                   by=c("Cell_Type")), stringsAsFactors=F)
  
  
  cat("MASTER_TABLE_STATUS_LINEAGE_ACTIVE_2\n")
  cat(str(MASTER_TABLE_STATUS_LINEAGE_ACTIVE))
  cat("\n")
  
  # quit(status = 1)
  
  
  step<-round((max(MASTER_TABLE_STATUS_LINEAGE_ACTIVE$Total_percentage_by_status) - 0)/10,0)
  
  breaks.Rank<-seq(0,max(MASTER_TABLE_STATUS_LINEAGE_ACTIVE$Total_percentage_by_status)+step,by=step)
  labels.Rank<-as.character(breaks.Rank)
  
  cat(sprintf(as.character(labels.Rank)))
  cat("\n")
  
  
  df.color_Cell_Type_sel<-df.color_Cell_Type[which(df.color_Cell_Type$Cell_Type%in%MASTER_TABLE_STATUS_LINEAGE_ACTIVE$Cell_Type),]
  
  
  df.color_Cell_Type_sel<-df.color_Cell_Type_sel[order(factor(df.color_Cell_Type_sel$Cell_Type,
                                                              levels=levels(MASTER_TABLE_STATUS_LINEAGE_ACTIVE$Cell_Type),
                                                              ordered=T)),]
  
  # MASTER_TABLE_STATUS_LINEAGE_ACTIVE<-as.data.frame(merge(MASTER_TABLE_STATUS_LINEAGE_ACTIVE,
  #                                                         df.color_Cell_Type_sel,
  #                                                         by=c("Cell_Type")), stringsAsFactors=F)
  # 
  # cat("MASTER_TABLE_STATUS_LINEAGE_ACTIVE_3\n")
  # cat(str(MASTER_TABLE_STATUS_LINEAGE_ACTIVE))
  # cat("\n")
 
  
  graph2<- MASTER_TABLE_STATUS_LINEAGE_ACTIVE %>%
    mutate(myaxis = paste0(Cell_Type, "\n", "n=", Total_associations_ACTIVE)) %>%
    mutate(myaxis=fct_reorder(myaxis,as.numeric(Cell_Type))) %>%
    ggplot(aes(x=Perc, 
                     y=myaxis, 
                     fill=LINEAGE_CLASSIF_DEF)) +
    geom_bar(stat="identity",colour='black')+
    theme_bw()+
    theme(axis.title.y=element_text(size=24, family="sans"),
          axis.title.x=element_text(size=24, family="sans"),
          axis.text.y=element_text(angle=0,size=18, color=df.color_Cell_Type_sel$color, family="sans"),
          axis.text.x=element_text(angle=0,size=18, color="black", family="sans"),
          legend.title=element_text(size=16,color="black", family="sans"),
          legend.text=element_text(size=12,color="black", family="sans"))+
    scale_x_continuous(name="Perc class active / total Elements",breaks=breaks.Rank,labels=labels.Rank, limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
    scale_fill_manual(values=c('#32A852','#6DB2EE','#553B68','#62D07F','#1877C9','#C9244B','#D45E85','#87447B','#D6B8E6'),drop=F)+
    scale_y_discrete(name=NULL, drop=F)+
    theme(legend.position="bottom")+
    guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
    theme(legend.title = element_blank())+
    ggeasy::easy_center_title()
  
  # scale_x_continuous(name=NULL,labels = scales::percent_format())+
  
  setwd(out)
  
  svgname<-paste("Lineage_plot_2_","INTRACELL_ACTIVE",".svg",sep='')
  makesvg = TRUE
  
  if (makesvg == TRUE)
  {
    ggsave(svgname, plot= graph2,
           device="svg",
           height=10, width=12)
  }
  
  
  # pdf(file="test.pdf")
  # print(graph)
  # print(graph2)
  # dev.off()
  
  
  
  ##### EXPECTED/OBSERVED ChiSq test --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  
  MASTER_TABLE_STATUS_LINEAGE<-as.data.frame(setDT(REP.LIN)[, .N, .(Cell_Type,LINEAGE_CLASSIF_DEF,Status)], stringsAsFactors=F)
  
  colnames(MASTER_TABLE_STATUS_LINEAGE)[which(colnames(MASTER_TABLE_STATUS_LINEAGE) == "N")]<-"Associations"
  
  cat("MASTER_TABLE_STATUS_LINEAGE_0\n")
  cat(str(MASTER_TABLE_STATUS_LINEAGE))
  cat("\n")
  
  MASTER_TABLE_STATUS<-as.data.frame(setDT(MASTER_TABLE_STATUS_LINEAGE)[, .(Total_associations=sum(Associations)), 
                                                                        .(Cell_Type,Status)], stringsAsFactors=F)
  
  
  MASTER_TABLE_STATUS$Status<-paste("Status_",MASTER_TABLE_STATUS$Status,sep='')
  
  cat("MASTER_TABLE_STATUS_\n")
  cat(str(MASTER_TABLE_STATUS))
  cat("\n")
  
  MASTER_TABLE_STATUS_wide<-as.data.frame(pivot_wider(MASTER_TABLE_STATUS,
                            id_cols=Cell_Type,
                            names_from=Status,
                            values_from=Total_associations))
  

  
  MASTER_TABLE_STATUS_wide$Frequency_ACTIVE<-MASTER_TABLE_STATUS_wide$Status_1/(MASTER_TABLE_STATUS_wide$Status_0+MASTER_TABLE_STATUS_wide$Status_1)
  
  cat("MASTER_TABLE_STATUS_wide_\n")
  cat(str(MASTER_TABLE_STATUS_wide))
  cat("\n")
  
 
  
  MASTER_TABLE_Lineage<-as.data.frame(setDT(MASTER_TABLE_STATUS_LINEAGE)[, .(Total_associations_per_Lineage=sum(Associations)),
                                                                         .(Cell_Type,LINEAGE_CLASSIF_DEF)], stringsAsFactors=F)




  cat("MASTER_TABLE_Lineage_\n")
  cat(str(MASTER_TABLE_Lineage))
  cat("\n")
  
  MASTER_TABLE_Lineage_wide<-as.data.frame(pivot_wider(MASTER_TABLE_Lineage,
                                                      id_cols=Cell_Type,
                                                      names_from=LINEAGE_CLASSIF_DEF,
                                                      values_from=Total_associations_per_Lineage))
  
  
  cat("MASTER_TABLE_Lineage_wide_0\n")
  cat(str(MASTER_TABLE_Lineage_wide))
  cat("\n")
  
  # MASTER_TABLE_Lineage_wide$Frequency_ACTIVE<-MASTER_TABLE_Lineage_wide$Status_1/(MASTER_TABLE_Lineage_wide$Status_0+MASTER_TABLE_Lineage_wide$Status_1)
  
  MASTER_TABLE_Lineage_wide$Total_associations<-MASTER_TABLE_Lineage_wide[,2]+MASTER_TABLE_Lineage_wide[,3]+MASTER_TABLE_Lineage_wide[,4]+MASTER_TABLE_Lineage_wide[,5]
  
  cat("MASTER_TABLE_Lineage_wide_\n")
  cat(str(MASTER_TABLE_Lineage_wide))
  cat("\n")
  
  MASTER_TABLE_Lineage_wide$Freq_gran_mono<-MASTER_TABLE_Lineage_wide[,2]/MASTER_TABLE_Lineage_wide$Total_associations
  MASTER_TABLE_Lineage_wide$Freq_pleio<-MASTER_TABLE_Lineage_wide[,3]/MASTER_TABLE_Lineage_wide$Total_associations
  MASTER_TABLE_Lineage_wide$Freq_erythro<-MASTER_TABLE_Lineage_wide[,4]/MASTER_TABLE_Lineage_wide$Total_associations
  MASTER_TABLE_Lineage_wide$Freq_mega<-MASTER_TABLE_Lineage_wide[,5]/MASTER_TABLE_Lineage_wide$Total_associations
  #MASTER_TABLE_Lineage_wide$Freq_lymph<-MASTER_TABLE_Lineage_wide[,6]/MASTER_TABLE_Lineage_wide$Total_associations
  
  cat("MASTER_TABLE_Lineage_wide_\n")
  cat(str(MASTER_TABLE_Lineage_wide))
  cat("\n")
  
  expectations_table<-merge(MASTER_TABLE_STATUS_wide,
                            MASTER_TABLE_Lineage_wide,
                            by="Cell_Type")
  
  cat("expectations_table_0\n")
  cat(str(expectations_table))
  cat("\n")
  
  
  expectations_table$expected_gran_mono_ACTIVE<-round(expectations_table$Freq_gran_mono*expectations_table$Frequency_ACTIVE*expectations_table$Total_associations,0)
  expectations_table$expected_pleio_ACTIVE<-round(expectations_table$Freq_pleio*expectations_table$Frequency_ACTIVE*expectations_table$Total_associations,0)
  expectations_table$expected_erythro_ACTIVE<-round(expectations_table$Freq_erythro*expectations_table$Frequency_ACTIVE*expectations_table$Total_associations,0)
  expectations_table$expected_mega_ACTIVE<-round(expectations_table$Freq_mega*expectations_table$Frequency_ACTIVE*expectations_table$Total_associations,0)
 # expectations_table$expected_lymph_ACTIVE<-round(expectations_table$Freq_lymph*expectations_table$Frequency_ACTIVE*expectations_table$Total_associations,0)
  
  
  cat("expectations_table_1\n")
  cat(str(expectations_table))
  cat("\n")
  
  # 
  # MASTER_TABLE_STATUS_LINEAGE_ACTIVE<-MASTER_TABLE_STATUS_LINEAGE[which(MASTER_TABLE_STATUS_LINEAGE$Status == "1"),]
  # 
  # cat("MASTER_TABLE_STATUS_LINEAGE_ACTIVE_1\n")
  # cat(str(MASTER_TABLE_STATUS_LINEAGE_ACTIVE))
  # cat("\n")
  
  MASTER_TABLE_STATUS_LINEAGE_wide<-as.data.frame(pivot_wider(MASTER_TABLE_STATUS_LINEAGE,
                                                       id_cols=c("Cell_Type","Status"),
                                                       names_from=LINEAGE_CLASSIF_DEF,
                                                       values_from=Associations), stringsAsFactors=F)
  
  cat("MASTER_TABLE_STATUS_LINEAGE_wide_1\n")
  cat(str(MASTER_TABLE_STATUS_LINEAGE_wide))
  cat("\n")
  
  colnames(MASTER_TABLE_STATUS_LINEAGE_wide)[-c(1:2)]<-c("observed_gran_mono_ACTIVE","observed_pleio_ACTIVE","observed_erythro_ACTIVE",
                                                         "observed_mega_ACTIVE")#,"observed_lymph_ACTIVE")
  
  
  cat("MASTER_TABLE_STATUS_LINEAGE_wide_2\n")
  cat(str(MASTER_TABLE_STATUS_LINEAGE_wide))
  cat("\n")
  
  
  MASTER_TABLE_STATUS_LINEAGE_wide_ACTIVE<-MASTER_TABLE_STATUS_LINEAGE_wide[which(MASTER_TABLE_STATUS_LINEAGE_wide$Status == "1"),]
  
  cat("MASTER_TABLE_STATUS_LINEAGE_wide_ACTIVE\n")
  cat(str(MASTER_TABLE_STATUS_LINEAGE_wide_ACTIVE))
  cat("\n")
  
  expectations_observations_table<-merge(expectations_table,
                   MASTER_TABLE_STATUS_LINEAGE_wide_ACTIVE,
                   by="Cell_Type",
                   all.x=T)
  
  cat("expectations_observations_table\n")
  cat(str(expectations_observations_table))
  cat("\n")
  
  expectations_observations_table[is.na(expectations_observations_table)]<-0
  
  
 
  
  expectations_observations_table$gran_mono_non_active<-expectations_observations_table[,8]-expectations_observations_table$observed_gran_mono_ACTIVE
  expectations_observations_table$gran_mono_non_active[(expectations_observations_table$gran_mono_non_active < 0)]<-0
  expectations_observations_table$no_gran_mono_non_active<-expectations_observations_table$Status_0-expectations_observations_table$gran_mono_non_active
  expectations_observations_table$no_gran_mono_non_active[(expectations_observations_table$no_gran_mono_non_active < 0)]<-0
  expectations_observations_table$no_gran_mono_active<-expectations_observations_table$Status_1-expectations_observations_table$observed_gran_mono_ACTIVE
  expectations_observations_table$no_gran_mono_active[(expectations_observations_table$no_gran_mono_active < 0)]<-0
  
  ##
  
 
  
  expectations_observations_table$pleio_non_active<-expectations_observations_table[,8]-expectations_observations_table$observed_pleio_ACTIVE
  expectations_observations_table$pleio_non_active[(expectations_observations_table$pleio_non_active < 0)]<-0
  expectations_observations_table$no_pleio_non_active<-expectations_observations_table$Status_0-expectations_observations_table$pleio_non_active
  expectations_observations_table$no_pleio_non_active[(expectations_observations_table$no_pleio_non_active < 0)]<-0
  expectations_observations_table$no_pleio_active<-expectations_observations_table$Status_1-expectations_observations_table$observed_pleio_ACTIVE
  expectations_observations_table$no_pleio_active[(expectations_observations_table$no_pleio_active < 0)]<-0
  
  ##
  
 
  expectations_observations_table$erythro_non_active<-expectations_observations_table[,8]-expectations_observations_table$observed_erythro_ACTIVE
  expectations_observations_table$erythro_non_active[(expectations_observations_table$erythro_non_active < 0)]<-0
  expectations_observations_table$no_erythro_non_active<-expectations_observations_table$Status_0-expectations_observations_table$erythro_non_active
  expectations_observations_table$no_erythro_non_active[(expectations_observations_table$no_erythro_non_active < 0)]<-0
  expectations_observations_table$no_erythro_active<-expectations_observations_table$Status_1-expectations_observations_table$observed_erythro_ACTIVE
  expectations_observations_table$no_erythro_active[(expectations_observations_table$no_erythro_active < 0)]<-0
  
  
  ##
  
  expectations_observations_table$mega_non_active<-expectations_observations_table[,8]-expectations_observations_table$observed_mega_ACTIVE
  expectations_observations_table$mega_non_active[(expectations_observations_table$mega_non_active < 0)]<-0
  expectations_observations_table$no_mega_non_active<-expectations_observations_table$Status_0-expectations_observations_table$mega_non_active
  expectations_observations_table$no_mega_non_active[(expectations_observations_table$no_mega_non_active < 0)]<-0
  expectations_observations_table$no_mega_active<-expectations_observations_table$Status_1-expectations_observations_table$observed_mega_ACTIVE
  expectations_observations_table$no_mega_active[(expectations_observations_table$no_mega_active < 0)]<-0
  
  
  ##
  
  # expectations_observations_table$lymph_non_active<-expectations_observations_table[,9]-expectations_observations_table$observed_lymph_ACTIVE
  # 
  # expectations_observations_table$no_lymph_non_active<-expectations_observations_table$Status_0-expectations_observations_table$lymph_non_active
  # 
  # expectations_observations_table$no_lymph_active<-expectations_observations_table$Status_1-expectations_observations_table$observed_lymph_ACTIVE
  # 
  ##
  
  
  
  
    
  cat("expectations_observations_table_1\n")
  cat(str(expectations_observations_table))
  cat("\n")
  
  # setwd(out)
  # write.table(file="test_3.tsv",expectations_observations_table, sep="\t",quote=F, row.names = F)
  
  
  
  indx.int<-c(colnames(expectations_observations_table)[which(colnames(expectations_observations_table) == "Cell_Type")],
              colnames(expectations_observations_table)[which(colnames(expectations_observations_table) == "gran_mono_non_active")],colnames(expectations_observations_table)[which(colnames(expectations_observations_table) == "no_gran_mono_non_active")],
              colnames(expectations_observations_table)[which(colnames(expectations_observations_table) == "observed_gran_mono_ACTIVE")],colnames(expectations_observations_table)[which(colnames(expectations_observations_table) == "no_gran_mono_active")],
              colnames(expectations_observations_table)[which(colnames(expectations_observations_table) == "pleio_non_active")],colnames(expectations_observations_table)[which(colnames(expectations_observations_table) == "no_pleio_non_active")],
              colnames(expectations_observations_table)[which(colnames(expectations_observations_table) == "observed_pleio_ACTIVE")],colnames(expectations_observations_table)[which(colnames(expectations_observations_table) == "no_pleio_active")],
              colnames(expectations_observations_table)[which(colnames(expectations_observations_table) == "erythro_non_active")],colnames(expectations_observations_table)[which(colnames(expectations_observations_table) == "no_erythro_non_active")],
              colnames(expectations_observations_table)[which(colnames(expectations_observations_table) == "observed_erythro_ACTIVE")],colnames(expectations_observations_table)[which(colnames(expectations_observations_table) == "no_erythro_active")],
              colnames(expectations_observations_table)[which(colnames(expectations_observations_table) == "mega_non_active")],colnames(expectations_observations_table)[which(colnames(expectations_observations_table) == "no_mega_non_active")],
              colnames(expectations_observations_table)[which(colnames(expectations_observations_table) == "observed_mega_ACTIVE")],colnames(expectations_observations_table)[which(colnames(expectations_observations_table) == "no_mega_active")])
  
  # colnames(expectations_observations_table)[which(colnames(expectations_observations_table) == "lymph_non_active")],colnames(expectations_observations_table)[which(colnames(expectations_observations_table) == "no_lymph_non_active")],
  # colnames(expectations_observations_table)[which(colnames(expectations_observations_table) == "observed_lymph_ACTIVE")],colnames(expectations_observations_table)[which(colnames(expectations_observations_table) == "no_lymph_active")]
  # 
  
  expectations_observations_table_subset<-unique(expectations_observations_table[,indx.int])
  
  
  cat("expectations_observations_table_0\n")
  cat(str(expectations_observations_table_subset))
  cat("\n")
  
    
  expectations_observations_table.m<-melt(expectations_observations_table_subset,
                                          id.variables=c("Cell_Type"),
                                          variable.name = "variable", 
                                          value.name = "n")
  
  cat("expectations_observations_table.m_1\n")
  cat(str(expectations_observations_table.m))
  cat("\n")
  
  
  expectations_observations_table.m$group<-"NA"
  
  expectations_observations_table.m$group[grep("gran_mono",expectations_observations_table.m$variable)]<-"gran_mono"
  expectations_observations_table.m$group[grep("erythro",expectations_observations_table.m$variable)]<-"erythro"
  expectations_observations_table.m$group[grep("pleio",expectations_observations_table.m$variable)]<-"pleio"
  # expectations_observations_table.m$group[grep("lymph",expectations_observations_table.m$variable)]<-"lymph"
  expectations_observations_table.m$group[grep("mega",expectations_observations_table.m$variable)]<-"mega"
  
 
  cat("expectations_observations_table.m_2\n")
  cat(str(expectations_observations_table.m))
  cat("\n")
  
  expectations_observations_table.m$status<-"NA"
  
  expectations_observations_table.m$status[grep("_ACTIVE$|_active$",expectations_observations_table.m$variable)]<-"active"
  expectations_observations_table.m$status[grep("non_active$",expectations_observations_table.m$variable)]<-"non_active"
  
  expectations_observations_table.m$belonging<-"belonging"
  
  expectations_observations_table.m$belonging[grep("^no_",expectations_observations_table.m$variable)]<-"not_belonging"

  # "lymph",
  
  expectations_observations_table.m$group<-factor(expectations_observations_table.m$group,
                                                  levels=c("erythro","mega","gran_mono","pleio"), 
                                                  ordered=T)
  expectations_observations_table.m$status<-factor(expectations_observations_table.m$status,
                                                  levels=c("active","non_active"),
                                                  ordered=T)
  expectations_observations_table.m$belonging<-factor(expectations_observations_table.m$belonging,
                                                  levels=c("belonging","not_belonging"),
                                                  ordered=T)
  
  
  cat(sprintf(as.character(names(summary(expectations_observations_table.m$group)))))
  cat("\n")
  cat(sprintf(as.character(summary(expectations_observations_table.m$group))))
  cat("\n")
  
  cat(sprintf(as.character(names(summary(expectations_observations_table.m$status)))))
  cat("\n")
  cat(sprintf(as.character(summary(expectations_observations_table.m$status))))
  cat("\n")
  
  cat(sprintf(as.character(names(summary(expectations_observations_table.m$belonging)))))
  cat("\n")
  cat(sprintf(as.character(summary(expectations_observations_table.m$belonging))))
  cat("\n")
  
  
  cat("expectations_observations_table.m_3\n")
  cat(str(expectations_observations_table.m))
  cat("\n")
  
  # setwd(out)
  # write.table(file="test_3.tsv",expectations_observations_table.m, sep="\t",quote=F, row.names = F)
  
  # quit(status = 1)
  
  
 
  
  ##### LOOP for ChiSq calculation INTRACELL_ER --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  
  list_ChiSq_DEF<-list()
  
  
  for(i in 1:length(levels(expectations_observations_table.m$group)))
  {
    
    level_sel<-levels(expectations_observations_table.m$group)[i]
    
    cat("--->\t")
    cat(sprintf(as.character(level_sel)))
    cat("\n")
    
    
    
    expectations_observations_table.m_sel<-expectations_observations_table.m[which(expectations_observations_table.m$group == level_sel),]
    
    cat("expectations_observations_table.m_sel\n")
    cat(str(expectations_observations_table.m_sel))
    cat("\n")
    
    
    array_Cell_Type<-levels(expectations_observations_table.m_sel$Cell_Type)
    
    
    cat("array_Cell_Type\n")
    cat(str(array_Cell_Type))
    cat("\n")
    
    list_ChiSq<-list()
    
    for(k in 1:length(array_Cell_Type))
    {
      Cell_Type_sel<-array_Cell_Type[k]
      
      cat(">\t")
      cat(sprintf(as.character(Cell_Type_sel)))
      cat("\n")
      
      expectations_observations_table.m_sel_Cell_Type_sel<-expectations_observations_table.m_sel[which(expectations_observations_table.m_sel$Cell_Type == Cell_Type_sel),]
      
      cat("expectations_observations_table.m_sel_Cell_Type_sel\n")
      cat(str(expectations_observations_table.m_sel_Cell_Type_sel))
      cat("\n")
      
      
      setwd(out)
      write.table(file="test_2.tsv",expectations_observations_table.m_sel_Cell_Type_sel, sep="\t",quote=F, row.names = F)
     
      
      tab <- xtabs(n ~ status+belonging, expectations_observations_table.m_sel_Cell_Type_sel)

      cat("tab\n")
      cat(str(tab))
      cat("\n")
      
      tab.chisq.test<-chisq.test(tab,correct = TRUE)
      
      cat("tab.chisq.test\n")
      cat(str(tab.chisq.test))
      cat("\n")
      
      log_pval<--1*log10(tab.chisq.test$p.value)
      
      cat("log_pval\n")
      cat(str(log_pval))
      cat("\n")
      
      a.df<-as.data.frame(cbind(Cell_Type_sel,log_pval), stringsAsFactors=F)
      
      colnames(a.df)<-c("Cell_Type","log_pval")
      
      cat("a.df\n")
      cat(str(a.df))
      cat("\n")
      
      
      list_ChiSq[[k]]<-a.df
      
      # quit(status = 1)
      
    
      
      
    }#k
    
    
    partial_tab = unique(as.data.frame(data.table::rbindlist(list_ChiSq, fill = T)))
    
    partial_tab$Traits<-level_sel
    
    cat("partial_tab\n")
    cat(str(partial_tab))
    cat("\n")
    
    
  
    
    
    list_ChiSq_DEF[[i]]<-partial_tab
    
    
    # quit(status = 1)
    
  }
  
  INTRACELL_ER = unique(as.data.frame(data.table::rbindlist(list_ChiSq_DEF, fill = T)))
  
  
  cat("INTRACELL_ER\n")
  cat(str(INTRACELL_ER))
  cat("\n")
  
  
  setwd(out)
  write.table(file="test_4.tsv",INTRACELL_ER, sep="\t",quote=F, row.names = F)
  # 
  
 
  
  #### Third plot only active breakdown by Cell Type -------------------------------------------------------------------------------------------------------------------------------------
  
  
  
  MASTER_TABLE_STATUS_LINEAGE<-as.data.frame(setDT(REP.LIN)[, .N, .(Cell_Type,LINEAGE_CLASSIF_DEF,Status)], stringsAsFactors=F)
  
  colnames(MASTER_TABLE_STATUS_LINEAGE)[which(colnames(MASTER_TABLE_STATUS_LINEAGE) == "N")]<-"Associations"
  
  cat("MASTER_TABLE_STATUS_LINEAGE_0\n")
  cat(str(MASTER_TABLE_STATUS_LINEAGE))
  cat("\n")
  
  
  
  MASTER_TABLE_STATUS_LINEAGE_NO_NA<-MASTER_TABLE_STATUS_LINEAGE[!is.na(MASTER_TABLE_STATUS_LINEAGE$Status),]
  
  cat("MASTER_TABLE_STATUS_LINEAGE_NO_NA_1\n")
  cat(str(MASTER_TABLE_STATUS_LINEAGE_NO_NA))
  cat("\n")
  

  
  MASTER_TABLE_Lineage<-as.data.frame(setDT(MASTER_TABLE_STATUS_LINEAGE_NO_NA)[, .(Total_associations_per_Lineage=sum(Associations)), 
                                                                         .(Status,LINEAGE_CLASSIF_DEF)], stringsAsFactors=F)
  
  
  
  MASTER_TABLE_Lineage$LINEAGE_CLASSIF_DEF<-factor(MASTER_TABLE_Lineage$LINEAGE_CLASSIF_DEF,
                                                   levels=c("Erythroid Traits","Megakaryocitic Traits","GRAN-MONO Traits","Lymphocitic Traits","Pleiotropic_variants"),
                                                   ordered=T)
  
  MASTER_TABLE_Lineage<-MASTER_TABLE_Lineage[order(MASTER_TABLE_Lineage$LINEAGE_CLASSIF_DEF,MASTER_TABLE_Lineage$Status),]
  
  cat("MASTER_TABLE_Lineage_\n")
  cat(str(MASTER_TABLE_Lineage))
  cat("\n")
  
  
  Plot3_table<-as.data.frame(merge(MASTER_TABLE_STATUS_LINEAGE_NO_NA,
                                   MASTER_TABLE_Lineage,
                                   by=c("Status","LINEAGE_CLASSIF_DEF"),all.x=T), stringsAsFactors=F)
  
  
  
  
  Plot3_table$Cell_Type<-factor(Plot3_table$Cell_Type,
                                levels=c("Kolf2","K562","HL60","THP1"),
                                ordered=T)
  Plot3_table$Status<-factor(Plot3_table$Status,
                             levels=c("0","1"),
                             ordered=T)
  
  Plot3_table$LINEAGE_CLASSIF_DEF<-factor(Plot3_table$LINEAGE_CLASSIF_DEF,
                                          levels=rev(c("Erythroid Traits","Megakaryocitic Traits","GRAN-MONO Traits","Lymphocitic Traits","Pleiotropic_variants")),
                                          ordered=T)
  
  Plot3_table<-Plot3_table[order(Plot3_table$LINEAGE_CLASSIF_DEF,Plot3_table$Status,Plot3_table$Cell_Type),]
  
  Plot3_table$Perc<-round(100*(Plot3_table$Associations/Plot3_table$Total_associations_per_Lineage),1)
  
  
  cat("Plot3_table_\n")
  cat(str(Plot3_table))
  cat("\n")
  
  Plot3_table_ACTIVE<-Plot3_table[which(Plot3_table$Status == "1"),]
  
  cat("Plot3_table_ACTIVE_\n")
  cat(str(Plot3_table_ACTIVE))
  cat("\n")
  
  
  local_MAX<-max(Plot3_table$Total_associations_per_Lineage)#[!is.na(Plot3_table$Total_associations_per_Lineage)]
  
  
  cat("local_MAX_\n")
  cat(str(local_MAX))
  cat("\n")
  
  
  step<-round((max(Plot3_table$Total_associations_per_Lineage) - 0)/10,2)
  
  
  cat("step_\n")
  cat(str(step))
  cat("\n")
  
  
  breaks.Rank<-seq(0,max(Plot3_table_ACTIVE$Total_associations_per_Lineage)+step,by=step)
  labels.Rank<-as.character(breaks.Rank)
  
  cat(sprintf(as.character(labels.Rank)))
  cat("\n")
  
  df.color_Cell_Type_sel<-df.color_Cell_Type[which(df.color_Cell_Type$Cell_Type%in%levels(Plot3_table_ACTIVE$Cell_Type)),]
  
  
  df.color_Cell_Type_sel<-df.color_Cell_Type_sel[order(factor(df.color_Cell_Type_sel$Cell_Type,
                                                              levels=levels(Plot3_table_ACTIVE$Cell_Type),
                                                              ordered=T)),]
  
  cat("df.color_Cell_Type_sel_\n")
  cat(str(df.color_Cell_Type_sel))
  cat("\n")
  
  graph3<- ggplot(data=Plot3_table_ACTIVE, 
                  aes(x=Associations, 
                      y=LINEAGE_CLASSIF_DEF, 
                      fill=Cell_Type)) +
    geom_bar(stat="identity",colour='black')+
    theme_bw()+
    theme(axis.title.y=element_text(size=24, family="sans"),
          axis.title.x=element_text(size=24, family="sans"),
          axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
          axis.text.x=element_text(angle=0,size=18, color="black", family="sans"),
          legend.title=element_text(size=16,color="black", family="sans"),
          legend.text=element_text(size=12,color="black", family="sans"))+
    scale_x_continuous(name="ACTIVE associations per Cell_Type/total ACTIVE associations",breaks=breaks.Rank,labels=labels.Rank, limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
    scale_fill_manual(values=df.color_Cell_Type_sel$color,drop=F)+
    scale_y_discrete(name=NULL, drop=F)+
    theme(legend.position="bottom")+
    guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
    theme(legend.title = element_blank())+
    ggeasy::easy_center_title()
  
  svgname<-paste("Lineage_plot_3_","INTERCELL_ACTIVE",".svg",sep='')
  makesvg = TRUE
  
  if (makesvg == TRUE)
  {
    ggsave(svgname, plot= graph3,
           device="svg",
           height=10, width=12)
  }
  
  
  
  # setwd(out)
  # pdf(file="test.pdf")
  # print(graph3)
  # dev.off()
  # quit(status = 1)
  
  #### ChiSq plot3----
  
  Plot3_table$LINEAGE_CLASSIF_DEF<-factor(Plot3_table$LINEAGE_CLASSIF_DEF,
                                          levels=c("Erythroid Traits","Megakaryocitic Traits","GRAN-MONO Traits","Lymphocitic Traits","Pleiotropic_variants"),
                                          ordered=T)
  
  
  cat("Plot3_table_\n")
  cat(str(Plot3_table))
  cat("\n")
  
  Plot3_table_wide<-as.data.frame(pivot_wider(Plot3_table,
                                              id_cols=c("LINEAGE_CLASSIF_DEF","Status","Total_associations_per_Lineage"),
                                              names_from=Cell_Type,
                                              values_from=Associations), stringsAsFactors=F)
  
  Plot3_table_wide<-Plot3_table_wide[order(Plot3_table_wide$LINEAGE_CLASSIF_DEF,Plot3_table_wide$Status),]
  
  Plot3_table_wide[is.na(Plot3_table_wide)]<-0
  
  cat("Plot3_table_wide_0\n")
  cat(str(Plot3_table_wide))
  cat("\n") 
  
  for(x in 1:length(levels(Plot3_table_ACTIVE$Cell_Type)))
  {
    
    Cell_Type_sel<-levels(Plot3_table_ACTIVE$Cell_Type)[x]
    
    
    cat("--->\t")
    cat(sprintf(as.character(Cell_Type_sel)))
    cat("\t")
    
    indx.int<-which(colnames(Plot3_table_wide) == Cell_Type_sel)
    
    
    cat(sprintf(as.character(indx.int)))
    cat("\n")
    
    if(length(indx.int)  ==0)
    {
      
      # quit(status = 1)
      
      indx.end<-dim(Plot3_table_wide)[2]
      
      
      Plot3_table_wide[,indx.end+1] <-0
      
      cat("Plot3_table_wide_0_5\n")
      cat(str(Plot3_table_wide))
      cat("\n") 
      
      colnames(Plot3_table_wide)[indx.end+1]<-Cell_Type_sel
      
      
      
    }
    
    
    
  }
  
  cat("Plot3_table_wide_1\n")
  cat(str(Plot3_table_wide))
  cat("\n") 
  
  # quit(status=1)
  
 
  
  
  Plot3_table_wide$non_Kolf2<-Plot3_table_wide$Total_associations_per_Lineage-Plot3_table_wide$Kolf2
  Plot3_table_wide$non_K562<-Plot3_table_wide$Total_associations_per_Lineage-Plot3_table_wide$K562
  Plot3_table_wide$non_HL60<-Plot3_table_wide$Total_associations_per_Lineage-Plot3_table_wide$HL60
  Plot3_table_wide$non_THP1<-Plot3_table_wide$Total_associations_per_Lineage-Plot3_table_wide$THP1
  
  cat("Plot3_table_wide_1\n")
  cat(str(Plot3_table_wide))
  cat("\n") 
  
  Plot3_table_wide_2<-as.data.frame(pivot_wider(Plot3_table_wide,
                                                id_cols=c("LINEAGE_CLASSIF_DEF"),
                                                names_from=Status,
                                                values_from=c("Total_associations_per_Lineage","K562","Kolf2","HL60","THP1",
                                                              "non_Kolf2","non_K562","non_HL60","non_THP1")), stringsAsFactors=F)
  
  Plot3_table_wide_2[is.na(Plot3_table_wide_2)]<-0
  
  
  cat("Plot3_table_wide_2\n")
  cat(str(Plot3_table_wide_2))
  cat("\n") 
  
  
  
  
  
  Plot3_table_wide_2.m<-melt(Plot3_table_wide_2[,-c(which(colnames(Plot3_table_wide_2) == "Total_associations_per_Lineage_0"),
                                                    which(colnames(Plot3_table_wide_2) == "Total_associations_per_Lineage_1"))],
                             id.variables=c("LINEAGE_CLASSIF_DEF"),
                             variable.name = "variable", 
                             value.name = "n")
  
  
  
  
  cat("Plot3_table_wide_2.m\n")
  cat(str(Plot3_table_wide_2.m))
  cat("\n") 
  
  
  Plot3_table_wide_2.m$Status<-"NA"
  
  Plot3_table_wide_2.m$Status[grep("_0$",Plot3_table_wide_2.m$variable)]<-"0"
  Plot3_table_wide_2.m$Status[grep("_1$",Plot3_table_wide_2.m$variable)]<-"1"
  
  
  Plot3_table_wide_2.m$belonging<-"belonging"
  
  Plot3_table_wide_2.m$belonging[grep("^non_",Plot3_table_wide_2.m$variable)]<-"not_belonging"
  
  Plot3_table_wide_2.m$Cell_Type <- "NA"
  
  Plot3_table_wide_2.m$Cell_Type[grep("K562",Plot3_table_wide_2.m$variable)]<-"K562"
  Plot3_table_wide_2.m$Cell_Type[grep("Kolf2",Plot3_table_wide_2.m$variable)]<-"Kolf2"
  Plot3_table_wide_2.m$Cell_Type[grep("HL60",Plot3_table_wide_2.m$variable)]<-"HL60"
  Plot3_table_wide_2.m$Cell_Type[grep("THP1",Plot3_table_wide_2.m$variable)]<-"THP1"
  
  
  # 
  # 
  cat("Plot3_table_wide_2.m_1\n")
  cat(str(Plot3_table_wide_2.m))
  cat("\n")
  
  Plot3_table_wide_2.m$Cell_Type<-factor(Plot3_table_wide_2.m$Cell_Type,
                                         levels=c("Kolf2","K562","HL60","THP1"),
                                         ordered=T)
  Plot3_table_wide_2.m$Status<-factor(Plot3_table_wide_2.m$Status,
                                      levels=c("0","1"),
                                      ordered=T)
  
  Plot3_table_wide_2.m$belonging<-factor(Plot3_table_wide_2.m$belonging,
                                         levels=c("not_belonging","belonging"),
                                         ordered=T)
  
  
  Plot3_table_wide_2.m<-Plot3_table_wide_2.m[order(Plot3_table_wide_2.m$LINEAGE_CLASSIF_DEF,Plot3_table_wide_2.m$Cell_Type,Plot3_table_wide_2.m$Status,Plot3_table_wide_2.m$belonging),]
  
  cat("Plot3_table_wide_2.m_2\n")
  cat(str(Plot3_table_wide_2.m))
  cat("\n")
  
  
  
  # quit(status = 1)
  

  
  ##### LOOP for ChiSq calculation INTERCELL_ER --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  
  list_ChiSq_DEF<-list()
  
  
  for(i in 1:length(levels(Plot3_table_wide_2.m$LINEAGE_CLASSIF_DEF)))
  {
    
    level_sel<-levels(Plot3_table_wide_2.m$LINEAGE_CLASSIF_DEF)[i]
    
    cat("--->\t")
    cat(sprintf(as.character(level_sel)))
    cat("\n")
    
    
    
    Plot3_table_wide_2.m_sel<-Plot3_table_wide_2.m[which(Plot3_table_wide_2.m$LINEAGE_CLASSIF_DEF == level_sel),]
    
    cat("Plot3_table_wide_2.m_sel\n")
    cat(str(Plot3_table_wide_2.m_sel))
    cat("\n")
    
    if(dim(Plot3_table_wide_2.m_sel)[1] > 0)
    {
      
      array_Cell_Type<-levels(Plot3_table_wide_2.m_sel$Cell_Type)
      
      
      cat("array_Cell_Type\n")
      cat(str(array_Cell_Type))
      cat("\n")
      
      list_ChiSq<-list()
      
      for(k in 1:length(array_Cell_Type))
      {
        Cell_Type_sel<-array_Cell_Type[k]
        
        cat(">\t")
        cat(sprintf(as.character(Cell_Type_sel)))
        cat("\n")
        
        Plot3_table_wide_2.m_sel_Cell_Type_sel<-Plot3_table_wide_2.m_sel[which(Plot3_table_wide_2.m_sel$Cell_Type == Cell_Type_sel),]
        
        cat("Plot3_table_wide_2.m_sel_Cell_Type_sel\n")
        cat(str(Plot3_table_wide_2.m_sel_Cell_Type_sel))
        cat("\n")
        
        
        setwd(out)
        write.table(file="test_2.tsv",Plot3_table_wide_2.m_sel_Cell_Type_sel, sep="\t",quote=F, row.names = F)
        
        
        tab <- xtabs(n ~ Status+belonging, Plot3_table_wide_2.m_sel_Cell_Type_sel)
        
        cat("tab\n")
        cat(str(tab))
        cat("\n")
        
        tab.chisq.test<-chisq.test(tab,correct = TRUE)
        
        cat("tab.chisq.test\n")
        cat(str(tab.chisq.test))
        cat("\n")
        
        log_pval<--1*log10(tab.chisq.test$p.value)
        
        cat("log_pval\n")
        cat(str(log_pval))
        cat("\n")
        
        a.df<-as.data.frame(cbind(Cell_Type_sel,log_pval), stringsAsFactors=F)
        
        colnames(a.df)<-c("Cell_Type","log_pval")
        
        cat("a.df\n")
        cat(str(a.df))
        cat("\n")
        
        
        list_ChiSq[[k]]<-a.df
        
        # quit(status = 1)
        
        
        
        
      }#k
      
      
      partial_tab = unique(as.data.frame(data.table::rbindlist(list_ChiSq, fill = T)))
      
      partial_tab$Traits<-level_sel
      
      cat("partial_tab\n")
      cat(str(partial_tab))
      cat("\n")
      
      
      list_ChiSq_DEF[[i]]<-partial_tab
    }
    
    
   
    
    
    # quit(status = 1)
    
  }
  
  INTERCELL_ER = unique(as.data.frame(data.table::rbindlist(list_ChiSq_DEF, fill = T)))
  
  
  cat("INTERCELL_ER\n")
  cat(str(INTERCELL_ER))
  cat("\n")
  
  # 
  # setwd(out)
  # write.table(file="test_4.tsv",INTERCELL_ER, sep="\t",quote=F, row.names = F)
  # 
  # 
  # 
  # 
  # 
  # quit(status = 1)
  
 
  
  
  ####    SAVE ----
  
  
  # quit(status = 1)
  
  setwd(out)

  filename_1<-paste("Lineages_INTRACELL_ER_",type,".tsv", sep='')
  write.table(file=filename_1,INTRACELL_ER, sep="\t",quote=F, row.names = F)
  
  filename_1<-paste("Lineages_INTERCELL_ER_",type,".tsv", sep='')
  write.table(file=filename_1,INTERCELL_ER, sep="\t",quote=F, row.names = F)

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
    make_option(c("--dB"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--genIE_RESULTS"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--TOME_correspondence"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--finemap_prob_Threshold"), type="numeric", default=NULL, 
                metavar="filename", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--genIE_Tier_1"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
  )
  parser = OptionParser(usage = "137_genIE_normalization_and_filtering_Rscript_v2.R
                        --regular_table FILE.txt
                        --replicas charac
                        --type1 type1
                        --type2 type2
                        --pvalThreshold integer
                        --FDRThreshold integer
                        --EquivalenceTable FILE.txt
                        --sharpr2Threshold charac",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  
  graph_function(opt)
  
}
  
  
  
 

###########################################################################

system.time( main() )
