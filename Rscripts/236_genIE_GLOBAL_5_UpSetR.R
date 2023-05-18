
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
library("reshape2",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")

library("ggeasy",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("viridisLite",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("ggrepel", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")


opt = NULL


opt = NULL

options(warn=1)


ALL_PURPOSE_function = function(option_list)
{
 
  
  packages_loaded<-sessionInfo()
  
  # cat("packages_loaded\n")
  # cat(str(packages_loaded))
  # cat("\n")
  
  
  
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
  
  #### READ INPUT FILES ----
  
  setwd(out)
  
  filename_1<-paste("genIE_results",".rds", sep='')
  
  cat("filename_1\n")
  cat(sprintf(as.character(filename_1)))
  cat("\n")
  
  
  Gather.m<-readRDS(file= filename_1)
  
  cat("Gather.m\n")
  str(Gather.m)
  cat("\n")
  cat(sprintf(as.character(levels(Gather.m$variable))))
  cat("\n")
  cat(sprintf(as.character(levels(as.factor(Gather.m$Rep)))))
  cat("\n")
  
  
  
  #### READ and transform DESIGN_DROPOUTS ----
  
  DESIGN_DROPOUTS = unlist(strsplit(split =",", opt$DESIGN_DROPOUTS))
  
  cat("DESIGN_DROPOUTS_\n")
  cat(sprintf(as.character(DESIGN_DROPOUTS)))
  cat("\n")
  # quit(status=1)
  
  #### READ and transform EDITION_DROPOUTS ----
  
  EDITION_DROPOUTS = unlist(strsplit(split =",", opt$EDITION_DROPOUTS))
  
  cat("EDITION_DROPOUTS_\n")
  cat(sprintf(as.character(EDITION_DROPOUTS)))
  cat("\n")
  # quit(status=1)
  
  #### READ and transform EXPRESSION_DROPOUTS ----
  
  EXPRESSION_DROPOUTS = unlist(strsplit(split =",", opt$EXPRESSION_DROPOUTS))
  
  cat("EXPRESSION_DROPOUTS_\n")
  cat(sprintf(as.character(EXPRESSION_DROPOUTS)))
  cat("\n")
  # quit(status=1)
  
  #### READ and transform ACCEPTED_REPLICATES ----
  
  ACCEPTED_REPLICATES = unlist(strsplit(split =",", opt$ACCEPTED_REPLICATES))
  
  cat("ACCEPTED_REPLICATES_\n")
  cat(sprintf(as.character(ACCEPTED_REPLICATES)))
  cat("\n")
  
  # quit(status=1)
  #### Gather.m accepted replicates ----
  
  Gather.m<-Gather.m[which(Gather.m$Rep%in%ACCEPTED_REPLICATES),]
  
  cat("Gather.m\n")
  str(Gather.m)
  cat("\n")
  
  #which(Gather.m$HGNC%in%EDITION_DROPOUTS),
  
  Gather.m<-Gather.m[-c(which(Gather.m$HGNC%in%DESIGN_DROPOUTS),
                        which(Gather.m$HGNC%in%EXPRESSION_DROPOUTS))
                     ,]
  
  cat("Gather.m\n")
  str(Gather.m)
  cat("\n")
  
  # quit(status=1)
  ### genIE_Rosetta ----
  
  filename_1<-paste("genIE_Rosetta",".rds", sep='')
  
  cat("filename_1\n")
  cat(sprintf(as.character(filename_1)))
  cat("\n")
  
  
  genIE_Rosetta<-readRDS(file= filename_1)
  
  cat("genIE_Rosetta\n")
  str(genIE_Rosetta)
  cat("\n")
  
  #### MERGE ----
  
  indx.int<-c(which(colnames(genIE_Rosetta) == "Cell_Type"),which(colnames(genIE_Rosetta) == "Rep"))
  
  
  genIE_Rosetta_subset<-unique(genIE_Rosetta[,indx.int])
  
  cat("genIE_Rosetta_subset\n")
  str(genIE_Rosetta_subset)
  cat("\n")
  
  
  Gather.m<-merge(Gather.m,
                  genIE_Rosetta_subset,
                  by="Rep",
                  all.x=T)
  
  cat("Gather.m\n")
  str(Gather.m)
  cat("\n")
  
  ##### Turn Inf values into something meaningful-----
  
  Gather.m_infinite<-Gather.m[is.infinite(Gather.m$value),]
  
  cat("Gather.m_infinite\n")
  str(Gather.m_infinite)
  cat("\n")
  
  
  Gather.m_finite<-Gather.m[is.finite(Gather.m$value),]
  
  cat("Gather.m_finite\n")
  str(Gather.m_finite)
  cat("\n")
  cat(sprintf(as.character(names(summary(Gather.m_finite$variable)))))
  cat("\n")
  cat(sprintf(as.character(summary(Gather.m_finite$variable))))
  cat("\n")
  
  #### Cell_Type colors ----
  
  Cell_Type_levels<-levels(Gather.m_finite$Cell_Type)
  colors_Cell_Type_levels<-c('firebrick2','#32A852','#553B68','#D45E85','#6DB2EE','#62D07F','#C9244B','#87447B','#D6B8E6')
  
  df.color_Cell_Type<-as.data.frame(cbind(Cell_Type_levels,colors_Cell_Type_levels[1:length(Cell_Type_levels)]), stringsAsFactors=F)
  colnames(df.color_Cell_Type)<-c("Cell_Type","colors")
  
  cat("df.color_Cell_Type_0\n")
  cat(str(df.color_Cell_Type))
  cat("\n")
  
  # quit(status = 1)
  
  #### READ INPUT FILES ----
  
  setwd(out)
  
  filename_1<-paste("genIE_results_per_metric",".rds", sep='')
  
  cat("filename_1\n")
  cat(sprintf(as.character(filename_1)))
  cat("\n")
  
  
  RESULTS<-readRDS(file= filename_1)
  
  cat("RESULTS\n")
  str(RESULTS)
  cat("\n")
  
  
  # del_logpval_df<-RESULTS$del_logpval
  # 
  # cat("del_logpval_df\n")
  # str(del_logpval_df)
  # cat("\n")
  # 
  # del_effect_df<-RESULTS$del_effect
  # 
  # cat("del_effect_df\n")
  # str(del_effect_df)
  # cat("\n")
  # 
  # Percentage_del_gDNA_df<-RESULTS$Percentage_del_gDNA
  # 
  # cat("Percentage_del_gDNA_df\n")
  # str(Percentage_del_gDNA_df)
  # cat("\n")
  # 
  # hdr_logpval_df<-RESULTS$hdr_logpval
  # 
  # cat("hdr_logpval_df\n")
  # str(hdr_logpval_df)
  # cat("\n")
  # 
  # hdr_effect_df<-RESULTS$hdr_effect
  # 
  # cat("hdr_effect_df\n")
  # str(hdr_effect_df)
  # cat("\n")
  # 
  # Percentage_hdr_gDNA_df<-RESULTS$Percentage_hdr_gDNA
  # 
  # cat("Percentage_hdr_gDNA_df\n")
  # str(Percentage_hdr_gDNA_df)
  # cat("\n")
  
  ##### LOOP volcano ----
  
  list_CellTypes<-list()
  
  
  for(i in 1:length(levels(Gather.m_finite$Cell_Type)))
  {
    Cell_type_sel<-levels(Gather.m_finite$Cell_Type)[i]
    
    cat("Cell_type_sel_\n")
    cat(sprintf(as.character(Cell_type_sel)))
    cat("\n")
    
    Gather_sel<-Gather.m_finite[which(Gather.m_finite$Cell_Type == Cell_type_sel),]
    
    cat("Gather_sel\n")
    str(Gather_sel)
    cat("\n")
    
    
    Gather_sel_wide<-as.data.frame(pivot_wider(Gather_sel,
                                               id_cols=c("Rep","VAR","rsid","HGNC","Cell_Type"),
                                               names_from=variable,
                                               values_from=value), stringsAsFactors=F)
    
    
    cat("Gather_sel_wide\n")
    str(Gather_sel_wide)
    cat("\n")
    
    Gather_sel_wide.dt<-data.table(Gather_sel_wide, key=c("VAR","rsid","HGNC","Cell_Type"))
    
    
    cat("Gather_sel_wide.dt\n")
    str(Gather_sel_wide.dt)
    cat("\n")    
    
    Gather_sel_wide.dt_median<-Gather_sel_wide.dt[,.(median_Percentage_hdr_gDNA=median(Percentage_hdr_gDNA, na.rm = T),
                                                     median_Percentage_del_gDNA=median(Percentage_del_gDNA, na.rm = T),
                                                     median_del_effect=median(del_effect, na.rm = T),
                                                     median_hdr_effect=median(hdr_effect, na.rm = T),
                                                     median_hdr_logpval=median(hdr_logpval, na.rm = T),
                                                     median_del_logpval=median(del_logpval, na.rm = T)),
                                                  .(VAR,rsid,HGNC,Cell_Type)]
    
   
    cat("Gather_sel_wide.dt_median\n")
    str(Gather_sel_wide.dt_median)
    cat("\n")    
    
    
    Gather_DEF<-as.data.frame(Gather_sel_wide.dt_median, stringsAsFactors=F)
    
    colnames(Gather_DEF)[which(colnames(Gather_DEF) == "median_Percentage_hdr_gDNA")]<-"Percentage_hdr_gDNA"
    colnames(Gather_DEF)[which(colnames(Gather_DEF) == "median_Percentage_del_gDNA")]<-"Percentage_del_gDNA"
    colnames(Gather_DEF)[which(colnames(Gather_DEF) == "median_hdr_logpval")]<-"hdr_logpval"
    colnames(Gather_DEF)[which(colnames(Gather_DEF) == "median_del_logpval")]<-"del_logpval"
    colnames(Gather_DEF)[which(colnames(Gather_DEF) == "median_hdr_effect")]<-"hdr_effect"
    colnames(Gather_DEF)[which(colnames(Gather_DEF) == "median_del_effect")]<-"del_effect"
    
    
    cat("Gather_DEF\n")
    str(Gather_DEF)
    cat("\n")
    
    
    ACTIVE_del<-subset(Gather_DEF, Gather_DEF$del_logpval >= 1.3 &
                                    Gather_DEF$Percentage_del_gDNA >= 2)
    
    cat("ACTIVE_del\n")
    str(ACTIVE_del)
    cat("\n")
    
    list_CellTypes$ACTIVE_del[[Cell_type_sel]]<-ACTIVE_del
    
    
    INACTIVE_del<-subset(Gather_DEF, Gather_DEF$del_logpval < 1.3 |
                         Gather_DEF$Percentage_del_gDNA < 2)
    
    cat("INACTIVE_del\n")
    str(INACTIVE_del)
    cat("\n")
    
    list_CellTypes$INACTIVE_del[[Cell_type_sel]]<-INACTIVE_del
    
    ACTIVE_hdr<-subset(Gather_DEF, Gather_DEF$hdr_logpval >= 1.3 &
                         Gather_DEF$Percentage_hdr_gDNA >= 2)
    
    cat("ACTIVE_hdr\n")
    str(ACTIVE_hdr)
    cat("\n")
    
    list_CellTypes$ACTIVE_hdr[[Cell_type_sel]]<-ACTIVE_hdr
    
    INACTIVE_hdr<-subset(Gather_DEF, Gather_DEF$hdr_logpval < 1.3 |
                           Gather_DEF$Percentage_hdr_gDNA < 2)
    
    cat("INACTIVE_hdr\n")
    str(INACTIVE_hdr)
    cat("\n")
    
    list_CellTypes$INACTIVE_hdr[[Cell_type_sel]]<-INACTIVE_hdr
    
    
    # quit(status = 1)
    
  } #i
  
  #### Graph 1 active barplot del ----
  
  
  SIGNIFICANT_del_DEF = unique(as.data.frame(data.table::rbindlist(list_CellTypes$ACTIVE_del, fill = T)))
  
  SIGNIFICANT_del_DEF$TAG_del<-"SIGNIFICANT"
  
  
  
  cat("SIGNIFICANT_del_DEF\n")
  cat(str(SIGNIFICANT_del_DEF))
  cat("\n")
  
  NOT_SIGNIFICANT_del_DEF = unique(as.data.frame(data.table::rbindlist(list_CellTypes$INACTIVE_del, fill = T)))
  
  NOT_SIGNIFICANT_del_DEF$TAG_del<-"NOT_SIGNIFICANT"
  
  
  
  cat("NOT_SIGNIFICANT_del_DEF\n")
  cat(str(NOT_SIGNIFICANT_del_DEF))
  cat("\n")
  
  
  del_DEF<-rbind(SIGNIFICANT_del_DEF,NOT_SIGNIFICANT_del_DEF)
  
  cat("del_DEF\n")
  cat(str(del_DEF))
  cat("\n")
  
  
  del_DEF$TAG_del[which(del_DEF$HGNC%in%DESIGN_DROPOUTS)]<-"DESIGN_DROPOUT"
  del_DEF$TAG_del[which(del_DEF$HGNC%in%EXPRESSION_DROPOUTS &
                          del_DEF$Cell_Type == "Kolf2")]<-"EXPRESSION_DROPOUT"
  del_DEF$TAG_del[which(del_DEF$Percentage_del_gDNA < 2)]<-"EDITION_DROPOUT"
  
  cat("as.factor(del_DEF$TAG_del)\n")
  cat(sprintf(as.character(names(summary(as.factor(interaction(del_DEF$TAG_del,del_DEF$Cell_Type)))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(interaction(del_DEF$TAG_del,del_DEF$Cell_Type))))))
  cat("\n")
  
  
  
  
  
  MASTER_del_by_Cell_Type<-as.data.frame(setDT(del_DEF)[, .N, by=.(Cell_Type,TAG_del)], stringsAsFactors=F)
  
  colnames(MASTER_del_by_Cell_Type)[which(colnames(MASTER_del_by_Cell_Type) == "N")]<-"Observations"
  
  
  cat("MASTER_del_by_Cell_Type_0\n")
  cat(str(MASTER_del_by_Cell_Type))
  cat("\n")
  
  
  MASTER_del_TOTAL_by_Cell_Type<-as.data.frame(setDT(MASTER_del_by_Cell_Type)[,.(Total_observations=sum(Observations)), by=.(Cell_Type)], stringsAsFactors=F)
  
  
  cat("MASTER_del_TOTAL_by_Cell_Type_\n")
  cat(str(MASTER_del_TOTAL_by_Cell_Type))
  cat("\n")
  
  MASTER_del_by_Cell_Type<-as.data.frame(merge(MASTER_del_by_Cell_Type,
                                               MASTER_del_TOTAL_by_Cell_Type,
                                               by="Cell_Type"), stringsAsFactors=F)
  
  
  
  cat("MASTER_del_by_Cell_Type_1\n")
  cat(str(MASTER_del_by_Cell_Type))
  cat("\n")
  
  
  
  MASTER_del_by_Cell_Type$Cell_Type<-factor(MASTER_del_by_Cell_Type$Cell_Type,
                                            levels=c("THP1","HL60","K562","Kolf2"),
                                            ordered=T)
  
  MASTER_del_by_Cell_Type$TAG_del<-factor(MASTER_del_by_Cell_Type$TAG_del,
                                          levels=c("NOT_SIGNIFICANT","EDITION_DROPOUT","EXPRESSION_DROPOUT","DESIGN_DROPOUT","SIGNIFICANT"),
                                          ordered=T)
  
  MASTER_del_by_Cell_Type$TAG_del<-droplevels(MASTER_del_by_Cell_Type$TAG_del)
  
  
  cat("MASTER_del_by_Cell_Type_3\n")
  cat(str(MASTER_del_by_Cell_Type))
  cat("\n")
  
  df.color_Cell_Type_del<-df.color_Cell_Type[order(factor(df.color_Cell_Type$Cell_Type,
                                                          levels=levels(MASTER_del_by_Cell_Type$Cell_Type),
                                                          ordered=T)),]
  
  
  
 
  
  
  
  graph_del<-MASTER_del_by_Cell_Type %>%
    mutate(myaxis = paste0(Cell_Type, "\n", "n=", Total_observations)) %>%
    mutate(myaxis=fct_reorder(myaxis,as.numeric(Cell_Type))) %>%
    ggplot(aes(x=myaxis, y=Observations, fill=TAG_del)) +
    geom_bar(position="fill",stat="identity",colour='black')+
    theme_bw()+
    scale_y_continuous(name=NULL,labels = scales::percent_format())+
    scale_x_discrete(name=NULL, drop=F)+
    theme(axis.title.y=element_text(size=24, family="sans"),
          axis.title.x=element_text(size=24, family="sans"),
          axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
          axis.text.x=element_text(angle=0,size=18, color="black", family="sans"),
          legend.text=element_text(size=12,color="black", family="sans"))+
    scale_fill_manual(values=c("white",'gray','#32A852','#6DB2EE','#62D07F','#1877C9','#C9244B','#D45E85'),drop=F)+
    theme(legend.position="bottom")+
    ggeasy::easy_center_title()
  
  graph_del<-graph_del+guides(fill=guide_legend(nrow=2,byrow=TRUE,title=paste("genIE","UDPs", sep="\n")))

  graph_del<-graph_del+coord_flip()+theme(axis.text.y=element_text(angle=0,size=18, color=df.color_Cell_Type_del$color, family="sans"))
  
  
  setwd(out)
  
  svglite(paste('del_SUMMARY','.svg',sep=''), width = 8, height = 8)
  print(graph_del)
  dev.off()
  
  # quit(status = 1)
  
 
  #### dotplot R del ----
  
  cat("del_DEF\n")
  cat(str(del_DEF))
  cat("\n")
  
  
  del_DEF<-del_DEF[order(del_DEF$Cell_Type),]
  
  
  del_DEF.dt<-data.table(del_DEF, key=c("HGNC","VAR"))
  
  del_DEF.summarised<-as.data.frame(del_DEF.dt[,.(ASSAY_COMPLETION=paste(Cell_Type,collapse="|")), by=.(HGNC,VAR)],stringsAsFactors=F)
  
  cat("del_DEF.summarised\n")
  cat(str(del_DEF.summarised))
  cat("\n")
  
  cat(sprintf(as.character((names(summary(as.factor(del_DEF.summarised$ASSAY_COMPLETION)))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(del_DEF.summarised$ASSAY_COMPLETION)))))
  cat("\n")
  # Subset_ALL_del<-
  
  del_DEF.summarised$ASSAY_COMPLETION<-factor(del_DEF.summarised$ASSAY_COMPLETION,
                                              levels = c("Kolf2|K562|HL60|THP1",
                                                         "Kolf2|K562|HL60","Kolf2|K562|THP1","Kolf2|HL60|THP1","K562|HL60|THP1",
                                                         "Kolf2|K562","Kolf2|HL60","Kolf2|THP1",
                                                         "K562|THP1","K562|HL60",
                                                         "HL60|THP1",
                                                         "Kolf2","K562","HL60","THP1"),
                                              ordered=T)
  
  del_DEF.summarised$ASSAY_COMPLETION<-droplevels(del_DEF.summarised$ASSAY_COMPLETION)
  
  cat(sprintf(as.character((names(summary(del_DEF.summarised$ASSAY_COMPLETION))))))
  cat("\n")
  cat(sprintf(as.character(summary(del_DEF.summarised$ASSAY_COMPLETION))))
  cat("\n")
  
  
  del_DEF<-merge(del_DEF,
                 del_DEF.summarised,
                 by=c("HGNC","VAR"),
                 all=T)
  
  
  cat("del_DEF_2\n")
  cat(str(del_DEF))
  cat("\n")
  
  ### subset del_DEF
  
  del_DEF_subset<-del_DEF[which(del_DEF$TAG_del != "EDITION_DROPOUT"),]
  
  del_DEF_subset<-del_DEF[which(del_DEF$TAG_del == "SIGNIFICANT"),]
  
  
  cat("del_DEF_subset_2\n")
  cat(str(del_DEF_subset))
  cat("\n")
  
  del_DEF_subset$TAG_del<-factor(del_DEF_subset$TAG_del,
                                 levels=c("NOT_SIGNIFICANT","EDITION_DROPOUT","EXPRESSION_DROPOUT","DESIGN_DROPOUT","SIGNIFICANT"),
                                 ordered=T)
  
  del_DEF_subset$TAG_del<-droplevels(del_DEF_subset$TAG_del)
  
  
  cat("del_DEF_subset_3\n")
  cat(str(del_DEF_subset))
  cat("\n")
  
  df.color_Cell_Type_del<-df.color_Cell_Type[order(factor(df.color_Cell_Type$Cell_Type,
                                                          levels=levels(del_DEF_subset$Cell_Type),
                                                          ordered=T)),]
  
  cat("df.color_Cell_Type_del_3\n")
  cat(str(df.color_Cell_Type_del))
  cat("\n")
  
  
  #### X and Y points
  
  ### del_effect LOG
  
  del_effect_vector<-del_DEF_subset$del_effect#[which(Gather_sel$variable == "del_effect")]
  
  cat(sprintf(as.character(summary(del_effect_vector))))
  cat("\n")
  
  
  breaks.del_effect<-as.numeric(summary(del_effect_vector[!is.na(del_effect_vector)]))
  
  cat("breaks.del_effect\n")
  cat(str(breaks.del_effect))
  cat("\n")
  
  breaks.del_effect[1]<-breaks.del_effect[1]-0.2
  
  
  labels.del_effect<-unique(as.character(sort(c(seq(breaks.del_effect[1],breaks.del_effect[length(breaks.del_effect)-1], by=0.2),
                                                breaks.del_effect[length(breaks.del_effect)]))))
  breaks.del_effect.log<-log10(as.numeric(labels.del_effect)+0.00001)
  
  labels.del_effect<-as.character(round(as.numeric(labels.del_effect),2))
  
  cat("labels.del_effect\n")
  cat(sprintf(as.character(labels.del_effect)))
  cat("\n")
  
  cat("breaks.del_effect.log\n")
  cat(sprintf(as.character(breaks.del_effect.log)))
  cat("\n")
  
  
  cat("graph\n")
  
  pos <- position_jitter(width=0.25, height=0.08, seed = 2)
  
  dotplot_del_effect<-ggplot(data=del_DEF_subset,
                             aes(x=ASSAY_COMPLETION, 
                                 y=log10(del_effect+0.00001),
                                 color=Cell_Type,
                                 label=HGNC)) +
    geom_jitter(size=5, position=pos)+
    scale_x_discrete(name=NULL)+
    scale_y_continuous(name="Del effect", 
                       breaks=breaks.del_effect.log,labels=labels.del_effect, 
                       limits=c(breaks.del_effect.log[1],breaks.del_effect.log[length(breaks.del_effect.log)]))+
    theme_bw()+
    theme(axis.title.y=element_text(size=24, family="sans"),
          axis.title.x=element_text(size=24, family="sans"),
          axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
          axis.text.x=element_text(angle=45,size=6, vjust=1,hjust=1,color="black", family="sans"))+
    scale_color_manual(values=df.color_Cell_Type$color, drop=F)+
    theme(legend.position="hidden",legend.title=element_blank(), legend.text = element_text(size=10))+
    geom_hline(yintercept=log10(1+0.00001), col="black")+
    geom_text_repel(box.padding = 0.8, 
                    max.overlaps = Inf,
                    segment.linetype = 6,
                    position=pos)+
    ggeasy::easy_center_title()
  
  dotplot_del_effect<-dotplot_del_effect+guides(color=guide_legend(nrow=1,byrow=TRUE))
  
  
  #### COWPLOT into single graphs del ----
  
  ### del
  
  df.color_Cell_Type_del<-df.color_Cell_Type[order(factor(df.color_Cell_Type$Cell_Type,
                                                          levels=levels(MASTER_del_by_Cell_Type$Cell_Type),
                                                          ordered=T)),]
  
  subgraph_del<-MASTER_del_by_Cell_Type %>%
    mutate(myaxis = paste0(Cell_Type, "\n", "n=", Total_observations)) %>%
    mutate(myaxis=fct_reorder(myaxis,as.numeric(Cell_Type))) %>%
    ggplot(aes(y=myaxis, x=Observations, fill=TAG_del)) +
    geom_bar(position="fill",stat="identity",colour='black')+
    theme_bw()+
    scale_x_reverse(name=NULL,labels = scales::percent_format())+
    scale_y_discrete(position="left",name=NULL, drop=F)+
    theme(axis.title.y=element_text(size=18, family="sans"),
          axis.title.x=element_text(size=18, family="sans"),
          axis.text.y=element_text(angle=0,size=12, color=df.color_Cell_Type_del$color, family="sans"),
          axis.text.x=element_text(angle=0,size=12, color="black", family="sans"),
          legend.text=element_text(size=12,color="black", family="sans"))+
    scale_fill_manual(values=c("white",'gray','#32A852','#6DB2EE','#62D07F','#1877C9','#C9244B','#D45E85'),drop=F)+
    theme(legend.position="bottom")+
    ggeasy::easy_center_title()
  
  subgraph_del<-subgraph_del+guides(fill=guide_legend(nrow=2,byrow=TRUE,title=paste("genIE","UDPs", sep="\n")))
  
  graph_del_FINAL<-plot_grid(NULL,dotplot_del_effect,subgraph_del,NULL,
                             nrow = 2,
                             ncol=2,
                             rel_heights = c(1, 0.5),
                             rel_widths=c(0.55,1))
  
  
  setwd(out)
  
  svglite(paste('del_SIGNIFICANT_RESULTS_UpSetR_like','.svg',sep=''), width = 8, height = 8)
  print(graph_del_FINAL)
  dev.off()
  
  
  
  #### Graph 1 active barplot hdr ----
  
  
  SIGNIFICANT_hdr_DEF = unique(as.data.frame(data.table::rbindlist(list_CellTypes$ACTIVE_hdr, fill = T)))
  
  SIGNIFICANT_hdr_DEF$TAG_hdr<-"SIGNIFICANT"
  
  
  
  cat("SIGNIFICANT_hdr_DEF\n")
  cat(str(SIGNIFICANT_hdr_DEF))
  cat("\n")
  
  NOT_SIGNIFICANT_hdr_DEF = unique(as.data.frame(data.table::rbindlist(list_CellTypes$INACTIVE_hdr, fill = T)))
  
  NOT_SIGNIFICANT_hdr_DEF$TAG_hdr<-"NOT_SIGNIFICANT"
  
  
  
  cat("NOT_SIGNIFICANT_hdr_DEF\n")
  cat(str(NOT_SIGNIFICANT_hdr_DEF))
  cat("\n")
  
  
  hdr_DEF<-rbind(SIGNIFICANT_hdr_DEF,NOT_SIGNIFICANT_hdr_DEF)
  
  cat("hdr_DEF\n")
  cat(str(hdr_DEF))
  cat("\n")
  
  
  hdr_DEF$TAG_hdr[which(hdr_DEF$HGNC%in%DESIGN_DROPOUTS)]<-"DESIGN_DROPOUT"
  hdr_DEF$TAG_hdr[which(hdr_DEF$HGNC%in%EXPRESSION_DROPOUTS &
                          hdr_DEF$Cell_Type == "Kolf2")]<-"EXPRESSION_DROPOUT"
  hdr_DEF$TAG_hdr[which(hdr_DEF$Percentage_hdr_gDNA < 2)]<-"EDITION_DROPOUT"
  
  cat("as.factor(hdr_DEF$TAG_hdr)\n")
  cat(sprintf(as.character(names(summary(as.factor(interaction(hdr_DEF$TAG_hdr,hdr_DEF$Cell_Type)))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(interaction(hdr_DEF$TAG_hdr,hdr_DEF$Cell_Type))))))
  cat("\n")
  
  
  
  
  
  MASTER_hdr_by_Cell_Type<-as.data.frame(setDT(hdr_DEF)[, .N, by=.(Cell_Type,TAG_hdr)], stringsAsFactors=F)
  
  colnames(MASTER_hdr_by_Cell_Type)[which(colnames(MASTER_hdr_by_Cell_Type) == "N")]<-"Observations"
  
  
  cat("MASTER_hdr_by_Cell_Type_0\n")
  cat(str(MASTER_hdr_by_Cell_Type))
  cat("\n")
  
  
  MASTER_hdr_TOTAL_by_Cell_Type<-as.data.frame(setDT(MASTER_hdr_by_Cell_Type)[,.(Total_observations=sum(Observations)), by=.(Cell_Type)], stringsAsFactors=F)
  
  
  cat("MASTER_hdr_TOTAL_by_Cell_Type_\n")
  cat(str(MASTER_hdr_TOTAL_by_Cell_Type))
  cat("\n")
  
  MASTER_hdr_by_Cell_Type<-as.data.frame(merge(MASTER_hdr_by_Cell_Type,
                                               MASTER_hdr_TOTAL_by_Cell_Type,
                                               by="Cell_Type"), stringsAsFactors=F)
  
  
  
  cat("MASTER_hdr_by_Cell_Type_1\n")
  cat(str(MASTER_hdr_by_Cell_Type))
  cat("\n")
  
  
  
  MASTER_hdr_by_Cell_Type$Cell_Type<-factor(MASTER_hdr_by_Cell_Type$Cell_Type,
                                            levels=c("THP1","HL60","K562","Kolf2"),
                                            ordered=T)
  
  MASTER_hdr_by_Cell_Type$TAG_hdr<-factor(MASTER_hdr_by_Cell_Type$TAG_hdr,
                                          levels=c("NOT_SIGNIFICANT","EDITION_DROPOUT","EXPRESSION_DROPOUT","DESIGN_DROPOUT","SIGNIFICANT"),
                                          ordered=T)
  
  MASTER_hdr_by_Cell_Type$TAG_hdr<-droplevels(MASTER_hdr_by_Cell_Type$TAG_hdr)
  
  
  cat("MASTER_hdr_by_Cell_Type_3\n")
  cat(str(MASTER_hdr_by_Cell_Type))
  cat("\n")
  
  df.color_Cell_Type_hdr<-df.color_Cell_Type[order(factor(df.color_Cell_Type$Cell_Type,
                                                          levels=levels(MASTER_hdr_by_Cell_Type$Cell_Type),
                                                          ordered=T)),]
  
  
  
  
  
  
  
  graph_hdr<-MASTER_hdr_by_Cell_Type %>%
    mutate(myaxis = paste0(Cell_Type, "\n", "n=", Total_observations)) %>%
    mutate(myaxis=fct_reorder(myaxis,as.numeric(Cell_Type))) %>%
    ggplot(aes(x=myaxis, y=Observations, fill=TAG_hdr)) +
    geom_bar(position="fill",stat="identity",colour='black')+
    theme_bw()+
    scale_y_continuous(name=NULL,labels = scales::percent_format())+
    scale_x_discrete(name=NULL, drop=F)+
    theme(axis.title.y=element_text(size=24, family="sans"),
          axis.title.x=element_text(size=24, family="sans"),
          axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
          axis.text.x=element_text(angle=0,size=18, color="black", family="sans"),
          legend.text=element_text(size=12,color="black", family="sans"))+
    scale_fill_manual(values=c("white",'gray','#32A852','#6DB2EE','#62D07F','#1877C9','#C9244B','#D45E85'),drop=F)+
    theme(legend.position="bottom")+
    ggeasy::easy_center_title()
  
  graph_hdr<-graph_hdr+guides(fill=guide_legend(nrow=2,byrow=TRUE,title=paste("genIE","UDPs", sep="\n")))
  
  graph_hdr<-graph_hdr+coord_flip()+theme(axis.text.y=element_text(angle=0,size=18, color=df.color_Cell_Type_hdr$color, family="sans"))
  
  
  setwd(out)
  
  svglite(paste('hdr_SUMMARY','.svg',sep=''), width = 8, height = 8)
  print(graph_hdr)
  dev.off()
  
  # quit(status = 1)
  
  
  #### dotplot R hdr ----
  
  cat("hdr_DEF\n")
  cat(str(hdr_DEF))
  cat("\n")
  
  
  hdr_DEF<-hdr_DEF[order(hdr_DEF$Cell_Type),]
  
  
  hdr_DEF.dt<-data.table(hdr_DEF, key=c("HGNC","VAR"))
  
  hdr_DEF.summarised<-as.data.frame(hdr_DEF.dt[,.(ASSAY_COMPLETION=paste(Cell_Type,collapse="|")), by=.(HGNC,VAR)],stringsAsFactors=F)
  
  cat("hdr_DEF.summarised\n")
  cat(str(hdr_DEF.summarised))
  cat("\n")
  
  cat(sprintf(as.character((names(summary(as.factor(hdr_DEF.summarised$ASSAY_COMPLETION)))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(hdr_DEF.summarised$ASSAY_COMPLETION)))))
  cat("\n")
  # Subset_ALL_hdr<-
  
  hdr_DEF.summarised$ASSAY_COMPLETION<-factor(hdr_DEF.summarised$ASSAY_COMPLETION,
                                              levels = c("Kolf2|K562|HL60|THP1",
                                                         "Kolf2|K562|HL60","Kolf2|K562|THP1","Kolf2|HL60|THP1","K562|HL60|THP1",
                                                         "Kolf2|K562","Kolf2|HL60","Kolf2|THP1",
                                                         "K562|THP1","K562|HL60",
                                                         "HL60|THP1",
                                                         "Kolf2","K562","HL60","THP1"),
                                              ordered=T)
  
  hdr_DEF.summarised$ASSAY_COMPLETION<-droplevels(hdr_DEF.summarised$ASSAY_COMPLETION)
  
  cat(sprintf(as.character((names(summary(hdr_DEF.summarised$ASSAY_COMPLETION))))))
  cat("\n")
  cat(sprintf(as.character(summary(hdr_DEF.summarised$ASSAY_COMPLETION))))
  cat("\n")
  
  
  hdr_DEF<-merge(hdr_DEF,
                 hdr_DEF.summarised,
                 by=c("HGNC","VAR"),
                 all=T)
  
  
  cat("hdr_DEF_2\n")
  cat(str(hdr_DEF))
  cat("\n")
  
  ### subset hdr_DEF
  
  hdr_DEF_subset<-hdr_DEF[which(hdr_DEF$TAG_hdr != "EDITION_DROPOUT"),]
  
  hdr_DEF_subset<-hdr_DEF[which(hdr_DEF$TAG_hdr == "SIGNIFICANT"),]
  
  
  cat("hdr_DEF_subset_2\n")
  cat(str(hdr_DEF_subset))
  cat("\n")
  
  hdr_DEF_subset$TAG_hdr<-factor(hdr_DEF_subset$TAG_hdr,
                                 levels=c("NOT_SIGNIFICANT","EDITION_DROPOUT","EXPRESSION_DROPOUT","DESIGN_DROPOUT","SIGNIFICANT"),
                                 ordered=T)
  
  hdr_DEF_subset$TAG_hdr<-droplevels(hdr_DEF_subset$TAG_hdr)
  
  
  cat("hdr_DEF_subset_3\n")
  cat(str(hdr_DEF_subset))
  cat("\n")
  
  df.color_Cell_Type_hdr<-df.color_Cell_Type[order(factor(df.color_Cell_Type$Cell_Type,
                                                          levels=levels(hdr_DEF_subset$Cell_Type),
                                                          ordered=T)),]
  
  cat("df.color_Cell_Type_hdr_3\n")
  cat(str(df.color_Cell_Type_hdr))
  cat("\n")
  
  
  #### X and Y points
  
  ### hdr_effect LOG
  
  hdr_effect_vector<-hdr_DEF_subset$hdr_effect#[which(Gather_sel$variable == "hdr_effect")]
  
  cat(sprintf(as.character(summary(hdr_effect_vector))))
  cat("\n")
  
  
  breaks.hdr_effect<-as.numeric(summary(hdr_effect_vector[!is.na(hdr_effect_vector)]))
  
  cat("breaks.hdr_effect\n")
  cat(str(breaks.hdr_effect))
  cat("\n")
  
  breaks.hdr_effect[1]<-breaks.hdr_effect[1]-0.2
  
  
  labels.hdr_effect<-unique(as.character(sort(c(seq(breaks.hdr_effect[1],breaks.hdr_effect[length(breaks.hdr_effect)-1], by=0.2),
                                                breaks.hdr_effect[length(breaks.hdr_effect)]))))
  breaks.hdr_effect.log<-log10(as.numeric(labels.hdr_effect)+0.00001)
  
  labels.hdr_effect<-as.character(round(as.numeric(labels.hdr_effect),2))
  
  cat("labels.hdr_effect\n")
  cat(sprintf(as.character(labels.hdr_effect)))
  cat("\n")
  
  cat("breaks.hdr_effect.log\n")
  cat(sprintf(as.character(breaks.hdr_effect.log)))
  cat("\n")
  
  
  cat("graph\n")
  
  pos <- position_jitter(width=0.25, height=0.08, seed = 2)
  
  dotplot_hdr_effect<-ggplot(data=hdr_DEF_subset,
                             aes(x=ASSAY_COMPLETION, 
                                 y=log10(hdr_effect+0.00001),
                                 color=Cell_Type,
                                 label=HGNC)) +
    geom_jitter(size=5, position=pos)+
    scale_x_discrete(name=NULL)+
    scale_y_continuous(name="HDR effect", 
                       breaks=breaks.hdr_effect.log,labels=labels.hdr_effect, 
                       limits=c(breaks.hdr_effect.log[1],breaks.hdr_effect.log[length(breaks.hdr_effect.log)]))+
    theme_bw()+
    theme(axis.title.y=element_text(size=24, family="sans"),
          axis.title.x=element_text(size=24, family="sans"),
          axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
          axis.text.x=element_text(angle=45,size=6, vjust=1,hjust=1,color="black", family="sans"))+
    scale_color_manual(values=df.color_Cell_Type$color, drop=F)+
    theme(legend.position="hidden",legend.title=element_blank(), legend.text = element_text(size=10))+
    geom_hline(yintercept=log10(1+0.00001), col="black")+
    geom_text_repel(box.padding = 0.8, 
                    max.overlaps = Inf,
                    segment.linetype = 6,
                    position=pos)+
    ggeasy::easy_center_title()
  
  dotplot_hdr_effect<-dotplot_hdr_effect+guides(color=guide_legend(nrow=1,byrow=TRUE))
  
  
  #### COWPLOT into single graphs hdr ----
  
  ### hdr
  
  df.color_Cell_Type_hdr<-df.color_Cell_Type[order(factor(df.color_Cell_Type$Cell_Type,
                                                          levels=levels(MASTER_hdr_by_Cell_Type$Cell_Type),
                                                          ordered=T)),]
  
  subgraph_hdr<-MASTER_hdr_by_Cell_Type %>%
    mutate(myaxis = paste0(Cell_Type, "\n", "n=", Total_observations)) %>%
    mutate(myaxis=fct_reorder(myaxis,as.numeric(Cell_Type))) %>%
    ggplot(aes(y=myaxis, x=Observations, fill=TAG_hdr)) +
    geom_bar(position="fill",stat="identity",colour='black')+
    theme_bw()+
    scale_x_reverse(name=NULL,labels = scales::percent_format())+
    scale_y_discrete(position="left",name=NULL, drop=F)+
    theme(axis.title.y=element_text(size=18, family="sans"),
          axis.title.x=element_text(size=18, family="sans"),
          axis.text.y=element_text(angle=0,size=12, color=df.color_Cell_Type_hdr$color, family="sans"),
          axis.text.x=element_text(angle=0,size=12, color="black", family="sans"),
          legend.text=element_text(size=12,color="black", family="sans"))+
    scale_fill_manual(values=c("white",'gray','#32A852','#6DB2EE','#62D07F','#1877C9','#C9244B','#D45E85'),drop=F)+
    theme(legend.position="bottom")+
    ggeasy::easy_center_title()
  
  subgraph_hdr<-subgraph_hdr+guides(fill=guide_legend(nrow=2,byrow=TRUE,title=paste("genIE","UDPs", sep="\n")))
  
  graph_hdr_FINAL<-plot_grid(NULL,dotplot_hdr_effect,subgraph_hdr,NULL,
                             nrow = 2,
                             ncol=2,
                             rel_heights = c(1, 0.5),
                             rel_widths=c(0.55,1))
  
  
  setwd(out)
  
  svglite(paste('hdr_SIGNIFICANT_RESULTS_UpSetR_like','.svg',sep=''), width = 8, height = 8)
  print(graph_hdr_FINAL)
  dev.off()
  
  
  
  #### SAVE THE TIERS & a dropout table---- 

  colnames(del_DEF)[which(colnames(del_DEF) == "TAG_del")]<-"STATUS_del"
  
  del_DEF<-as.data.frame(del_DEF, stringsAsFactors=F)
  
  cat("del_DEF_3\n")
  cat(str(del_DEF))
  cat("\n")
  
  colnames(hdr_DEF)[which(colnames(hdr_DEF) == "TAG_hdr")]<-"STATUS_hdr"
  
  hdr_DEF<-as.data.frame(hdr_DEF, stringsAsFactors=F)
  
  
  cat("hdr_DEF_3\n")
  cat(str(hdr_DEF))
  cat("\n")
  
  
  DEF_TABLE<-merge(del_DEF,
                   hdr_DEF,
                   by=colnames(del_DEF)[which(colnames(del_DEF)%in%colnames(hdr_DEF))],
                   all=T)
  
  cat("DEF_TABLE_3\n")
  cat(str(DEF_TABLE))
  cat("\n")
  
  
  DEF_TABLE_DESIGN_DROPOUTS<-DEF_TABLE[which(DEF_TABLE$HGNC%in%DESIGN_DROPOUTS),]
  
  cat("DEF_TABLE_DESIGN_DROPOUTS\n")
  cat(str(DEF_TABLE_DESIGN_DROPOUTS))
  cat("\n")
  
  if(dim(DEF_TABLE_DESIGN_DROPOUTS)[1] ==0)
  {
    
    temp<- data.frame(matrix(vector(), length(DESIGN_DROPOUTS), 
                             dim(DEF_TABLE_DESIGN_DROPOUTS)[2],
                             dimnames=list(c(),
                                           colnames(DEF_TABLE_DESIGN_DROPOUTS)
                             )),
                      stringsAsFactors=F)
    
    # colnames(temp)<-colnames(DEF_TABLE_DESIGN_DROPOUTS)
    
    cat("temp_0\n")
    cat(str(temp))
    cat("\n")
    
    temp$HGNC<-DESIGN_DROPOUTS
    temp$STATUS_del<-"DESIGN_DROPOUT"
    temp$STATUS_hdr<-"DESIGN_DROPOUT"
    temp$Cell_Type<-"Kolf2"
    
    cat("temp_1\n")
    cat(str(temp))
    cat("\n")
    
    DEF_TABLE<-rbind(DEF_TABLE,temp)
    
  }
  
 
  
  
  DEF_TABLE_EXPRESSION_DROPOUTS<-DEF_TABLE[which(DEF_TABLE$HGNC%in%EXPRESSION_DROPOUTS),]
  
  cat("DEF_TABLE_EXPRESSION_DROPOUTS\n")
  cat(str(DEF_TABLE_EXPRESSION_DROPOUTS))
  cat("\n")
  
  if(dim(DEF_TABLE_EXPRESSION_DROPOUTS)[1] ==0)
  {
    
  
    
    temp<- data.frame(matrix(vector(), length(EXPRESSION_DROPOUTS), 
                                        dim(DEF_TABLE_EXPRESSION_DROPOUTS)[2],
                                 dimnames=list(c(),
                                               colnames(DEF_TABLE_EXPRESSION_DROPOUTS)
                                               )),
                          stringsAsFactors=F)
    
   # colnames(temp)<-colnames(DEF_TABLE_EXPRESSION_DROPOUTS)
    
    cat("temp_0\n")
    cat(str(temp))
    cat("\n")
    
    temp$HGNC<-EXPRESSION_DROPOUTS
    temp$STATUS_del<-"EXPRESSION_DROPOUT"
    temp$STATUS_hdr<-"EXPRESSION_DROPOUT"
    temp$Cell_Type<-"Kolf2"
    
    cat("temp_1\n")
    cat(str(temp))
    cat("\n")
    
    DEF_TABLE<-rbind(DEF_TABLE,temp)
  }
  
  cat("DEF_TABLE\n")
  cat(str(DEF_TABLE))
  cat("\n")
  
  Tier_del<-DEF_TABLE[which(DEF_TABLE$STATUS_del == "SIGNIFICANT"),]
  
  cat("Tier_del\n")
  cat(str(Tier_del))
  cat("\n")
  
  Tier_hdr<-DEF_TABLE[which(DEF_TABLE$STATUS_hdr == "SIGNIFICANT"),]
  
  cat("Tier_hdr\n")
  cat(str(Tier_hdr))
  cat("\n")
  
  
  Tier_del_and_hdr<-DEF_TABLE[which(DEF_TABLE$STATUS_del == "SIGNIFICANT" &
                                      DEF_TABLE$STATUS_hdr == "SIGNIFICANT"),]
  
  cat("Tier_del_and_hdr\n")
  cat(str(Tier_del_and_hdr))
  cat("\n")
  
  
  #### SAVE THE ACTIVE TILES ---- 
  
 
  setwd(out)
  
  write.table(file="genIE_RESULTS_TABLE.tsv",
              DEF_TABLE,
              sep="\t",
              quote = F,
              row.names = F)
  
  saveRDS(file="genIE_RESULTS_TABLE.rds",DEF_TABLE)
  
  
  
  write.table(file="genIE_Tier_del_and_hdr.tsv",
              Tier_del_and_hdr,
              sep="\t",
              quote = F,
              row.names = F)
  
  saveRDS(Tier_del_and_hdr,file="genIE_Tier_del_and_hdr.rds")
  
  write.table(file="genIE_Tier_del.tsv",
              Tier_del,
              sep="\t",
              quote = F,
              row.names = F)
  
  saveRDS(Tier_del,file="genIE_Tier_del.rds")
  
  write.table(file="genIE_Tier_hdr.tsv",
              Tier_hdr,
              sep="\t",
              quote = F,
              row.names = F)
  
  saveRDS(Tier_hdr,file="genIE_Tier_hdr.rds")

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
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="out", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
    
  )
  parser = OptionParser(usage = "137_MPRA_normalization_and_filtering_Rscript_v2.R
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
  
  
  ALL_PURPOSE_function(opt)
  
}
  
  
  
 

###########################################################################

system.time( main() )
