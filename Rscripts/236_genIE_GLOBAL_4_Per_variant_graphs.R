

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



opt = NULL

options(warn=1)



per_variant_graphs = function(option_list)
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
  
  filename_1<-paste("genIE_results_per_metric",".rds", sep='')
  
  cat("filename_1\n")
  cat(sprintf(as.character(filename_1)))
  cat("\n")
  
  
  RESULTS<-readRDS(file= filename_1)
  
  cat("RESULTS\n")
  str(RESULTS)
  cat("\n")
  
  
  del_logpval_df<-RESULTS$del_logpval
  
  cat("del_logpval_df\n")
  str(del_logpval_df)
  cat("\n")
  
  del_effect_df<-RESULTS$del_effect
  
  cat("del_effect_df\n")
  str(del_effect_df)
  cat("\n")
  
  Percentage_del_gDNA_df<-RESULTS$Percentage_del_gDNA
  
  cat("Percentage_del_gDNA_df\n")
  str(Percentage_del_gDNA_df)
  cat("\n")
  
  hdr_logpval_df<-RESULTS$hdr_logpval
  
  cat("hdr_logpval_df\n")
  str(hdr_logpval_df)
  cat("\n")
  
  hdr_effect_df<-RESULTS$hdr_effect
  
  cat("hdr_effect_df\n")
  str(hdr_effect_df)
  cat("\n")
  
  Percentage_hdr_gDNA_df<-RESULTS$Percentage_hdr_gDNA
  
  cat("Percentage_hdr_gDNA_df\n")
  str(Percentage_hdr_gDNA_df)
  cat("\n")
  
  #### Cell_Type colors ----
  
  Cell_Type_levels<-levels(Percentage_del_gDNA_df$Cell_Type)
  colors_Cell_Type_levels<-c('#ff1493','#32A852','#553B68','#D45E85','#6DB2EE','#62D07F','#C9244B','#87447B','#D6B8E6')
  
  df.color_Cell_Type<-as.data.frame(cbind(Cell_Type_levels,colors_Cell_Type_levels[1:length(Cell_Type_levels)]), stringsAsFactors=F)
  colnames(df.color_Cell_Type)<-c("Cell_Type","colors")
  
  cat("df.color_Cell_Type_0\n")
  cat(str(df.color_Cell_Type))
  cat("\n")
  
  
  Percentage_del_gDNA_df_subset<-unique(Percentage_del_gDNA_df[,c(which(colnames(Percentage_del_gDNA_df) == "Cell_Type"),
                                                           which(colnames(Percentage_del_gDNA_df) == "Rep"))])
  
  cat("Percentage_del_gDNA_df_subset_\n")
  cat(str(Percentage_del_gDNA_df_subset))
  cat("\n")
  
  df.color_Cell_Type<-merge(df.color_Cell_Type,
                            Percentage_del_gDNA_df_subset,
                            by="Cell_Type")
  
  
  cat("df.color_Cell_Type_1\n")
  cat(str(df.color_Cell_Type))
  cat("\n")
  # quit(status = 1)
  
  #### X and Y points ----
  
  
  breaks.del_logpval<-as.numeric(summary(del_logpval_df$value[!is.na(del_logpval_df$value)]))
  breaks.del_logpval[1]<-0
  breaks.del_logpval[length(breaks.del_logpval)]<-breaks.del_logpval[length(breaks.del_logpval)]+2
  labels.del_logpval<-as.character(seq(breaks.del_logpval[1],breaks.del_logpval[length(breaks.del_logpval)], by=2))
  breaks.del_logpval<-as.numeric(labels.del_logpval)
  
  cat("labels.del_logpval\n")
  cat(str(labels.del_logpval))
  cat("\n")
  
  breaks.hdr_logpval<-as.numeric(summary(hdr_logpval_df$value[!is.na(hdr_logpval_df$value)]))
  breaks.hdr_logpval[1]<-0
  breaks.hdr_logpval[length(breaks.hdr_logpval)]<-breaks.hdr_logpval[length(breaks.hdr_logpval)]+2
  labels.hdr_logpval<-as.character(seq(breaks.hdr_logpval[1],breaks.hdr_logpval[length(breaks.hdr_logpval)], by=2))
  breaks.hdr_logpval<-as.numeric(labels.hdr_logpval)
  
  cat("labels.hdr_logpval\n")
  cat(str(labels.hdr_logpval))
  cat("\n")
  
  
  
  breaks.del_effect<-as.numeric(summary(del_effect_df$value[!is.na(del_effect_df$value)]))
  breaks.del_effect<-unique(c(breaks.del_effect[1],breaks.del_effect[1]*2,breaks.del_effect[1]*5,breaks.del_effect[2],breaks.del_effect[5],breaks.del_effect[5]*2,breaks.del_effect[length(breaks.del_effect)]))
  breaks.del_effect.log<-log10(breaks.del_effect + 0.001)
  labels.del_effect<-as.character(round(breaks.del_effect,2))
  
  cat("labels.del_effect\n")
  cat(str(labels.del_effect))
  cat("\n")
  cat(str(breaks.del_effect.log))
  cat("\n")
  
  
  breaks.hdr_effect<-as.numeric(summary(hdr_effect_df$value[!is.na(hdr_effect_df$value)]))
  breaks.hdr_effect<-unique(c(breaks.hdr_effect[1],breaks.hdr_effect[1]*2,breaks.hdr_effect[1]*5,breaks.hdr_effect[2],breaks.hdr_effect[5],breaks.hdr_effect[5]*2,breaks.hdr_effect[length(breaks.hdr_effect)]))
  breaks.hdr_effect.log<-log10(breaks.hdr_effect + 0.001)
  labels.hdr_effect<-as.character(round(breaks.hdr_effect,2))
  
  cat("labels.hdr_effect\n")
  cat(str(labels.hdr_effect))
  cat("\n")
  cat(str(breaks.hdr_effect.log))
  cat("\n")
  
  
  breaks.Percentage_del_gDNA<-as.numeric(summary(Percentage_del_gDNA_df$value[!is.na(Percentage_del_gDNA_df$value)]))
  breaks.Percentage_del_gDNA<-unique(c(breaks.Percentage_del_gDNA[1],breaks.Percentage_del_gDNA[1]*2,breaks.Percentage_del_gDNA[1]*5,breaks.Percentage_del_gDNA[2],breaks.Percentage_del_gDNA[5],breaks.Percentage_del_gDNA[5]*2,breaks.Percentage_del_gDNA[length(breaks.Percentage_del_gDNA)]))
  breaks.Percentage_del_gDNA.log<-log10(breaks.Percentage_del_gDNA + 0.001)
  labels.Percentage_del_gDNA<-as.character(round(breaks.Percentage_del_gDNA,2))
  
  cat("labels.Percentage_del_gDNA\n")
  cat(str(labels.Percentage_del_gDNA))
  cat("\n")
  cat(str(breaks.Percentage_del_gDNA))
  cat("\n")
  cat(str(breaks.Percentage_del_gDNA.log))
  cat("\n")
  
  # quit(status=1)
  
  
  breaks.Percentage_hdr_gDNA<-as.numeric(summary(Percentage_hdr_gDNA_df$value[!is.na(Percentage_hdr_gDNA_df$value)]))
  breaks.Percentage_hdr_gDNA<-unique(c(breaks.Percentage_hdr_gDNA[1],breaks.Percentage_hdr_gDNA[1]*2,breaks.Percentage_hdr_gDNA[1]*5,breaks.Percentage_hdr_gDNA[2],breaks.Percentage_hdr_gDNA[5],breaks.Percentage_hdr_gDNA[5]*2,breaks.Percentage_hdr_gDNA[length(breaks.Percentage_hdr_gDNA)]))
  breaks.Percentage_hdr_gDNA.log<-log10(breaks.Percentage_hdr_gDNA + 0.001)
  labels.Percentage_hdr_gDNA<-as.character(round(breaks.Percentage_hdr_gDNA,2))
  
  cat("labels.Percentage_hdr_gDNA\n")
  cat(str(labels.Percentage_hdr_gDNA))
  cat("\n")
  cat(str(breaks.Percentage_hdr_gDNA.log))
  cat("\n")
  
  #### LOOP HGNC array ----
  
  HGNC_array<-unique(del_logpval_df$HGNC)
  
  cat("HGNC_array\n")
  cat(str(HGNC_array))
  cat("\n")
  
  for(i in 1:length(HGNC_array))
  {
    HGNC_array_sel<-HGNC_array[i]
    
    cat("--->\t")
    cat(sprintf(as.character(HGNC_array_sel)))
    cat("\t")
    
    del_logpval_df_sel<-del_logpval_df[which(del_logpval_df$HGNC == HGNC_array_sel),]
    
    # cat("del_logpval_df_sel\n")
    # cat(str(del_logpval_df_sel))
    # cat("\n")
    
    VAR_sel<-unique(del_logpval_df_sel$VAR)
    cat("--->\t")
    cat(sprintf(as.character(VAR_sel)))
    cat("\t")
    
    rsid_sel<-unique(del_logpval_df_sel$rsid)
    cat("--->\t")
    cat(sprintf(as.character(rsid_sel)))
    cat("\n")
    
    del_effect_df_sel<-del_effect_df[which(del_effect_df$HGNC == HGNC_array_sel),]
    
    cat("del_effect_df_sel\n")
    cat(str(del_effect_df_sel))
    cat("\n")
    
    Percentage_del_gDNA_df_sel<-Percentage_del_gDNA_df[which(Percentage_del_gDNA_df$HGNC == HGNC_array_sel),]
    
    cat("Percentage_del_gDNA_df_sel\n")
    cat(str(Percentage_del_gDNA_df_sel))
    cat("\n")
    
    hdr_logpval_df_sel<-hdr_logpval_df[which(hdr_logpval_df$HGNC == HGNC_array_sel),]
    
    cat("hdr_logpval_df_sel\n")
    cat(str(hdr_logpval_df_sel))
    cat("\n")
    
    hdr_effect_df_sel<-hdr_effect_df[which(hdr_effect_df$HGNC == HGNC_array_sel),]
    
    cat("hdr_effect_df_sel\n")
    cat(str(hdr_effect_df_sel))
    cat("\n")
    
    Percentage_hdr_gDNA_df_sel<-Percentage_hdr_gDNA_df[which(Percentage_hdr_gDNA_df$HGNC == HGNC_array_sel),]
    
    cat("Percentage_hdr_gDNA_df_sel\n")
    cat(str(Percentage_hdr_gDNA_df_sel))
    cat("\n")
    
    
    #### G1 Percentage_del_gDNA_df_sel vs del_logpval_df_sel ----
    
    G1.df<-merge(del_logpval_df_sel,
                 Percentage_del_gDNA_df_sel,
                 by=c("VAR","rsid","HGNC","MOCK_COORD","Cell_Type","Rep"),
                 all.x=T)
    
    colnames(G1.df)[which(colnames(G1.df) == "value.x")]<-"del_logpval"
    colnames(G1.df)[which(colnames(G1.df) == "value.y")]<-"Percentage_del_gDNA"
    
    
    cat("G1.df\n")
    cat(str(G1.df))
    cat("\n")
    
    fill_color<-df.color_Cell_Type$colors[which(df.color_Cell_Type$Cell_Type%in%G1.df$Cell_Type &
                                                  df.color_Cell_Type$Rep%in%G1.df$Rep)]
    
    
    
    cat("fill_color_\n")
    cat(str(fill_color))
    cat("\n")
    
    
    # quit(status = 1)
    
    graph_del_logpval_Percentage_del_gDNA<-ggplot()+
      geom_point(data=G1.df,
                 aes(x=log10(na.omit(Percentage_del_gDNA)), 
                     y=del_logpval, 
                     color=Cell_Type),size=5)+
      theme_bw()+
      theme(axis.title.y=element_text(size=18, family="sans"),
            axis.title.x=element_text(size=18, family="sans"),
            axis.text.y=element_text(angle=0,size=12, color="black", family="sans"),
            axis.text.x=element_text(angle=0,size=12, color="black", family="sans"),
            legend.title=element_text(size=16,color="black", family="sans"),
            legend.text=element_text(size=12,color="black", family="sans"))+
      scale_x_continuous(name="genIE Percentage_del_gDNA", breaks=breaks.Percentage_del_gDNA.log,
                         labels=labels.Percentage_del_gDNA, 
                         limits=c(breaks.Percentage_del_gDNA.log[1],breaks.Percentage_del_gDNA.log[length(breaks.Percentage_del_gDNA.log)]))+
      scale_y_continuous(name="del_logpval",breaks=breaks.del_logpval,
                         labels=labels.del_logpval, 
                         limits=c(breaks.del_logpval[1],breaks.del_logpval[length(breaks.del_logpval)]))+
      scale_color_manual(values=fill_color, drop=T)+
      theme(legend.position="bottom",legend.title=element_blank(), legend.text = element_text(size=10))+
      guides(color=guide_legend(nrow=2,byrow=TRUE))+
      geom_hline(yintercept=1.3, linetype='dotted', col = 'red')+
      geom_vline(xintercept=log10(4), linetype='dotted', col = 'red')+
      ggeasy::easy_center_title()
    
   graph_del_logpval_Percentage_del_gDNA<-graph_del_logpval_Percentage_del_gDNA+
      geom_text_repel(data=G1.df,
                      aes(x=log10(na.omit(Percentage_del_gDNA)),
                          y=del_logpval,
                          label=Rep),
                      nudge_x = .15,
                      box.padding = 0.5,
                      nudge_y = 1,
                      segment.curvature = -0.1,
                      segment.ncp = 3,
                      segment.angle = 20,
                      max.overlaps = Inf)+
      geom_label_repel(size = 8)
    
    
   #  setwd(out)
   # svglite(paste('del_logpval_Percentage_del_gDNA_',HGNC_array_sel,'.svg',sep=''), width = 8, height = 8)
   # print(graph_del_logpval_Percentage_del_gDNA)
   # dev.off()
    
    #### G1 del_effect_df_sel vs del_logpval_df_sel ----
    
    G1.df<-merge(del_logpval_df_sel,
                 del_effect_df_sel,
                 by=c("VAR","rsid","HGNC","MOCK_COORD","Cell_Type","Rep"),
                 all.x=T)
    
    colnames(G1.df)[which(colnames(G1.df) == "value.x")]<-"del_logpval"
    colnames(G1.df)[which(colnames(G1.df) == "value.y")]<-"del_effect"
    
    
    cat("G1.df\n")
    cat(str(G1.df))
    cat("\n")
    
    fill_color<-df.color_Cell_Type$colors[which(df.color_Cell_Type$Cell_Type%in%G1.df$Cell_Type &
                                                  df.color_Cell_Type$Rep%in%G1.df$Rep)]
    
    
    
    cat("fill_color_\n")
    cat(str(fill_color))
    cat("\n")
    
    
    graph_del_logpval_del_effect<-ggplot()+
      geom_point(data=G1.df,
                 aes(x=log10(na.omit(del_effect)), 
                     y=del_logpval, 
                     color=Cell_Type),size=5)+
      theme_bw()+
      theme(axis.title.y=element_text(size=18, family="sans"),
            axis.title.x=element_text(size=18, family="sans"),
            axis.text.y=element_text(angle=0,size=12, color="black", family="sans"),
            axis.text.x=element_text(angle=0,size=12, color="black", family="sans"),
            legend.title=element_text(size=16,color="black", family="sans"),
            legend.text=element_text(size=12,color="black", family="sans"))+
      scale_x_continuous(name="genIE del_effect", breaks=breaks.del_effect.log,
                         labels=labels.del_effect, 
                         limits=c(breaks.del_effect.log[1],breaks.del_effect.log[length(breaks.del_effect.log)]))+
      scale_y_continuous(name="del_logpval",breaks=breaks.del_logpval,
                         labels=labels.del_logpval, 
                         limits=c(breaks.del_logpval[1],breaks.del_logpval[length(breaks.del_logpval)]))+
      scale_color_manual(values=fill_color, drop=T)+
      theme(legend.position="bottom",legend.title=element_blank(), legend.text = element_text(size=10))+
      guides(color=guide_legend(nrow=2,byrow=TRUE))+
      geom_hline(yintercept=1.3, linetype='dotted', col = 'red')+
      geom_vline(xintercept=log10(1), linetype='dotted', col = 'red')+
      ggeasy::easy_center_title()
    
    graph_del_logpval_del_effect<-graph_del_logpval_del_effect+
      geom_text_repel(data=G1.df,
                      aes(x=log10(na.omit(del_effect)),
                          y=del_logpval,
                          label=Rep),
                      nudge_x = .15,
                      box.padding = 0.5,
                      nudge_y = 1,
                      segment.curvature = -0.1,
                      segment.ncp = 3,
                      segment.angle = 20,
                      max.overlaps = Inf)+
      geom_label_repel(size = 8)
    
    #### G1 Percentage_hdr_gDNA_df_sel vs hdr_logpval_df_sel ----
    
    G1.df<-merge(hdr_logpval_df_sel,
                 Percentage_hdr_gDNA_df_sel,
                 by=c("VAR","rsid","HGNC","MOCK_COORD","Cell_Type","Rep"),
                 all.x=T)
    
    colnames(G1.df)[which(colnames(G1.df) == "value.x")]<-"hdr_logpval"
    colnames(G1.df)[which(colnames(G1.df) == "value.y")]<-"Percentage_hdr_gDNA"
    
    
    cat("G1.df\n")
    cat(str(G1.df))
    cat("\n")
    
    fill_color<-df.color_Cell_Type$colors[which(df.color_Cell_Type$Cell_Type%in%G1.df$Cell_Type &
                                                  df.color_Cell_Type$Rep%in%G1.df$Rep)]
    
    
    
    cat("fill_color_\n")
    cat(str(fill_color))
    cat("\n")
    
    
    # quit(status = 1)
    
    graph_hdr_logpval_Percentage_hdr_gDNA<-ggplot()+
      geom_point(data=G1.df,
                 aes(x=log10(na.omit(Percentage_hdr_gDNA)), 
                     y=hdr_logpval, 
                     color=Cell_Type),size=5)+
      theme_bw()+
      theme(axis.title.y=element_text(size=18, family="sans"),
            axis.title.x=element_text(size=18, family="sans"),
            axis.text.y=element_text(angle=0,size=12, color="black", family="sans"),
            axis.text.x=element_text(angle=0,size=12, color="black", family="sans"),
            legend.title=element_text(size=16,color="black", family="sans"),
            legend.text=element_text(size=12,color="black", family="sans"))+
      scale_x_continuous(name="genIE Percentage_hdr_gDNA", breaks=breaks.Percentage_hdr_gDNA.log,
                         labels=labels.Percentage_hdr_gDNA, 
                         limits=c(breaks.Percentage_hdr_gDNA.log[1],breaks.Percentage_hdr_gDNA.log[length(breaks.Percentage_hdr_gDNA.log)]))+
      scale_y_continuous(name="hdr_logpval",breaks=breaks.hdr_logpval,
                         labels=labels.hdr_logpval, 
                         limits=c(breaks.hdr_logpval[1],breaks.hdr_logpval[length(breaks.hdr_logpval)]))+
      scale_color_manual(values=fill_color, drop=T)+
      theme(legend.position="bottom",legend.title=element_blank(), legend.text = element_text(size=10))+
      guides(color=guide_legend(nrow=2,byrow=TRUE))+
      geom_hline(yintercept=1.3, linetype='dotted', col = 'red')+
      geom_vline(xintercept=log10(4), linetype='dotted', col = 'red')+
      ggeasy::easy_center_title()
    
    graph_hdr_logpval_Percentage_hdr_gDNA<-graph_hdr_logpval_Percentage_hdr_gDNA+
      geom_text_repel(data=G1.df,
                      aes(x=log10(na.omit(Percentage_hdr_gDNA)),
                          y=hdr_logpval,
                          label=Rep),
                      nudge_x = .15,
                      box.padding = 0.5,
                      nudge_y = 1,
                      segment.curvature = -0.1,
                      segment.ncp = 3,
                      segment.angle = 20,
                      max.overlaps = Inf)+
      geom_label_repel(size = 8)
    
    
    
    
    
    #### G1 hdr_effect_df_sel vs hdr_logpval_df_sel ----
    
    G1.df<-merge(hdr_logpval_df_sel,
                 hdr_effect_df_sel,
                 by=c("VAR","rsid","HGNC","MOCK_COORD","Cell_Type","Rep"),
                 all.x=T)
    
    colnames(G1.df)[which(colnames(G1.df) == "value.x")]<-"hdr_logpval"
    colnames(G1.df)[which(colnames(G1.df) == "value.y")]<-"hdr_effect"
    
    
    cat("G1.df\n")
    cat(str(G1.df))
    cat("\n")
    
    fill_color<-df.color_Cell_Type$colors[which(df.color_Cell_Type$Cell_Type%in%G1.df$Cell_Type &
                                                  df.color_Cell_Type$Rep%in%G1.df$Rep)]
    
    
    
    cat("fill_color_\n")
    cat(str(fill_color))
    cat("\n")
    
    
    graph_hdr_logpval_hdr_effect<-ggplot()+
      geom_point(data=G1.df,
                 aes(x=log10(na.omit(hdr_effect)), 
                     y=hdr_logpval, 
                     color=Cell_Type),size=5)+
      theme_bw()+
      theme(axis.title.y=element_text(size=18, family="sans"),
            axis.title.x=element_text(size=18, family="sans"),
            axis.text.y=element_text(angle=0,size=12, color="black", family="sans"),
            axis.text.x=element_text(angle=0,size=12, color="black", family="sans"),
            legend.title=element_text(size=16,color="black", family="sans"),
            legend.text=element_text(size=12,color="black", family="sans"))+
      scale_x_continuous(name="genIE hdr_effect", breaks=breaks.hdr_effect.log,
                         labels=labels.hdr_effect, 
                         limits=c(breaks.hdr_effect.log[1],breaks.hdr_effect.log[length(breaks.hdr_effect.log)]))+
      scale_y_continuous(name="hdr_logpval",breaks=breaks.hdr_logpval,
                         labels=labels.hdr_logpval, 
                         limits=c(breaks.hdr_logpval[1],breaks.hdr_logpval[length(breaks.hdr_logpval)]))+
      scale_color_manual(values=fill_color, drop=T)+
      theme(legend.position="bottom",legend.title=element_blank(), legend.text = element_text(size=10))+
      guides(color=guide_legend(nrow=2,byrow=TRUE))+
      geom_hline(yintercept=1.3, linetype='dotted', col = 'red')+
      geom_vline(xintercept=log10(1), linetype='dotted', col = 'red')+
      ggeasy::easy_center_title()
    
    graph_hdr_logpval_hdr_effect<-graph_hdr_logpval_hdr_effect+
      geom_text_repel(data=G1.df,
                      aes(x=log10(na.omit(hdr_effect)),
                          y=hdr_logpval,
                          label=Rep),
                      nudge_x = .15,
                      box.padding = 0.5,
                      nudge_y = 1,
                      segment.curvature = -0.1,
                      segment.ncp = 3,
                      segment.angle = 20,
                      max.overlaps = Inf)+
      geom_label_repel(size = 8)
    

    #### SAVE path2 ----
    
    source_name<-paste(HGNC_array_sel,
                       VAR_sel,
                       rsid_sel,sep='__')
    
    cat("--->\t")
    cat(sprintf(as.character(source_name)))
    cat("\n")
    
    path2<-paste(out,'Per_variant_graphs','/',source_name, sep='')
    
    cat("path2\n")
    cat(sprintf(as.character(path2)))
    cat("\n")
    
    if (file.exists(path2)){
      
      setwd(path2)
    } else {
      dir.create(file.path(path2))
      setwd(path2)
      
    }
    
    svglite(paste('del_logpval_Percentage_del_gDNA_',HGNC_array_sel,'.svg',sep=''), width = 8, height = 8)
    print(graph_del_logpval_Percentage_del_gDNA)
    dev.off()
    
    svglite(paste('del_logpval_del_effect_',HGNC_array_sel,'.svg',sep=''), width = 8, height = 8)
    print(graph_del_logpval_del_effect)
    dev.off()
    
    svglite(paste('hdr_logpval_Percentage_hdr_gDNA_',HGNC_array_sel,'.svg',sep=''), width = 8, height = 8)
    print(graph_hdr_logpval_Percentage_hdr_gDNA)
    dev.off()
    
    svglite(paste('hdr_logpval_hdr_effect_',HGNC_array_sel,'.svg',sep=''), width = 8, height = 8)
    print(graph_hdr_logpval_hdr_effect)
    dev.off()
   
    # quit(status=1)
    
  }#i HGNC
  
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
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
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
  
 

per_variant_graphs(opt)
  
  
  
}


###########################################################################

system.time( main() )
