

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

volcano_plots_per_CellType = function(option_list)
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
  
  
  
  

 
  

  # quit(status=1)
  
  #### Volcano per Cell Type LOOP ----
  
  
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
    
    ##### X and Y points -----
    
    Percentage_hdr_gDNA_vector<-Gather_sel$value[which(Gather_sel$variable == "Percentage_hdr_gDNA")]
    
    cat(sprintf(as.character(summary(Percentage_hdr_gDNA_vector))))
    cat("\n")
    
    
    breaks.Percentage_hdr_gDNA<-as.numeric(summary(Percentage_hdr_gDNA_vector[!is.na(Percentage_hdr_gDNA_vector)]))
    
    cat("breaks.Percentage_hdr_gDNA\n")
    cat(str(breaks.Percentage_hdr_gDNA))
    cat("\n")
    
    breaks.Percentage_hdr_gDNA[1]<-0
    breaks.Percentage_hdr_gDNA[length(breaks.Percentage_hdr_gDNA)]<-breaks.Percentage_hdr_gDNA[length(breaks.Percentage_hdr_gDNA)]+5
    labels.Percentage_hdr_gDNA<-as.character(seq(breaks.Percentage_hdr_gDNA[1],breaks.Percentage_hdr_gDNA[length(breaks.Percentage_hdr_gDNA)], by=5))
    breaks.Percentage_hdr_gDNA<-as.numeric(labels.Percentage_hdr_gDNA)
    
    cat("labels.Percentage_hdr_gDNA\n")
    cat(sprintf(as.character(labels.Percentage_hdr_gDNA)))
    cat("\n")
    
    ###
    
    Percentage_del_gDNA_vector<-Gather_sel$value[which(Gather_sel$variable == "Percentage_del_gDNA")]
    
    cat(sprintf(as.character(summary(Percentage_del_gDNA_vector))))
    cat("\n")
    
    
    breaks.Percentage_del_gDNA<-as.numeric(summary(Percentage_del_gDNA_vector[!is.na(Percentage_del_gDNA_vector)]))
    
    cat("breaks.Percentage_del_gDNA\n")
    cat(str(breaks.Percentage_del_gDNA))
    cat("\n")
    
    breaks.Percentage_del_gDNA[1]<-0
    breaks.Percentage_del_gDNA[length(breaks.Percentage_del_gDNA)]<-breaks.Percentage_del_gDNA[length(breaks.Percentage_del_gDNA)]+5
    labels.Percentage_del_gDNA<-as.character(seq(breaks.Percentage_del_gDNA[1],breaks.Percentage_del_gDNA[length(breaks.Percentage_del_gDNA)], by=5))
    breaks.Percentage_del_gDNA<-as.numeric(labels.Percentage_del_gDNA)
    
    cat("labels.Percentage_del_gDNA\n")
    cat(sprintf(as.character(labels.Percentage_del_gDNA)))
    cat("\n")
    
    ### hdr_effect LOG
    
    hdr_effect_vector<-Gather_sel$value[which(Gather_sel$variable == "hdr_effect")]
    
    cat(sprintf(as.character(summary(hdr_effect_vector))))
    cat("\n")
    
    
    breaks.hdr_effect<-as.numeric(summary(hdr_effect_vector[!is.na(hdr_effect_vector)]))
    
    cat("breaks.hdr_effect\n")
    cat(str(breaks.hdr_effect))
    cat("\n")
    
    breaks.hdr_effect[1]<-breaks.hdr_effect[1]
    
    
    labels.hdr_effect<-unique(as.character(sort(c(seq(breaks.hdr_effect[1],breaks.hdr_effect[length(breaks.hdr_effect)]+0.5, by=0.5),
                                                  breaks.hdr_effect[length(breaks.hdr_effect)]))))
    breaks.hdr_effect.log<-log10(as.numeric(labels.hdr_effect)+0.00001)
    
    cat("labels.hdr_effect\n")
    cat(sprintf(as.character(labels.hdr_effect)))
    cat("\n")
    
    cat("breaks.hdr_effect.log\n")
    cat(sprintf(as.character(breaks.hdr_effect.log)))
    cat("\n")
    
    ### del_effect LOG
    
    del_effect_vector<-Gather_sel$value[which(Gather_sel$variable == "del_effect")]
    
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
    
    cat("labels.del_effect\n")
    cat(sprintf(as.character(labels.del_effect)))
    cat("\n")
    
    cat("breaks.del_effect.log\n")
    cat(sprintf(as.character(breaks.del_effect.log)))
    cat("\n")
    
    
    #### hdr_logpval_vector
    
    hdr_logpval_vector<-Gather_sel$value[which(Gather_sel$variable == "hdr_logpval")]
    
    cat(sprintf(as.character(summary(hdr_logpval_vector))))
    cat("\n")
    
    
    breaks.hdr_logpval<-as.numeric(summary(hdr_logpval_vector[!is.na(hdr_logpval_vector)]))
    
    cat("breaks.hdr_logpval\n")
    cat(str(breaks.hdr_logpval))
    cat("\n")
    
    
    
    breaks.hdr_logpval[1]<-0
    breaks.hdr_logpval[length(breaks.hdr_logpval)]<-breaks.hdr_logpval[length(breaks.hdr_logpval)]+2
    labels.hdr_logpval<-as.character(seq(breaks.hdr_logpval[1],breaks.hdr_logpval[length(breaks.hdr_logpval)], by=2))
    breaks.hdr_logpval<-as.numeric(labels.hdr_logpval)
    
    cat("labels.hdr_logpval\n")
    cat(str(labels.hdr_logpval))
    cat("\n")
    
    #### del_logpval_vector
    
    del_logpval_vector<-Gather_sel$value[which(Gather_sel$variable == "del_logpval")]
    
    cat(sprintf(as.character(summary(del_logpval_vector))))
    cat("\n")
    
    
    breaks.del_logpval<-as.numeric(summary(del_logpval_vector[!is.na(del_logpval_vector)]))
    
    cat("breaks.del_logpval\n")
    cat(str(breaks.del_logpval))
    cat("\n")
    
    
    
    breaks.del_logpval[1]<-0
    breaks.del_logpval[length(breaks.del_logpval)]<-breaks.del_logpval[length(breaks.del_logpval)]+2
    labels.del_logpval<-as.character(seq(breaks.del_logpval[1],breaks.del_logpval[length(breaks.del_logpval)], by=2))
    breaks.del_logpval<-as.numeric(labels.del_logpval)
    
    cat("labels.del_logpval\n")
    cat(str(labels.del_logpval))
    cat("\n")
    
    #### pivot to wider by variable ----
    
    Gather_sel_wide<-as.data.frame(pivot_wider(Gather_sel,
                                               id_cols=c("Rep","VAR","rsid","HGNC","Cell_Type"),
                                               names_from=variable,
                                               values_from=value), stringsAsFactors=F)
    
    
    cat("Gather_sel_wide\n")
    str(Gather_sel_wide)
    cat("\n")                                                    
                                                                
    
    #### for kolf2 calculate median of values ----
    
    # if(Cell_type_sel == "Kolf2")
    # {
    #   
    #  
    #   
    #   
    #   
    #   
    # }else{
    #   
    #   
    #   Gather_DEF<-Gather_sel_wide
    # }
    
    
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
    
    Gather_sel_wide.dt_median$Rep<-"DUMMY"
    
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
    
    ###### graphs ----
   
    
    
    Volcano_del_effect<-ggplot(data=Gather_DEF,
                               aes(x=log10(del_effect+0.00001), 
                                  y=del_logpval,
                                  label=HGNC)) +
      geom_point(color="gray",
                 size=10)+
      geom_point(data=subset(Gather_DEF, Gather_DEF$del_logpval >= 1.3 &
                               Gather_DEF$Percentage_del_gDNA >= 2),
                 aes(x=log10(del_effect+0.00001), 
                     y=del_logpval), color="red",size=10)+
      theme_bw()+
      scale_x_continuous(name="Del effect", 
                         breaks=breaks.del_effect.log,labels=labels.del_effect, 
                         limits=c(breaks.del_effect.log[1],breaks.del_effect.log[length(breaks.del_effect.log)]))+
      scale_y_continuous(name="Del logpval",
                         breaks=breaks.del_logpval,labels=labels.del_logpval, 
                         limits=c(breaks.del_logpval[1],breaks.del_logpval[length(breaks.del_logpval)]))+
      theme(axis.title.y=element_text(size=24, family="sans"),
            axis.title.x=element_text(size=24, family="sans"),
            axis.text.y=element_text(angle=0,size=24, color="black", family="sans"),
            axis.text.x=element_text(angle=45,vjust=1,hjust=1,size=24, color="black", family="sans"))+
      theme(legend.position="hidden",legend.title=element_blank(), legend.text = element_text(size=20))+
      guides(color=guide_legend(nrow=3,byrow=TRUE))+
      geom_vline(xintercept=log10(1+0.00001), col="gray")+
      geom_label_repel(size = 8)+
      ggeasy::easy_center_title()
    
    # Volcano_del_effect<-Volcano_del_effect+
    #   geom_text_repel(data=Gather_DEF,
    #                   aes(x=log10(del_effect+0.00001), 
    #                       y=del_logpval, label=HGNC),
    #                   box.padding = 0.8, max.overlaps = Inf,segment.linetype = 6)
      #geom_text_repel(box.padding = 0.8, max.overlaps = Inf,segment.linetype = 6)+
    
    
    cat("del_effect_graph\n")
    
    Volcano_hdr_effect<-ggplot(data=Gather_DEF,
                               aes(x=log10(hdr_effect+0.00001), 
                                   y=hdr_logpval,
                                   label=HGNC)) +
      geom_point(color="gray",
                 size=10)+
      geom_point(data=subset(Gather_DEF, Gather_DEF$hdr_logpval >= 1.3 &
                               Gather_DEF$Percentage_hdr_gDNA >= 2),
                 aes(x=log10(hdr_effect+0.00001), 
                     y=hdr_logpval), color="red",size=10)+
      theme_bw()+
      scale_x_continuous(name="hdr effect", 
                         breaks=breaks.hdr_effect.log,labels=labels.hdr_effect, 
                         limits=c(breaks.hdr_effect.log[1],breaks.hdr_effect.log[length(breaks.hdr_effect.log)]))+
      scale_y_continuous(name="hdr logpval",
                         breaks=breaks.hdr_logpval,labels=labels.hdr_logpval, 
                         limits=c(breaks.hdr_logpval[1],breaks.hdr_logpval[length(breaks.hdr_logpval)]))+
      theme(axis.title.y=element_text(size=24, family="sans"),
            axis.title.x=element_text(size=24, family="sans"),
            axis.text.y=element_text(angle=0,size=24, color="black", family="sans"),
            axis.text.x=element_text(angle=45,vjust=1,hjust=1,size=24, color="black", family="sans"))+
      theme(legend.position="bottom",legend.title=element_blank(), legend.text = element_text(size=20))+
      geom_label_repel(size = 8)+
      guides(color=guide_legend(nrow=3,byrow=TRUE))+
      geom_vline(xintercept=log10(1+0.00001), col="gray")+
      ggeasy::easy_center_title()
    
    
    cat("hdr_effect_graph\n")
    
    
    Volcano_Percentage_hdr_gDNA<-ggplot(data=Gather_DEF,
                               aes(x=Percentage_hdr_gDNA, 
                                   y=hdr_logpval,
                                   label=HGNC)) +
      geom_point(color="gray",
                 size=10)+
      geom_point(data=subset(Gather_DEF, Gather_DEF$hdr_logpval >= 1.3  &
                               Gather_DEF$Percentage_hdr_gDNA >= 2),
                 aes(x=Percentage_hdr_gDNA, 
                     y=hdr_logpval), color="red",size=10)+
      theme_bw()+
      scale_x_continuous(name="Percentage_hdr_gDNA", 
                         breaks=breaks.Percentage_hdr_gDNA,labels=labels.Percentage_hdr_gDNA, 
                         limits=c(breaks.Percentage_hdr_gDNA[1],breaks.Percentage_hdr_gDNA[length(breaks.Percentage_hdr_gDNA)]))+
      scale_y_continuous(name="HDR logpval",
                         breaks=breaks.hdr_logpval,labels=labels.hdr_logpval, 
                         limits=c(breaks.hdr_logpval[1],breaks.hdr_logpval[length(breaks.hdr_logpval)]))+
      theme(axis.title.y=element_text(size=24, family="sans"),
            axis.title.x=element_text(size=24, family="sans"),
            axis.text.y=element_text(angle=0,size=24, color="black", family="sans"),
            axis.text.x=element_text(angle=45,vjust=1,hjust=1,size=24, color="black", family="sans"))+
      theme(legend.position="hidden",legend.title=element_blank(), legend.text = element_text(size=10))+
      guides(color=guide_legend(nrow=3,byrow=TRUE))+
      geom_vline(xintercept=4, col="gray")+
      geom_label_repel(size = 8)+
      ggeasy::easy_center_title()
    
    # Volcano_Percentage_hdr_gDNA<-Volcano_Percentage_hdr_gDNA+
    #   geom_text_repel(data=Gather_DEF,
    #                   aes(x=Percentage_hdr_gDNA, 
    #                       y=hdr_logpval, label=HGNC),
    #                   box.padding = 0.8, max.overlaps = Inf,segment.linetype = 6)
    
    
    cat("Percentage_hdr_gDNA_graph\n")
    
    
    Volcano_Percentage_del_gDNA<-ggplot(data=Gather_DEF,
                                        aes(x=Percentage_del_gDNA, 
                                            y=del_logpval,
                                            label=HGNC)) +
      geom_point(color="gray",
                 size=10)+
      geom_point(data=subset(Gather_DEF, Gather_DEF$del_logpval >= 1.3  &
                               Gather_DEF$Percentage_del_gDNA >= 2),
                 aes(x=Percentage_del_gDNA, 
                     y=del_logpval), color="red",size=10)+
      theme_bw()+
      scale_x_continuous(name="Percentage_del_gDNA", 
                         breaks=breaks.Percentage_del_gDNA,labels=labels.Percentage_del_gDNA, 
                         limits=c(breaks.Percentage_del_gDNA[1],breaks.Percentage_del_gDNA[length(breaks.Percentage_del_gDNA)]))+
      scale_y_continuous(name="del logpval",
                         breaks=breaks.del_logpval,labels=labels.del_logpval, 
                         limits=c(breaks.del_logpval[1],breaks.del_logpval[length(breaks.del_logpval)]))+
      theme(axis.title.y=element_text(size=24, family="sans"),
            axis.title.x=element_text(size=24, family="sans"),
            axis.text.y=element_text(angle=0,size=24, color="black", family="sans"),
            axis.text.x=element_text(angle=45,vjust=1,hjust=1,size=24, color="black", family="sans"))+
      theme(legend.position="hidden",legend.title=element_blank(), legend.text = element_text(size=10))+
      guides(color=guide_legend(nrow=3,byrow=TRUE))+
      geom_vline(xintercept=4, col="gray")+
      geom_label_repel(size = 8)+
      ggeasy::easy_center_title()
    
    
    cat("Percentage_del_gDNA_graph\n")
    
    
    cat("Volcano_del_effect_svg_graph\n")
    
    setwd(out)
    
    
    svgname<-paste("Volcano_del_effect_",Cell_type_sel,".svg", sep='')
    makesvg = TRUE
    
    if (makesvg == TRUE)
    {
      
      ggsave(svgname, plot= Volcano_del_effect,
             device="svg",
             height=10, width=12)
    }
    
    cat("Volcano_hdr_effect_svg_graph\n")
    
    
    svgname<-paste("Volcano_hdr_effect_",Cell_type_sel,".svg", sep='')
    makesvg = TRUE
    
    if (makesvg == TRUE)
    {
      
      ggsave(svgname, plot= Volcano_hdr_effect,
             device="svg",
             height=10, width=12)
    }
    
    cat("Volcano_Percentage_hdr_gDNA_svg_graph\n")
    
    
    svgname<-paste("Volcano_Percentage_hdr_gDNA_",Cell_type_sel,".svg", sep='')
    makesvg = TRUE
    
    if (makesvg == TRUE)
    {
      
      ggsave(svgname, plot= Volcano_Percentage_hdr_gDNA,
             device="svg",
             height=10, width=12)
    }
    
    
    cat("Volcano_Percentage_del_gDNA_svg_graph\n")
    
    
    svgname<-paste("Volcano_Percentage_del_gDNA_",Cell_type_sel,".svg", sep='')
    makesvg = TRUE
    
    if (makesvg == TRUE)
    {
      
      ggsave(svgname, plot= Volcano_Percentage_del_gDNA,
             device="svg",
             height=10, width=12)
    }
    
    
    # setwd(out)
    # 
    # pdf(paste("test","_",Cell_type_sel,".pdf",sep=''))
    # print(Volcano_del_effect)
    # print(Volcano_Percentage_del_gDNA)
    # print(Volcano_hdr_effect)
    # print(Volcano_Percentage_hdr_gDNA)
    # dev.off()
    
    
    
    # scale_color_manual(values=c('#32A852','#6DB2EE','#62D07F','#1877C9','#C9244B','#D45E85','#87447B','#553B68','#D6B8E6',
    #                             "#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C","#D7C1B1", "#AD6F3B", "#689030", "red"), drop=F,
    #                    name="Totals",breaks=Freq.total_sel$Label3,
    #                    labels=paste(Freq.total_sel$Label3,
    #                                 Freq.total_sel$Total, sep =' n= '))+
    
    
    
    # quit(status=1)
    
    
   
    
    
  }
  
  
  
}



printer = function(option_list)
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
  
  
  #### READ INPUT FILES ----
  
  setwd(out)
  
  
  filename_1<-paste("corr_PLOTS_QC_PASS",".rds", sep='')
  
  plots_corr<-readRDS(file=filename_1)
  
  
  corr_Percentage_hdr_gDNA<-plots_corr$Percentage_hdr_gDNA
  corr_Percentage_del_gDNA<-plots_corr$Percentage_del_gDNA
  corr_hdr_effect<-plots_corr$hdr_effect
  corr_hdr_logpval<-plots_corr$hdr_logpval
  corr_del_effect<-plots_corr$del_effect
  corr_del_logpval<-plots_corr$del_logpval
  
  
 
  #### PRINTING ----
  
  setwd(out)
  
  cat("PDF REPORT\n")
  
  pdfname<-paste(type,".pdf", sep='')
  
  makepdf = TRUE
  
  if (makepdf == TRUE)
  {
    pdf ( pdfname , height=10, width=12)
  }
  
  layout(matrix(c(1,1,2,2,3,3),1,1, byrow = TRUE))
  par(mar=c(3,4,1,1))
  
  
  
  print(corr_Percentage_hdr_gDNA)
  print(corr_Percentage_del_gDNA)
  print(corr_hdr_effect)
  print(corr_hdr_logpval)
  print(corr_del_effect)
  print(corr_del_logpval)
  
  
  
  if (makepdf == TRUE)
  {
    dev.off()
  }
  
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
  
 

volcano_plots_per_CellType(opt)
# printer(opt)

  
  
  
}


###########################################################################

system.time( main() )
