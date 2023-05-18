
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
  
  cat("Gather.m_0\n")
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
  
  cat("Gather.m_1\n")
  str(Gather.m)
  cat("\n")
  
  
  #### Consolidate Replicates, 1 replicate per CT ----
  
  Gather.m_subset<-unique(Gather.m[,-c(which(colnames(Gather.m) == "variable"),
                                which(colnames(Gather.m) == "value"))])
  
  cat("Gather.m_subset\n")
  str(Gather.m_subset)
  cat("\n")
  
  
  Gather.m_subset.dt<-data.table(Gather.m_subset, key=c("VAR","rsid","HGNC","Cell_Type"))
  
  
  cat("Gather.m_subset.dt\n")
  str(Gather.m_subset.dt)
  cat("\n")
  
  
  Gather.m_subset_collapsed<-as.data.frame(Gather.m_subset.dt[, .(string_Rep=paste(Rep, collapse="|")), by=key(Gather.m_subset.dt)], stringsAsFactors=F)
  
  cat("Gather.m_subset_collapsed_0\n")
  str(Gather.m_subset_collapsed)
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Gather.m_subset_collapsed$string_Rep))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Gather.m_subset_collapsed$string_Rep)))))
  cat("\n")
  
  
  Gather.m_subset_collapsed$Rep<-"NA"
  
  Gather.m_subset_collapsed$Rep[grep("Rep_3",Gather.m_subset_collapsed$string_Rep)]<-"Rep_3"
  Gather.m_subset_collapsed$Rep[grep("Rep_1",Gather.m_subset_collapsed$string_Rep)]<-"Rep_1"
  Gather.m_subset_collapsed$Rep[grep("Rep_4",Gather.m_subset_collapsed$string_Rep)]<-"Rep_4"
  
  Gather.m_subset_collapsed$Rep[grep("Rep_5",Gather.m_subset_collapsed$string_Rep)]<-"Rep_5"
  Gather.m_subset_collapsed$Rep[grep("Rep_6",Gather.m_subset_collapsed$string_Rep)]<-"Rep_6"
  Gather.m_subset_collapsed$Rep[grep("Rep_7",Gather.m_subset_collapsed$string_Rep)]<-"Rep_7"
  
  cat("Gather.m_subset_collapsed_1\n")
  str(Gather.m_subset_collapsed)
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Gather.m_subset_collapsed$Rep))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Gather.m_subset_collapsed$Rep)))))
  cat("\n")
  
  Gather.m_consolidated<-merge(Gather.m,
                               Gather.m_subset_collapsed,
                               merge=c("VAR","rsid","HGNC","Cell_Type","Rep"))
  
  cat("Gather.m_consolidated\n")
  str(Gather.m_consolidated)
  cat("\n")
  
 
  
  ##### Turn Inf values into something meaningful-----
  
  # Gather.m_consolidated_infinite<-Gather.m_consolidated[is.infinite(Gather.m_consolidated$del_effect),]
  # 
  # cat("Gather.m_consolidated_infinite\n")
  # str(Gather.m_consolidated_infinite)
  # cat("\n")
  
  
  Gather.m_consolidated_finite<-Gather.m_consolidated#[is.finite(Gather.m_consolidated$del_effect),]
  
  cat("Gather.m_consolidated_finite_0\n")
  str(Gather.m_consolidated_finite)
  cat("\n")
  cat(sprintf(as.character(names(summary(Gather.m_consolidated_finite$variable)))))
  cat("\n")
  cat(sprintf(as.character(summary(Gather.m_consolidated_finite$variable))))
  cat("\n")
  
  
  
  #### Cell_Type colors ----
  
  Cell_Type_levels<-levels(Gather.m_consolidated_finite$Cell_Type)
  colors_Cell_Type_levels<-c('firebrick2','#32A852','#553B68','#D45E85','#6DB2EE','#62D07F','#C9244B','#87447B','#D6B8E6')
  
  df.color_Cell_Type<-as.data.frame(cbind(Cell_Type_levels,colors_Cell_Type_levels[1:length(Cell_Type_levels)]), stringsAsFactors=F)
  colnames(df.color_Cell_Type)<-c("Cell_Type","colors")
  
  cat("df.color_Cell_Type_0\n")
  cat(str(df.color_Cell_Type))
  cat("\n")
  
  #### RMV CTRL GENE ----
  
  Gather.m_consolidated_finite<-Gather.m_consolidated_finite[which(Gather.m_consolidated_finite$HGNC != "EROS"),]
  
  
  cat("Gather.m_consolidated_finite_1\n")
  str(Gather.m_consolidated_finite)
  cat("\n")
  
  #### string Rsid|HGNC ----
  
  Gather.m_consolidated_finite$string_ID<-paste(Gather.m_consolidated_finite$rsid, Gather.m_consolidated_finite$HGNC, sep='|')
  
  #### mycols ----
  
  mycols <- c("red","black","blue","dark cyan","goldenrod1","yellow2")
  my_HGNC<- c("C2CD5","SH2B3","FOXP1","CUX1","UGCG","BID")
  
  
  #### Del Analysis ----
  
  path7<-paste(out,'/','Deletion_analysis','/', sep='')
  
  cat("path7\n")
  cat(sprintf(as.character(path7)))
  cat("\n")
  
  
  if (file.exists(path7)){
    
    
    
    
  } else {
    dir.create(file.path(path7))
    
  }
  
  
  
  Gather.m_consolidated_finite_wide<-as.data.frame(pivot_wider(Gather.m_consolidated_finite,
                                        names_from=variable,
                                        values_from=value), stringsAsFactors=F)
  
  cat("Gather.m_consolidated_finite_wide_0\n")
  cat(str(Gather.m_consolidated_finite_wide))
  cat("\n")
  
  Gather.m_consolidated_finite_wide<-Gather.m_consolidated_finite_wide[order(Gather.m_consolidated_finite_wide$Percentage_del_gDNA, decreasing =F),]
  
  levels_VAR<-unique(as.character(Gather.m_consolidated_finite_wide$VAR))
  
  cat("levels_VAR_0\n")
  cat(str(levels_VAR))
  cat("\n")
  
  Gather.m_consolidated_finite_wide$VAR<-factor(Gather.m_consolidated_finite_wide$VAR,
                         levels=levels_VAR,
                         ordered=T)
  
  Gather.m_consolidated_finite_wide$QC<-"NA"
  
  Gather.m_consolidated_finite_wide$QC[which(Gather.m_consolidated_finite_wide$Percentage_del_gDNA >=4)]<-"PASS"
  Gather.m_consolidated_finite_wide$QC[which(Gather.m_consolidated_finite_wide$Percentage_del_gDNA < 4)]<-"EDITION_DROPOUT"
  
  
  Gather.m_consolidated_finite_wide$QC<-factor(Gather.m_consolidated_finite_wide$QC,
                        levels=c("EDITION_DROPOUT","PASS"),
                        ordered=T)
  
  
  
  cat("Gather.m_consolidated_finite_wide_1\n")
  cat(str(Gather.m_consolidated_finite_wide))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Gather.m_consolidated_finite_wide$QC))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Gather.m_consolidated_finite_wide$QC)))))
  cat("\n")
  
  
  
  Gather.m_consolidated_finite_wide_PASS<-Gather.m_consolidated_finite_wide[which(Gather.m_consolidated_finite_wide$QC != "EDITION_DROPOUT"),]
  
  
  cat("Gather.m_consolidated_finite_wide_PASS_0\n")
  cat(str(Gather.m_consolidated_finite_wide_PASS))
  cat("\n")
  

  
  
  
  
  if(dim(Gather.m_consolidated_finite_wide_PASS)[1] >0)
  {
    
    A<-summary(Gather.m_consolidated_finite_wide_PASS$del_effect)
    
    
    
    cat("A\n")
    cat(sprintf(as.character(names(A))))
    cat("\n")
    cat(sprintf(as.character(A)))
    cat("\n")
    
    
    
    max_value<-A[6]
    min_value<-A[1]
    # min_value<-0.4
    # max_value<-3
    
    step<-(max_value-min_value)/5
    
    if(step == 0)
    {
      step<-0.1
      
    }
    
    breaks.value<-sort(c(seq(min_value-step,max_value+step,by=step)))
    labels_value<-as.character(round(breaks.value,1))
    
    
    cat("labels_value\n")
    cat(sprintf(as.character(labels_value)))
    cat("\n")
    
    
    A<-summary(Gather.m_consolidated_finite_wide_PASS$del_logpval)
    
    
    
    cat("A\n")
    cat(sprintf(as.character(names(A))))
    cat("\n")
    cat(sprintf(as.character(A)))
    cat("\n")
    
    max_value<-A[6]
    min_value<-0
    
    
    step<-round((max_value-min_value)/5,0)
    
    if(step == 0)
    {
      step<-0.5
      
    }
    
    
    cat("step_ASE_logpval\n")
    cat(sprintf(as.character(step)))
    cat("\n")
    # if(step)
    
    breaks.del_logpval<-unique(seq(min_value,max_value+step, by=step))
    labels.del_logpval<-as.character(round(breaks.del_logpval,1))
    
    
    cat("labels.del_logpval\n")
    cat(sprintf(as.character(labels.del_logpval)))
    cat("\n")
    
    
    
    
    
    #### CT loop del Analysis ----
    
    
    Cell_Type_array<-unique(Gather.m_consolidated_finite_wide_PASS$Cell_Type)
    
    
    cat("Cell_Type_array\n")
    cat(str(Cell_Type_array))
    cat("\n")
    
    for(i in 1:length(Cell_Type_array))
    {
      Cell_Type_array_sel<-Cell_Type_array[i]
      
      cat("--------------------------------------------------->\t")
      cat(sprintf(as.character(Cell_Type_array_sel)))
      cat("\n")
      
      
      df_CT<-Gather.m_consolidated_finite_wide_PASS[which(Gather.m_consolidated_finite_wide_PASS$Cell_Type == Cell_Type_array_sel),]
      
      cat("df_CT_0\n")
      cat(str(df_CT))
      cat("\n")
      
     
      
      #### PDF Klaudia
   
      setwd(path7)
      #setwd(out)
      
      
      pdfname<-paste("volcano_","del_logpval","_","del_effect_size","_",Cell_Type_array_sel,".pdf",sep='')
      pdf(file=pdfname, width=5, height=4, pointsize=12)
      
      par(mai=c(0.9,0.9,0.3,0.2))
      
      
      
      
      
      ind <- which(df_CT$del_logpval < 1.3); length(ind)
      
      
      plot(df_CT$del_effect, df_CT$del_logpval, ty="n", xlab="del effect", ylab="-log10pval", 
           axes=F, cex.lab=1.2, cex.lab=1.3, xlim=c(breaks.value[1],breaks.value[length(breaks.value)]),
                                             ylim=c(breaks.del_logpval[1],breaks.del_logpval[length(breaks.del_logpval)]))
      abline(v=1, col="black",lty=1,lwd=2)
      abline(h=1.3, col="black",lty=2,lwd=2)
      
      
      points(df_CT$del_effect[ind], df_CT$del_logpval[ind], col="darkgrey", pch=19)
      
      sdf_CT <- df_CT[-ind,]
      
      lab <- as.character(unique(sdf_CT$string_ID))
      
      
      cat("lab\n")
      cat(sprintf(as.character(lab)))
      cat("\n")
      
      sel_cols<-NULL
      
      for (iteration_length_lab in 1:length(lab))  {
        
        cat(sprintf(as.character(lab[iteration_length_lab])))
        cat("\n")
        
        ind <- which(sdf_CT$string_ID==lab[iteration_length_lab])
        
        cat("------------------->ind\n")
        cat(str(ind))
        cat("\n")
        
        
        my_HGNC_index<-which(my_HGNC%in%sdf_CT$HGNC[ind])
        
        # cat("------------------->my_HGNC_index\n")
        # cat(str(sdf_CT$HGNC[ind]))
        # cat("\n")
        # cat(str(my_HGNC_index))
        # cat("\n")
        # cat(str(mycols[my_HGNC_index]))
        # cat("\n")
        
        sel_cols[iteration_length_lab]<-mycols[my_HGNC_index]
        
        
        points(sdf_CT$del_effect[ind], sdf_CT$del_logpval[ind], pch=19, col=mycols[my_HGNC_index])
      }
      legend("topright", legend=lab, fill=sel_cols, border=sel_cols, bty="n")
      axis(1, las=1)
      axis(2, las=1)
      
      # axis(1, at=seq(breaks.del_logpval[1],breaks.del_logpval[length(breaks.del_logpval)]))
      # axis(2, at=seq(breaks.value[1],breaks.value[length(breaks.value)]))
      # axis(2, las=1)
      
      dev.off()
      
      
      # #######################################
      # quit(status = 1)
      
    }# i Cell_Type_array
    
    
    
    
    
    
    ##### G1 percentage del gDNA ----

    breaks.Rank<-seq(0,100,by=10)
    labels.Rank<-as.character(breaks.Rank)
    

    dot_plot_Percentage_del_gDNA<-Gather.m_consolidated_finite_wide_PASS %>%
      mutate(myaxis = paste0(rsid, "|", HGNC)) %>%
      mutate(myaxis=fct_reorder(myaxis,as.numeric(VAR))) %>%
      ggplot(aes(x=myaxis, y=Percentage_del_gDNA, color=QC)) +
      geom_point(size=5)+
      theme_classic()+
      scale_y_continuous(name="Percentage del gDNA",breaks=breaks.Rank,labels=labels.Rank, limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]+1))+
      scale_x_discrete(name=NULL, drop=F)+
      theme(axis.title.y=element_text(size=24, family="sans"),
            axis.title.x=element_text(size=24, family="sans"),
            axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
            axis.text.x=element_text(angle=90,vjust=1,hjust=1,size=14, color="black", family="sans"),
            legend.text=element_text(size=12,color="black", family="sans"))+
      scale_color_manual(values=c("green","gray"))+
      geom_hline(yintercept=4, col="red", linetype="dashed")+
      facet_grid(~Cell_Type)+
      theme(legend.position="hidden")+
      
      
      ggeasy::easy_center_title()
     


    setwd(path7)

    svglite(paste('Percentage_del_gDNA','.svg',sep=''), width = 8, height = 8)
    print(dot_plot_Percentage_del_gDNA)
    dev.off()

    cat("dot_plot_Percentage_del_gDNA\n")

    
    
    
    
    
    
    
  }# dim(Gather.m_consolidated_finite_wide_PASS)[1] >0
  
  Del_analysis = Gather.m_consolidated_finite_wide
  
  Del_analysis$QC<-"NA"
  
  Del_analysis$QC[which(Del_analysis$Percentage_del_gDNA >=4 &
                          Del_analysis$del_logpval < 1.3 )]<-"NON_SIGNIFICANT"
  Del_analysis$QC[which(Del_analysis$Percentage_del_gDNA >=4 &
                          Del_analysis$del_logpval >= 1.3 )]<-"SIGNIFICANT"
  Del_analysis$QC[which(Del_analysis$Percentage_del_gDNA < 4)]<-"EDITION_DROPOUT"
  
  
  Del_analysis$QC<-factor(Del_analysis$QC,
                          levels=c("EDITION_DROPOUT","NON_SIGNIFICANT","SIGNIFICANT"),
                          ordered=T)
  
  
  
  cat("Del_analysis_0\n")
  cat(str(Del_analysis))
  cat("\n")
  
  
  colnames(Del_analysis)[which(colnames(Del_analysis) == "QC")]<-"Del_status"
  
  Del_analysis_subset<-unique(Del_analysis[,c(which(colnames(Del_analysis) == "HGNC"),
                                              which(colnames(Del_analysis) == "VAR"),
                                              which(colnames(Del_analysis) == "rsid"),
                                              which(colnames(Del_analysis) == "Del_status"),
                                              which(colnames(Del_analysis) == "Cell_Type"))])
  
  cat("Del_analysis_subset\n")
  str(Del_analysis_subset)
  cat("\n")
  
  Del_analysis_subset_wide<-as.data.frame(pivot_wider(Del_analysis_subset,
                                                      names_from=Cell_Type,
                                                      values_from=Del_status,
                                                      names_prefix="Del_status_"), stringsAsFactors=F)
  
  cat("Del_analysis_subset_wide\n")
  str(Del_analysis_subset_wide)
  cat("\n")
    
  #### HDR Analysis ----
  
  path7<-paste(out,'/','HDR_analysis','/', sep='')
  
  cat("path7\n")
  cat(sprintf(as.character(path7)))
  cat("\n")
  
  
  if (file.exists(path7)){
    
    
    
    
  } else {
    dir.create(file.path(path7))
    
  }
  
  
  
  Gather.m_consolidated_finite_wide<-as.data.frame(pivot_wider(Gather.m_consolidated_finite,
                                                               names_from=variable,
                                                               values_from=value), stringsAsFactors=F)
  
  cat("Gather.m_consolidated_finite_wide_0\n")
  cat(str(Gather.m_consolidated_finite_wide))
  cat("\n")
  
  Gather.m_consolidated_finite_wide<-Gather.m_consolidated_finite_wide[order(Gather.m_consolidated_finite_wide$Percentage_hdr_gDNA, decreasing =F),]
  
  levels_VAR<-unique(as.character(Gather.m_consolidated_finite_wide$VAR))
  
  cat("levels_VAR_0\n")
  cat(str(levels_VAR))
  cat("\n")
  
  Gather.m_consolidated_finite_wide$VAR<-factor(Gather.m_consolidated_finite_wide$VAR,
                                                levels=levels_VAR,
                                                ordered=T)
  
  Gather.m_consolidated_finite_wide$QC<-"NA"
  
  Gather.m_consolidated_finite_wide$QC[which(Gather.m_consolidated_finite_wide$Percentage_hdr_gDNA >=1)]<-"PASS"
  Gather.m_consolidated_finite_wide$QC[which(Gather.m_consolidated_finite_wide$Percentage_hdr_gDNA < 1)]<-"EDITION_DROPOUT"
  
  
  Gather.m_consolidated_finite_wide$QC<-factor(Gather.m_consolidated_finite_wide$QC,
                                               levels=c("EDITION_DROPOUT","PASS"),
                                               ordered=T)
  
  
  
  cat("Gather.m_consolidated_finite_wide_1\n")
  cat(str(Gather.m_consolidated_finite_wide))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Gather.m_consolidated_finite_wide$QC))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Gather.m_consolidated_finite_wide$QC)))))
  cat("\n")
  
  
  
  Gather.m_consolidated_finite_wide_PASS<-Gather.m_consolidated_finite_wide[which(Gather.m_consolidated_finite_wide$QC != "EDITION_DROPOUT"),]
  
  
  cat("Gather.m_consolidated_finite_wide_PASS_0\n")
  cat(str(Gather.m_consolidated_finite_wide_PASS))
  cat("\n")
  
  
  
  if(dim(Gather.m_consolidated_finite_wide_PASS)[1] >0)
  {
    
    A<-summary(Gather.m_consolidated_finite_wide_PASS$hdr_effect)
    
    
    
    cat("A\n")
    cat(sprintf(as.character(names(A))))
    cat("\n")
    cat(sprintf(as.character(A)))
    cat("\n")
    
    
    
    max_value<-A[6]
    min_value<-A[1]
    # min_value<-0.4
    # max_value<-3
    
    step<-(max_value-min_value)/5
    
    if(step == 0)
    {
      step<-0.1
      
    }
    
    breaks.value<-sort(c(seq(min_value-step,max_value+step,by=step)))
    labels_value<-as.character(round(breaks.value,1))
    
    
    cat("labels_value\n")
    cat(sprintf(as.character(labels_value)))
    cat("\n")
    
    
    A<-summary(Gather.m_consolidated_finite_wide_PASS$hdr_logpval)
    
    
    
    cat("A\n")
    cat(sprintf(as.character(names(A))))
    cat("\n")
    cat(sprintf(as.character(A)))
    cat("\n")
    
    max_value<-A[6]
    min_value<-0
    
    
    step<-round((max_value-min_value)/5,0)
    
    if(step == 0)
    {
      step<-0.5
      
    }
    
    
    cat("step_ASE_logpval\n")
    cat(sprintf(as.character(step)))
    cat("\n")
    # if(step)
    
    breaks.hdr_logpval<-unique(seq(min_value,max_value+step, by=step))
    labels.hdr_logpval<-as.character(round(breaks.hdr_logpval,1))
    
    
    cat("labels.hdr_logpval\n")
    cat(sprintf(as.character(labels.hdr_logpval)))
    cat("\n")
    
    
    
    
    
    
    #### CT loop del Analysis ----
    
    
    Cell_Type_array<-unique(Gather.m_consolidated_finite_wide_PASS$Cell_Type)
    
    
    cat("Cell_Type_array\n")
    cat(str(Cell_Type_array))
    cat("\n")
    
    for(i in 1:length(Cell_Type_array))
    {
      Cell_Type_array_sel<-Cell_Type_array[i]
      
      cat("--------------------------------------------------->\t")
      cat(sprintf(as.character(Cell_Type_array_sel)))
      cat("\n")
      
      
      df_CT<-Gather.m_consolidated_finite_wide_PASS[which(Gather.m_consolidated_finite_wide_PASS$Cell_Type == Cell_Type_array_sel),]
      
      cat("df_CT_0\n")
      cat(str(df_CT))
      cat("\n")
      
      
      
      #### PDF Klaudia
      
      setwd(path7)
      #setwd(out)
      
      
      pdfname<-paste("volcano_","hdr_logpval","_","hdr_effect_size","_",Cell_Type_array_sel,".pdf",sep='')
      pdf(file=pdfname, width=5, height=4, pointsize=12)
      
      par(mai=c(0.9,0.9,0.3,0.2))
      
      
      
      
      
      ind <- which(df_CT$hdr_logpval < 1.3); length(ind)
      
      
      plot(df_CT$hdr_effect, df_CT$hdr_logpval, ty="n", xlab="hdr effect", ylab="-log10pval", 
           axes=F, cex.lab=1.2, cex.lab=1.3, xlim=c(breaks.value[1],breaks.value[length(breaks.value)]),
           ylim=c(breaks.hdr_logpval[1],breaks.hdr_logpval[length(breaks.hdr_logpval)]))
      abline(v=1, col="black",lty=1,lwd=2)
      abline(h=1.3, col="black",lty=2,lwd=2)
      
      
      points(df_CT$hdr_effect[ind], df_CT$hdr_logpval[ind], col="darkgrey", pch=19)
      
      sdf_CT <- df_CT[-ind,]
      
      lab <- as.character(unique(sdf_CT$string_ID))
      
      
      cat("lab\n")
      cat(sprintf(as.character(lab)))
      cat("\n")
      
      sel_cols<-NULL
      
      if(length(lab) >0)
      {
        for (iteration_length_lab in 1:length(lab))  {
          
          cat(sprintf(as.character(lab[iteration_length_lab])))
          cat("\n")
          
          ind <- which(sdf_CT$string_ID==lab[iteration_length_lab])
          
          cat("------------------->ind\n")
          cat(str(ind))
          cat("\n")
          
          my_HGNC_index<-which(my_HGNC%in%sdf_CT$HGNC[ind])
          
          sel_cols[iteration_length_lab]<-mycols[my_HGNC_index]
          
          
          points(sdf_CT$hdr_effect[ind], sdf_CT$hdr_logpval[ind], pch=19, col=mycols[my_HGNC_index])
        }
        legend("topright", legend=lab, fill=sel_cols, border=sel_cols, bty="n")
        axis(1, las=1)
        axis(2, las=1)
      }else{
        
        
        axis(1, las=1)
        axis(2, las=1)
      }
      
     
      
      # axis(1, at=seq(breaks.hdr_logpval[1],breaks.hdr_logpval[length(breaks.hdr_logpval)]))
      # axis(2, at=seq(breaks.value[1],breaks.value[length(breaks.value)]))
      # axis(2, las=1)
      
      dev.off()
    }# i Cell_Type_array
    
    
    
    
    
    
    ##### G1 percentage del gDNA ----
    
    breaks.Rank<-seq(0,100,by=10)
    labels.Rank<-as.character(breaks.Rank)
    
    
    dot_plot_Percentage_hdr_gDNA<-Gather.m_consolidated_finite_wide_PASS %>%
      mutate(myaxis = paste0(rsid, "|", HGNC)) %>%
      mutate(myaxis=fct_reorder(myaxis,as.numeric(VAR))) %>%
      ggplot(aes(x=myaxis, y=Percentage_hdr_gDNA, color=QC)) +
      geom_point(size=5)+
      theme_classic()+
      scale_y_continuous(name="Percentage del gDNA",breaks=breaks.Rank,labels=labels.Rank, limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]+1))+
      scale_x_discrete(name=NULL, drop=F)+
      theme(axis.title.y=element_text(size=24, family="sans"),
            axis.title.x=element_text(size=24, family="sans"),
            axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
            axis.text.x=element_text(angle=90,vjust=1,hjust=1,size=14, color="black", family="sans"),
            legend.text=element_text(size=12,color="black", family="sans"))+
      scale_color_manual(values=c("green","gray"))+
      geom_hline(yintercept=1, col="red", linetype="dashed")+
      facet_grid(~Cell_Type)+
      theme(legend.position="hidden")+
      
      
      ggeasy::easy_center_title()
    
    
    
    setwd(path7)
    
    svglite(paste('Percentage_hdr_gDNA','.svg',sep=''), width = 8, height = 8)
    print(dot_plot_Percentage_hdr_gDNA)
    dev.off()
    
    cat("dot_plot_Percentage_hdr_gDNA\n")
    
    
    
    
    
    
    
    
  }# dim(Gather.m_consolidated_finite_wide_PASS)[1] >0
  
  HDR_analysis = Gather.m_consolidated_finite_wide
  
  HDR_analysis$QC<-"NA"
  
  HDR_analysis$QC[which(HDR_analysis$Percentage_del_gDNA >=1 &
                          HDR_analysis$del_logpval < 1.3 )]<-"NON_SIGNIFICANT"
  HDR_analysis$QC[which(HDR_analysis$Percentage_del_gDNA >=1 &
                          HDR_analysis$del_logpval >= 1.3 )]<-"SIGNIFICANT"
  HDR_analysis$QC[which(HDR_analysis$Percentage_del_gDNA < 1)]<-"EDITION_DROPOUT"
  
  
  HDR_analysis$QC<-factor(HDR_analysis$QC,
                          levels=c("EDITION_DROPOUT","NON_SIGNIFICANT","SIGNIFICANT"),
                          ordered=T)
  
  
  
  cat("HDR_analysis_0\n")
  cat(str(HDR_analysis))
  cat("\n")
  
  
  colnames(HDR_analysis)[which(colnames(HDR_analysis) == "QC")]<-"HDR_status"
  
  HDR_analysis_subset<-unique(HDR_analysis[,c(which(colnames(HDR_analysis) == "HGNC"),
                                              which(colnames(HDR_analysis) == "VAR"),
                                              which(colnames(HDR_analysis) == "rsid"),
                                              which(colnames(HDR_analysis) == "HDR_status"),
                                              which(colnames(HDR_analysis) == "Cell_Type"))])
  
  cat("HDR_analysis_subset\n")
  str(HDR_analysis_subset)
  cat("\n")
  
  HDR_analysis_subset_wide<-as.data.frame(pivot_wider(HDR_analysis_subset,
                                                      names_from=Cell_Type,
                                                      values_from=HDR_status,
                                                      names_prefix="HDR_status_"), stringsAsFactors=F)
  
  cat("HDR_analysis_subset_wide\n")
  str(HDR_analysis_subset_wide)
  cat("\n")
  
  ############################### genIE Tiers --------------
  
  DEF<-merge(Del_analysis_subset_wide,
             HDR_analysis_subset_wide,
             by=c("HGNC","VAR","rsid"),#,"Rep","string_Rep"),
             all=T)
  
  cat("DEF_0\n")
  str(DEF)
  cat("\n")
  
  DEF$DEF_CLASS<-"EDITION_DROPOUT"
  
  DEF$DEF_CLASS[grep("BID|CUX1|C2CD5|FOXP1|FUT8|SH2B3|UGCG",DEF$HGNC)]<-"SIGNIFICANT"
  DEF$DEF_CLASS[grep("SLC9A3R1|BRAP",DEF$HGNC)]<-"NON_SIGNIFICANT"
  
 
  
  DEF$DEF_CLASS<-factor(DEF$DEF_CLASS,
                        levels = c("EDITION_DROPOUT","NON_SIGNIFICANT","SIGNIFICANT"))
  
  setwd(out)
 # write.table(DEF,file="test.tsv",sep="\t",quote=F,row.names = F)
  
  saveRDS(DEF, file="genIE_Tiers_DEF.rds")
  
 
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
