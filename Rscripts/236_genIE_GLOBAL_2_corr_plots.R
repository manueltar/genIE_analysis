

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


corr.graphs = function(option_list)
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
  
  
  # quit(status=1)
  
  #### list_results ----
  
  list_results<-list()
  
  list_heatmaps<-list()
  
  #### levels_variable_array ----
  
  levels_variable_array<-levels(Gather.m_finite$variable)
  
  
  cat("levels_variable_array_\n")
  cat(sprintf(as.character(levels_variable_array)))
  cat("\n")
  
  #### LOOP variables array ----
  
  
  for(i in 1:length(levels_variable_array))
  {
    levels_variable_array_sel<-levels_variable_array[i]
    
    cat("---->\t")
    cat(sprintf(as.character(levels_variable_array_sel)))
    cat("\n")
    
    Gather.m_finite_sel<-Gather.m_finite[which(Gather.m_finite$variable == levels_variable_array_sel),]
    
    cat("Gather.m_finite_sel\n")
    str(Gather.m_finite_sel)
    cat("\n")
    
    
    
    
    
    Gather.m_finite_sel_wide<-as.data.frame(pivot_wider(Gather.m_finite_sel,
                                                                        id_cols=c("VAR","variable","rsid","HGNC"),
                                                                        names_from=Rep,
                                                                        values_from=value), stringsAsFactors=F)
    
    cat("Gather.m_finite_sel_wide\n")
    cat(str(Gather.m_finite_sel_wide))
    cat("\n")
    
    Gather.m_finite_sel_wide_ACCEPTED<-Gather.m_finite_sel_wide[-which(Gather.m_finite_sel_wide$HGNC%in%DESIGN_DROPOUTS |
                                                                                                         Gather.m_finite_sel_wide$HGNC%in%EDITION_DROPOUTS |
                                                                                                         Gather.m_finite_sel_wide$HGNC%in%EXPRESSION_DROPOUTS),]
    
    cat("Gather.m_finite_sel_wide_ACCEPTED\n")
    str(Gather.m_finite_sel_wide_ACCEPTED)
    cat("\n")
    
    
    # quit(status=1)
    
    Gather.m_finite_sel_wide_ACCEPTED_matrix<-as.matrix(Gather.m_finite_sel_wide_ACCEPTED[,-c(which(colnames(Gather.m_finite_sel_wide_ACCEPTED) == "VAR"),which(colnames(Gather.m_finite_sel_wide_ACCEPTED) == "variable"),
                                                                                                                              which(colnames(Gather.m_finite_sel_wide_ACCEPTED) == "rsid"),which(colnames(Gather.m_finite_sel_wide_ACCEPTED) == "HGNC"),
                                                                                                                              which(colnames(Gather.m_finite_sel_wide_ACCEPTED) == "Cell_Type"))])
    row.names(Gather.m_finite_sel_wide_ACCEPTED_matrix)<-Gather.m_finite_sel_wide_ACCEPTED$VAR
    
    cat("Gather.m_finite_sel_wide_ACCEPTED_matrix_2\n")
    cat(str(Gather.m_finite_sel_wide_ACCEPTED_matrix))
    cat("\n")
    
    
    # quit(status=1)
    
    Corr.Matrix<-cor(Gather.m_finite_sel_wide_ACCEPTED_matrix, use="pairwise.complete.obs",method = c("pearson"))
    
    cat("Corr.Matrix\n")
    str(Corr.Matrix)
    cat("\n")
    
    # quit(status=1)
    
    # Reorder the correlation matrix according to clusters
    
    # cormat <- reorder_cormat(Corr.Matrix)
    cormat<-Corr.Matrix
    # 
    # 
    # check<-colnames(cormat)
    # 
    # cat("check_reorder\n")
    # cat(sprintf(as.character(check)))
    # cat("\n")
    
    
    
    
    # cormat<-Corr.Matrix
    lower_tri <- get_lower_tri(cormat)
    
    cat("lower_tri\n")
    str(lower_tri)
    cat("\n")
    
    new_order<-row.names(lower_tri)
    
    cat("new_order\n")
    cat(sprintf(as.character(new_order)))
    cat("\n")
    
    
    #### Graphs ----
    
    
    Gather.m_finite_sel$Rep<-factor(Gather.m_finite_sel$Rep,
                                                    levels=new_order,
                                                    ordered = T)
    
    Gather.m_finite_sel<-Gather.m_finite_sel[order(Gather.m_finite_sel$Rep),]
    
    
    Gather.m_finite_sel$MOCK_COORD<-"MOCK"
    
    
    
    
    
    
    cat("Gather.m_finite_sel\n")
    str(Gather.m_finite_sel)
    cat("\n")
    cat(sprintf(as.character(Gather.m_finite_sel$Cell_Type)))
    cat("\n")
    cat(sprintf(as.character(levels(Gather.m_finite_sel$Cell_Type))))
    cat("\n")
    cat(sprintf(as.character(Gather.m_finite_sel$Rep)))
    cat("\n")
    cat(sprintf(as.character(levels(Gather.m_finite_sel$Rep))))
    cat("\n")
    
    
    # quit(status=1)
    
    h2 <- ggplot(data=Gather.m_finite_sel, aes(x=MOCK_COORD,
                                                               y=Rep,
                                                               fill = Cell_Type))+
      geom_tile(color = "white")+
      scale_fill_manual(values = c('#D6B8E6','#32A852','#553B68','#D45E85','#6DB2EE','#62D07F','#C9244B','#87447B'),
                        drop=F)+
      ggeasy::easy_center_title()+
      theme_minimal()+ # minimal theme
      scale_x_discrete(name=element_blank(), drop =F) +
      scale_y_discrete(name=NULL, drop=F)+
      theme(legend.position = "hidden", legend.title = element_blank())+
      theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
      theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
    
    h3 <- ggplot(data=Gather.m_finite_sel, aes(x=Rep,
                                                               y=MOCK_COORD,
                                                               fill = Cell_Type))+
      geom_tile(color = "white")+
      scale_fill_manual(values = c('#D6B8E6','#32A852','#553B68','#D45E85','#6DB2EE','#62D07F','#C9244B','#87447B'),
                        drop=F)+
      ggeasy::easy_center_title()+
      theme_minimal()+ # minimal theme
      scale_x_discrete(name=element_blank(), drop =F) +
      scale_y_discrete(name=NULL, drop=F)+
      theme(legend.position = "hidden", legend.title = element_blank())+
      theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
      theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
    
    
    # quit(status=1)
    
    
    # Melt the correlation matrix
    
    cat("lower_tri\n")
    str(lower_tri)
    cat("\n")
    
    melted_cormat <- melt(lower_tri, na.rm = TRUE)
    
    colnames(melted_cormat)[which(colnames(melted_cormat) == "Var2")]<-"Rep"
    
    cat("melted_cormat_0\n")
    str(melted_cormat)
    cat("\n")
    
    
    Gather.m_finite_sel_subset<-unique(Gather.m_finite_sel[,c(which(colnames(Gather.m_finite_sel) == "Rep"),
                                                       which(colnames(Gather.m_finite_sel) == "Cell_Type"))])
    
    cat("Gather.m_finite_sel_subset\n")
    str(Gather.m_finite_sel_subset)
    cat("\n")
    
    melted_cormat<-merge(melted_cormat,
                         Gather.m_finite_sel_subset,
                         by="Rep")
    
    cat("melted_cormat_1\n")
    str(melted_cormat)
    cat("\n")
                         
    
    
    levels_Cell_Type<-levels(melted_cormat$Cell_Type)
    
    vector_COORDS<-0
    
    FINAL_COORDS<-NULL
    
    for(i in 1:length(levels_Cell_Type))
    {
      
      Cell_Type_sel<-levels_Cell_Type[i]
      
      cat("--->\t")
      cat(sprintf(as.character(Cell_Type_sel)))
      cat("\n")
      
      melted_cormat_sel<-melted_cormat[which(melted_cormat$Cell_Type == Cell_Type_sel),]
      
      cat("melted_cormat_sel\n")
      str(melted_cormat_sel)
      cat("\n")
      
      melted_cormat_sel_Reps<-unique(as.character(melted_cormat_sel$Rep))
      
      cat("melted_cormat_sel_Reps\n")
      str(melted_cormat_sel_Reps)
      cat("\n")
      
      COORD_X<-length(melted_cormat_sel_Reps)
      
      cat("COORD_X--->\t")
      cat(sprintf(as.character(COORD_X)))
      cat("\n")
      
      COORD_DEF<-sum(vector_COORDS,COORD_X)
      
      cat("COORD_DEF--->\t")
      cat(sprintf(as.character(COORD_DEF)))
      cat("\n")
      
      vector_COORDS[i]<-COORD_X
      
      
      # COORD_Y<-length(melted_cormat_sel_Reps)
      # 
      # cat("COORD_Y--->\t")
      # cat(sprintf(as.character(COORD_Y)))
      # cat("\n")
      
      FINAL_COORDS[i]<-as.numeric(COORD_DEF)
      
    }
    
    FINAL_COORDS<-FINAL_COORDS+0.5
    
    FINAL_COORDS<-FINAL_COORDS[-length(FINAL_COORDS)]
    
    cat("FINAL_COORDS--->\t")
    cat(sprintf(as.character(FINAL_COORDS)))
    cat("\n")
    
    # Create a ggheatmap
    
    breaks.Rank<-c(-1.1,seq(-1,1,by=0.5),1.1)
    labels.Rank<-as.character(breaks.Rank)
    
    
    ggheatmap <- ggplot(melted_cormat, aes(Rep, Var1, fill = value))+
      geom_tile(color = "white")+
      scale_fill_gradient2(name = paste("Pearson","Correlation","Coefficient", sep="\n"),
                           low = "red", high = "green",mid="white",midpoint=0,
                           na.value = NA,
                           breaks=breaks.Rank,
                           labels=labels.Rank,
                           limits=c(breaks.Rank[1],
                                    breaks.Rank[length(breaks.Rank)]))+
      theme_minimal()+ # minimal theme
      scale_x_discrete(name=element_blank(), drop =F) +
      scale_y_discrete(name=NULL, drop=F)+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                       size = 12, hjust = 1))+
      geom_hline(yintercept = FINAL_COORDS, linetype=3, color="black")+
      geom_vline(xintercept = FINAL_COORDS, linetype=3, color="black")+
      coord_fixed()
    
    
    
    ggheatmap2<-ggheatmap+
      theme(axis.text.y  = element_blank())+
      theme(axis.text.x  = element_blank())+
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(0.8, 0.2),
        legend.direction = "horizontal")+
      guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                   title.position = "top", title.hjust = 0.5))+
      theme(text=element_text(size=24,  family="sans"))
    
    
    h3<-plot_grid(NULL,h3,NULL,
                  nrow = 1,
                  ncol = 3,
                  align = "hv",
                  rel_widths = c(0.18, 1,0.1),
                  rel_heights = 1)
    
    
    graph_DEF<-plot_grid(h2,ggheatmap2,
                         nrow = 1,
                         ncol = 2,
                         align = "hv",
                         rel_widths = c(0.065, 1),
                         rel_heights = 1)
    
    graph_DEF<-plot_grid(h3,graph_DEF,
                         nrow = 2,
                         ncol = 1,
                         align = "hv",
                         rel_widths = 1,
                         rel_heights = c(0.065,1))
    
    
    title <- ggdraw() + draw_label(paste("sel",type,sep="  "))
    
    
    graph_DEF2<-plot_grid(h2,ggheatmap,
                          nrow = 2,
                          ncol = 2,
                          align = "hv",
                          rel_widths = c(0.065, 1),
                          rel_heights = c(1,0.05))
    
    graph_DEF3<-plot_grid(title,graph_DEF2,
                          ncol=1,
                          rel_heights = c(0.1,1))
    
    
    
    
    
    
    
    cat("svg_graph\n")
    
    setwd(out)
    
    
    svgname<-paste("Graph_corr_",levels_variable_array_sel,".svg", sep='')
    makesvg = TRUE
    
    if (makesvg == TRUE)
    {
      
      ggsave(svgname, plot= graph_DEF,
             device="svg",
             height=10, width=12)
    }
    
    list_heatmaps[[levels_variable_array_sel]]<-graph_DEF3
    list_results[[levels_variable_array_sel]]<-Gather.m_finite_sel
    
    
    # quit(status = 1)
    
    
  }#i levels_variable_array
  
  
  #### SAVE ----
  
  setwd(out)
  
  filename_1<-paste("corr_PLOTS_QC_PASS",".rds", sep='')
  
  cat("filename_1_\n")
  str(filename_1)
  cat("\n")
  
  saveRDS(list_heatmaps,
          file=filename_1)
  
  
  filename_1<-paste("genIE_results_per_metric",".rds", sep='')
  
  cat("filename_1_\n")
  str(filename_1)
  cat("\n")
  
  saveRDS(list_results,
          file=filename_1)
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
  
 

corr.graphs(opt)
# printer(opt)

  
  
  
}


###########################################################################

system.time( main() )
