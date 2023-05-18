
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
library("gtools", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("reshape2", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")



opt = NULL

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
  
  genIE_Tier_1$abs_shift_del_effect<-abs(1-genIE_Tier_1$del_effect)
  
  genIE_Tier_1$abs_shift_hdr_effect<-abs(1-genIE_Tier_1$hdr_effect)
  
  
  cat("genIE_Tier_1\n")
  str(genIE_Tier_1)
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
                                                 values_from=c("STATUS_del",
                                                               "abs_shift_del_effect",
                                                               "abs_shift_hdr_effect")))
    
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
                                                 values_from=c("STATUS_hdr",
                                                               "abs_shift_del_effect",
                                                               "abs_shift_hdr_effect")))
    
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
                                                 values_from=c("STATUS_del",
                                                               "abs_shift_del_effect",
                                                               "abs_shift_hdr_effect")))
    
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
  dB$absolute_finemap_z<-abs(dB$finemap_z)
  
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
  
  # quit(status=1)
  
 
  
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
  
 
 
  
  
  
  ##### LINEAGE specificty 0 vs 1 per Cell Type ----
  
  indx.int<-c(which(colnames(traits_merge) == "VAR"),which(colnames(traits_merge) == "LINEAGE_CLASSIF_DEF"),which(colnames(traits_merge) == "Lineage_CLASS"),which(colnames(traits_merge) == "absolute_finemap_z"))
  
  traits_merge_subset<-unique(traits_merge[,indx.int])
  
  cat("traits_merge_subset_1\n")
  cat(str(traits_merge_subset))
  cat("\n")
  
  
  # quit(status = 1)
  
  traits_merge_subset<-unique(traits_merge_subset[which(traits_merge_subset$LINEAGE_CLASSIF_DEF != ""),])
  
  traits_merge_subset$LINEAGE_CLASSIF_DEF<-factor(traits_merge_subset$LINEAGE_CLASSIF_DEF,
                                                  levels=c("Erythroid Traits","Megakaryocitic Traits","GRAN-MONO Traits","Lymphocitic Traits","Pleiotropic_variants"),
                                                  ordered=T)
  
  
 
    
    
  
  cat("traits_merge_subset_2\n")
  cat(str(traits_merge_subset))
  cat("\n")
  
  #### KEY merge ----
  
  
  REP<-merge(genIE_Tier_1,
             traits_merge_subset,
             by="VAR")
  
 
  
  cat("REP\n")
  str(REP)
  cat("\n")
  
  
  # quit(status = 1)
  
#### LOOP CELL TYPES ----
  
  LINEAGE_CLASSIF_DEF_array<-levels(REP$LINEAGE_CLASSIF_DEF)
  
  cat("LINEAGE_CLASSIF_DEF_array_\n")
  cat(str(LINEAGE_CLASSIF_DEF_array))
  cat("\n")
  
  Cell_Type_array<-levels(REP$Cell_Type)
  
  cat("Cell_Type_array_\n")
  cat(str(Cell_Type_array))
  cat("\n")
  
  list_RESULTS_DEF<-list()
  
  # list_stats<-list()
  # list_stats_Vockley_REF<-list()
  
  
for(i in 1:length(LINEAGE_CLASSIF_DEF_array))
{
    
    LINEAGE_CLASSIF_DEF_array_sel<-LINEAGE_CLASSIF_DEF_array[i]
    
    cat("----------------------------------------------------------------->LINEAGE_CLASSIF_DEF_array_sel\n")
    cat(sprintf(as.character(LINEAGE_CLASSIF_DEF_array_sel)))
    cat("\n")
    
    
    REP_sel<-REP[which(REP$LINEAGE_CLASSIF_DEF == LINEAGE_CLASSIF_DEF_array_sel),]
    
    cat("REP_sel_\n")
    cat(str(REP_sel))
    cat("\n")
    
    
    
    list_RESULTS<-list()
    
    for(k in 1:length(Cell_Type_array))
    {
      
      Cell_Type_array_sel<-Cell_Type_array[k]
      
      cat("----------------------------------------------------------------->Cell_Type_array_sel\n")
      cat(sprintf(as.character(Cell_Type_array_sel)))
      cat("\n")
      
      REP_sel_Cell_Type_sel<-REP_sel[which(REP_sel$Cell_Type == Cell_Type_array_sel),]
      
      cat("REP_sel_Cell_Type_sel_\n")
      cat(str(REP_sel_Cell_Type_sel))
      cat("\n")
      
      
     
      
      
      if(dim(REP_sel_Cell_Type_sel)[1] >0)
      {
     
        #### abs_shift_del_effect and GWAS ----
        
      
        Max.abs_shift_del_effect_GWAS_effect_size<-as.data.frame(setDT(REP_sel_Cell_Type_sel)[, .(Max_abs_shift_del_effect=max(abs_shift_del_effect),
                                                                                                   Max_absolute_finemap_z=max(absolute_finemap_z)),
                                                                                               .(VAR,LINEAGE_CLASSIF_DEF)], stringsAsFactors=F)
        
        
        Max.abs_shift_del_effect_GWAS_effect_size$Quad_Max_abs_shift_del_effect<-(Max.abs_shift_del_effect_GWAS_effect_size$Max_abs_shift_del_effect)^2
        
        
        
        cat("Max.abs_shift_del_effect_GWAS_effect_size_\n")
        cat(str(Max.abs_shift_del_effect_GWAS_effect_size))
        cat("\n")
        
        
        
        #### linear model 
        
        linearModel_Max_abs_shift_del_effect <- summary(lm(Max_absolute_finemap_z ~ Max_abs_shift_del_effect, 
                                            data=Max.abs_shift_del_effect_GWAS_effect_size))
        
        # cat("linearModel_Max_abs_shift_del_effect\n")
        # cat(str(linearModel_Max_abs_shift_del_effect))
        # cat("\n")
        
        
        r.squared_linearModel_Max_abs_shift_del_effect<-round(linearModel_Max_abs_shift_del_effect$r.squared,3)
        adj.r.squared_linearModel_Max_abs_shift_del_effect<-round(linearModel_Max_abs_shift_del_effect$adj.r.squared,3)
        
        
        
        linearModel_Max_abs_shift_del_effect_coeffcient_df.m<-melt(linearModel_Max_abs_shift_del_effect$coefficients)
        
        colnames(linearModel_Max_abs_shift_del_effect_coeffcient_df.m)[which(colnames(linearModel_Max_abs_shift_del_effect_coeffcient_df.m)=="Var1")]<-"Terms"
        colnames(linearModel_Max_abs_shift_del_effect_coeffcient_df.m)[which(colnames(linearModel_Max_abs_shift_del_effect_coeffcient_df.m)=="Var2")]<-"Parameters"
        
        linearModel_Max_abs_shift_del_effect_coeffcient_df.m$Terms<-as.character(linearModel_Max_abs_shift_del_effect_coeffcient_df.m$Terms)
        linearModel_Max_abs_shift_del_effect_coeffcient_df.m$Parameters<-as.character(linearModel_Max_abs_shift_del_effect_coeffcient_df.m$Parameters)
        linearModel_Max_abs_shift_del_effect_coeffcient_df.m$Cell_Type<-Cell_Type_array_sel
        linearModel_Max_abs_shift_del_effect_coeffcient_df.m$LINEAGE_CLASSIF_DEF<-LINEAGE_CLASSIF_DEF_array_sel
        linearModel_Max_abs_shift_del_effect_coeffcient_df.m$adj.r.squared_linearModel_Max_abs_shift_del_effect<-adj.r.squared_linearModel_Max_abs_shift_del_effect
        
        
        
        cat("linearModel_Max_abs_shift_del_effect_coeffcient_df.m\n")
        cat(str(linearModel_Max_abs_shift_del_effect_coeffcient_df.m))
        cat("\n")
        
        #### quadratic model 
        
        quadraticModel_Max_abs_shift_del_effect <- summary(lm(Max_absolute_finemap_z ~ Max_abs_shift_del_effect+Quad_Max_abs_shift_del_effect, 
                                               data=Max.abs_shift_del_effect_GWAS_effect_size))
        
        # cat("quadraticModel_Max_abs_shift_del_effect\n")
        # cat(str(quadraticModel_Max_abs_shift_del_effect))
        # cat("\n")
        
        
        r.squared_quadraticModel_Max_abs_shift_del_effect<-round(quadraticModel_Max_abs_shift_del_effect$r.squared,3)
        adj.r.squared_quadraticModel_Max_abs_shift_del_effect<-round(quadraticModel_Max_abs_shift_del_effect$adj.r.squared,3)
        
        cat("adj.r.squared_quadraticModel_Max_abs_shift_del_effect_0\n")
        cat(str(adj.r.squared_quadraticModel_Max_abs_shift_del_effect))
        cat("\n")
        
        if(is.na(adj.r.squared_quadraticModel_Max_abs_shift_del_effect))
        {
          
          adj.r.squared_quadraticModel_Max_abs_shift_del_effect<-r.squared_quadraticModel_Max_abs_shift_del_effect
        }
        
        cat("adj.r.squared_quadraticModel_Max_abs_shift_del_effect_0\n")
        cat(str(adj.r.squared_quadraticModel_Max_abs_shift_del_effect))
        cat("\n")
        
        quadraticModel_Max_abs_shift_del_effect_coeffcient_df.m<-melt(quadraticModel_Max_abs_shift_del_effect$coefficients)
        
        # cat("quadraticModel_Max_abs_shift_del_effect_coeffcient_df.m_0\n")
        # cat(str(quadraticModel_Max_abs_shift_del_effect_coeffcient_df.m))
        # cat("\n")
        
        colnames(quadraticModel_Max_abs_shift_del_effect_coeffcient_df.m)[which(colnames(quadraticModel_Max_abs_shift_del_effect_coeffcient_df.m)=="Var1")]<-"Terms"
        colnames(quadraticModel_Max_abs_shift_del_effect_coeffcient_df.m)[which(colnames(quadraticModel_Max_abs_shift_del_effect_coeffcient_df.m)=="Var2")]<-"Parameters"
        
        quadraticModel_Max_abs_shift_del_effect_coeffcient_df.m$Terms<-as.character(quadraticModel_Max_abs_shift_del_effect_coeffcient_df.m$Terms)
        quadraticModel_Max_abs_shift_del_effect_coeffcient_df.m$Parameters<-as.character(quadraticModel_Max_abs_shift_del_effect_coeffcient_df.m$Parameters)
        quadraticModel_Max_abs_shift_del_effect_coeffcient_df.m$Cell_Type<-Cell_Type_array_sel
        quadraticModel_Max_abs_shift_del_effect_coeffcient_df.m$LINEAGE_CLASSIF_DEF<-LINEAGE_CLASSIF_DEF_array_sel
        quadraticModel_Max_abs_shift_del_effect_coeffcient_df.m$adj.r.squared_quadraticModel_Max_abs_shift_del_effect<-adj.r.squared_quadraticModel_Max_abs_shift_del_effect
        
        
        
        cat("quadraticModel_Max_abs_shift_del_effect_coeffcient_df.m\n")
        cat(str(quadraticModel_Max_abs_shift_del_effect_coeffcient_df.m))
        cat("\n")
        
        #### FLAG linear vs quadratic ----
        
        FLAG_linear_vs_quadratic<-sum(abs(unique(quadraticModel_Max_abs_shift_del_effect_coeffcient_df.m$adj.r.squared_quadraticModel_Max_abs_shift_del_effect)) > abs(unique(linearModel_Max_abs_shift_del_effect_coeffcient_df.m$adj.r.squared_linearModel_Max_abs_shift_del_effect)))
        
        cat("FLAG_linear_vs_quadratic\n")
        cat(sprintf(as.character(FLAG_linear_vs_quadratic)))
        cat("\n")
        
        
        if(is.na(FLAG_linear_vs_quadratic) == TRUE)
        {
          FLAG_NA<-1
          # quit(status = 1)
        }else{
          FLAG_NA<-0
          
        }
        # quit(status = 1)
        
        # if(LINEAGE_CLASSIF_DEF_array_sel == "Pleiotropic_variants" & Cell_Type_array_sel == "Kolf2")
        # {
        #   quit(status = 1)
        #   
        # }
        
        
        #### X and Y points
        
        step<-round((max(Max.abs_shift_del_effect_GWAS_effect_size$Max_abs_shift_del_effect) - min(Max.abs_shift_del_effect_GWAS_effect_size$Max_abs_shift_del_effect))/10,4)
        
        cat("step\n")
        cat(str(step))
        cat("\n")
        
        breaks.x<-round(seq(sum(min(Max.abs_shift_del_effect_GWAS_effect_size$Max_abs_shift_del_effect),-1*step),sum(max(Max.abs_shift_del_effect_GWAS_effect_size$Max_abs_shift_del_effect),step), by=step),4)
        
        cat("breaks.x\n")
        cat(str(breaks.x))
        cat("\n")
        
        
        labels.x<-as.character(round(breaks.x,3))
        
        cat("labels.x\n")
        cat(sprintf(as.character(labels.x)))
        cat("\n")
        
        
        # quit(status = 1)
        
        step<-round((max(Max.abs_shift_del_effect_GWAS_effect_size$Max_absolute_finemap_z) -0 )/10,4)
        
        cat("step\n")
        cat(str(step))
        cat("\n")
        
        breaks.y<-seq(0,sum(max(Max.abs_shift_del_effect_GWAS_effect_size$Max_absolute_finemap_z),step), by=step)
        labels.y<-as.character(round(breaks.y,0))
        
        cat("labels.y\n")
        cat(sprintf(as.character(labels.y)))
        cat("\n")
        
        if(FLAG_linear_vs_quadratic == 0 & FLAG_NA == 0)
        {
          
          graph_abs_Vockley_REF<-ggplot(Max.abs_shift_del_effect_GWAS_effect_size, 
                                        aes(x=Max_abs_shift_del_effect, 
                                            y=Max_absolute_finemap_z, 
                                            color=LINEAGE_CLASSIF_DEF)) +
            geom_point(size=3) +
            theme_bw()+
            theme(axis.title.y=element_text(size=24, family="sans"),
                  axis.title.x=element_text(size=24, family="sans"),
                  axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
                  axis.text.x=element_text(angle=0,size=18, color="black", family="sans"),
                  legend.title=element_text(size=16,color="black", family="sans"),
                  legend.text=element_text(size=12,color="black", family="sans"))+
            scale_x_continuous(name="Max_abs_shift_del_effect", breaks=breaks.x,labels=labels.x, limits=c(breaks.x[1],breaks.x[length(breaks.x)]))+
            scale_y_continuous(name="Max_absolute_finemap_z",breaks=breaks.y,labels=labels.y, limits=c(breaks.y[1],breaks.y[length(breaks.y)]))+
            scale_color_manual(values = c('#32A852','#6DB2EE','#553B68','#62D07F','#1877C9','#C9244B','#D45E85','#87447B','#D6B8E6'),
                               drop=F)+
            theme(legend.position="hidden",legend.title=element_blank(), legend.text = element_text(size=10))+
            guides(fill=guide_legend(nrow=2,byrow=TRUE))+
            stat_smooth(aes(y=Max_absolute_finemap_z, x=Max_abs_shift_del_effect), method = "lm", formula = y ~ x, se=F)+
            ggeasy::easy_center_title()
          
          linearModel_Max_abs_shift_del_effect_coeffcient_df.m$type<-"linear"
          
          
          list_RESULTS$Max_abs_shift_del_effect[[k]]<-linearModel_Max_abs_shift_del_effect_coeffcient_df.m
          
          
        }
        if(FLAG_linear_vs_quadratic > 0 & FLAG_NA == 0)
        {
          
          graph_abs_Vockley_REF<-ggplot(Max.abs_shift_del_effect_GWAS_effect_size, 
                                        aes(x=Max_abs_shift_del_effect, 
                                            y=Max_absolute_finemap_z, 
                                            color=LINEAGE_CLASSIF_DEF)) +
            geom_point(size=3) +
            theme_bw()+
            theme(axis.title.y=element_text(size=24, family="sans"),
                  axis.title.x=element_text(size=24, family="sans"),
                  axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
                  axis.text.x=element_text(angle=0,size=18, color="black", family="sans"),
                  legend.title=element_text(size=16,color="black", family="sans"),
                  legend.text=element_text(size=12,color="black", family="sans"))+
            scale_x_continuous(name="Max_abs_shift_del_effect", breaks=breaks.x,labels=labels.x, limits=c(breaks.x[1],breaks.x[length(breaks.x)]))+
            scale_y_continuous(name="Max_absolute_finemap_z",breaks=breaks.y,labels=labels.y, limits=c(breaks.y[1],breaks.y[length(breaks.y)]))+
            scale_color_manual(values = c('#32A852','#6DB2EE','#553B68','#62D07F','#1877C9','#C9244B','#D45E85','#87447B','#D6B8E6'),
                               drop=F)+
            theme(legend.position="hidden",legend.title=element_blank(), legend.text = element_text(size=10))+
            guides(fill=guide_legend(nrow=2,byrow=TRUE))+
            stat_smooth(aes(y=Max_absolute_finemap_z, x=Max_abs_shift_del_effect), method = "lm", formula = y ~ x+ I(x^2), se=F)+
            ggeasy::easy_center_title()
          
          quadraticModel_Max_abs_shift_del_effect_coeffcient_df.m$type<-"quadratic"
          
          
          list_RESULTS$Max_abs_shift_del_effect[[k]]<-quadraticModel_Max_abs_shift_del_effect_coeffcient_df.m
          
        }
        if(FLAG_NA > 0)
        {
          graph_abs_Vockley_REF<-ggplot(Max.abs_shift_del_effect_GWAS_effect_size, 
                                        aes(x=Max_abs_shift_del_effect, 
                                            y=Max_absolute_finemap_z, 
                                            color=LINEAGE_CLASSIF_DEF)) +
            geom_point(size=3) +
            theme_bw()+
            theme(axis.title.y=element_text(size=24, family="sans"),
                  axis.title.x=element_text(size=24, family="sans"),
                  axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
                  axis.text.x=element_text(angle=0,size=18, color="black", family="sans"),
                  legend.title=element_text(size=16,color="black", family="sans"),
                  legend.text=element_text(size=12,color="black", family="sans"))+
            scale_x_continuous(name="Max_abs_shift_del_effect", breaks=breaks.x,labels=labels.x, limits=c(breaks.x[1],breaks.x[length(breaks.x)]))+
            scale_y_continuous(name="Max_absolute_finemap_z",breaks=breaks.y,labels=labels.y, limits=c(breaks.y[1],breaks.y[length(breaks.y)]))+
            scale_color_manual(values = c('#32A852','#6DB2EE','#553B68','#62D07F','#1877C9','#C9244B','#D45E85','#87447B','#D6B8E6'),
                               drop=F)+
            theme(legend.position="hidden",legend.title=element_blank(), legend.text = element_text(size=10))+
            guides(fill=guide_legend(nrow=2,byrow=TRUE))+
            ggeasy::easy_center_title()
        }
        
        
        setwd(out)
        
        cat("\t")
        cat(sprintf(as.character(LINEAGE_CLASSIF_DEF_array_sel)))
        cat("\t")
        cat(sprintf(as.character(Cell_Type_array_sel)))
        cat("\t")
        
        
        tag_LINEAGE<-gsub(" ","_",LINEAGE_CLASSIF_DEF_array_sel)
        
        cat(sprintf(as.character(tag_LINEAGE)))
        cat("\n")
        
        svgname<-paste("GWAS_effect_size_shift_del_effect_",tag_LINEAGE,"_",Cell_Type_array_sel,".svg", sep='')
        makesvg = TRUE
        
        if (makesvg == TRUE)
        {
          ggsave(svgname, plot= graph_abs_Vockley_REF,
                 device="svg",
                 height=10, width=12)
        }
        
        
       
        
        
        # quit(status = 1)
        
        # if(is.na(FLAG_linear_vs_quadratic) == TRUE)
        # {
        #   
        #    quit(status = 1)
        # }
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        #### abs_shift_hdr_effect and GWAS ----
        
        
        Max.abs_shift_hdr_effect_GWAS_effect_size<-as.data.frame(setDT(REP_sel_Cell_Type_sel)[, .(Max_abs_shift_hdr_effect=max(abs_shift_hdr_effect),
                                                                                                  Max_absolute_finemap_z=max(absolute_finemap_z)),
                                                                                              .(VAR,LINEAGE_CLASSIF_DEF)], stringsAsFactors=F)
        
        
        Max.abs_shift_hdr_effect_GWAS_effect_size$Quad_Max_abs_shift_hdr_effect<-(Max.abs_shift_hdr_effect_GWAS_effect_size$Max_abs_shift_hdr_effect)^2
        
        
        
        cat("Max.abs_shift_hdr_effect_GWAS_effect_size_\n")
        cat(str(Max.abs_shift_hdr_effect_GWAS_effect_size))
        cat("\n")
        
        
        
        #### linear mohdr 
        
        linearMohdr_Max_abs_shift_hdr_effect <- summary(lm(Max_absolute_finemap_z ~ Max_abs_shift_hdr_effect, 
                                                           data=Max.abs_shift_hdr_effect_GWAS_effect_size))
        
        # cat("linearMohdr_Max_abs_shift_hdr_effect\n")
        # cat(str(linearMohdr_Max_abs_shift_hdr_effect))
        # cat("\n")
        
        
        r.squared_linearMohdr_Max_abs_shift_hdr_effect<-round(linearMohdr_Max_abs_shift_hdr_effect$r.squared,3)
        adj.r.squared_linearMohdr_Max_abs_shift_hdr_effect<-round(linearMohdr_Max_abs_shift_hdr_effect$adj.r.squared,3)
        
        
        
        linearMohdr_Max_abs_shift_hdr_effect_coeffcient_df.m<-melt(linearMohdr_Max_abs_shift_hdr_effect$coefficients)
        
        colnames(linearMohdr_Max_abs_shift_hdr_effect_coeffcient_df.m)[which(colnames(linearMohdr_Max_abs_shift_hdr_effect_coeffcient_df.m)=="Var1")]<-"Terms"
        colnames(linearMohdr_Max_abs_shift_hdr_effect_coeffcient_df.m)[which(colnames(linearMohdr_Max_abs_shift_hdr_effect_coeffcient_df.m)=="Var2")]<-"Parameters"
        
        linearMohdr_Max_abs_shift_hdr_effect_coeffcient_df.m$Terms<-as.character(linearMohdr_Max_abs_shift_hdr_effect_coeffcient_df.m$Terms)
        linearMohdr_Max_abs_shift_hdr_effect_coeffcient_df.m$Parameters<-as.character(linearMohdr_Max_abs_shift_hdr_effect_coeffcient_df.m$Parameters)
        linearMohdr_Max_abs_shift_hdr_effect_coeffcient_df.m$Cell_Type<-Cell_Type_array_sel
        linearMohdr_Max_abs_shift_hdr_effect_coeffcient_df.m$LINEAGE_CLASSIF_DEF<-LINEAGE_CLASSIF_DEF_array_sel
        linearMohdr_Max_abs_shift_hdr_effect_coeffcient_df.m$adj.r.squared_linearMohdr_Max_abs_shift_hdr_effect<-adj.r.squared_linearMohdr_Max_abs_shift_hdr_effect
        
        
        
        cat("linearMohdr_Max_abs_shift_hdr_effect_coeffcient_df.m\n")
        cat(str(linearMohdr_Max_abs_shift_hdr_effect_coeffcient_df.m))
        cat("\n")
        
        #### quadratic mohdr 
        
        quadraticMohdr_Max_abs_shift_hdr_effect <- summary(lm(Max_absolute_finemap_z ~ Max_abs_shift_hdr_effect+Quad_Max_abs_shift_hdr_effect, 
                                                              data=Max.abs_shift_hdr_effect_GWAS_effect_size))
        
        # cat("quadraticMohdr_Max_abs_shift_hdr_effect\n")
        # cat(str(quadraticMohdr_Max_abs_shift_hdr_effect))
        # cat("\n")
        
        
        r.squared_quadraticMohdr_Max_abs_shift_hdr_effect<-round(quadraticMohdr_Max_abs_shift_hdr_effect$r.squared,3)
        adj.r.squared_quadraticMohdr_Max_abs_shift_hdr_effect<-round(quadraticMohdr_Max_abs_shift_hdr_effect$adj.r.squared,3)
        
        cat("adj.r.squared_quadraticMohdr_Max_abs_shift_hdr_effect_0\n")
        cat(str(adj.r.squared_quadraticMohdr_Max_abs_shift_hdr_effect))
        cat("\n")
        
        if(is.na(adj.r.squared_quadraticMohdr_Max_abs_shift_hdr_effect))
        {
          
          adj.r.squared_quadraticMohdr_Max_abs_shift_hdr_effect<-r.squared_quadraticMohdr_Max_abs_shift_hdr_effect
        }
        
        cat("adj.r.squared_quadraticMohdr_Max_abs_shift_hdr_effect_0\n")
        cat(str(adj.r.squared_quadraticMohdr_Max_abs_shift_hdr_effect))
        cat("\n")
        
        quadraticMohdr_Max_abs_shift_hdr_effect_coeffcient_df.m<-melt(quadraticMohdr_Max_abs_shift_hdr_effect$coefficients)
        
        # cat("quadraticMohdr_Max_abs_shift_hdr_effect_coeffcient_df.m_0\n")
        # cat(str(quadraticMohdr_Max_abs_shift_hdr_effect_coeffcient_df.m))
        # cat("\n")
        
        colnames(quadraticMohdr_Max_abs_shift_hdr_effect_coeffcient_df.m)[which(colnames(quadraticMohdr_Max_abs_shift_hdr_effect_coeffcient_df.m)=="Var1")]<-"Terms"
        colnames(quadraticMohdr_Max_abs_shift_hdr_effect_coeffcient_df.m)[which(colnames(quadraticMohdr_Max_abs_shift_hdr_effect_coeffcient_df.m)=="Var2")]<-"Parameters"
        
        quadraticMohdr_Max_abs_shift_hdr_effect_coeffcient_df.m$Terms<-as.character(quadraticMohdr_Max_abs_shift_hdr_effect_coeffcient_df.m$Terms)
        quadraticMohdr_Max_abs_shift_hdr_effect_coeffcient_df.m$Parameters<-as.character(quadraticMohdr_Max_abs_shift_hdr_effect_coeffcient_df.m$Parameters)
        quadraticMohdr_Max_abs_shift_hdr_effect_coeffcient_df.m$Cell_Type<-Cell_Type_array_sel
        quadraticMohdr_Max_abs_shift_hdr_effect_coeffcient_df.m$LINEAGE_CLASSIF_DEF<-LINEAGE_CLASSIF_DEF_array_sel
        quadraticMohdr_Max_abs_shift_hdr_effect_coeffcient_df.m$adj.r.squared_quadraticMohdr_Max_abs_shift_hdr_effect<-adj.r.squared_quadraticMohdr_Max_abs_shift_hdr_effect
        
        
        
        cat("quadraticMohdr_Max_abs_shift_hdr_effect_coeffcient_df.m\n")
        cat(str(quadraticMohdr_Max_abs_shift_hdr_effect_coeffcient_df.m))
        cat("\n")
        
        #### FLAG linear vs quadratic ----
        
        FLAG_linear_vs_quadratic<-sum(abs(unique(quadraticMohdr_Max_abs_shift_hdr_effect_coeffcient_df.m$adj.r.squared_quadraticMohdr_Max_abs_shift_hdr_effect)) > abs(unique(linearMohdr_Max_abs_shift_hdr_effect_coeffcient_df.m$adj.r.squared_linearMohdr_Max_abs_shift_hdr_effect)))
        
        cat("FLAG_linear_vs_quadratic\n")
        cat(sprintf(as.character(FLAG_linear_vs_quadratic)))
        cat("\n")
        
        
        if(is.na(FLAG_linear_vs_quadratic) == TRUE)
        {
          FLAG_NA<-1
          # quit(status = 1)
        }else{
          FLAG_NA<-0
          
        }
        # quit(status = 1)
        
        # if(LINEAGE_CLASSIF_DEF_array_sel == "Pleiotropic_variants" & Cell_Type_array_sel == "Kolf2")
        # {
        #   quit(status = 1)
        #   
        # }
        
        
        #### X and Y points
        
        step<-round((max(Max.abs_shift_hdr_effect_GWAS_effect_size$Max_abs_shift_hdr_effect) - min(Max.abs_shift_hdr_effect_GWAS_effect_size$Max_abs_shift_hdr_effect))/10,4)
        
        cat("step\n")
        cat(str(step))
        cat("\n")
        
        breaks.x<-round(seq(sum(min(Max.abs_shift_hdr_effect_GWAS_effect_size$Max_abs_shift_hdr_effect),-1*step),sum(max(Max.abs_shift_hdr_effect_GWAS_effect_size$Max_abs_shift_hdr_effect),step), by=step),4)
        
        cat("breaks.x\n")
        cat(str(breaks.x))
        cat("\n")
        
        
        labels.x<-as.character(round(breaks.x,3))
        
        cat("labels.x\n")
        cat(sprintf(as.character(labels.x)))
        cat("\n")
        
        
        # quit(status = 1)
        
        step<-round((max(Max.abs_shift_hdr_effect_GWAS_effect_size$Max_absolute_finemap_z) -0 )/10,4)
        
        cat("step\n")
        cat(str(step))
        cat("\n")
        
        breaks.y<-seq(0,sum(max(Max.abs_shift_hdr_effect_GWAS_effect_size$Max_absolute_finemap_z),step), by=step)
        labels.y<-as.character(round(breaks.y,0))
        
        cat("labels.y\n")
        cat(sprintf(as.character(labels.y)))
        cat("\n")
        
        if(FLAG_linear_vs_quadratic == 0 & FLAG_NA == 0)
        {
          
          graph_abs_Vockley_REF<-ggplot(Max.abs_shift_hdr_effect_GWAS_effect_size, 
                                        aes(x=Max_abs_shift_hdr_effect, 
                                            y=Max_absolute_finemap_z, 
                                            color=LINEAGE_CLASSIF_DEF)) +
            geom_point(size=3) +
            theme_bw()+
            theme(axis.title.y=element_text(size=24, family="sans"),
                  axis.title.x=element_text(size=24, family="sans"),
                  axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
                  axis.text.x=element_text(angle=0,size=18, color="black", family="sans"),
                  legend.title=element_text(size=16,color="black", family="sans"),
                  legend.text=element_text(size=12,color="black", family="sans"))+
            scale_x_continuous(name="Max_abs_shift_hdr_effect", breaks=breaks.x,labels=labels.x, limits=c(breaks.x[1],breaks.x[length(breaks.x)]))+
            scale_y_continuous(name="Max_absolute_finemap_z",breaks=breaks.y,labels=labels.y, limits=c(breaks.y[1],breaks.y[length(breaks.y)]))+
            scale_color_manual(values = c('#32A852','#6DB2EE','#553B68','#62D07F','#1877C9','#C9244B','#D45E85','#87447B','#D6B8E6'),
                               drop=F)+
            theme(legend.position="hidden",legend.title=element_blank(), legend.text = element_text(size=10))+
            guides(fill=guide_legend(nrow=2,byrow=TRUE))+
            stat_smooth(aes(y=Max_absolute_finemap_z, x=Max_abs_shift_hdr_effect), method = "lm", formula = y ~ x, se=F)+
            ggeasy::easy_center_title()
          
          linearMohdr_Max_abs_shift_hdr_effect_coeffcient_df.m$type<-"linear"
          
          
          list_RESULTS$Max_abs_shift_hdr_effect[[k]]<-linearMohdr_Max_abs_shift_hdr_effect_coeffcient_df.m
          
          
        }
        if(FLAG_linear_vs_quadratic > 0 & FLAG_NA == 0)
        {
          
          graph_abs_Vockley_REF<-ggplot(Max.abs_shift_hdr_effect_GWAS_effect_size, 
                                        aes(x=Max_abs_shift_hdr_effect, 
                                            y=Max_absolute_finemap_z, 
                                            color=LINEAGE_CLASSIF_DEF)) +
            geom_point(size=3) +
            theme_bw()+
            theme(axis.title.y=element_text(size=24, family="sans"),
                  axis.title.x=element_text(size=24, family="sans"),
                  axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
                  axis.text.x=element_text(angle=0,size=18, color="black", family="sans"),
                  legend.title=element_text(size=16,color="black", family="sans"),
                  legend.text=element_text(size=12,color="black", family="sans"))+
            scale_x_continuous(name="Max_abs_shift_hdr_effect", breaks=breaks.x,labels=labels.x, limits=c(breaks.x[1],breaks.x[length(breaks.x)]))+
            scale_y_continuous(name="Max_absolute_finemap_z",breaks=breaks.y,labels=labels.y, limits=c(breaks.y[1],breaks.y[length(breaks.y)]))+
            scale_color_manual(values = c('#32A852','#6DB2EE','#553B68','#62D07F','#1877C9','#C9244B','#D45E85','#87447B','#D6B8E6'),
                               drop=F)+
            theme(legend.position="hidden",legend.title=element_blank(), legend.text = element_text(size=10))+
            guides(fill=guide_legend(nrow=2,byrow=TRUE))+
            stat_smooth(aes(y=Max_absolute_finemap_z, x=Max_abs_shift_hdr_effect), method = "lm", formula = y ~ x+ I(x^2), se=F)+
            ggeasy::easy_center_title()
          
          quadraticMohdr_Max_abs_shift_hdr_effect_coeffcient_df.m$type<-"quadratic"
          
          
          list_RESULTS$Max_abs_shift_hdr_effect[[k]]<-quadraticMohdr_Max_abs_shift_hdr_effect_coeffcient_df.m
          
        }
        if(FLAG_NA > 0)
        {
          graph_abs_Vockley_REF<-ggplot(Max.abs_shift_hdr_effect_GWAS_effect_size, 
                                        aes(x=Max_abs_shift_hdr_effect, 
                                            y=Max_absolute_finemap_z, 
                                            color=LINEAGE_CLASSIF_DEF)) +
            geom_point(size=3) +
            theme_bw()+
            theme(axis.title.y=element_text(size=24, family="sans"),
                  axis.title.x=element_text(size=24, family="sans"),
                  axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
                  axis.text.x=element_text(angle=0,size=18, color="black", family="sans"),
                  legend.title=element_text(size=16,color="black", family="sans"),
                  legend.text=element_text(size=12,color="black", family="sans"))+
            scale_x_continuous(name="Max_abs_shift_hdr_effect", breaks=breaks.x,labels=labels.x, limits=c(breaks.x[1],breaks.x[length(breaks.x)]))+
            scale_y_continuous(name="Max_absolute_finemap_z",breaks=breaks.y,labels=labels.y, limits=c(breaks.y[1],breaks.y[length(breaks.y)]))+
            scale_color_manual(values = c('#32A852','#6DB2EE','#553B68','#62D07F','#1877C9','#C9244B','#D45E85','#87447B','#D6B8E6'),
                               drop=F)+
            theme(legend.position="hidden",legend.title=element_blank(), legend.text = element_text(size=10))+
            guides(fill=guide_legend(nrow=2,byrow=TRUE))+
            ggeasy::easy_center_title()
        }
        
        
        setwd(out)
        
        cat("\t")
        cat(sprintf(as.character(LINEAGE_CLASSIF_DEF_array_sel)))
        cat("\t")
        cat(sprintf(as.character(Cell_Type_array_sel)))
        cat("\t")
        
        
        tag_LINEAGE<-gsub(" ","_",LINEAGE_CLASSIF_DEF_array_sel)
        
        cat(sprintf(as.character(tag_LINEAGE)))
        cat("\n")
        
        svgname<-paste("GWAS_effect_size_shift_hdr_effect_",tag_LINEAGE,"_",Cell_Type_array_sel,".svg", sep='')
        makesvg = TRUE
        
        if (makesvg == TRUE)
        {
          ggsave(svgname, plot= graph_abs_Vockley_REF,
                 device="svg",
                 height=10, width=12)
        }
        
        
        
        
        
        # quit(status = 1)
        
        # if(is.na(FLAG_linear_vs_quadratic) == TRUE)
        # {
        #   
        #    quit(status = 1)
        # }
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
      }#dim(REP_sel_Cell_Type_sel)[1] >0
      
      }#k  for(k in 1:length(LINEAGE_CLASSIF_DEF_array))
    
   
    
   
    
    Max_abs_shift_del_effect_stats = unique(as.data.frame(data.table::rbindlist(list_RESULTS$Max_abs_shift_del_effect, fill = T)))
    
    if(dim(Max_abs_shift_del_effect_stats)[1] >0)
    {
      Max_abs_shift_del_effect_stats$LINEAGE_CLASSIF_DEF<-LINEAGE_CLASSIF_DEF_array_sel
      
      
      cat("Max_abs_shift_del_effect_stats\n")
      cat(str(Max_abs_shift_del_effect_stats))
      cat("\n")
      
      list_RESULTS_DEF$Max_abs_shift_del_effect[[i]]<-Max_abs_shift_del_effect_stats
      
    }
    
    Max_abs_shift_hdr_effect_stats = unique(as.data.frame(data.table::rbindlist(list_RESULTS$Max_abs_shift_hdr_effect, fill = T)))
    
    if(dim(Max_abs_shift_hdr_effect_stats)[1] >0)
    {
      Max_abs_shift_hdr_effect_stats$LINEAGE_CLASSIF_DEF<-LINEAGE_CLASSIF_DEF_array_sel
      
      
      cat("Max_abs_shift_hdr_effect_stats\n")
      cat(str(Max_abs_shift_hdr_effect_stats))
      cat("\n")
      
      list_RESULTS_DEF$Max_abs_shift_hdr_effect[[i]]<-Max_abs_shift_hdr_effect_stats
      
    }
    
    
   
    
}#i LINEAGE_CLASSIF_DEF_array
    
 
  
  Max_abs_shift_del_effect_stats_DEF = unique(as.data.frame(data.table::rbindlist(list_RESULTS_DEF$Max_abs_shift_del_effect, fill = T)))
  
  cat("Max_abs_shift_del_effect_stats_DEF\n")
  cat(str(Max_abs_shift_del_effect_stats_DEF))
  cat("\n")
  
  Max_abs_shift_hdr_effect_stats_DEF = unique(as.data.frame(data.table::rbindlist(list_RESULTS_DEF$Max_abs_shift_hdr_effect, fill = T)))
  
  cat("Max_abs_shift_hdr_effect_stats_DEF\n")
  cat(str(Max_abs_shift_hdr_effect_stats_DEF))
  cat("\n")
  
  # ###########################################################################################################################################################################################################################################################
  # quit(status=1)

  #### SAVE ----
  
  setwd(out)
  
 
  
  filename_1<-paste("GWAS_Max_abs_shift_del_effect_effect_size_plots_LM_",type,".tsv", sep='')
  
  write.table(file=filename_1,
              Max_abs_shift_del_effect_stats_DEF,
              sep="\t",
              quote = F,
              row.names = F)
 
  filename_1<-paste("GWAS_Max_abs_shift_hdr_effect_effect_size_plots_LM_",type,".tsv", sep='')
  
  write.table(file=filename_1,
              Max_abs_shift_hdr_effect_stats_DEF,
              sep="\t",
              quote = F,
              row.names = F)
  
  #quit(status=1)
  
 
  
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
  
  
  graph_function(opt)
 
  
}
  
  
  
 

###########################################################################

system.time( main() )
