
# library("desc",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
# library("ps",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
# library("usethis",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
# library("withr",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
# library("devtools",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")



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
library("vroom", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
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
library("FNN", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("gridExtra", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("egg", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")


library("rgenie", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
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
suppressMessages(library("ggdendro", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))

opt = NULL

options(warn=1)

rgenie_script = function(option_list)
{
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("type\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("out\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### READ regions ----
  
  regions = readr::read_tsv(opt$regions)
  
  cat("regions\n")
  cat(str(regions))
  cat("\n")
  
 
  #### READ QC_PASS ----
  
  QC_PASS = opt$QC_PASS
  QC_PASS<-unlist(strsplit(QC_PASS, split=","))
  
  cat("------------------>QC_PASS:\t")
  cat(sprintf(as.character(QC_PASS)))
  cat("\n")
  
  #### READ dropouts_ALL_PLOTS ----
  
  dropouts_ALL_PLOTS = opt$dropouts_ALL_PLOTS
  dropouts_ALL_PLOTS<-unlist(strsplit(dropouts_ALL_PLOTS, split=","))
  
  cat("------------------>dropouts_ALL_PLOTS:\t")
  cat(sprintf(as.character(dropouts_ALL_PLOTS)))
  cat("\n")
  
  #### READ and transform required_match_left_THRESHOLD ----
  
  required_match_left_THRESHOLD = opt$required_match_left_THRESHOLD
  
  cat("required_match_left_THRESHOLD\n")
  cat(sprintf(as.character(required_match_left_THRESHOLD)))
  cat("\n")
  
  
  #### READ and transform required_match_right_THRESHOLD ----
  
  required_match_right_THRESHOLD = opt$required_match_right_THRESHOLD
  
  cat("required_match_right_THRESHOLD\n")
  cat(sprintf(as.character(required_match_right_THRESHOLD)))
  cat("\n")
  
  
  #### Read input_TABLE_subset & TRANSFORM ----
  
  setwd(out)
  
  input_TABLE<-read.table(opt$input, header=T, sep="\t", stringsAsFactors = F)
  
  cat("input_TABLE_\n")
  str(input_TABLE)
  cat("\n")
  
  input_TABLE_subset<-input_TABLE[which(input_TABLE$HGNC%in%QC_PASS |
                                          input_TABLE$HGNC%in%dropouts_ALL_PLOTS),]
  
  cat("input_TABLE_subset_\n")
  str(input_TABLE_subset)
  cat("\n")
  
  
  
 
  #### READ replicates ----
  
  replicates = readr::read_tsv(opt$replicates)
  
  cat("replicates\n")
  cat(str(replicates))
  cat("\n")
  
  # replicates$type<-"NA"
  # 
  # #### PATCH
  # 
  # replicates$type[grep("_C[0-9]+",replicates$replicate)]<-"cDNA"
  # replicates$type[grep("_G[0-9]+",replicates$replicate)]<-"gDNA"
  # 
  # cat("replicates\n")
  # cat(str(replicates))
  # cat("\n")
  
 
  replicates_PRE<-as.data.frame(replicates)
  
  cat("replicates_PRE\n")
  cat(str(replicates_PRE))
  cat("\n")
  
  
 levels_type<-levels(as.factor(replicates$type)) 
 
 cat("levels_type\n")
 cat(sprintf(as.character(levels_type)))
 cat("\n")
 
 
  ACCEPTED_types<-c('cDNA','gDNA','ATAC')
  
  check_levels_type<-which(levels_type%in%ACCEPTED_types)
  
  cat("check_levels_type\n")
  cat(sprintf(as.character(check_levels_type)))
  cat("\n")
  
  if(length(check_levels_type) == 0)
  {
    cat("unaccepted_format\n")
    
    indx.cDNA<-grep("[cC]", replicates$type)
    
    check.cDNA<-unique(replicates$type[indx.cDNA])
    
    cat("check.cDNA\n")
    cat(sprintf(as.character(check.cDNA)))
    cat("\n")
    
    replicates$type[indx.cDNA]<-'cDNA'
    
    
    indx.gDNA<-grep("[gG]", replicates$type)
    
    check.gDNA<-unique(replicates$type[indx.gDNA])
    
    cat("check.gDNA\n")
    cat(sprintf(as.character(check.gDNA)))
    cat("\n")
    
    replicates$type[indx.gDNA]<-'gDNA'
    
    # quit(status=1)
    
  }else{
    
    
    
    
  }
  
  
 
  
  
  replicates_OUT<-replicates_PRE[-which(replicates_PRE$name%in%QC_PASS |
                                          replicates_PRE$name%in%dropouts_ALL_PLOTS),]
  
  cat("replicates_OUT\n")
  cat(str(replicates_OUT))
  cat("\n")
  
  
  ####  subset replicates in QC_PASS ----
  
  replicates = replicates %>% dplyr::filter(name %in% QC_PASS |
                                              name %in% dropouts_ALL_PLOTS)
  
  cat("replicates_after_QC_PASS\n")
  cat(str(replicates))
  cat("\n")
  
  
  
  regions = regions[which(regions$name%in%replicates$name),]
  
  # regions$name<-regions$sequence_name
  
  cat("regions_AFTER\n")
  cat(str(regions))
  cat("\n")
  
  # quit(status = 1)
  
  #### GREP analysis ----
  
  cat("grep_analysis\40\40\n")
  cat("\n")
  
  grep_results = grep_analysis(regions,
                               replicates,
                               required_match_left = required_match_left_THRESHOLD,
                               required_match_right = required_match_right_THRESHOLD,
                               min_mapq = 0,
                               quiet = F)
  
  
  cat("THE END\n")
  cat("\n")
  # 
  # cat("grep_results\n")
  # cat(str(grep_results))
  # 
  # 
  # 
  
  #### del analysis ----
  # 
  cat("del_analysis\40\40\n")
  cat("\n")
 
 
  
  del_results = rgenie::alignment_analysis(regions,
                                           replicates,
                                           required_match_left = required_match_left_THRESHOLD,
                                           required_match_right = required_match_right_THRESHOLD,
                                           crispr_del_window = 50,
                                           min_mapq = 0,
                                           max_mismatch_frac = 0.05,
                                           min_aligned_bases = 50,
                                           exclude_multiple_deletions = FALSE,
                                           exclude_nonspanning_reads = TRUE,
                                           allele_profile = FALSE,
                                           del_span_start = -2*required_match_left_THRESHOLD,
                                           del_span_end = 2*required_match_right_THRESHOLD,
                                           quiet = F)
  
  cat("\n")
  
  cat("THE END\n")
  cat("\n")

  # cat("del_results\n")
  # cat(str(del_results))
  # cat("\n")

  
  # quit(status=1)
    
 # 
 
  #### LOOP ----
  
  myplots_TOTAL<-list()
  
  my_stats<-list()
  
  array_extension<-length(grep_results)
  
  cat("array_extension\n")
  cat(str(array_extension))
  cat("\n")
  
  
  VARS<-unique(input_TABLE_subset$VAR)
  
  cat("VARS\n")
  cat(str(VARS))
  cat("\n")

  # for(i in 1:array_extension)
  # {
  for(i in 1:length(VARS))
  {
    cat("------->\t")
    VAR_sel<-VARS[i]
    
    cat(sprintf(as.character(VAR_sel)))
    cat("\t")
    
    HGNC<-input_TABLE_subset$HGNC[input_TABLE_subset$VAR == VAR_sel]
    
    cat(sprintf(as.character(HGNC)))
    cat("\t")
    
    amplicon_reconstructed<-paste(VAR_sel,paste(HGNC,"haplotype","REF",sep="_"),sep="__")
    
    cat(sprintf(as.character(amplicon_reconstructed)))
    cat("\n")
    
    pivot_point<-QC_PASS[which(QC_PASS == HGNC)]
    
    if(length(pivot_point) >0)
    {
      for(k in 1:array_extension)
      {
        HGNC_object<-as.character(grep_results[[k]]$region$name)
        
        
        
        sequence_name_object<-as.character(grep_results[[k]]$region$sequence_name)
        
        
        
        if(sequence_name_object == amplicon_reconstructed)
        {
          cat("--Hello_world-->\n")
          
          # cat("----sequence_name_object--->\t")
          
          cat(sprintf(as.character(HGNC_object)))
          cat("\t")
          cat(sprintf(as.character(sequence_name_object)))
          cat("\n")
          
          grep_summary_plot_instance<-rgenie::grep_summary_plot(grep_results[[k]])
          
          # cat("grep_summary_plot_instance\n")
          # cat(str(grep_summary_plot_instance))
          # cat("\n")
          
          object_grep<-grep_results[[k]]
          
          # cat("object_grep\n")
          # cat(str(object_grep))
          # cat("\n")
          
          myplots_TOTAL$grep_analysis[[VAR_sel]][[HGNC]]<-grep_summary_plot_instance
          
          setwd(out)
          pdf(paste("test_grep_",VAR_sel,"_",HGNC, '.pdf',sep=""))
          print(grep_summary_plot_instance)
          dev.off()
          
          object_del<-del_results[[k]]
          
          # cat("object_del\n")
          # cat(str(object_del))
          cat("\n")
          
          
          all_plots = rgenie::alignment_analysis_plots(del_results[[k]],
                                                       opts = genie_plot_options(),
                                                       variance_components_plot = FALSE,
                                                       power_plots = FALSE)
          
          
          setwd(out)
          pdf(paste("test_alignment_",VAR_sel,"_",HGNC, '.pdf',sep=""))
          print(all_plots)
          dev.off()
          
          
          # cat("all_plots\n")
          # cat(str(all_plots))
          # cat("\n")
          
          # myplots_TOTAL$experiment_summary_plot[[HGNC]]<-experiment_summary_plot_instance
          
          
          
          
          myplots_TOTAL$deletion_analysis[[VAR_sel]][[HGNC]]<-all_plots
          
          
          a.df<-as.data.frame(del_results[[k]]$region_stats, stringsAsFactors = F)
          
          a.df$VAR<-VAR_sel
          a.df$HGNC<-HGNC
          
          
          cat("a.df\n")
          cat(str(a.df))
          cat("\n")
          
          my_stats$region_stats[[k]]<-a.df
          
          b.df<-as.data.frame(del_results[[k]]$replicate_stats, stringsAsFactors = F)
          
          b.df$VAR<-VAR_sel
          b.df$HGNC<-HGNC
          
          cat("b.df\n")
          cat(str(b.df))
          cat("\n")
          
          my_stats$replicate_stats[[k]]<-b.df
          
          
          cat("THE END\n")
          cat("\n")
          
          # quit(status=1)
        }
        
      }# k
      
    }#length(pivot_point) >0
    else{
      
      pivot_point2<-dropouts_ALL_PLOTS[which(dropouts_ALL_PLOTS == HGNC)]
      
      if(length(pivot_point2) >0)
      {
        for(k in 1:array_extension)
        {
          HGNC_object<-as.character(grep_results[[k]]$region$name)
          
          
          
          sequence_name_object<-as.character(grep_results[[k]]$region$sequence_name)
          
          
          
          if(sequence_name_object == amplicon_reconstructed)
          {
            cat("--Hello_world-2->\n")
            
            # cat("----sequence_name_object--->\t")
            
            cat(sprintf(as.character(HGNC_object)))
            cat("\t")
            cat(sprintf(as.character(sequence_name_object)))
            cat("\n")
            
            grep_summary_plot_instance<-rgenie::grep_summary_plot(grep_results[[k]])
            
            # cat("grep_summary_plot_instance\n")
            # cat(str(grep_summary_plot_instance))
            # cat("\n")
            
            object_grep<-grep_results[[k]]
            
            # cat("object_grep\n")
            # cat(str(object_grep))
            # cat("\n")
            
            myplots_TOTAL$grep_analysis[[VAR_sel]][[HGNC]]<-grep_summary_plot_instance
            
            setwd(out)
            pdf(paste("test_grep_",VAR_sel,"_",HGNC, '.pdf',sep=""))
            print(grep_summary_plot_instance)
            dev.off()
            
            object_del<-del_results[[k]]
            
            # cat("object_del\n")
            # cat(str(object_del))
            cat("\n")
            
            
            
            del_profile_plot_instance<-rgenie::deletion_profile_plot(object_del,
                                          viewing_window = 40)
            
            deletion_alleles_plot_instance<-rgenie::deletion_alleles_plot(object_del,
                                          viewing_window = 40,
                                          color_by="window")
            
            ####
            
            cat("alignment_summary_plot\n")
            
            alignment_summary_plot_instance<-rgenie::alignment_summary_plot(object_del)
            
            
            cat("replicate_summary_plot\n")
            
            replicate_summary_plot_instance<-rgenie::replicate_summary_plot(object_del,
                                           outlier_threshold = NA)
            cat("replicate_qc_plot\n")
            
            
            replicate_qc_plot_instance<-rgenie::replicate_qc_plot(object_del,
                                      outlier_threshold = NA)
            
            cat("allele_effect_plot\n")
            
            
            allele_effect_plot_instance<-rgenie::allele_effect_plot(object_del,
                                       viewing_window = 40,
                                       max_alleles = 40)
            
            list_graphs<-list()
            
            list_graphs$alignment_summary<-alignment_summary_plot_instance
            list_graphs$del_profile<-del_profile_plot_instance
            list_graphs$deletion_alleles<-deletion_alleles_plot_instance
            list_graphs$replicate_summary<-replicate_summary_plot_instance
            list_graphs$replicate_qc<-replicate_qc_plot_instance
            list_graphs$allele_effect<-allele_effect_plot_instance
            
            
         #   $alignment_summary
            
         
            
            
            
            # all_plots = rgenie::alignment_analysis_plots(del_results[[k]],
            #                                              opts = genie_plot_options(),
            #                                              variance_components_plot = FALSE,
            #                                              power_plots = FALSE)
            # 
            
            # cat("PRINT\n")
            # 
            # 
            # setwd(out)
            # pdf(paste("test_alignment_",VAR_sel,"_",HGNC, '.pdf',sep=""))
            # print(del_profile_plot_instance)
            # print(deletion_alleles_plot_instance)
            # print(replicate_summary_plot_instance)
            # print(replicate_qc_plot_instance)
            # print(allele_effect_plot_instance)
            # dev.off()
            
            
            # cat("all_plots\n")
            # cat(str(all_plots))
            # cat("\n")
            
            # myplots_TOTAL$experiment_summary_plot[[HGNC]]<-experiment_summary_plot_instance
            
            
            
            
            myplots_TOTAL$deletion_analysis[[VAR_sel]][[HGNC]]<-list_graphs
            
            
            a.df<-as.data.frame(del_results[[k]]$region_stats, stringsAsFactors = F)
            
            a.df$VAR<-VAR_sel
            a.df$HGNC<-HGNC
            
            
            cat("a.df\n")
            cat(str(a.df))
            cat("\n")
            
            my_stats$region_stats[[k]]<-a.df
            
            b.df<-as.data.frame(del_results[[k]]$replicate_stats, stringsAsFactors = F)
            
            b.df$VAR<-VAR_sel
            b.df$HGNC<-HGNC
            
            cat("b.df\n")
            cat(str(b.df))
            cat("\n")
            
            my_stats$replicate_stats[[k]]<-b.df
            
            
            cat("THE END\n")
            cat("\n")
            
            # quit(status=1)
          }
          
        }# k
        
      }# length(pivot_point2) >0
      else{
        
        cat("HARD DROPOUT\t")
        cat(sprintf(as.character(VAR_sel)))
        cat("\t")
        cat(sprintf(as.character(HGNC)))
        cat("\n")
      }
      
    }
   

  }#i VARS LOOP
  
  
  #### Gathering ----
  
  region_stats_df = unique(as.data.frame(data.table::rbindlist(my_stats$region_stats, fill = T)))
  
  
  cat("region_stats_df\n")
  cat(str(region_stats_df))
  cat("\n")
  
  replicate_stats_df = unique(as.data.frame(data.table::rbindlist(my_stats$replicate_stats, fill = T)))
  
  
  cat("replicate_stats_df\n")
  cat(str(replicate_stats_df))
  cat("\n")
  
  
  
  #### SAVE RDS ----
  
  setwd(out)
  
  type<-paste(type,required_match_left_THRESHOLD,required_match_right_THRESHOLD, sep="_")
  
  
  filename=paste(type,".rds", sep='')

  saveRDS(myplots_TOTAL, file = filename)
  
  
  filename=paste(type,"_","region_stats",".tsv", sep='')
  
  write.table(region_stats_df, file=filename,sep="\t", row.names=F, quote=F)
  
  filename=paste(type,"_","replicate_stats",".tsv", sep='')
  
  write.table(replicate_stats_df, file=filename,sep="\t", row.names=F, quote=F)
  
  
  #### Print PDF ----
  
  # pdfname<-paste(type,"_Graphs",".pdf", sep='')
  # makepdf = TRUE
  # 
  # if (makepdf == TRUE)
  # {
  #   pdf ( pdfname , height=10, width=12)
  # }
  # 
  # 
  # print(plots)
  # 
  # 
  # 
  # if (makepdf == TRUE)
  # {
  #   dev.off()
  # }
  
}

Report_printer = function(option_list)
{
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("type\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("out\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### READ QC_PASS ----
  
  QC_PASS = opt$QC_PASS
  QC_PASS<-unlist(strsplit(QC_PASS, split=","))
  
  cat("------------------>QC_PASS:\t")
  cat(sprintf(as.character(QC_PASS)))
  cat("\n")
  
  #### READ dropouts_ALL_PLOTS ----
  
  dropouts_ALL_PLOTS = opt$dropouts_ALL_PLOTS
  dropouts_ALL_PLOTS<-unlist(strsplit(dropouts_ALL_PLOTS, split=","))
  
  cat("------------------>dropouts_ALL_PLOTS:\t")
  cat(sprintf(as.character(dropouts_ALL_PLOTS)))
  cat("\n")
  
  #### READ and transform required_match_left_THRESHOLD ----
  
  required_match_left_THRESHOLD = opt$required_match_left_THRESHOLD
  
  cat("required_match_left_THRESHOLD\n")
  cat(sprintf(as.character(required_match_left_THRESHOLD)))
  cat("\n")
  
  
  #### READ and transform required_match_right_THRESHOLD ----
  
  required_match_right_THRESHOLD = opt$required_match_right_THRESHOLD
  
  cat("required_match_right_THRESHOLD\n")
  cat(sprintf(as.character(required_match_right_THRESHOLD)))
  cat("\n")
  
  #### Read input_TABLE_subset & TRANSFORM ----
  
  setwd(out)
  
  input_TABLE<-read.table(opt$input, header=T, sep="\t", stringsAsFactors = F)
  
  cat("input_TABLE_\n")
  str(input_TABLE)
  cat("\n")
  
  input_TABLE_subset<-input_TABLE[which(input_TABLE$HGNC%in%QC_PASS |
                                          input_TABLE$HGNC%in%dropouts_ALL_PLOTS),]
  
  cat("input_TABLE_subset_\n")
  str(input_TABLE_subset)
  cat("\n")
  
  #### read TOTAL ----
  
  type<-paste(type,required_match_left_THRESHOLD,required_match_right_THRESHOLD, sep="_")
  
  
  filename=paste(type,".rds", sep='')
  List_TOTAL<-readRDS(file = filename)
  
  
 
  #### PRINTING LOOP ----
  
  
  
  VARS<-unique(input_TABLE_subset$VAR)
  
  cat("VARS\n")
  cat(str(VARS))
  cat("\n")
  
  setwd(out)
  
  pdfname<-paste(type,"_Graphs",".pdf", sep='')
  makepdf = TRUE
  
  if (makepdf == TRUE)
  {
    pdf ( pdfname , height=10, width=12)
  }
  
  
  for(i in 1:length(VARS))
  {
    cat("------->\n")
    VAR_sel<-VARS[i]
    
    cat(sprintf(as.character(VAR_sel)))
    cat("\t")
    
    HGNC<-input_TABLE_subset$HGNC[input_TABLE_subset$VAR == VAR_sel]
    
    cat(sprintf(as.character(HGNC)))
    cat("\t")
    
    amplicon_reconstructed<-paste(VAR_sel,paste(HGNC,"haplotype","REF",sep="_"),sep="__")
    
    cat(sprintf(as.character(amplicon_reconstructed)))
    cat("\n")
    
   
    graph_deletion_analysis<-List_TOTAL$deletion_analysis[[VAR_sel]][[HGNC]]
    
    graph_grep_analysis<-List_TOTAL$grep_analysis[[VAR_sel]][[HGNC]]
    
    
    print(graph_grep_analysis)
    
    print(graph_deletion_analysis)
    
   
    
    
    
    
  }#i VAR_sel
  
  
  
  if (makepdf == TRUE)
  {
    dev.off()
  }
  
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
    make_option(c("--regions"), type="character", default=NULL,
                metavar="FILE.tsv",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--QC_PASS"), type="character", default=NULL, 
                metavar="FILE.tsv", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--dropouts_ALL_PLOTS"), type="character", default=NULL, 
                metavar="FILE.tsv", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--required_match_left_THRESHOLD"), type="numeric", default=NULL, 
                metavar="FILE.tsv", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--required_match_right_THRESHOLD"), type="numeric", default=NULL, 
                metavar="FILE.tsv", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--input"), type="character", default=NULL, 
                metavar="FILE.txt", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--replicates"), type="character", default=NULL, 
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
                        --regions FILE.tsv
                        --mismatch_CONDENSED FILE.tsv 
                        --TOTAL_READS FILE.tsv 
                        --DEMULTIPLEX_RESULT FILE.tsv
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)

  rgenie_script(opt)
  Report_printer(opt)
 
  
}


###########################################################################

system.time( main() )
