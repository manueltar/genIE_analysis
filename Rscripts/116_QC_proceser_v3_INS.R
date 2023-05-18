

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
suppressMessages(library("svglite", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("ggeasy", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("sandwich", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
library("digest", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")

library("ggforce", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
suppressMessages(library("cowplot", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))


##### Taken from https://www.biostars.org/p/9335/


opt = NULL



Graph_INS = function(option_list)
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
  
  #### READ and transform COUNT_MATRIX ----
  
  setwd(out)
  
  
  filename_3<-paste(type,"_COUNT_matrix",".rds", sep='')
  
  
  COUNT_matrix = readRDS(file=filename_3)
  
  COUNT_matrix$name<-gsub("-.+","",COUNT_matrix$sample)
  
  levels_replica_order_factor<-levels(COUNT_matrix$replica_order_factor)
  
  
  cat("COUNT_matrix\n")
  cat(str(COUNT_matrix))
  cat("\n")
  #quit(status = 1)
  
  

  #### READ and transform mismatch_CONDENSED ----
  
  mismatch_CONDENSED = read.table(opt$mismatch_CONDENSED, sep="\t", stringsAsFactors = F, header = T)
  
  colnames(mismatch_CONDENSED)[which(colnames(mismatch_CONDENSED) == "Sample_Name")]<-"sample"
  
  
  mismatch_CONDENSED$name<-gsub("-.+","",mismatch_CONDENSED$sample)
  mismatch_CONDENSED$sample<-gsub("-","_",mismatch_CONDENSED$sample)
  
  
  cat("mismatch_CONDENSED\n")
  cat(str(mismatch_CONDENSED))
  cat("\n")
  
  mismatch_CONDENSED$Distance_editing<-as.integer(gsub("^NM:[^:]+:","",mismatch_CONDENSED$NM))
  
  cat("mismatch_CONDENSED\n")
  cat(str(mismatch_CONDENSED))
  cat("\n")
  
  
  #mismatch_CONDENSED_subset<-mismatch_CONDENSED[which(mismatch_CONDENSED$Distance_editing <= 5),]
  
  mismatch_CONDENSED_subset<-mismatch_CONDENSED
  cat("mismatch_CONDENSED_subset\n")
  cat(str(mismatch_CONDENSED_subset))
  cat("\n")
  
  ### Discard Soft & Hard clipped
  
  indx.dep<-grep("S[0-9]+",mismatch_CONDENSED_subset$CIGAR)
  
  cat("indx.dep\n")
  cat(str(indx.dep))
  cat("\n")
  
  indx.dep2<-grep("[0-9]+S",mismatch_CONDENSED_subset$CIGAR)
  
  cat("indx.dep2\n")
  cat(str(indx.dep2))
  cat("\n")
  
  indx.dep3<-grep("H[0-9]+",mismatch_CONDENSED_subset$CIGAR)
  
  cat("indx.dep3\n")
  cat(str(indx.dep3))
  cat("\n")
  
  indx.dep4<-grep("[0-9]+H",mismatch_CONDENSED_subset$CIGAR)
  
  cat("indx.dep4\n")
  cat(str(indx.dep4))
  cat("\n")
  
  indx.dep5<-grep("D[0-9]+",mismatch_CONDENSED_subset$CIGAR)
  
  cat("indx.dep5\n")
  cat(str(indx.dep5))
  cat("\n")
  
  indx.dep6<-grep("[0-9]+D",mismatch_CONDENSED_subset$CIGAR)
  
  cat("indx.dep6\n")
  cat(str(indx.dep6))
  cat("\n")
  
  indx.dep7<-grep("I[0-9]+",mismatch_CONDENSED_subset$CIGAR)
  
  cat("indx.dep7\n")
  cat(str(indx.dep7))
  cat("\n")
  
  indx.dep8<-grep("[0-9]+I",mismatch_CONDENSED_subset$CIGAR)
  
  cat("indx.dep8\n")
  cat(str(indx.dep8))
  cat("\n")
  
  indx.dep.DEF<-unique(c(indx.dep7,indx.dep8))
  
  mismatch_CONDENSED_subset<-mismatch_CONDENSED_subset[indx.dep.DEF,]
  
  cat("mismatch_CONDENSED_subset\n")
  cat(str(mismatch_CONDENSED_subset))
  cat("\n")
  
  if(dim(mismatch_CONDENSED_subset)[1] >0)
  {
    cat("mismatch_CONDENSED_subset$CIGAR\n")
    cat(sprintf(as.character(names(summary(as.factor(mismatch_CONDENSED_subset$CIGAR))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(mismatch_CONDENSED_subset$CIGAR)))))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(mismatch_CONDENSED_subset$Distance_editing))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(mismatch_CONDENSED_subset$Distance_editing)))))
    cat("\n")
    
    ACCEPTED_CIGAR<-names(summary(as.factor(mismatch_CONDENSED_subset$CIGAR)))
    
    cat("ACCEPTED_CIGAR\n")
    cat(sprintf(as.character(ACCEPTED_CIGAR)))
    cat("\n")
    
    
    #### READ and transform TOTAL_READS ----
    
    TOTAL_READS = read.table(opt$TOTAL_READS, sep="\t", stringsAsFactors = F, header = T)
    
    colnames(TOTAL_READS)[which(colnames(TOTAL_READS) == "Sample_Name")]<-"sample"
    TOTAL_READS$name<-gsub("-.+","",TOTAL_READS$sample)
    
    TOTAL_READS$sample<-gsub("-","_",TOTAL_READS$sample)
    
    
    
    cat("TOTAL_READS\n")
    cat(str(TOTAL_READS))
    cat("\n")
    
    
    #### READ and transform DEMULTIPLEX_RESULT ---- 
    
    DEMULTIPLEX_RESULT = read.table(opt$DEMULTIPLEX_RESULT, sep="\t", stringsAsFactors = F, header = T)
    
    DEMULTIPLEX_RESULT$sample<-paste( DEMULTIPLEX_RESULT$HGNC,
                                      DEMULTIPLEX_RESULT$Sample_Name,
                                      sep="_")
    
    DEMULTIPLEX_RESULT$name<-gsub("-.+","",DEMULTIPLEX_RESULT$sample)
    
    DEMULTIPLEX_RESULT$replica_order_factor<-gsub("[^_]+_","",DEMULTIPLEX_RESULT$Sample_Name)
    DEMULTIPLEX_RESULT$replica_order_factor<-tolower(DEMULTIPLEX_RESULT$replica_order_factor)
    
    
    
    
    check.vector<-DEMULTIPLEX_RESULT$replica_order_factor[which(DEMULTIPLEX_RESULT$replica_order_factor%in%levels_replica_order_factor)]
    
    cat("check.vector\n")
    cat(sprintf(as.character(check.vector)))
    cat("\n")
    
    cat("levels_replica_order_factor\n")
    cat(sprintf(as.character(levels_replica_order_factor)))
    cat("\n")
    
    
    
    if(length(check.vector) == 0)
    {
      COUNT_matrix$replica_order_factor<-gsub("_","",COUNT_matrix$replica_order_factor)
      
      levels_replica_order_factor<-gsub("_","",levels_replica_order_factor)
      
      COUNT_matrix$replica_order_factor<-factor(COUNT_matrix$replica_order_factor,
                                                levels=levels_replica_order_factor,
                                                ordered = T)
      
      cat("COUNT_matrix_PRE\n")
      cat(str(COUNT_matrix))
      cat("\n")
      
      DEMULTIPLEX_RESULT$replica_order_factor<-factor(DEMULTIPLEX_RESULT$replica_order_factor,
                                                      levels=levels_replica_order_factor,
                                                      ordered = T)
      
      cat("DEMULTIPLEX_RESULT_PRE\n")
      cat(str(DEMULTIPLEX_RESULT))
      cat("\n")
      
      
      
      
    }else{
      
      DEMULTIPLEX_RESULT$replica_order_factor<-factor(DEMULTIPLEX_RESULT$replica_order_factor,
                                                      levels=levels_replica_order_factor,
                                                      ordered = T)
    }
    
    
    
    cat("DEMULTIPLEX_RESULT_PRE\n")
    cat(str(DEMULTIPLEX_RESULT))
    cat("\n")
    
    #quit(status = 1)
    #### READ decompose CIGAR string ----
    
    setwd(out)
    
    filename_3<-paste(type,"_DECOMPOSED_CIGAR",".rds", sep='')
    
    
    Gather_decomposed_cigar<-readRDS(file=filename_3)
    
    
    cat("Gather_decomposed_cigar\n")
    cat(str(Gather_decomposed_cigar))
    cat("\n")
    # cat(sprintf(as.character(Gather_decomposed_cigar$CIGAR[1])))
    # cat("\n")
    # cat(sprintf(as.character(Gather_decomposed_cigar$reads[1])))
    # cat("\n")
    
    #### ACCEPTED CIGAR ----
    
    Gather_decomposed_cigar_subset<-Gather_decomposed_cigar[which(Gather_decomposed_cigar$CIGAR%in%ACCEPTED_CIGAR),]
    
    
    cat("Gather_decomposed_cigar_subset\n")
    cat(str(Gather_decomposed_cigar_subset))
    cat("\n")
    
    cat("Gather_decomposed_cigar_subset$CIGAR\n")
    cat(sprintf(as.character(names(summary(as.factor(Gather_decomposed_cigar_subset$CIGAR))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(Gather_decomposed_cigar_subset$CIGAR)))))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(Gather_decomposed_cigar_subset$mismatches))))))
    cat("\n")
    # cat(sprintf(as.character(summary(as.factor(Gather_decomposed_cigar_subset$mismatches)))))
    # cat("\n")
    
    
    Brackets<-unique(as.character(c(0,1,2,5,as.character(names(summary(as.factor(Gather_decomposed_cigar_subset$mismatches))))[length(as.character(names(summary(as.factor(Gather_decomposed_cigar_subset$mismatches)))))])))
    
    cat("Brackets\n")
    cat(sprintf(Brackets))
    cat("\n")
    
    # quit(status=1)
    
    #### Discretize mismatches----
    
    
    Gather_decomposed_cigar_subset$BIN_mismatches<-cut(Gather_decomposed_cigar_subset$mismatches,breaks = Brackets,right = FALSE)
    
    
    cat("Gather_decomposed_cigar_subset\n")
    cat(str(Gather_decomposed_cigar_subset))
    cat("\n")
    
    levels_BIN_mismatches<-levels(Gather_decomposed_cigar_subset$BIN_mismatches)
    
    
    cat("levels_BIN_mismatches\n")
    cat(sprintf(as.character(levels_BIN_mismatches)))
    cat("\n")
    
    
    
    Gather_decomposed_cigar_subset$BIN_order_factor<-factor(Gather_decomposed_cigar_subset$BIN_mismatches,
                                                            levels=levels_BIN_mismatches,
                                                            ordered = T)
    
    cat("Gather_decomposed_cigar_subset_POST\n")
    cat(str(Gather_decomposed_cigar_subset))
    cat("\n")
    
    # quit(status = 1)
    
    
    #### select INS ----
    
    Gather_decomposed_cigar_INS<-Gather_decomposed_cigar_subset#Gather_decomposed_cigar[which(Gather_decomposed_cigar$S > 0),]
    
    cat("Gather_decomposed_cigar_INS\n")
    cat(str(Gather_decomposed_cigar_INS))
    cat("\n")
    
    
    #### Add up the CIGARS by Sample_Name and BIN_order_factor----
    
    Gather_decomposed_cigar_INS.dt <- data.table(Gather_decomposed_cigar_INS) 
    
    DT.reads<-Gather_decomposed_cigar_INS.dt[ , .(Totalcount = sum(reads)), 
                                              by = .(Sample_Name,BIN_order_factor,replica_order_factor)]
    
    #### Add up the reads by Sample_Name TOTAL reads that belonng to INS in every replica ----
    
    DT.reads2<-Gather_decomposed_cigar_INS.dt[ , .(Totalcount = sum(reads)), 
                                               by = .(Sample_Name,replica_order_factor)]
    
    colnames(DT.reads2)[which(colnames(DT.reads2) == "Totalcount")]<-"Count_of_event"
    
    
    
    #TOTAL_reads_sel<-TOTAL_reads[which(TOTAL_reads$name == name_sel),]
    
    #### merge with total INS reads and with total Sample_Name reads  ----
    
    DEF0<-merge(as.data.frame(DT.reads),
                as.data.frame(DT.reads2),
                by=c("Sample_Name","replica_order_factor"))
    
    colnames(DEF0)[which(colnames(DEF0) == "Sample_Name")]<-"sample"
    
    
    cat("DEF0\n")
    cat(str(DEF0))
    cat("\n")
    
    DEF_INS<-merge(DEF0,
                   TOTAL_READS,
                   by=c("sample"),
                   all.y = T)
    
    DEF_INS$Perc<-round(100*DEF_INS$Totalcount/DEF_INS$TOTAL_reads,2)
    
    
    cat("DEF_INS\n")
    cat(str(DEF_INS))
    cat("\n")
    
    check<-DEF_INS[is.na(DEF_INS$Totalcount),]
    
    cat("check\n")
    cat(str(check))
    cat("\n")
    
    
    # quit(status = 1)
    
    
    #### Calculate Perc over TOTAL_reads in sample ----
    
    
    cat("DEF_INS_POST\n")
    cat(str(DEF_INS))
    cat("\n")
    
    ##### LOOP TO GET Get the distribution of the size INSs ----
    
    Gather_sub_decomposed_cigar<- data.frame(matrix(vector(), 0, 3,
                                                    dimnames=list(c(), 
                                                                  c("name","INS_size",
                                                                    "replica_order_factor"
                                                                  ))),
                                             stringsAsFactors=F)
    
    CIGARS_INS<-unique(as.character(Gather_decomposed_cigar_INS$CIGAR))
    
    cat("CIGARS_INS\n")
    cat(sprintf(as.character(CIGARS_INS)))
    cat("\n")
    
    list_S_DEF<-list()
    
    
    for(k in 1:length(CIGARS_INS))
    {
      CIGAR_sel<-CIGARS_INS[k]
      
      cat("------------------------------------------------>CIGAR_sel\t")
      cat(sprintf(as.character(CIGAR_sel)))
      cat("\n")
      
      Gather_decomposed_cigar_INS_CIGAR<-Gather_decomposed_cigar_INS[which(Gather_decomposed_cigar_INS$CIGAR == CIGAR_sel),]
      
      colnames(Gather_decomposed_cigar_INS_CIGAR)[which(colnames(Gather_decomposed_cigar_INS_CIGAR) == "Sample_Name")]<-"sample"
      
      # cat("Gather_decomposed_cigar_INS_CIGAR_POST\n")
      # cat(str(Gather_decomposed_cigar_INS_CIGAR))
      # cat("\n")
      
      Gather_decomposed_cigar_INS_CIGAR.dt<-data.table(Gather_decomposed_cigar_INS_CIGAR)
      collapse.reads<-Gather_decomposed_cigar_INS_CIGAR.dt[ , .(Totalcount = sum(reads)), 
                                                            by = .(name,I,replica_order_factor,
                                                                   sample)]
      # cat("collapse.reads_POST\n")
      # cat(str(collapse.reads))
      # cat("\n")
      # 
      list_S<-list()
      
      for(l in 1:length(collapse.reads$sample))
      {
        sample_sel<-as.character(collapse.reads$sample[l])
        
        # cat("------------------------------------------------>sample_sel\n")
        # cat(sprintf(as.character(sample_sel)))
        # cat("\n")
        
        Gather_decomposed_cigar_INS_CIGAR_sample<-Gather_decomposed_cigar_INS_CIGAR[which(Gather_decomposed_cigar_INS_CIGAR$sample == sample_sel),]
        
        # cat("Gather_decomposed_cigar_INS_CIGAR_sample_POST\n")
        # cat(str(Gather_decomposed_cigar_INS_CIGAR_sample))
        # cat("\n")
        
        
        collapse.reads_sel<-collapse.reads[which(collapse.reads$sample ==sample_sel),]
        
        # cat("collapse.reads_sel_POST\n")
        # cat(str(collapse.reads_sel))
        # cat("\n")
        
        # quit(status = 1)
        
        
        if(collapse.reads_sel$Totalcount == 0)
        {
          collapse.reads_sel$Totalcount = 1
        }
        
        list_S[[l]]<-collapse.reads_sel
        
        
      }# l sample
      
      TABLE_CIGAR_SIZES = unique(as.data.frame(data.table::rbindlist(list_S, fill = T)))
      
      # cat("TABLE_CIGAR_SIZES\n")
      # str(TABLE_CIGAR_SIZES)
      # cat("\n")
      
      list_S_DEF[[k]]<-TABLE_CIGAR_SIZES
      
      # quit(status = 1)
      
    }# k CIGARS
    
    TABLE_CIGAR_SIZES_DEF = unique(as.data.frame(data.table::rbindlist(list_S_DEF, fill = T)))
    
    colnames(TABLE_CIGAR_SIZES_DEF)[which(colnames(TABLE_CIGAR_SIZES_DEF) == "I")]<-"INS_size"
    
    TABLE_CIGAR_SIZES_DEF$type<-"NA"
    
    TABLE_CIGAR_SIZES_DEF$type[grep("g",TABLE_CIGAR_SIZES_DEF$replica_order_factor)]<-"gDNA"
    TABLE_CIGAR_SIZES_DEF$type[grep("c",TABLE_CIGAR_SIZES_DEF$replica_order_factor)]<-"cDNA"
    
    TABLE_CIGAR_SIZES_DEF$type<-factor(TABLE_CIGAR_SIZES_DEF$type,
                                       levels=c("gDNA","cDNA"),
                                       ordered=T)
    
    
    cat("TABLE_CIGAR_SIZES_DEF\n")
    str(TABLE_CIGAR_SIZES_DEF)
    cat("\n")
    
    
    TABLE_CIGAR_SIZES_DEF.dt<-data.table(TABLE_CIGAR_SIZES_DEF, key=c("name","replica_order_factor"))
    
    TABLE_CIGAR_SIZES_DEF.dt.agg<-as.data.frame(TABLE_CIGAR_SIZES_DEF.dt[,.(Absolutecount = sum(Totalcount)),by=.(name,replica_order_factor)], strinsAsFactors=F)
    
    cat("TABLE_CIGAR_SIZES_DEF.dt.agg\n")
    str(TABLE_CIGAR_SIZES_DEF.dt.agg)
    cat("\n")
    
    
    # quit(status = 1)
    
    
    
    #### LOOP GRAPH INS ----
    
    VARS<-unique(COUNT_matrix$VAR)
    
    cat("VARS\n")
    cat(str(VARS))
    cat("\n")
    
    cat("COUNT_matrix\n")
    cat(str(COUNT_matrix))
    cat("\n")
    
    # quit(status = 1)
    
    
    # HGNC.vector<-unique(COUNT_matrix$HGNC)
    # 
    # cat("HGNC.vector\n")
    # cat(str(HGNC.vector))
    # cat("\n")
    
    
    breaks.INS<-seq(0,100, by=10)
    labels.INS<-as.character(breaks.INS)
    
    
    
    
    # quit(status=1)
    
    
    
    myplots_INS<-list()
    
    # quit(status=1)
    
    for(i in 1:length(VARS))
    {
      VAR_sel<-VARS[i]
      
      cat("--------------->\t")
      cat(sprintf(as.character(VAR_sel)))
      cat("\n")
      
      COUNT_matrix_sel<-COUNT_matrix[which(COUNT_matrix$VAR ==VAR_sel),]
      
      cat("COUNT_matrix_sel\n")
      cat(str(COUNT_matrix_sel))
      cat("\n")
      
      if(dim(COUNT_matrix_sel)[1]>0)
      {
        HGNC.vector<-unique(COUNT_matrix_sel$HGNC)
        
        cat("HGNC.vector\n")
        cat(str(HGNC.vector))
        cat("\n")
        
        for(k in 1:length(HGNC.vector))
        {
          HGNC_sel<-HGNC.vector[k]
          
          cat("--------------->\t")
          cat(sprintf(as.character(HGNC_sel)))
          cat("\n")
          
          
          
          COUNT_matrix_sel_HGNC_sel<-COUNT_matrix_sel[which(COUNT_matrix_sel$HGNC ==HGNC_sel),]
          
          cat("COUNT_matrix_sel_HGNC_sel\n")
          cat(str(COUNT_matrix_sel_HGNC_sel))
          cat("\n")
          
          DEF_INS_sel<-DEF_INS[which(DEF_INS$sample%in%COUNT_matrix_sel_HGNC_sel$sample),]
          
          cat("DEF_INS_sel\n")
          cat(str(DEF_INS_sel))
          cat("\n")
          
          if(dim(DEF_INS_sel)[1] > 0)
          {
            
            cat("----------------------------->HELLO_WORLD\n")
            cat("\n")
            
            
            g.INS<-DEF_INS_sel %>%
              mutate(myaxis = paste0(replica_order_factor, "\n", "n=", TOTAL_reads)) %>%
              mutate(myaxis=fct_reorder(myaxis,as.numeric(replica_order_factor))) %>%
              ggplot(aes(x=myaxis, y=Perc, fill=BIN_order_factor)) +
              geom_bar(stat="identity",colour='black')+
              theme_bw()+
              theme(plot.title=element_text(size=11))+
              scale_x_discrete(name=NULL, drop=F)+
              theme(axis.text.x=element_text(angle=45,size=12,colour="black",vjust=1,hjust=1,face="bold"), 
                    axis.title.y=element_text(size=12,face="bold"),
                    legend.title=element_text(size=12,face="bold"),
                    title=element_text(size=18,face="bold"),
                    legend.text=element_text(size=12,face="bold",colour="black"),
                    axis.text.y=element_text(colour="black",size=12,face="bold"))+
              scale_y_continuous(name ="Percentage of TOTAL READS", breaks =breaks.INS, 
                                 labels=labels.INS, 
                                 limits = c(breaks.INS[1],breaks.INS[length(breaks.INS)]))+
              labs(x="",y="Relative abundance of each class(%)", fill=paste("mismatches", sep="\n"))+
              scale_fill_manual(values=c('#32A852','#6DB2EE','#553B68','red','#1877C9','#C9244B','#D45E85','#87447B','#D6B8E6'),drop=F)+
              ggeasy::easy_center_title()
            
            #ggtitle(paste("Reads INS",VAR_sel, HGNC_sel, sep= " ")) +
            
            # myplots_INS[[VAR_sel]][[HGNC_sel]]<-g.INS
            # 
            # setwd(out)
            # 
            # pdf(file=paste("test","_",VAR_sel,"_",HGNC_sel,".pdf",sep=''))
            # print(g.INS)
            # dev.off()
            
            cat("----------------------------->THE END\n")
            cat("\n")
            
            # quit(status = 1)
            
          } #dim(DEF_INS_sel)[1] > 0
          
          
          #### sizes
          
          TABLE_CIGAR_SIZES_DEF_sel<-TABLE_CIGAR_SIZES_DEF[which(TABLE_CIGAR_SIZES_DEF$sample%in%COUNT_matrix_sel_HGNC_sel$sample),]
          
          cat("TABLE_CIGAR_SIZES_DEF_sel\n")
          cat(str(TABLE_CIGAR_SIZES_DEF_sel))
          cat("\n")
          
          # quit(status=1)
          
          if(dim(TABLE_CIGAR_SIZES_DEF_sel)[1] > 0)
          {
            
            cat("----------------------------->HELLO_WORLD_2\n")
            cat("\n")
            
            max.log.INS_size<-max(log10(TABLE_CIGAR_SIZES_DEF_sel$INS_size+0.001))
            
            cat("max.log.INS_size\n")
            cat(sprintf(as.character(max.log.INS_size)))
            cat("\n")
            
            min.log.INS_size<-min(log10(TABLE_CIGAR_SIZES_DEF_sel$INS_size+0.001))
            
            cat("min.log.INS_size\n")
            cat(sprintf(as.character(min.log.INS_size)))
            cat("\n")
            
            breaks.INS_size.log<-sort(unique(c(log10(1),min.log.INS_size,seq(min.log.INS_size,max.log.INS_size, by=0.5),max.log.INS_size+0.01)), decreasing = F)
            
            cat("breaks.INS_size.log\n")
            cat(sprintf(as.character(breaks.INS_size.log)))
            cat("\n")
            
            labels.INS_size<-as.character(round(10^breaks.INS_size.log,1))
            
            cat("labels.INS_size\n")
            cat(sprintf(as.character(labels.INS_size)))
            cat("\n")
            
            
            max.log.Totalcount<-max(log10(TABLE_CIGAR_SIZES_DEF_sel$Totalcount+0.001))
            
            cat("max.log.Totalcount\n")
            cat(sprintf(as.character(max.log.Totalcount)))
            cat("\n")
            
            min.log.Totalcount<-min(log10(TABLE_CIGAR_SIZES_DEF_sel$Totalcount+0.001))
            
            cat("min.log.Totalcount\n")
            cat(sprintf(as.character(min.log.Totalcount)))
            cat("\n")
            
            breaks.Totalcount.log<-sort(unique(c(log10(1),min.log.Totalcount,seq(min.log.Totalcount,max.log.Totalcount, by=0.5),max.log.Totalcount+0.01)), decreasing = F)
            
            cat("breaks.Totalcount.log\n")
            cat(sprintf(as.character(breaks.Totalcount.log)))
            cat("\n")
            
            labels.Totalcount<-as.character(round(10^breaks.Totalcount.log,1))
            
            cat("labels.Totalcount\n")
            cat(sprintf(as.character(labels.Totalcount)))
            cat("\n")
            
            
            g.INS.sizes<-ggplot(TABLE_CIGAR_SIZES_DEF_sel, aes(x=log10(INS_size+0.001), 
                                                               y=log10(Totalcount+0.001), 
                                                               color=replica_order_factor,
                                                               shape=type)) +
              geom_point() +
              theme_bw()+
              theme(legend.position="bottom",legend.title=element_blank(), legend.text = element_text(size=10))+
              guides(fill=guide_legend(nrow=2,byrow=TRUE))+
              theme(plot.title = element_text(size=11)) +
              theme(plot.title=element_text(size=11))+
              scale_x_continuous(name="INS size", 
                                 breaks=breaks.INS_size.log,
                                 labels=labels.INS_size, 
                                 limits=c(breaks.INS_size.log[1],breaks.INS_size.log[length(breaks.INS_size.log)]))+
              scale_y_continuous(name="Total reads",
                                 breaks=breaks.Totalcount.log,
                                 labels=labels.Totalcount, 
                                 limits=c(breaks.Totalcount.log[1],breaks.Totalcount.log[length(breaks.Totalcount.log)]))+
              theme(axis.title.y=element_text(size=14,face="bold"),
                    axis.title.x=element_text(size=14,face="bold"),
                    legend.text=element_text(size=8,colour="black"),
                    axis.text.y=element_text(angle=45, vjust=1,hjust=1,colour="black",size=8),
                    axis.text.x=element_text(angle=45, vjust=1,hjust=1,colour="black",size=8))+
              scale_color_manual(values=c('#32A852','#6DB2EE','#62D07F','#1877C9','#C9244B','#D45E85','#87447B','#553B68','#D6B8E6',
                                          "#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", 
                                          "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
                                          "#D14285", "#6DDE88", "#652926", "#D1A33D", "#C84248","#89C5DA", "#DA5724", "#74D944", "#CE50CA"),
                                 drop=F,
                                 name="Totals",breaks=TABLE_CIGAR_SIZES_DEF.dt.agg$replica_order_factor,
                                 labels=paste(TABLE_CIGAR_SIZES_DEF.dt.agg$replica_order_factor,
                                              TABLE_CIGAR_SIZES_DEF.dt.agg$Absolutecount, sep =' n= '))+
              ggeasy::easy_center_title()
            
            #ggtitle(paste("Reads INS sizes",VAR_sel, HGNC_sel, sep= " ")) +
            
            
            title <- ggdraw() + draw_label(paste("INS",VAR_sel,HGNC_sel,sep=" "))
            
            graph_INS_DEF<-plot_grid(title,g.INS, g.INS.sizes,
                                     nrow = 3,
                                     align = "v",
                                     rel_heights = c(0.1,1,1))
            
            
            
            myplots_INS[[VAR_sel]][[HGNC_sel]]<-graph_INS_DEF
            
            
            setwd(out)
            
            pdf(file=paste("test","_",VAR_sel,"_",HGNC_sel,".pdf",sep=''))
            print(graph_INS_DEF)
            dev.off()
            
            
            cat("----------------------------->THE END_2\n")
            cat("\n")
            
            # quit(status = 1)
            
            
            
            
          } #dim(TABLE_CIGAR_SIZES_DEF_sel)[1] > 0
          
        }#k HGNC
        
        
      }# dim(COUNT_matrix_sel)[1]>0
      
    }#i VARS
    
    
    #### SAVE RDS ----
    
    filename=paste(type, "_CIGAR_INS",".rds", sep='')
    
    saveRDS(myplots_INS, file = filename)
  }# dim(mismatch_CONDENSED_subset)[1] >0)
  
  # quit(status=1)
  
  
 
  # quit(status = 1)
  
 
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
    make_option(c("--COUNT_MATRIX"), type="character", default=NULL,
                metavar="FILE.tsv",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--mismatch_CONDENSED"), type="character", default=NULL, 
                metavar="FILE.tsv", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--TOTAL_READS"), type="character", default=NULL, 
                metavar="FILE.tsv", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--DEMULTIPLEX_RESULT"), type="character", default=NULL, 
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
                        --COUNT_MATRIX FILE.tsv
                        --mismatch_CONDENSED FILE.tsv 
                        --TOTAL_READS FILE.tsv 
                        --DEMULTIPLEX_RESULT FILE.tsv
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)


 
  Graph_INS(opt)
 
  
}


###########################################################################

system.time( main() )
