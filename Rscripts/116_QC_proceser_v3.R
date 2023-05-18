

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
suppressMessages(library("BiocGenerics", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("S4Vectors", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("IRanges", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("GenomeInfoDb", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("GenomicRanges", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("Biobase", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("AnnotationDbi", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("GenomicFeatures", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("OrganismDbi", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
library("XVector", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
suppressMessages(library("Biostrings", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("seqinr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("seqRFLP", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("splitstackshape", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("cowplot", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))




#suppressMessages(library("udunits2", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
#suppressMessages(library("ggforce", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))

##### Taken from https://www.biostars.org/p/9335/

matcher <- function(pattern, x) {
  
  ind = gregexpr(pattern, x)[[1]]
  start = as.numeric(ind)
  end = start + attr(ind, "match.length")- 2
  apply(cbind(start,end), 1, function(y) substr(x, start=y[1], stop=y[2]));
}
doone <- function(c, cigar) {
  pat <- paste("\\d+", c , sep="")
  sum(as.numeric(matcher(pat, cigar)), na.rm=T)
}


cigarsums <- function(cigar, chars=c("M","N","D","I","S","H", "P", "X", "=")) {
  sapply (chars, doone, cigar)
}


opt = NULL

TOTAL_and_HMAP_graphs = function(option_list)
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
  
  COUNT_matrix = read.table(opt$COUNT_MATRIX, sep="\t", stringsAsFactors = F, header = T)
  
  COUNT_matrix$name<-gsub("-.+","",COUNT_matrix$sample)
  
  COUNT_matrix$replica_order_factor<-gsub("DNA_","",COUNT_matrix$Replicate)
  
  cat("COUNT_matrix_0\n")
  cat(str(COUNT_matrix))
  cat("\n")

  
  
  
  levels_replica_order_factor<-sort(unique(COUNT_matrix$replica_order_factor))
  
  cat("levels_replica_order_factor\n")
  cat(sprintf(as.character(levels_replica_order_factor)))
  cat("\n")
  
  gDNA_part<-levels_replica_order_factor[grep("g",levels_replica_order_factor)]
  
  cat("gDNA_part\n")
  cat(sprintf(as.character(gDNA_part)))
  cat("\n")
  
  
  cDNA_part<-levels_replica_order_factor[grep("c",levels_replica_order_factor)]
  
  cat("cDNA_part\n")
  cat(sprintf(as.character(cDNA_part)))
  cat("\n")
  
  levels_replica_order_factor<-c(gDNA_part,cDNA_part)
  
  cat("levels_replica_order_factor\n")
  cat(sprintf(as.character(levels_replica_order_factor)))
  cat("\n")
  
  # HERE HERE 
  
  
  
  
  
  
  COUNT_matrix$replica_order_factor<-factor(COUNT_matrix$replica_order_factor,
                                            levels=levels_replica_order_factor,
                                            ordered = T)
  
  cat("COUNT_matrix\n")
  cat(str(COUNT_matrix))
  cat("\n")
  
  # quit(status = 1)
 
  #### READ and transform mismatch_CONDENSED ----
  
  mismatch_CONDENSED = read.table(opt$mismatch_CONDENSED, sep="\t", stringsAsFactors = F, header = T)
  
  
  
  colnames(mismatch_CONDENSED)[which(colnames(mismatch_CONDENSED) == "Sample_Name")]<-"sample"
  
  mismatch_CONDENSED$name<-gsub("-.+","",mismatch_CONDENSED$sample)
  mismatch_CONDENSED$sample<-gsub("-","_",mismatch_CONDENSED$sample)
  
  
  cat("mismatch_CONDENSED\n")
  cat(str(mismatch_CONDENSED))
  cat("\n")
  
    
  
  #
  
  #### READ and transform TOTAL_READS ----
  
  TOTAL_READS = read.table(opt$TOTAL_READS, sep="\t", stringsAsFactors = F, header = T)
  
  colnames(TOTAL_READS)[which(colnames(TOTAL_READS) == "Sample_Name")]<-"sample"
  TOTAL_READS$name<-gsub("-.+","",TOTAL_READS$sample)
  
  TOTAL_READS$sample<-gsub("-","_",TOTAL_READS$sample)
  
  
  
  cat("TOTAL_READS\n")
  cat(str(TOTAL_READS))
  cat("\n")
  
  # quit(status = 1)
  
  #### READ and transform DEMULTIPLEX_RESULT ----
  
  DEMULTIPLEX_RESULT = read.table(opt$DEMULTIPLEX_RESULT, sep="\t", stringsAsFactors = F, header = T)
  
  cat("DEMULTIPLEX_RESULT_0\n")
  cat(str(DEMULTIPLEX_RESULT))
  cat("\n")
  
  DEMULTIPLEX_RESULT$sample<-paste( DEMULTIPLEX_RESULT$HGNC,
                                    DEMULTIPLEX_RESULT$Sample_Name,
                                    sep="_")
  
  DEMULTIPLEX_RESULT$name<-gsub("[-_][^-_]+$","",DEMULTIPLEX_RESULT$sample)
  
  

  
  cat("DEMULTIPLEX_RESULT_PRE\n")
  cat(str(DEMULTIPLEX_RESULT))
  cat("\n")
 
  # quit(status = 1)
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
  
  ####  FIRST get TOTAL amount of mapped reads in COUNT_MATRIX ----
  
  indx.column<-c(grep("haplotype",colnames(COUNT_matrix)),
                grep("UNDETERMINED",colnames(COUNT_matrix)))
  
  cat("indx.column\n")
  cat(sprintf(as.character(indx.column)))
  cat("\n")
  
  COUNT_matrix.m<-as.matrix(COUNT_matrix[,indx.column])
  
  cat("COUNT_matrix.m\n")
  cat(str(COUNT_matrix.m))
  cat("\n")
  
  
  row.names(COUNT_matrix.m)<-COUNT_matrix$Sample_Name
  indexTOTAL<-dim(COUNT_matrix.m)[2]+1
  
  cat("indexTOTAL\n")
  cat(sprintf(as.character(indexTOTAL)))
  cat("\n")
  
  COUNT_matrix$TOTAL<-apply(COUNT_matrix.m,1,sum)
  
  cat("COUNT_matrix_TOTAL\n")
  cat(str(COUNT_matrix))
  cat("\n")
  
  
 
  
  #### Second merge COUNT_matrix and DEMULTIPLEX_RESULT----
  
  cat("COUNT_matrix.m\n")
  cat(str(COUNT_matrix.m))
  cat("\n")
  
  cat("DEMULTIPLEX_RESULT\n")
  cat(str(DEMULTIPLEX_RESULT))
  cat("\n")
  
  DEF_TOT_vs_CLASS<-merge(COUNT_matrix,
                          DEMULTIPLEX_RESULT,
                          by=c("sample","name",
                               "Sample_Name","HGNC"),
                          all.x=T)
  
  DEF_TOT_vs_CLASS$TOTAL[is.na(DEF_TOT_vs_CLASS$TOTAL)]<-0
  
  
  
  # quit(status=1)
  
  #### SAVE DEMULTIPLEX_and_COUNT ----
  
  setwd(out)
  
  cat("DEF_TOT_vs_CLASS--------------->SAVE\n")
  cat(str(DEF_TOT_vs_CLASS))
  cat("\n")
  
  
  filename_1<-paste(type,"_DEMULTIPLEX_and_COUNT",".rds", sep='')

  #saveRDS(DEF_TOT_vs_CLASS, file=filename_1, sep="\t", quote=F, row.names = F)
  saveRDS(DEF_TOT_vs_CLASS, file = filename_1)
  
  
  # # quit(status = 1)
  
  
  #### LOOP TO CALCULATE RELATIVE AMOUNTS OF alignment in COUNT MATRIX LIABILITY !!! ----
  
  #### initialize variables
  
  indx.column<-c(grep("haplotype",colnames(COUNT_matrix)),
                grep("UNDETERMINED",colnames(COUNT_matrix)))
  
  COUNT_name_vector<-unique(COUNT_matrix$name)
  
  
  Gather_Heat_map<- data.frame(matrix(vector(), 0, 5,
                                      dimnames=list(c(), 
                                                    c("name","Sample_Name",
                                                      "replica_order_factor","Amplicon",
                                                      "Percentage_aligned"))),
                               stringsAsFactors=F)
  
  #### loop
  for(i in 1:length(COUNT_name_vector))
  {
    name_sel<-COUNT_name_vector[i]
    
    COUNT_matrix_sel<-COUNT_matrix[which(COUNT_matrix$name == name_sel),]
        
    Sample_Name_sel<-COUNT_matrix_sel$Sample_Name
    replica_order_factor_sel<-COUNT_matrix_sel$replica_order_factor
    
    for(k in 1:length(indx.column))
    {
      indx.column_sel<-indx.column[k]
      colnames_sel<-colnames(COUNT_matrix_sel)[indx.column_sel]
      
      Percentage_aligned_sel<-round(100*COUNT_matrix_sel[,indx.column_sel]/COUNT_matrix_sel$TOTAL,2)
      
      tmp<-as.data.frame(cbind(rep(name_sel,length(Percentage_aligned_sel)),
                               Sample_Name_sel,as.character(replica_order_factor_sel),
                               rep(colnames_sel,length(Percentage_aligned_sel)),
                               Percentage_aligned_sel), stringsAsFactors=F)
      
      # cat("tmp\n")
      # cat(str(tmp))
      # cat("\n")
      
      
      tmp$Percentage_aligned_sel<-as.numeric(tmp$Percentage_aligned_sel)
      
      colnames(tmp)<-colnames(Gather_Heat_map)
      
      # cat("tmp\n")
      # cat(str(tmp))
      # cat("\n")
      # 
      # quit(status = 1)
      
      
      Gather_Heat_map<-rbind(Gather_Heat_map,tmp)
      
    } #k indexes of amplicons
    
  } # names
  
  
  
  #### Add and order the result of the loop
  
  Gather_Heat_map$Y1 <- cut(Gather_Heat_map$Percentage_aligned,
                            breaks = c(0,0.5,20,50,70,80,90,95,Inf),
                            right = FALSE)
  
  Gather_Heat_map$replica_order_factor<-factor(Gather_Heat_map$replica_order_factor,
                                               levels=levels_replica_order_factor,
                                               ordered = T)
  
  Gather_Heat_map$HGNC<-gsub(".+__","",Gather_Heat_map$Amplicon)
  Gather_Heat_map$HGNC<-gsub("_.+$","",Gather_Heat_map$HGNC)
  
  cat("Gather_Heat_map\n")
  cat(str(Gather_Heat_map))
  cat("\n")
  
  
  HGNC_order_vector<-levels(as.factor(Gather_Heat_map$HGNC))
  
  HGNC_order_vector_UNT<-HGNC_order_vector[which(HGNC_order_vector == "UNDETERMINED")]
  HGNC_order_vector_REST<-HGNC_order_vector[-which(HGNC_order_vector == "UNDETERMINED")]
  
  HGNC_order_vector_DEF<-c(HGNC_order_vector_REST,
                           HGNC_order_vector_UNT)
  
  Gather_Heat_map$HGNC<-factor(Gather_Heat_map$HGNC,
                            levels=HGNC_order_vector_DEF,
                            ordered = T)
  
  
 
  
  # Gather_Heat_map$Y1<-factor(Gather_Heat_map$Y1,
  #                    levels=c("[0,0.5)","[0.5,20)","[20,50)",
  #                             "[50,70)","[70,80)","[80,90)",
  #                             "[90,95)","[95,Inf)"),
  #                    ordered = T)
  
  cat("Gather_Heat_map\n")
  cat(str(Gather_Heat_map))
  cat("\n")
  cat("HGNC\n")
  cat(sprintf(as.character(unique(Gather_Heat_map$HGNC))))
  cat("\n")
  cat("name\n")
  cat(sprintf(as.character(unique(Gather_Heat_map$name))))
  cat("\n")
  
  # quit(status=1)
  
  #### SAVE ----
  setwd(out)
  
  filename_2<-paste(type,"_Gather_Heat_map",".txt", sep='')
  
  write.table(Gather_Heat_map, file=filename_2, sep="\t", quote=F, row.names = F)
  
  # quit(status = 1)
  
  
 
  
  #### LOOP GRAPH DEL ----
  
  VARS<-unique(COUNT_matrix$VAR)
  
  cat("VARS\n")
  cat(str(VARS))
  cat("\n")
  
  myplots_TOTAL<-list()
  myplots_HEATMAP<-list()
  
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
        
        Gather_Heat_map_sel<-Gather_Heat_map[which(Gather_Heat_map$name%in%COUNT_matrix_sel_HGNC_sel$sample),]
        
        cat("Gather_Heat_map_sel\n")
        cat(str(Gather_Heat_map_sel))
        cat("\n")
        
        DEF_TOT_vs_CLASS_sel<-DEF_TOT_vs_CLASS[which(DEF_TOT_vs_CLASS$sample%in%COUNT_matrix_sel_HGNC_sel$sample),]
        
        cat("DEF_TOT_vs_CLASS_sel\n")
        cat(str(DEF_TOT_vs_CLASS_sel))
        cat("\n")
        
        # quit(status=1)
        
        #### GRAPH 1: TOT_vs_CLASS ----
        
        MAX.TOT<-round(max(DEF_TOT_vs_CLASS_sel$TOTAL)+max(DEF_TOT_vs_CLASS_sel$TOTAL)/10,0)
        break.TOT<-round(seq(0,MAX.TOT,by=MAX.TOT/10),0)
        labels.TOT<-as.character(break.TOT)
        
        # str(DEF_TOT_vs_CLASS_sel)
        # str(MAX.TOT)
        
        
        
        g.TOT<-ggplot(DEF_TOT_vs_CLASS_sel, aes(fill=CLASS, y=TOTAL,
                                                x=replica_order_factor)) +
          geom_bar(stat="identity")+
          scale_fill_manual(values = c("green","#DA5724","#89C5DA","#C0717C" , "#CE50CA", "#74D944", "#CBD588", "#5F7FC7",
                                       "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD",
                                       "#D14285", "#6DDE88", "#652926", "#D1A33D", "#C84248","#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588","#673770", "#D3D93E", "#38333E"),
                            drop=F)+
          scale_y_continuous(name ="Total reads", breaks =break.TOT,
                             labels=labels.TOT,
                             limits = c(break.TOT[1],break.TOT[length(break.TOT)]))+
          scale_x_discrete(name=NULL)+
          ggtitle(paste("TOTAL READS vs DEMULTIPLEX CLASS", name_sel, sep= " ")) +
          theme(legend.position="right")+ labs(fill = "DEMULTIPLEX_CLASS")
        
        
        myplots_TOTAL[[VAR_sel]][[HGNC_sel]]<-g.TOT
        
        #### graph heat map ----
        
        g.Gather_Heat_map_sel<-ggplot(data = Gather_Heat_map_sel, 
                                      aes(x=replica_order_factor, 
                                          y=Amplicon, 
                                          fill=Y1)) + 
          geom_tile()+
          scale_fill_manual(breaks=c("[0,0.5)","[0.5,20)","[20,50)",
                                     "[50,70)","[70,80)","[80,90)",
                                     "[90,95)","[95,Inf)"),
                            values=c("white","gray","orange","yellow",
                                     "yellow3","yellowgreen","green",
                                     "darkblue"),
                            drop =FALSE )+
          scale_y_discrete(name ="Amplicons")+
          scale_x_discrete(name=NULL)+
          ggtitle(paste("Alignment distribution", name_sel, sep= " ")) +
          theme(legend.position="right")+ labs(fill = "Percentages of reads")+
          geom_hline(yintercept =c(3,6,9,12,15,18,21), linetype="dotted", color="gray")
        
        myplots_HEATMAP[[VAR_sel]][[HGNC_sel]]<-g.Gather_Heat_map_sel
        
        # quit(status = 1)
        

        
      }# k HGNC.vector
    }# dim(COUNT_matrix_sel)[1]>0
  }#i VAR_sel
  
  
  

  
  #### SAVE RDS ----
  
  setwd(out)
  
  filename=paste(type, "_TOTAL",".rds", sep='')
  
  saveRDS(myplots_TOTAL, file = filename)
  
  filename=paste(type, "_HEATMAP",".rds", sep='')
  
  saveRDS(myplots_HEATMAP, file = filename)
  
  filename_3<-paste(type,"_COUNT_matrix",".rds", sep='')
  
  #write.table(Gather_decomposed_cigar, file=filename_3, sep="\t", quote=F, row.names = F)
  saveRDS(COUNT_matrix, file=filename_3)
  
  # quit(status = 1)
  
  
}

Decompose_CIGAR_String = function(option_list)
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
  
  
  COUNT_matrix = readRDS(file=filename_3)#read.table(opt$COUNT_MATRIX, sep="\t", stringsAsFactors = F, header = T)
   
  levels_replica_order_factor<-levels(COUNT_matrix$replica_order_factor)
  
  # COUNT_matrix$name<-gsub("-.+","",COUNT_matrix$sample)
  # 
  # COUNT_matrix$replica_order_factor<-gsub("DNA_","",COUNT_matrix$Replicate)
  # 
  # 
  # levels_replica_order_factor<-sort(unique(COUNT_matrix$replica_order_factor))
  # 
  # cat("levels_replica_order_factor\n")
  # cat(sprintf(as.character(levels_replica_order_factor)))
  # cat("\n")
  # 
  # gDNA_part<-levels_replica_order_factor[grep("g",levels_replica_order_factor)]
  # 
  # cat("gDNA_part\n")
  # cat(sprintf(as.character(gDNA_part)))
  # cat("\n")
  # 
  # 
  # cDNA_part<-levels_replica_order_factor[grep("c",levels_replica_order_factor)]
  # 
  # cat("cDNA_part\n")
  # cat(sprintf(as.character(cDNA_part)))
  # cat("\n")
  # 
 
  # 
  # cat("levels_replica_order_factor\n")
  # cat(sprintf(as.character(levels_replica_order_factor)))
  # cat("\n")
  # 
  # # HERE HERE 
  # 
  # 
  # COUNT_matrix$replica_order_factor<-factor(COUNT_matrix$replica_order_factor,
  #                                           levels=levels_replica_order_factor,
  #                                           ordered = T)
  # 
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
  
  
  
  #
  
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
  
  cat("DEMULTIPLEX_RESULT_PRE\n")
  cat(str(DEMULTIPLEX_RESULT))
  cat("\n")
  
  DEMULTIPLEX_RESULT$name<-gsub("[-_][^-_]+$","",DEMULTIPLEX_RESULT$sample)
  
  DEMULTIPLEX_RESULT$replica_order_factor<-gsub("[^_]+_","",DEMULTIPLEX_RESULT$Sample_Name)
  DEMULTIPLEX_RESULT$replica_order_factor<-tolower(DEMULTIPLEX_RESULT$replica_order_factor)
  
  cat("DEMULTIPLEX_RESULT_PRE\n")
  cat(str(DEMULTIPLEX_RESULT))
  cat("\n")
  
  
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
  
  
  # quit(status=1)
  
  # 
  
  
  
 
  
  
  
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
  
  #### LOOP to decompose CIGAR string ----
  
  
  Gather_decomposed_cigar<- data.frame(matrix(vector(), 0, 14,
                                              dimnames=list(c(), 
                                                            c("name","Sample_Name","CIGAR",
                                                              "reads","mismatches",
                                                              "M","N","D",
                                                              "I","S","H",
                                                              "P","X","SAME"))),
                                       stringsAsFactors=F)
  
  for(i in 1:length(mismatch_CONDENSED$reads))
  {
    name_sel<-mismatch_CONDENSED$name[i]
    Sample_Name_sel<-mismatch_CONDENSED$sample[i]
    CIGAR_sel<-mismatch_CONDENSED$CIGAR[i]
    reads_sel<-mismatch_CONDENSED$reads[i]
    NM_sel<-mismatch_CONDENSED$NM[i]
    NM_sel<-as.numeric(gsub("NM:i:","",NM_sel))
    
    cat("-------------------->name_sel\n")
    cat(sprintf(as.character(name_sel)))
    cat("\n")
    cat("Sample_Name_sel\n")
    cat(sprintf(as.character(Sample_Name_sel)))
    cat("\n")
    cat("CIGAR_sel\n")
    cat(sprintf(as.character(CIGAR_sel)))
    cat("\n")
    cat("reads_sel\n")
    cat(sprintf(as.character(reads_sel)))
    cat("\n")
    cat("NM_sel\n")
    cat(sprintf(as.character(NM_sel)))
    cat("\n")
   
    #quit(status = 1)
    
    
    A<-cigarsums(CIGAR_sel)
    
    cat("A\n")
    cat(str(A))
    cat("\n")
    
    B<-data.frame(matrix(as.numeric(A), 1, 9,
                         dimnames=list(c(), 
                                       c("M","N","D",
                                         "I","S","H",
                                         "P","X","SAME"))),
                  stringsAsFactors=F)
    cat("B\n")
    cat(str(B))
    cat("\n")
    
    D<-as.data.frame(cbind(name_sel,Sample_Name_sel,CIGAR_sel,reads_sel,
                           NM_sel,B))
    
    cat("D_PRE\n")
    cat(str(D))
    cat("\n")
    
    colnames(D)<-colnames(Gather_decomposed_cigar)
    
    cat("D_POST\n")
    cat(str(D))
    cat("\n")
    
    Gather_decomposed_cigar<-rbind(D,Gather_decomposed_cigar)
    
    cat("Gather_decomposed_cigar\n")
    cat(str(Gather_decomposed_cigar))
    cat("\n")
    
    
  }
  
  #### SAVE DECOMPOSED CIGARS ----
  
  setwd(out)
  
  cat("Gather_decomposed_cigar\n")
  cat(str(Gather_decomposed_cigar))
  cat("\n")
  
  
  # quit(status=1)
  
  Gather_decomposed_cigar$replica_order_factor<-gsub("[^_]+_","",Gather_decomposed_cigar$Sample_Name)
  Gather_decomposed_cigar$replica_order_factor<-tolower(Gather_decomposed_cigar$replica_order_factor)
  
  cat("Gather_decomposed_cigar\n")
  cat(str(Gather_decomposed_cigar))
  cat("\n")
  
  Gather_decomposed_cigar$replica_order_factor<-factor(Gather_decomposed_cigar$replica_order_factor,
                                                  levels=levels_replica_order_factor,
                                                  ordered = T)
  
  
  cat("Gather_decomposed_cigar\n")
  cat(str(Gather_decomposed_cigar))
  cat("\n")

  filename_3<-paste(type,"_DECOMPOSED_CIGAR",".rds", sep='')
  
  #write.table(Gather_decomposed_cigar, file=filename_3, sep="\t", quote=F, row.names = F)
  saveRDS(Gather_decomposed_cigar, file=filename_3)
  
  
  
  
  
  
  # quit(status=1)
  
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

   TOTAL_and_HMAP_graphs(opt)
   Decompose_CIGAR_String(opt)
  
}


###########################################################################

system.time( main() )
