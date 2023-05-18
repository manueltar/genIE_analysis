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
suppressMessages(library("ggdendro", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("sandwich", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
library("digest", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("ggforce", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
suppressMessages(library("BiocGenerics", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("S4Vectors", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("IRanges", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("GenomeInfoDb", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("GenomicRanges", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("Biobase", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
library("NOISeq", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("ggrepel", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")


opt = NULL


Normalize_STAR_alignment_Kolf2 = function(option_list)
{
  library("NOISeq")
  
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
  
  #### Read STAR_kolf2_1 ----
  
  STAR_kolf2_1<-as.data.frame(fread(opt$STAR_kolf2_1,
                              sep="\t", header=T), stringsAsFactors = F)
  
  
  cat("STAR_kolf2_1:\n")
  cat(str(STAR_kolf2_1))
  cat("\n")
  
  indx.int<-c(which(colnames(STAR_kolf2_1) == "gene_id"),
              which(colnames(STAR_kolf2_1) == "transcript_id"),
              which(colnames(STAR_kolf2_1) == "pme_FPKM"))
  
  STAR_kolf2_1_subset<-unique(STAR_kolf2_1[,indx.int])
  colnames(STAR_kolf2_1_subset)<-c("ensembl_gene_id","ensembl_transcript_id","Rep_1")
  
  
  cat("STAR_kolf2_1_subset:\n")
  cat(str(STAR_kolf2_1_subset))
  cat("\n")
  
  
  #### Read STAR_kolf2_2 ----
  
  STAR_kolf2_2<-as.data.frame(fread(opt$STAR_kolf2_2,
                                    sep="\t", header=T), stringsAsFactors = F)
  
  
  cat("STAR_kolf2_2:\n")
  cat(str(STAR_kolf2_2))
  cat("\n")
  
  indx.int<-c(which(colnames(STAR_kolf2_2) == "gene_id"),
              which(colnames(STAR_kolf2_2) == "transcript_id"),
              which(colnames(STAR_kolf2_2) == "pme_FPKM"))
  
  STAR_kolf2_2_subset<-unique(STAR_kolf2_2[,indx.int])
  colnames(STAR_kolf2_2_subset)<-c("ensembl_gene_id","ensembl_transcript_id","Rep_2")
  
  
  cat("STAR_kolf2_2_subset:\n")
  cat(str(STAR_kolf2_2_subset))
  cat("\n")
  
  # quit(status=1)
  
  #### Read gene_lengths ----
  
  indx.int<-c(which(colnames(STAR_kolf2_2) == "gene_id"),
              which(colnames(STAR_kolf2_2) == "transcript_id"),
              which(colnames(STAR_kolf2_2) == "length"))
  
  transcript_lengths<-STAR_kolf2_2[,indx.int]
  
  colnames(transcript_lengths)<-c("ensembl_gene_id","ensembl_transcript_id","length")
  
  
  
  cat("transcript_lengths:\n")
  cat(str(transcript_lengths))
  cat("\n")
  

  #### merge both ----
  
  STAR_kolf2<-merge(STAR_kolf2_1_subset,STAR_kolf2_2_subset,
                    by=c("ensembl_gene_id","ensembl_transcript_id"),
                    all=T)
  
  cat("STAR_kolf2:\n")
  cat(str(STAR_kolf2))
  cat("\n")
  
  # quit(status = 1)
  
  transcript_lengths_subset<-transcript_lengths[which(transcript_lengths$ensembl_gene_id%in%STAR_kolf2$ensembl_gene_id &
                                                        transcript_lengths$ensembl_transcript_id%in%STAR_kolf2$ensembl_transcript_id),]
  
  cat("transcript_lengths_subset:\n")
  cat(str(transcript_lengths_subset))
  cat("\n")
  
  mylength<-as.vector(transcript_lengths_subset$length)
  names(mylength) = transcript_lengths_subset$ensembl_transcript_id
  
  cat("mylength:\n")
  cat(str(mylength))
  cat("\n")
  
  #### NOISeq part ----
  
  mycounts<-as.matrix(STAR_kolf2[,-c(which(colnames(STAR_kolf2) == "ensembl_gene_id"),
                                     which(colnames(STAR_kolf2) == "ensembl_transcript_id"))])
  
  row.names(mycounts)<-STAR_kolf2[,which(colnames(STAR_kolf2) == "ensembl_transcript_id")]
  
  cat("mycounts:\n")
  cat(str(mycounts))
  cat("\n")
  
  myfactors = data.frame("group" = as.factor(c(rep(c("DNA.seq.iPSCs"),2))))
  
  cat("myfactors:\n")
  cat(str(myfactors))
  cat("\n")
  
  # creation of a NOISeq object storing all the information: counts, gene length and factors
  mydata = readData(data = mycounts, length = mylength, factors = myfactors)
  
  cat("mydata:\n")
  cat(str(mydata))
  cat("\n")
  
  #### normalization ----
  

  ## RPKM --> sequencing depth and length correction
  rpkm.counts = rpkm(mycounts, lc = 1)
  
  ## TMM --> sequencing depth and DNA composition correction
  tmm.counts = tmm(mycounts, lc = 1)
  
  ## TMM + length correction 
  tmm.counts.length = tmm(mycounts, long = as.numeric(mylength), lc = 1)
  
  cat("rpkm.counts:\n")
  cat(str(rpkm.counts))
  cat("\n")
  
  #### QC plots ----
  
 
  
  path2<-paste(out,type,'/',sep='')
  
  cat("path2\n")
  cat(sprintf(as.character(path2)))
  cat("\n")
  
  if (file.exists(path2)){
    
    setwd(path2)
  } else {
    dir.create(file.path(path2))
    setwd(path2)
    
  }
  
  #General distribution of reads
  
  pdfname<-paste(type,"_QC_plots_normalization.pdf", sep='')
  makepdf = TRUE
  
  if (makepdf == TRUE)
  {
    pdf ( pdfname , height=10, width=12)
  }
  
  
  
  par(mfrow = c(1,2))
  boxplot(log(as.matrix(mycounts+1)) ~ col(mycounts), main = "Before normalization", names=colnames(mycounts))
  boxplot(log(rpkm.counts+1) ~ col(rpkm.counts), main = "After rpkm normalization", names=colnames(mycounts))
  boxplot(log(tmm.counts+1) ~ col(tmm.counts), main = "After tmm normalization", names=colnames(mycounts))
  boxplot(log(tmm.counts.length+1) ~ col(tmm.counts.length), main = "After tmm+length normalization", names=colnames(mycounts))
  #dev.off()
  
  
  # Enough sequencing depth to retrieve all expressed genes?
  
  #png(file="00000_detected_features_DNA.png")
  # par(mfrow = c(1,1))
  # mysaturation = dat(mydata,  type="saturation", ndepth=5)
  # explo.plot(mysaturation, toplot=1, Idxs=1:5 ,yleftlim = NULL, yrightlim = NULL)
  #dev.off()
  
  if (makepdf == TRUE)
  {
    dev.off()
  }
  
  
  #### SAVE ----
  
  
  
  filename18<-paste(type,"_rpkm_counts.tsv",sep='')
  
  cat("\n")
  cat(sprintf(as.character(filename18)))
  cat("\n")
  
  write.table(rpkm.counts, file=filename18, sep="\t", quote =F)
  
  
  filename18<-paste(type,"_tmm_counts.tsv",sep='')
  
  cat("\n")
  cat(sprintf(as.character(filename18)))
  cat("\n")
  
  write.table(tmm.counts, file=filename18, sep="\t", quote =F)
  
  
  
  filename18<-paste(type,"_tmm_plus_length_counts.tsv",sep='')
  
  cat("\n")
  cat(sprintf(as.character(filename18)))
  cat("\n")
  
  write.table(tmm.counts.length, file=filename18, sep="\t", quote =F)
  
}

GENE_EXP_and_TRatios = function(option_list)
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
  
  path2<-paste(out,type,'/',sep='')
  
  cat("path2\n")
  cat(sprintf(as.character(path2)))
  cat("\n")
  
  if (file.exists(path2)){
    
    setwd(path2)
  } else {
    dir.create(file.path(path2))
    setwd(path2)
    
  }
  
  #### Read STAR_kolf2_1 ----
  
  Equiv.table<-as.data.frame(fread(opt$Equiv.table,
                                    sep="\t", header=F), stringsAsFactors = F)
  
  colnames(Equiv.table)<-c("ensembl_gene_id","HGNC","ensembl_transcript_id","Transcript_name")
  
  
  cat("Equiv.table:\n")
  cat(str(Equiv.table))
  cat("\n")
  
   
 
  #### Read expression matrix ----
  
  filename18<-paste(type,"_tmm_plus_length_counts.tsv",sep='')
  
  cat("\n")
  cat(sprintf(as.character(filename18)))
  cat("\n")
  
  EXP.df<-as.data.frame(fread(file=filename18, sep="\t", 
                              header =F), stringsAsFactors=F)
  
  colnames(EXP.df)<-c("ensembl_transcript_id","Rep_1","Rep_2")
  
  
  cat("EXP.df:0\n")
  cat(str(EXP.df))
  cat("\n")
  
  #### Merge with Equiv.table ----
  
  EXP.df<-merge(EXP.df,
                Equiv.table,
        by="ensembl_transcript_id",
        all.x=T)
  
  cat("EXP.df:1\n")
  cat(str(EXP.df))
  cat("\n")
  
  #### Mean Transcript EXP ----
  mean.matrix.Transcript_EXP<-EXP.df[,c(which(colnames(EXP.df) == "Rep_1"),
                                        which(colnames(EXP.df) == "Rep_2"))]
  
  cat("--------------------------------->mean.matrix.Transcript_EXP\n")
  cat(str(mean.matrix.Transcript_EXP))
  cat("\n")
  
  mean.matrix.Transcript_EXP.mean<-apply(mean.matrix.Transcript_EXP,1,mean)
  
  # cat("--------------------------------->mean.matrix.Transcript_EXP.mean\n")
  # cat(str(mean.matrix.Transcript_EXP.mean))
  # cat("\n")
  
  mean.matrix.Transcript_EXP.sd<-apply(mean.matrix.Transcript_EXP,1,sd)
  
  # cat("--------------------------------->mean.matrix.Transcript_EXP.sd\n")
  # cat(str(mean.matrix.Transcript_EXP.sd))
  # cat("\n")
  
  EXP.df$mean_Transcript_EXP<-mean.matrix.Transcript_EXP.mean
  EXP.df$sd_Transcript_EXP<-mean.matrix.Transcript_EXP.sd
  EXP.df$cv_Transcript_EXP<-100*(EXP.df$sd_Transcript_EXP/EXP.df$mean_Transcript_EXP)
  
  cat("EXP.df:1.5\n")
  cat(str(EXP.df))
  cat("\n")
  
  # quit(status=1)
  #### GENE_EXP ----
  
  Matrix_file.dt<-data.table(EXP.df, key="ensembl_gene_id")
  
  
  
  GENE_EXP<-as.data.frame(Matrix_file.dt[,.(GENE_EXP_Rep_1=sum(Rep_1),
                                            GENE_EXP_Rep_2=sum(Rep_2))
                                         ,by=.(ensembl_gene_id)])
  
  cat("GENE_EXP:0\n")
  cat(str(GENE_EXP))
  cat("\n")
  
  mean.matrix.GENE_EXP<-GENE_EXP[,-1]
  
  cat("--------------------------------->mean.matrix.GENE_EXP\n")
  cat(str(mean.matrix.GENE_EXP))
  cat("\n")
  
  mean.matrix.GENE_EXP.mean<-apply(mean.matrix.GENE_EXP,1,mean)
  
  # cat("--------------------------------->mean.matrix.GENE_EXP.mean\n")
  # cat(str(mean.matrix.GENE_EXP.mean))
  # cat("\n")
  
  mean.matrix.GENE_EXP.sd<-apply(mean.matrix.GENE_EXP,1,sd)
  
  # cat("--------------------------------->mean.matrix.GENE_EXP.sd\n")
  # cat(str(mean.matrix.GENE_EXP.sd))
  # cat("\n")
  
  GENE_EXP$mean_GENE_EXP<-mean.matrix.GENE_EXP.mean
  GENE_EXP$sd_GENE_EXP<-mean.matrix.GENE_EXP.sd
  GENE_EXP$cv_GENE_EXP<-100*(GENE_EXP$sd/GENE_EXP$mean)
  
  
  
  EXP.df<-merge(EXP.df,
                GENE_EXP,
                by="ensembl_gene_id",
                all.x=T)
  
  cat("EXP.df:2\n")
  cat(str(EXP.df))
  cat("\n")
  
  #### Transcript Ratio ----
  
  Matrix_file.dt<-data.table(EXP.df, key="ensembl_gene_id")
  
  
  
  Ratios<-Matrix_file.dt[,.(ensembl_transcript_id,
                            Transcript_Ratio_Rep_1=Rep_1/GENE_EXP_Rep_1,
                            Transcript_Ratio_Rep_2=Rep_2/GENE_EXP_Rep_2)
                         ,by=.(ensembl_gene_id)]
  
  cat("Ratios:0\n")
  cat(str(Ratios))
  cat("\n")
  
  mean.matrix.Ratios<-Ratios[,-c(1:2)]
  
  cat("--------------------------------->mean.matrix.Ratios\n")
  cat(str(mean.matrix.Ratios))
  cat("\n")
  
  mean.matrix.Ratios.mean<-apply(mean.matrix.Ratios,1,mean)
  
  # cat("--------------------------------->mean.matrix.Ratios.mean\n")
  # cat(str(mean.matrix.Ratios.mean))
  # cat("\n")
  
  mean.matrix.Ratios.sd<-apply(mean.matrix.Ratios,1,sd)
  
  # cat("--------------------------------->mean.matrix.Ratios.sd\n")
  # cat(str(mean.matrix.Ratios.sd))
  # cat("\n")
  
  Ratios$mean_Transcript_Ratio<-mean.matrix.Ratios.mean
  Ratios$sd_Transcript_Ratio<-mean.matrix.Ratios.sd
  Ratios$cv_Transcript_Ratio<-100*(Ratios$sd/Ratios$mean)
  
  EXP.df<-merge(EXP.df,
                Ratios,
                by=c("ensembl_gene_id","ensembl_transcript_id"),
                all.x=T)
  
  cat("EXP.df:3\n")
  cat(str(EXP.df))
  cat("\n")
  
  # quit(status=1)
  
  #### SAVE ----
  
  setwd(path2)
  
  
  filename18<-paste(type,"_object.tsv",sep='')
  
  cat("\n")
  cat(sprintf(as.character(filename18)))
  cat("\n")
  
  write.table(EXP.df, file=filename18, sep="\t", row.names=F,quote =F)
 
}

ploting_GENE_EXP = function(option_list)
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
  
  ###### Input files -----
  
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
  cat(sprintf(as.character(unique(Gather.m$HGNC))))
  cat("\n")
  
  Gather.m_subset<-Gather.m[which(Gather.m$HGNC%in%c("CUX1","C2CD5","FOXP1","SH2B3","UGCG","BID","BRAP","TNRC6A","NBN","EPB41")),]
  
  cat("Gather.m_subset\n")
  str(Gather.m_subset)
  cat("\n")
  cat(sprintf(as.character(levels(Gather.m_subset$variable))))
  cat("\n")
  cat(sprintf(as.character(levels(as.factor(Gather.m_subset$Rep)))))
  cat("\n")
  cat(sprintf(as.character(unique(Gather.m_subset$HGNC))))
  cat("\n")
  
  
  
  
  #### READ and transform SELECTED_GENES ----

  SELECTED_GENES = unlist(strsplit(split =",", opt$SELECTED_GENES))

  cat("SELECTED_GENES_\n")
  cat(sprintf(as.character(SELECTED_GENES)))
  cat("\n")
  # quit(status=1)
  # 
  # #### READ and transform EDITION_DROPOUTS ----
  # 
  # EDITION_DROPOUTS = unlist(strsplit(split =",", opt$EDITION_DROPOUTS))
  # 
  # cat("EDITION_DROPOUTS_\n")
  # cat(sprintf(as.character(EDITION_DROPOUTS)))
  # cat("\n")
  # # quit(status=1)
  # 
  # #### READ and transform EXPRESSION_DROPOUTS ----
  # 
  # EXPRESSION_DROPOUTS = unlist(strsplit(split =",", opt$EXPRESSION_DROPOUTS))
  # 
  # cat("EXPRESSION_DROPOUTS_\n")
  # cat(sprintf(as.character(EXPRESSION_DROPOUTS)))
  # cat("\n")
  
  ###### Path 2 -----
  
  path2<-paste(out,type,'/',sep='')
  
  cat("path2\n")
  cat(sprintf(as.character(path2)))
  cat("\n")
  
  if (file.exists(path2)){
    
    setwd(path2)
  } else {
    dir.create(file.path(path2))
    setwd(path2)
    
  }
  
  #### Read expression matrix ----
  
  filename18<-paste(type,"_object.tsv",sep='')
  
  cat("\n")
  cat(sprintf(as.character(filename18)))
  cat("\n")
  
  EXP.df<-as.data.frame(fread(file=filename18, sep="\t", 
                              header =T), stringsAsFactors=F)
  
  cat("EXP.df\n")
  cat(str(EXP.df))
  cat("\n")
  
  EXP.df<-EXP.df[which(EXP.df$mean_GENE_EXP > 0),]
  
  cat("EXP.df\n")
  cat(str(EXP.df))
  cat("\n")
  
  # quit(status = 1)
  
  #### GENE_EXP plot ----
  
  indx.int<-c(which(colnames(EXP.df) == "ensembl_gene_id"),
              which(colnames(EXP.df) == "HGNC"),
              which(colnames(EXP.df) == "mean_GENE_EXP"))
  
  GENE_EXP<-unique(EXP.df[,indx.int])
  
  GENE_EXP$MOCK_COORD <- "ALL GENES"
  
  GENE_EXP$MOCK_COORD[which(GENE_EXP$HGNC%in%Gather.m_subset$HGNC)] <- "genIE host genes"

  GENE_EXP$MOCK_COORD<-factor(GENE_EXP$MOCK_COORD,
                              levels=c("ALL GENES", "genIE host genes"),
                              ordered=T)
  
  cat("GENE_EXP\n")
  cat(str(GENE_EXP))
  cat("\n")
  
  cat(sprintf(as.character(names(summary(GENE_EXP$MOCK_COORD)))))
  cat("\n")
  cat(sprintf(as.character(summary(GENE_EXP$MOCK_COORD))))
  cat("\n")
  
  # #######################################################
  # quit(status = 1)
  
  ### x and y plot ----
  
  Max_CS<-max(GENE_EXP$mean_GENE_EXP[!is.na(GENE_EXP$mean_GENE_EXP)])
  log_10_Max_CS<-log10(Max_CS+0.0001)
  Min_CS<-min(GENE_EXP$mean_GENE_EXP[!is.na(GENE_EXP$mean_GENE_EXP)])
  log_10_Min_CS<-log10(Min_CS+0.0001)
  
  
  cat("Max_CS_\n")
  cat(sprintf(as.character(Max_CS)))
  cat("\n")
  
  cat("log_10_Max_CS_\n")
  cat(sprintf(as.character(log_10_Max_CS)))
  cat("\n")
  
  cat("Min_CS_\n")
  cat(sprintf(as.character(Min_CS)))
  cat("\n")
  
  cat("log_10_Min_CS_\n")
  cat(sprintf(as.character(log_10_Min_CS)))
  cat("\n")
  
   # quit(status = 1)
  
 # breaks.GENE_EXP<-round(seq(Min_CS,(Max_CS+1),by=1),2)
  breaks.GENE_EXP.log<-unique(c(seq(log_10_Min_CS,log_10_Max_CS,by=1),log_10_Max_CS,log_10_Max_CS+0.2))
  breaks.GENE_EXP<-10^(breaks.GENE_EXP.log)
  
  breaks.GENE_EXP[(breaks.GENE_EXP >= 1)]<-round(breaks.GENE_EXP[(breaks.GENE_EXP >= 1)],0)
  breaks.GENE_EXP[(breaks.GENE_EXP < 1)]<-round(breaks.GENE_EXP[(breaks.GENE_EXP < 1)],3)
  
  labels.GENE_EXP<-as.character(breaks.GENE_EXP)
  
  
  cat("breaks.GENE_EXP.log_\n")
  cat(sprintf(as.character(breaks.GENE_EXP.log)))
  cat("\n")
  
  cat("breaks.GENE_EXP_\n")
  cat(sprintf(as.character(breaks.GENE_EXP)))
  cat("\n")
  cat("labels.GENE_EXP_\n")
  cat(sprintf(as.character(labels.GENE_EXP)))
  cat("\n")
  
 
  
  ### GENE_EXP_ALL ----
  
  cat("Hello_world1\n")
  
  # GENE_EXP  %>%
  # mutate(myaxis = paste0(Fig3_Annot_Category, "\n", "n=", Total_Fig3_Annot_Category), drop=F) %>%
  #   mutate(myaxis=fct_reorder(myaxis,as.numeric(Fig3_Annot_Category)), drop=F) %>%
  
  
  graph_GENE_EXP<-ggplot() +
    geom_violin(data=GENE_EXP[which(GENE_EXP$MOCK_COORD == "ALL GENES"),],
                aes(x=MOCK_COORD, y=log10(mean_GENE_EXP + 0.0001)), color="gray", fill="gray")+
    geom_boxplot(data=GENE_EXP[which(GENE_EXP$MOCK_COORD == "ALL GENES"),], 
                 aes(x=MOCK_COORD, y=log10(mean_GENE_EXP + 0.0001)), color="gray", width = 0.2)+
    geom_point(data=GENE_EXP[which(GENE_EXP$MOCK_COORD == "genIE host genes"),], 
              aes(x=MOCK_COORD, y=log10(mean_GENE_EXP + 0.0001), color=MOCK_COORD),
              size=4)+
    theme(plot.title = element_text(size=11)) +
    theme(plot.title=element_text(size=11))+
    theme_bw()+
    scale_x_discrete(name=NULL, drop=F)+
    theme(axis.text.x=element_text(angle=45,size=14,colour="black",vjust=1,hjust=1,face="bold"), 
          axis.title.y=element_text(size=16),
          legend.title=element_text(size=16),
          legend.text=element_text(size=12,colour="black"),
          axis.text.y=element_text(colour="black",size=16))+
    scale_y_continuous(name="mean Gene Expression (rpkm)",breaks=breaks.GENE_EXP.log,labels=labels.GENE_EXP, 
                       limits=c(breaks.GENE_EXP.log[1],breaks.GENE_EXP.log[length(breaks.GENE_EXP.log)]))+
    labs(color=paste("Gene","Expression","(rpkm)", sep="\n"))+
    scale_color_manual(values = c('#32A852','#1877C9','#553B68','#D45E85','#6DB2EE','#62D07F','#C9244B','#87447B','#D6B8E6'),
                       drop=F)+
    geom_hline(yintercept=log10(10+0.0001), linetype='dashed', col = 'black')+
    theme(legend.key=element_blank(), legend.key.size=unit(1,"point"), 
          legend.position="none")+
    ggeasy::easy_center_title()
  
  graph_GENE_EXP<-graph_GENE_EXP+
                  geom_text_repel(data=GENE_EXP[which(GENE_EXP$MOCK_COORD == "genIE host genes"),],
                                  aes(x=MOCK_COORD,y=log10(mean_GENE_EXP + 0.0001),label=HGNC),
                                  box.padding = 0.8, max.overlaps = Inf,segment.linetype = 1)
  
  setwd(path2)
  
  # geom_mark_ellipse(data=GENE_EXP[which(GENE_EXP$MOCK_COORD == "genIE host genes"),], 
  #                   aes(x=MOCK_COORD,y=log10(mean_GENE_EXP + 0.0001),label=HGNC))
  cat("Hello_world2\n")
  
  svglite(paste('Kolf2_GENE_EXP','.svg',sep=''), width = 8, height = 8)
  print(graph_GENE_EXP)
  dev.off()
  
  
  # quit(status=1)
  
 
  ### GENE_EXP_DNMT1andS1PRD ----
  
  GENE_EXP$MOCK_COORD_2 <- "ALL GENES"
  
  GENE_EXP$MOCK_COORD_2[which(GENE_EXP$HGNC%in%SELECTED_GENES)] <- "DNMT1 overlapping"
  
  GENE_EXP$MOCK_COORD_2<-factor(GENE_EXP$MOCK_COORD_2,
                              levels=c("ALL GENES", "DNMT1 overlapping"),
                              ordered=T)
  
  cat("GENE_EXP\n")
  cat(str(GENE_EXP))
  cat("\n")
  
  cat(sprintf(as.character(names(summary(GENE_EXP$MOCK_COORD_2)))))
  cat("\n")
  cat(sprintf(as.character(summary(GENE_EXP$MOCK_COORD_2))))
  cat("\n")
  
  
   # quit(status = 1)
  
  cat("Hello_world1\n")
  

  
  graph_GENE_EXP<-ggplot() +
    geom_violin(data=GENE_EXP[which(GENE_EXP$MOCK_COORD_2 == "ALL GENES"),],
              aes(x=MOCK_COORD_2, y=log10(mean_GENE_EXP + 0.0001)), color="gray", fill="gray")+
    geom_boxplot(data=GENE_EXP[which(GENE_EXP$MOCK_COORD_2 == "ALL GENES"),],
                 aes(x=MOCK_COORD_2, y=log10(mean_GENE_EXP + 0.0001), color=MOCK_COORD_2), width = 0.2)+
    geom_point(data=GENE_EXP[which(GENE_EXP$MOCK_COORD_2 == "DNMT1 overlapping"),],
               aes(x=MOCK_COORD_2, y=log10(mean_GENE_EXP + 0.0001), color=MOCK_COORD_2), size=5)+
    theme_bw()+
    theme(plot.title = element_text(size=11)) +
    theme(plot.title=element_text(size=11))+
    scale_x_discrete(name=NULL, drop=F)+
    theme(axis.text.x=element_text(angle=45,size=14,colour="black",vjust=1,hjust=1,face="bold"),
          axis.title.y=element_text(size=16,face="bold"),
          legend.title=element_text(size=16,face="bold"),
          legend.text=element_text(size=12,face="bold",colour="black"),
          axis.text.y=element_text(colour="black",size=12,face="bold"))+
    scale_y_continuous(name="mean Gene Expression (rpkm)",breaks=breaks.GENE_EXP.log,labels=labels.GENE_EXP,
                       limits=c(breaks.GENE_EXP.log[1],breaks.GENE_EXP.log[length(breaks.GENE_EXP.log)]))+
    labs(color=paste("Gene","Expression","(rpkm)", sep="\n"))+
    scale_color_manual(values = c('#32A852','#1877C9','#553B68','#D45E85','#6DB2EE','#62D07F','#C9244B','#87447B','#D6B8E6'),
                       drop=F)+
    geom_hline(yintercept=log10(10+0.0001), linetype='dashed', col = 'black')+
    theme(legend.key=element_blank(), legend.key.size=unit(1,"point"),
          legend.position="none")+
    ggeasy::easy_center_title()

  graph_GENE_EXP<-graph_GENE_EXP+
    geom_text_repel(data=GENE_EXP[which(GENE_EXP$MOCK_COORD_2 == "DNMT1 overlapping"),],
                    aes(x=MOCK_COORD_2,y=log10(mean_GENE_EXP + 0.0001),label=HGNC),
                    box.padding = 0.8, max.overlaps = Inf,segment.linetype = 6)

  setwd(path2)

  # geom_mark_ellipse(data=GENE_EXP[which(GENE_EXP$MOCK_COORD_2 == "DNMT1 overlapping"),],
  #                   aes(x=MOCK_COORD_2,y=log10(mean_GENE_EXP + 0.0001),label=HGNC))
  cat("Hello_world2\n")

  svglite(paste('Kolf2_GENE_EXP_DNMT1_overlaping','.svg',sep=''), width = 8, height = 8)
  print(graph_GENE_EXP)
  dev.off()
  # 
  # 
  # # quit(status=1)
  
  
}

ploting_Transcript_Ratio = function(option_list)
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
  
  ###### Input files -----
  
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
  
  ###### Path 2 -----
  
  path2<-paste(out,type,'/',sep='')
  
  cat("path2\n")
  cat(sprintf(as.character(path2)))
  cat("\n")
  
  if (file.exists(path2)){
    
    setwd(path2)
  } else {
    dir.create(file.path(path2))
    setwd(path2)
    
  }
  
  #### Read expression matrix ----
  
  filename18<-paste(type,"_object.tsv",sep='')
  
  cat("\n")
  cat(sprintf(as.character(filename18)))
  cat("\n")
  
  EXP.df<-as.data.frame(fread(file=filename18, sep="\t", 
                              header =T), stringsAsFactors=F)
  
  cat("EXP.df\n")
  cat(str(EXP.df))
  cat("\n")
  
  EXP.df<-EXP.df[which(EXP.df$mean_GENE_EXP > 0),]
  
  cat("EXP.df\n")
  cat(str(EXP.df))
  cat("\n")
  
  #### READ and transform VEP_route ----
  
  VEP_route = opt$VEP_route
  
  cat("VEP_route_\n")
  cat(sprintf(as.character(VEP_route)))
  cat("\n")
  
  setwd(VEP_route)
  
  
  file_list <- list.files(path=VEP_route, include.dirs = FALSE)
  
  cat("file_list\n")
  cat(str(file_list))
  cat("\n")
  
  indexes_sel <- grep("VEP_consequence_Plus_Transcript_ID\\.csv$",file_list)
  
  cat("indexes_sel\n")
  cat(sprintf(as.character(indexes_sel)))
  cat("\n")
  
  file_list_sel <- as.data.frame(file_list[indexes_sel], stringsAsFactors=F)
  
  colnames(file_list_sel)<-"file"
  
  cat("file_list_sel\n")
  cat(str(file_list_sel))
  cat("\n")
  
  
  VEP_list<-list()
  
  if(dim(file_list_sel)[1] > 0)
  {
    for(i in 1:dim(file_list_sel)[1])
    {
      candidate_file <- file_list_sel$file[i]
      
      # cat("candidate_file\n")
      # cat(sprintf(as.character(candidate_file)))
      # cat("\n")
      
      output_vep <- as.data.frame(fread(file=candidate_file,
                          sep=",", header=T), stringsAsFactors = F)
      
      # cat("output_vep\n")
      # cat(str(output_vep))
      # cat("\n")
      
      indx.dep<-c(which(colnames(output_vep) == "Rank_EXP_in_Cell_Type"),
                  which(colnames(output_vep) == "value"),
                  which(colnames(output_vep) == "BP_CELL_LABELS"),
                  which(colnames(output_vep) == "pos37"))
      
      output_vep_subset<-unique(output_vep[,-indx.dep])
      
      # cat("output_vep_subset\n")
      # cat(str(output_vep_subset))
      # cat("\n")
      
      VEP_list[[i]]<-output_vep_subset
      
      
      # quit(status=1)
      
      
    }
    
    
  }
  
  
  VEP_df = as.data.frame(unique(data.table::rbindlist(VEP_list, fill=T), stringsAsFactors=F))
  
  colnames(VEP_df)[which(colnames(VEP_df) == "transcript_id")]<-"ensembl_transcript_id"
  
  cat("VEP_df\n")
  cat(str(VEP_df))
  cat("\n")
  
  VEP_df_subset<-VEP_df[which(VEP_df$VAR%in%Gather.m$VAR),]
  
  
  cat("VEP_df_subset\n")
  cat(str(VEP_df_subset))
  cat("\n")
  
  
  
  HGNC_intronic<-unique(VEP_df_subset[which(VEP_df_subset$VEP_DEF_LABELS == "INTRON"),-c(which(colnames(VEP_df_subset) == "ensembl_transcript_id"),
                                                                    which(colnames(VEP_df_subset) == "VAR"))])
  
  cat("HGNC_intronic\n")
  cat(str(HGNC_intronic))
  cat("\n")
  
  
  VEP_df_subset_selected<-VEP_df_subset[which(VEP_df_subset$ensembl_gene_id%in%HGNC_intronic$ensembl_gene_id),]
  
  cat("VEP_df_subset_selected\n")
  cat(str(VEP_df_subset_selected))
  cat("\n")
  
  ##### Merge with Transcript Ratio object -----
  
  DEF<-merge(VEP_df_subset_selected,
             EXP.df,
             by=c("ensembl_transcript_id","HGNC","ensembl_gene_id"),
             all.x=T)
  
  DEF$VEP_DEF_LABELS<-factor(DEF$VEP_DEF_LABELS,
                                       levels=c("LOF","MISS","SYN","UTR5","UTR3",
                                                "INTRON","INTERGENIC","UPSTREAM","DOWNSTREAM","REGULATORY",
                                                "TFBS","SPLICE","OTHER","NMD","NCT","PCHIC_Relevant_link"),ordered = T)
  
  
  # check<-unique(DEF$VEP_DEF_LABELS[is.na(DEF$VEP_DEF_LABELS_2)])
  # 
  cat("DEF\n")
  cat(str(DEF))
  cat("\n")
  # cat("check_\n")
  # cat(sprintf(as.character(check)))
  # cat("\n")
  # 
  # quit(status=1)
  
  breaks.Rank<-seq(0,1,by=0.1)
  labels.Rank<-as.character(breaks.Rank)
  
  
  cat("labels.Rank_\n")
  cat(sprintf(as.character(labels.Rank)))
  cat("\n")
  
  
  
  p <- ggplot(data=DEF,
              aes(x=HGNC,y=mean_Transcript_Ratio,
                  color=VEP_DEF_LABELS)) +
    geom_jitter(size = 3,width = 0.05)+
    scale_color_manual(values=c('#32A852','#6DB2EE','#62D07F','#1877C9','#C9244B','#D45E85','#87447B','#553B68','#D6B8E6',
                               "#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "gray", "#CBD588", "#5F7FC7", 
                               "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
                               "#D14285", "#6DDE88", "#652926", "#D1A33D", "#C84248","#89C5DA", "#DA5724", "#74D944", "#CE50CA"),drop=F)+
    scale_y_continuous(name="mean Transcript Ratio",breaks=breaks.Rank,labels=labels.Rank, 
                       limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
    scale_x_discrete(name=element_blank(), drop =FALSE)+
    theme_bw()+
    theme(axis.text.x=element_text(angle=45,size=14,colour="black",vjust=1,hjust=1,face="bold"), 
          axis.title.y=element_text(size=16,face="bold"),
          legend.title=element_text(size=16,face="bold"),
          legend.text=element_text(size=12,face="bold",colour="black"),
          axis.text.y=element_text(colour="black",size=12,face="bold"))+
    labs(color=paste("Transcript","consequences", sep="\n"))+
    theme(legend.key=element_blank(), legend.key.size=unit(1,"point"), 
          legend.position="right")+
    ggeasy::easy_center_title()
    
  setwd(path2)
  
  # geom_mark_ellipse(data=GENE_EXP[which(GENE_EXP$MOCK_COORD == "genIE host genes"),], 
  #                   aes(x=MOCK_COORD,y=log10(mean_GENE_EXP + 0.0001),label=HGNC))
  cat("Hello_world2\n")
  
  svglite(paste('Kolf2_Transcript_Ratio','.svg',sep=''), width = 8, height = 8)
  print(p)
  dev.off()
  
  # quit(status=1)
  
  ##### Transcript EXP -----
  
  
  
  Max_Transcript_EXP<-max(DEF$mean_Transcript_EXP[!is.na(DEF$mean_Transcript_EXP)])
  log_10_Max_Transcript_EXP<-log10(Max_Transcript_EXP+0.0001)
  Min_Transcript_EXP<-min(DEF$mean_Transcript_EXP[!is.na(DEF$mean_Transcript_EXP)])
  log_10_Min_Transcript_EXP<-log10(Min_Transcript_EXP+0.0001)
  
  
  cat("Max_Transcript_EXP_\n")
  cat(sprintf(as.character(Max_Transcript_EXP)))
  cat("\n")
  
  cat("log_10_Max_Transcript_EXP_\n")
  cat(sprintf(as.character(log_10_Max_Transcript_EXP)))
  cat("\n")
  
  cat("Min_Transcript_EXP_\n")
  cat(sprintf(as.character(Min_Transcript_EXP)))
  cat("\n")
  
  cat("log_10_Min_Transcript_EXP_\n")
  cat(sprintf(as.character(log_10_Min_Transcript_EXP)))
  cat("\n")
  
  # quit(status = 1)
  
  # breaks.Transcript_EXP<-round(seq(Min_Transcript_EXP,(Max_Transcript_EXP+1),by=1),2)
  breaks.Transcript_EXP.log<-unique(c(seq(log_10_Min_Transcript_EXP,log_10_Max_Transcript_EXP,by=1),log_10_Max_Transcript_EXP,log_10_Max_Transcript_EXP+0.2))
  breaks.Transcript_EXP<-10^(breaks.Transcript_EXP.log)
  labels.Transcript_EXP<-as.character(round(breaks.Transcript_EXP,3))
  
  
  cat("breaks.Transcript_EXP.log_\n")
  cat(sprintf(as.character(breaks.Transcript_EXP.log)))
  cat("\n")
  
  cat("breaks.Transcript_EXP_\n")
  cat(sprintf(as.character(breaks.Transcript_EXP)))
  cat("\n")
  cat("labels.Transcript_EXP_\n")
  cat(sprintf(as.character(labels.Transcript_EXP)))
  cat("\n")
  
 
  # quit(status = 1)
  
  
  
  p <- ggplot(data=DEF,
              aes(x=HGNC,y=log10(mean_Transcript_EXP+0.0001),
                  color=VEP_DEF_LABELS)) +
    geom_jitter(size = 3,width = 0.05)+
    scale_color_manual(values=c('#32A852','#6DB2EE','#62D07F','#1877C9','#C9244B','#D45E85','#87447B','#553B68','#D6B8E6',
                                "#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "gray", "#CBD588", "#5F7FC7", 
                                "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
                                "#D14285", "#6DDE88", "#652926", "#D1A33D", "#C84248","#89C5DA", "#DA5724", "#74D944", "#CE50CA"),drop=F)+
    scale_y_continuous(name="mean Transcript Expression (rpkm)",
                       breaks=breaks.Transcript_EXP.log,
                       labels=labels.Transcript_EXP, 
                       limits=c(breaks.Transcript_EXP.log[1],breaks.Transcript_EXP.log[length(breaks.Transcript_EXP.log)]))+
    scale_x_discrete(name=element_blank(), drop =FALSE)+
    theme_bw()+
    theme(axis.text.x=element_text(angle=45,size=14,colour="black",vjust=1,hjust=1,face="bold"), 
          axis.title.y=element_text(size=16,face="bold"),
          legend.title=element_text(size=16,face="bold"),
          legend.text=element_text(size=12,face="bold",colour="black"),
          axis.text.y=element_text(colour="black",size=12,face="bold"))+
    labs(color=paste("Transcript","consequences", sep="\n"))+
    theme(legend.key=element_blank(), legend.key.size=unit(1,"point"), 
          legend.position="right")+
    geom_hline(yintercept=log10(10+0.0001), linetype='dotted', col = 'red')+
    ggeasy::easy_center_title()
  
  setwd(path2)
  
  # geom_mark_ellipse(data=GENE_EXP[which(GENE_EXP$MOCK_COORD == "genIE host genes"),], 
  #                   aes(x=MOCK_COORD,y=log10(mean_GENE_EXP + 0.0001),label=HGNC))
  cat("Hello_world2\n")
  
  svglite(paste('Kolf2_Transcript_EXP','.svg',sep=''), width = 8, height = 8)
  print(p)
  dev.off()
  
  ##### ADD BY CSQ ------
 
 
  Matrix_file.dt<-data.table(DEF, key=c("ensembl_gene_id","HGNC","VEP_DEF_LABELS","VAR"))
  
  
  
  CSQ_SUM<-as.data.frame(Matrix_file.dt[,.(CSQ_SUM_Rep_1=sum(Rep_1),
                                            CSQ_SUM_Rep_2=sum(Rep_2))
                                         ,by=.(ensembl_gene_id,HGNC,VEP_DEF_LABELS,VAR)])
  
  cat("CSQ_SUM:0\n")
  cat(str(CSQ_SUM))
  cat("\n")
  
  mean.matrix.CSQ_SUM<-CSQ_SUM[,-c(1:4)]
  
  cat("--------------------------------->mean.matrix.CSQ_SUM\n")
  cat(str(mean.matrix.CSQ_SUM))
  cat("\n")
  
  mean.matrix.CSQ_SUM.mean<-apply(mean.matrix.CSQ_SUM,1,mean)
  
  cat("mean.matrix.CSQ_SUM.mean:0\n")
  cat(str(mean.matrix.CSQ_SUM.mean))
  cat("\n")
  
  CSQ_SUM$mean_CSQ_SUM<-mean.matrix.CSQ_SUM.mean
  
  
  cat("CSQ_SUM:2\n")
  cat(str(CSQ_SUM))
  cat("\n")
  
  #### graph
  
  Max_CSQ_SUM<-max(CSQ_SUM$mean_CSQ_SUM[!is.na(CSQ_SUM$mean_CSQ_SUM)])
  log_10_Max_CSQ_SUM<-log10(Max_CSQ_SUM+0.0001)
  Min_CSQ_SUM<-min(CSQ_SUM$mean_CSQ_SUM[!is.na(CSQ_SUM$mean_CSQ_SUM)])
  log_10_Min_CSQ_SUM<-log10(Min_CSQ_SUM+0.0001)
  
  
  cat("Max_CSQ_SUM_\n")
  cat(sprintf(as.character(Max_CSQ_SUM)))
  cat("\n")
  
  cat("log_10_Max_CSQ_SUM_\n")
  cat(sprintf(as.character(log_10_Max_CSQ_SUM)))
  cat("\n")
  
  cat("Min_CSQ_SUM_\n")
  cat(sprintf(as.character(Min_CSQ_SUM)))
  cat("\n")
  
  cat("log_10_Min_CSQ_SUM_\n")
  cat(sprintf(as.character(log_10_Min_CSQ_SUM)))
  cat("\n")
  
  # quit(status = 1)
  
  # breaks.CSQ_SUM<-round(seq(Min_CSQ_SUM,(Max_CSQ_SUM+1),by=1),2)
  breaks.CSQ_SUM.log<-unique(c(seq(log_10_Min_CSQ_SUM,log_10_Max_CSQ_SUM,by=1),log_10_Max_CSQ_SUM,log_10_Max_CSQ_SUM+0.2))
  breaks.CSQ_SUM<-10^(breaks.CSQ_SUM.log)
  labels.CSQ_SUM<-as.character(round(breaks.CSQ_SUM,3))
  
  
  cat("breaks.CSQ_SUM.log_\n")
  cat(sprintf(as.character(breaks.CSQ_SUM.log)))
  cat("\n")
  
  cat("breaks.CSQ_SUM_\n")
  cat(sprintf(as.character(breaks.CSQ_SUM)))
  cat("\n")
  cat("labels.CSQ_SUM_\n")
  cat(sprintf(as.character(labels.CSQ_SUM)))
  cat("\n")
  
  
  # quit(status = 1)
  
  
  
  p <- ggplot(data=CSQ_SUM,
              aes(x=HGNC,y=log10(mean_CSQ_SUM+0.0001),
                  color=VEP_DEF_LABELS)) +
    geom_jitter(size = 3,width = 0.05)+
    scale_color_manual(values=c('#32A852','#6DB2EE','#62D07F','#1877C9','#C9244B','#D45E85','#87447B','#553B68','#D6B8E6',
                                "#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "gray", "#CBD588", "#5F7FC7", 
                                "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
                                "#D14285", "#6DDE88", "#652926", "#D1A33D", "#C84248","#89C5DA", "#DA5724", "#74D944", "#CE50CA"),drop=F)+
    scale_y_continuous(name="Aggregate transcript expression per consequence (rpkm)",
                       breaks=breaks.CSQ_SUM.log,
                       labels=labels.CSQ_SUM, 
                       limits=c(breaks.CSQ_SUM.log[1],breaks.CSQ_SUM.log[length(breaks.CSQ_SUM.log)]))+
    scale_x_discrete(name=element_blank(), drop =FALSE)+
    theme_bw()+
    theme(axis.text.x=element_text(angle=45,size=14,colour="black",vjust=1,hjust=1,face="bold"), 
          axis.title.y=element_text(size=16,face="bold"),
          legend.title=element_text(size=16,face="bold"),
          legend.text=element_text(size=12,face="bold",colour="black"),
          axis.text.y=element_text(colour="black",size=12,face="bold"))+
    labs(color=paste("Transcript","consequences", sep="\n"))+
    theme(legend.key=element_blank(), legend.key.size=unit(1,"point"), 
          legend.position="right")+
    geom_hline(yintercept=log10(10+0.0001), linetype='dotted', col = 'red')+
    ggeasy::easy_center_title()
  
  setwd(path2)
  
  # geom_mark_ellipse(data=GENE_EXP[which(GENE_EXP$MOCK_COORD == "genIE host genes"),], 
  #                   aes(x=MOCK_COORD,y=log10(mean_GENE_EXP + 0.0001),label=HGNC))
  cat("Hello_world2\n")
  
  svglite(paste('Kolf2_CSQ_SUM','.svg',sep=''), width = 8, height = 8)
  print(p)
  dev.off()
  
  
  
 # quit(status = 1)
  ##### SAVE -----
  
  write.table(DEF,file="Expression_Kolf2_genIE.tsv", sep="\t", row.names = F, quote = F)
  
  write.table(CSQ_SUM,file="Expression_Kolf2_genIE_aggregate_by_csq.tsv", sep="\t", row.names = F, quote = F)
  
  

}

ploting_Transcript_Ratio_DNMT1 = function(option_list)
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
  
  ###### Input files -----
  
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
  
  ###### Path 2 -----
  
  path2<-paste(out,type,'/',sep='')
  
  cat("path2\n")
  cat(sprintf(as.character(path2)))
  cat("\n")
  
  if (file.exists(path2)){
    
    setwd(path2)
  } else {
    dir.create(file.path(path2))
    setwd(path2)
    
  }
  
  #### Read expression matrix ----
  
  filename18<-paste(type,"_object.tsv",sep='')
  
  cat("\n")
  cat(sprintf(as.character(filename18)))
  cat("\n")
  
  EXP.df<-as.data.frame(fread(file=filename18, sep="\t", 
                              header =T), stringsAsFactors=F)
  
  cat("EXP.df\n")
  cat(str(EXP.df))
  cat("\n")
  
  EXP.df<-EXP.df[which(EXP.df$mean_GENE_EXP > 0),]
  
  cat("EXP.df\n")
  cat(str(EXP.df))
  cat("\n")
  
  #### READ and transform SELECTED_GENES ----
  
  SELECTED_GENES = unlist(strsplit(split =",", opt$SELECTED_GENES))
  
  cat("SELECTED_GENES_\n")
  cat(sprintf(as.character(SELECTED_GENES)))
  cat("\n")
  
  #### READ and transform VEP_route ----
  
  VEP_route = opt$VEP_route
  
  cat("VEP_route_\n")
  cat(sprintf(as.character(VEP_route)))
  cat("\n")
  
  setwd(VEP_route)
  
  
  file_list <- list.files(path=VEP_route, include.dirs = FALSE)
  
  cat("file_list\n")
  cat(str(file_list))
  cat("\n")
  
  indexes_sel <- grep("VEP_consequence_Plus_Transcript_ID\\.csv$",file_list)
  
  cat("indexes_sel\n")
  cat(sprintf(as.character(indexes_sel)))
  cat("\n")
  
  file_list_sel <- as.data.frame(file_list[indexes_sel], stringsAsFactors=F)
  
  colnames(file_list_sel)<-"file"
  
  cat("file_list_sel\n")
  cat(str(file_list_sel))
  cat("\n")
  
  
  VEP_list<-list()
  
  if(dim(file_list_sel)[1] > 0)
  {
    for(i in 1:dim(file_list_sel)[1])
    {
      candidate_file <- file_list_sel$file[i]
      
      # cat("candidate_file\n")
      # cat(sprintf(as.character(candidate_file)))
      # cat("\n")
      
      output_vep <- as.data.frame(fread(file=candidate_file,
                                        sep=",", header=T), stringsAsFactors = F)
      
      # cat("output_vep\n")
      # cat(str(output_vep))
      # cat("\n")
      
      indx.dep<-c(which(colnames(output_vep) == "Rank_EXP_in_Cell_Type"),
                  which(colnames(output_vep) == "value"),
                  which(colnames(output_vep) == "BP_CELL_LABELS"),
                  which(colnames(output_vep) == "pos37"))
      
      output_vep_subset<-unique(output_vep[,-indx.dep])
      
      # cat("output_vep_subset\n")
      # cat(str(output_vep_subset))
      # cat("\n")
      
      VEP_list[[i]]<-output_vep_subset
      
      
      # quit(status=1)
      
      
    }
    
    
  }
  
  
  VEP_df = as.data.frame(unique(data.table::rbindlist(VEP_list, fill=T), stringsAsFactors=F))
  
  colnames(VEP_df)[which(colnames(VEP_df) == "transcript_id")]<-"ensembl_transcript_id"
  
  cat("VEP_df\n")
  cat(str(VEP_df))
  cat("\n")
  
  VEP_df_subset<-VEP_df[which(VEP_df$VAR%in%Gather.m$VAR),]
  
  
  cat("VEP_df_subset\n")
  cat(str(VEP_df_subset))
  cat("\n")
  
  
  
  # HGNC_intronic<-unique(VEP_df_subset[which(VEP_df_subset$VEP_DEF_LABELS == "INTRON"),-c(which(colnames(VEP_df_subset) == "ensembl_transcript_id"),
  #                                                                                        which(colnames(VEP_df_subset) == "VAR"))])
  # 
  # cat("HGNC_intronic\n")
  # cat(str(HGNC_intronic))
  # cat("\n")
  
  
  VEP_df_subset_selected<-VEP_df_subset[which(VEP_df_subset$HGNC%in%SELECTED_GENES),]
  
  cat("VEP_df_subset_selected\n")
  cat(str(VEP_df_subset_selected))
  cat("\n")
  
  # quit(status=1)
  
  ##### Merge with Transcript Ratio object -----
  
  DEF<-merge(VEP_df_subset_selected,
             EXP.df,
             by=c("ensembl_transcript_id","HGNC","ensembl_gene_id"),
             all.x=T)
  
  DEF$VEP_DEF_LABELS<-factor(DEF$VEP_DEF_LABELS,
                             levels=c("LOF","MISS","SYN","UTR5","UTR3",
                                      "INTRON","INTERGENIC","UPSTREAM","DOWNSTREAM","REGULATORY",
                                      "TFBS","SPLICE","OTHER","NMD","NCT","PCHIC_Relevant_link"),ordered = T)
  
  
  # check<-unique(DEF$VEP_DEF_LABELS[is.na(DEF$VEP_DEF_LABELS_2)])
  # 
  cat("DEF\n")
  cat(str(DEF))
  cat("\n")
  # cat("check_\n")
  # cat(sprintf(as.character(check)))
  # cat("\n")
  # 
  # quit(status=1)
  
  breaks.Rank<-seq(0,1,by=0.1)
  labels.Rank<-as.character(breaks.Rank)
  
  
  cat("labels.Rank_\n")
  cat(sprintf(as.character(labels.Rank)))
  cat("\n")
  
  
  
  p <- ggplot(data=DEF,
              aes(x=HGNC,y=mean_Transcript_Ratio,
                  color=VEP_DEF_LABELS)) +
    geom_jitter(size = 3,width = 0.05)+
    scale_color_manual(values=c('#32A852','#6DB2EE','#62D07F','#1877C9','#C9244B','#D45E85','#87447B','#553B68','#D6B8E6',
                                "#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "gray", "#CBD588", "#5F7FC7", 
                                "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
                                "#D14285", "#6DDE88", "#652926", "#D1A33D", "#C84248","#89C5DA", "#DA5724", "#74D944", "#CE50CA"),drop=F)+
    scale_y_continuous(name="mean Transcript Ratio",breaks=breaks.Rank,labels=labels.Rank, 
                       limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
    scale_x_discrete(name=element_blank(), drop =FALSE)+
    theme_bw()+
    theme(axis.text.x=element_text(angle=45,size=14,colour="black",vjust=1,hjust=1,face="bold"), 
          axis.title.y=element_text(size=16,face="bold"),
          legend.title=element_text(size=16,face="bold"),
          legend.text=element_text(size=12,face="bold",colour="black"),
          axis.text.y=element_text(colour="black",size=12,face="bold"))+
    labs(color=paste("Transcript","consequences", sep="\n"))+
    theme(legend.key=element_blank(), legend.key.size=unit(1,"point"), 
          legend.position="right")+
    ggeasy::easy_center_title()
  
  setwd(path2)
  
  # geom_mark_ellipse(data=GENE_EXP[which(GENE_EXP$MOCK_COORD == "genIE host genes"),], 
  #                   aes(x=MOCK_COORD,y=log10(mean_GENE_EXP + 0.0001),label=HGNC))
  cat("Hello_world2\n")
  
  svglite(paste('Kolf2_Transcript_Ratio_DNMT1_S1PR2','.svg',sep=''), width = 8, height = 8)
  print(p)
  dev.off()
  
  # quit(status=1)
  
  ##### Transcript EXP -----
  
  
  
  Max_Transcript_EXP<-max(DEF$mean_Transcript_EXP[!is.na(DEF$mean_Transcript_EXP)])
  log_10_Max_Transcript_EXP<-log10(Max_Transcript_EXP+0.0001)
  Min_Transcript_EXP<-min(DEF$mean_Transcript_EXP[!is.na(DEF$mean_Transcript_EXP)])
  log_10_Min_Transcript_EXP<-log10(Min_Transcript_EXP+0.0001)
  
  
  cat("Max_Transcript_EXP_\n")
  cat(sprintf(as.character(Max_Transcript_EXP)))
  cat("\n")
  
  cat("log_10_Max_Transcript_EXP_\n")
  cat(sprintf(as.character(log_10_Max_Transcript_EXP)))
  cat("\n")
  
  cat("Min_Transcript_EXP_\n")
  cat(sprintf(as.character(Min_Transcript_EXP)))
  cat("\n")
  
  cat("log_10_Min_Transcript_EXP_\n")
  cat(sprintf(as.character(log_10_Min_Transcript_EXP)))
  cat("\n")
  
  # quit(status = 1)
  
  # breaks.Transcript_EXP<-round(seq(Min_Transcript_EXP,(Max_Transcript_EXP+1),by=1),2)
  breaks.Transcript_EXP.log<-unique(c(seq(log_10_Min_Transcript_EXP,log_10_Max_Transcript_EXP,by=1),log_10_Max_Transcript_EXP,log_10_Max_Transcript_EXP+0.2))
  breaks.Transcript_EXP<-10^(breaks.Transcript_EXP.log)
  labels.Transcript_EXP<-as.character(round(breaks.Transcript_EXP,3))
  
  
  cat("breaks.Transcript_EXP.log_\n")
  cat(sprintf(as.character(breaks.Transcript_EXP.log)))
  cat("\n")
  
  cat("breaks.Transcript_EXP_\n")
  cat(sprintf(as.character(breaks.Transcript_EXP)))
  cat("\n")
  cat("labels.Transcript_EXP_\n")
  cat(sprintf(as.character(labels.Transcript_EXP)))
  cat("\n")
  
  
  # quit(status = 1)
  
  
  
  p <- ggplot(data=DEF,
              aes(x=HGNC,y=log10(mean_Transcript_EXP+0.0001),
                  color=VEP_DEF_LABELS)) +
    geom_jitter(size = 3,width = 0.05)+
    scale_color_manual(values=c('#32A852','#6DB2EE','#62D07F','#1877C9','#C9244B','#D45E85','#87447B','#553B68','#D6B8E6',
                                "#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "gray", "#CBD588", "#5F7FC7", 
                                "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
                                "#D14285", "#6DDE88", "#652926", "#D1A33D", "#C84248","#89C5DA", "#DA5724", "#74D944", "#CE50CA"),drop=F)+
    scale_y_continuous(name="mean Transcript Expression (rpkm)",
                       breaks=breaks.Transcript_EXP.log,
                       labels=labels.Transcript_EXP, 
                       limits=c(breaks.Transcript_EXP.log[1],breaks.Transcript_EXP.log[length(breaks.Transcript_EXP.log)]))+
    scale_x_discrete(name=element_blank(), drop =FALSE)+
    theme_bw()+
    theme(axis.text.x=element_text(angle=45,size=14,colour="black",vjust=1,hjust=1,face="bold"), 
          axis.title.y=element_text(size=16,face="bold"),
          legend.title=element_text(size=16,face="bold"),
          legend.text=element_text(size=12,face="bold",colour="black"),
          axis.text.y=element_text(colour="black",size=12,face="bold"))+
    labs(color=paste("Transcript","consequences", sep="\n"))+
    theme(legend.key=element_blank(), legend.key.size=unit(1,"point"), 
          legend.position="right")+
    geom_hline(yintercept=log10(10+0.0001), linetype='dotted', col = 'red')+
    ggeasy::easy_center_title()
  
  setwd(path2)
  
  # geom_mark_ellipse(data=GENE_EXP[which(GENE_EXP$MOCK_COORD == "genIE host genes"),], 
  #                   aes(x=MOCK_COORD,y=log10(mean_GENE_EXP + 0.0001),label=HGNC))
  cat("Hello_world2\n")
  
  svglite(paste('Kolf2_Transcript_EXP_DNMT1_S1PR2','.svg',sep=''), width = 8, height = 8)
  print(p)
  dev.off()
  
  ##### ADD BY CSQ ------
  
  
  Matrix_file.dt<-data.table(DEF, key=c("ensembl_gene_id","HGNC","VEP_DEF_LABELS","VAR"))
  
  
  
  CSQ_SUM<-as.data.frame(Matrix_file.dt[,.(CSQ_SUM_Rep_1=sum(Rep_1),
                                           CSQ_SUM_Rep_2=sum(Rep_2))
                                        ,by=.(ensembl_gene_id,HGNC,VEP_DEF_LABELS,VAR)])
  
  cat("CSQ_SUM:0\n")
  cat(str(CSQ_SUM))
  cat("\n")
  
  mean.matrix.CSQ_SUM<-CSQ_SUM[,-c(1:4)]
  
  cat("--------------------------------->mean.matrix.CSQ_SUM\n")
  cat(str(mean.matrix.CSQ_SUM))
  cat("\n")
  
  mean.matrix.CSQ_SUM.mean<-apply(mean.matrix.CSQ_SUM,1,mean)
  
  cat("mean.matrix.CSQ_SUM.mean:0\n")
  cat(str(mean.matrix.CSQ_SUM.mean))
  cat("\n")
  
  CSQ_SUM$mean_CSQ_SUM<-mean.matrix.CSQ_SUM.mean
  
  
  cat("CSQ_SUM:2\n")
  cat(str(CSQ_SUM))
  cat("\n")
  
  #### graph
  
  Max_CSQ_SUM<-max(CSQ_SUM$mean_CSQ_SUM[!is.na(CSQ_SUM$mean_CSQ_SUM)])
  log_10_Max_CSQ_SUM<-log10(Max_CSQ_SUM+0.0001)
  Min_CSQ_SUM<-min(CSQ_SUM$mean_CSQ_SUM[!is.na(CSQ_SUM$mean_CSQ_SUM)])
  log_10_Min_CSQ_SUM<-log10(Min_CSQ_SUM+0.0001)
  
  
  cat("Max_CSQ_SUM_\n")
  cat(sprintf(as.character(Max_CSQ_SUM)))
  cat("\n")
  
  cat("log_10_Max_CSQ_SUM_\n")
  cat(sprintf(as.character(log_10_Max_CSQ_SUM)))
  cat("\n")
  
  cat("Min_CSQ_SUM_\n")
  cat(sprintf(as.character(Min_CSQ_SUM)))
  cat("\n")
  
  cat("log_10_Min_CSQ_SUM_\n")
  cat(sprintf(as.character(log_10_Min_CSQ_SUM)))
  cat("\n")
  
  # quit(status = 1)
  
  # breaks.CSQ_SUM<-round(seq(Min_CSQ_SUM,(Max_CSQ_SUM+1),by=1),2)
  breaks.CSQ_SUM.log<-unique(c(seq(log_10_Min_CSQ_SUM,log_10_Max_CSQ_SUM,by=1),log_10_Max_CSQ_SUM,log_10_Max_CSQ_SUM+0.2))
  breaks.CSQ_SUM<-10^(breaks.CSQ_SUM.log)
  labels.CSQ_SUM<-as.character(round(breaks.CSQ_SUM,3))
  
  
  cat("breaks.CSQ_SUM.log_\n")
  cat(sprintf(as.character(breaks.CSQ_SUM.log)))
  cat("\n")
  
  cat("breaks.CSQ_SUM_\n")
  cat(sprintf(as.character(breaks.CSQ_SUM)))
  cat("\n")
  cat("labels.CSQ_SUM_\n")
  cat(sprintf(as.character(labels.CSQ_SUM)))
  cat("\n")
  
  
  # quit(status = 1)
  
  
  
  p <- ggplot(data=CSQ_SUM,
              aes(x=HGNC,y=log10(mean_CSQ_SUM+0.0001),
                  color=VEP_DEF_LABELS)) +
    geom_jitter(size = 3,width = 0.05)+
    scale_color_manual(values=c('#32A852','#6DB2EE','#62D07F','#1877C9','#C9244B','#D45E85','#87447B','#553B68','#D6B8E6',
                                "#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "gray", "#CBD588", "#5F7FC7", 
                                "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
                                "#D14285", "#6DDE88", "#652926", "#D1A33D", "#C84248","#89C5DA", "#DA5724", "#74D944", "#CE50CA"),drop=F)+
    scale_y_continuous(name="Aggregate transcript expression per consequence (rpkm)",
                       breaks=breaks.CSQ_SUM.log,
                       labels=labels.CSQ_SUM, 
                       limits=c(breaks.CSQ_SUM.log[1],breaks.CSQ_SUM.log[length(breaks.CSQ_SUM.log)]))+
    scale_x_discrete(name=element_blank(), drop =FALSE)+
    theme_bw()+
    theme(axis.text.x=element_text(angle=45,size=14,colour="black",vjust=1,hjust=1,face="bold"), 
          axis.title.y=element_text(size=16,face="bold"),
          legend.title=element_text(size=16,face="bold"),
          legend.text=element_text(size=12,face="bold",colour="black"),
          axis.text.y=element_text(colour="black",size=12,face="bold"))+
    labs(color=paste("Transcript","consequences", sep="\n"))+
    theme(legend.key=element_blank(), legend.key.size=unit(1,"point"), 
          legend.position="right")+
    geom_hline(yintercept=log10(10+0.0001), linetype='dotted', col = 'red')+
    ggeasy::easy_center_title()
  
  setwd(path2)
  
  # geom_mark_ellipse(data=GENE_EXP[which(GENE_EXP$MOCK_COORD == "genIE host genes"),], 
  #                   aes(x=MOCK_COORD,y=log10(mean_GENE_EXP + 0.0001),label=HGNC))
  cat("Hello_world2\n")
  
  svglite(paste('Kolf2_CSQ_SUM_DNMT1_S1PR2','.svg',sep=''), width = 8, height = 8)
  print(p)
  dev.off()
  
  
  
  # quit(status = 1)
  ##### SAVE -----
  
  write.table(DEF,file="Expression_Kolf2_genIE_DNMT1_S1PR2.tsv", sep="\t", row.names = F, quote = F)
  
  write.table(CSQ_SUM,file="Expression_Kolf2_genIE_aggregate_by_csq_DNMT1_S1PR2.tsv", sep="\t", row.names = F, quote = F)
  
  
  
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
    make_option(c("--STAR_kolf2_1"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--STAR_kolf2_2"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Equiv.table"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--gene_lengths"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--VEP_route"), type="numeric", default=NULL, 
                metavar="filename", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--genIE_threshold_replicate"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--SELECTED_GENES"), type="character", default=NULL, 
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
                metavar="filename", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
  )
  parser = OptionParser(usage = "140__Rscript_v106.R
                        --subset type
                        --TranscriptEXP FILE.txt
                        --cadd FILE.txt
                        --ncboost FILE.txt
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  
 # Normalize_STAR_alignment_Kolf2(opt)
  #GENE_EXP_and_TRatios(opt)
  ploting_GENE_EXP(opt)
  # ploting_Transcript_Ratio(opt)
  # ploting_Transcript_Ratio_DNMT1(opt)
  
}


###########################################################################

system.time( main() )

