
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

opt = NULL


PREFIXES_Table_EXPANDED_generator = function(option_list)
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
  
  #### Read master file ----
  
  master_file<-as.data.frame(fread(opt$master_file, header=T, sep="\t"), stringsAsFactors = F)
  
    
  cat("master_file_\n")
  str(master_file)
  cat("\n")
  
  indx.int<-c(which(colnames(master_file) == "HGNC"),which(colnames(master_file) == "VAR"),which(colnames(master_file) == "id"))
  
  master_file_subset<-unique(master_file[,indx.int])
  
  cat("master_file_subset_\n")
  str(master_file_subset)
  cat("\n")
  
  #### Read manifest ----
  
  manifest<-read.csv(opt$manifest, header=T, stringsAsFactors = F, skip=14)
  
  manifest$sample<-gsub("_.+$","",manifest$Sample_Name)
  colnames(manifest)[which(colnames(manifest) == "index")]<-"index.P7"
  colnames(manifest)[which(colnames(manifest) == "index2")]<-"index.P5"
  manifest$type<-gsub("^[^_]+_","",manifest$Sample_Name)
  manifest$number<-gsub("^[^0-9+]+","",manifest$type)
  manifest$type<-gsub("[0-9]+","",manifest$type)
  
  
  
  
  manifest$type[manifest$type == "G"]<-"gDNA"
  manifest$type[manifest$type == "C"]<-"cDNA"
  manifest$Replicate<-paste(manifest$type,manifest$number,sep="_")
  
  
  cat("manifest_\n")
  str(manifest)
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(manifest$type))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(manifest$type)))))
  cat("\n")
  
  
  setwd(out)
  
  write.table(manifest,file="test_manifest.tsv",sep="\t",quote=F,row.names=F)
  #### Read equivalence ----
  
  equivalence<-read.csv(opt$equivalence, header=F, stringsAsFactors = F)
  
  cat("equivalence_\n")
  str(equivalence)
  cat("\n")
  
  
  check_A1_1<-equivalence[1,1]
  
  cat("check_A1_1\n")
  cat(sprintf(as.character(check_A1_1)))
  cat("\n")
  
  if(check_A1_1 == "GenIE ID")
  {
    
    equivalence<-equivalence[-1,]
    
    colnames(equivalence)[which(colnames(equivalence) == "V1")]<-"sample"
    colnames(equivalence)[which(colnames(equivalence) == "V2")]<-"chr"
    colnames(equivalence)[which(colnames(equivalence) == "V3")]<-"VAR"
    colnames(equivalence)[which(colnames(equivalence) == "V13")]<-"HGNC"
    
    
    
    cat("equivalence_\n")
    str(equivalence)
    cat("\n")
    
    
    equivalence_subset<-unique(equivalence[which(equivalence$HGNC%in%master_file$HGNC),c(which(colnames(equivalence) == "sample"),which(colnames(equivalence) == "HGNC"))])
    
    cat("equivalence_subset_\n")
    str(equivalence_subset)
    cat("\n")
    
    
    # quit(status=1)
    
    
  }else{
    
    if(check_A1_1 == "GenERA round")
    {
      
      equivalence<-equivalence[-1,]
      
      colnames(equivalence)[which(colnames(equivalence) == "V1")]<-"sample"
      colnames(equivalence)[which(colnames(equivalence) == "V3")]<-"chr"
      colnames(equivalence)[which(colnames(equivalence) == "V4")]<-"VAR"
      colnames(equivalence)[which(colnames(equivalence) == "V14")]<-"HGNC"
      
      
      
      cat("equivalence_\n")
      str(equivalence)
      cat("\n")
      
      
      equivalence_subset<-unique(equivalence[which(equivalence$HGNC%in%master_file$HGNC),c(which(colnames(equivalence) == "sample"),which(colnames(equivalence) == "HGNC"))])
      
      cat("equivalence_subset_\n")
      str(equivalence_subset)
      cat("\n")
      
      
      # quit(status=1)
      
      
    }else{
      
      colnames(equivalence)[which(colnames(equivalence) == "V1")]<-"sample"
      colnames(equivalence)[which(colnames(equivalence) == "V2")]<-"HGNC"
      
      
      
      equivalence_subset<-unique(equivalence[which(equivalence$HGNC%in%master_file$HGNC),c(which(colnames(equivalence) == "sample"),which(colnames(equivalence) == "HGNC"))])
      
      cat("equivalence_subset_\n")
      str(equivalence_subset)
      cat("\n")
    }
   
    
  }
    

  
  # quit(status = 1)
  
  
  
  setwd(out)
  
  write.table(equivalence_subset,file="test_equivalence.tsv",sep="\t",quote=F,row.names=F)
  
  #### merge manifest and equivalence ----
  
  manifest2<-merge(manifest,
                  equivalence_subset,
                  by="sample")
  
  cat("manifest2_\n")
  str(manifest2)
  cat("\n")
  
 
  if(dim(manifest2)[1] ==0)
  {
    manifest$sample<-gsub("^[^_]+_[^_]+_","",manifest$Sample_Name)
    
    
    
    manifest$type<-gsub("^[^_]+_","",manifest$sample)
    manifest$number<-gsub("^[^0-9+]+","",manifest$type)
    manifest$type<-gsub("[0-9]+","",manifest$type)
    
    gDNA_vector<-c("g","G")
    cDNA_vector<-c("c","C")
    
    
    manifest$type[manifest$type%in%gDNA_vector]<-"gDNA"
    manifest$type[manifest$type%in%cDNA_vector]<-"cDNA"
    manifest$Replicate<-paste(manifest$type,manifest$number,sep="_")
    
    manifest$sample<-gsub("_.+$","",manifest$sample)
    
    
    cat("manifest_NEW\n")
    str(manifest)
    cat("\n")
    
    manifest<-merge(manifest,
                     equivalence_subset,
                     by="sample")
    
    cat("manifest_\n")
    str(manifest)
    cat("\n")
    
    
  }else{
    
    manifest<-manifest2
  }
  
 # quit(status = 1)
  
  index.interest<-c(which(colnames(manifest) == "sample"),
                    which(colnames(manifest) == "HGNC"),
                    which(colnames(manifest) == "Replicate"),
                    which(colnames(manifest) == "Sample_Name"),
                    which(colnames(manifest) == "index.P5"),
                    which(colnames(manifest) == "index.P7"))
  
  manifest<-unique(as.data.frame(manifest)[,index.interest])
  
  cat("manifest_\n")
  str(manifest)
  cat("\n")
  
  manifest<-merge(manifest,
                   master_file_subset,
                   by="HGNC")
  
  cat("manifest_\n")
  str(manifest)
  cat("\n")
  
  
  cat("HGNC\n")
  cat(sprintf(as.character(levels(as.factor(manifest$HGNC)))))
  cat("\n")
  cat("Replicate\n")
  cat(sprintf(as.character(levels(as.factor(manifest$Replicate)))))
  cat("\n")
  
  manifest$Sample_Name<-gsub("^Genie2_|^Genie1_","",manifest$Sample_Name)
 
  # quit(status = 1)
  

  #### SAVE  PREFIXES_Table_EXPANDED.txt ----
  
  filename<-paste(type,"_PREFIXES_Table_EXPANDED.txt",sep="")
  
  write.table(manifest, file=filename,
              sep="\t", quote=F,
              row.names = F)
}

Table_of_files_generator = function(option_list)
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
  
  #### READ and transform FILES ----
  
  FILES<-read.table(opt$FILES, sep="\t", header=T, stringsAsFactors = F)
  
  cat("FILES_0\n")
  str(FILES)
  cat("\n")
  
  FILES$Prefix<-gsub("_.+$","",FILES$File)
  FILES$Prefix<-gsub("^Genie2[_-]|^Genie1[_-]","",FILES$Prefix)
  
  cat("FILES_\n")
  str(FILES)
  cat("\n")
  
  
  # quit(status=1)
  
  #### READ and transform master_path ----
  
  master_path = opt$master_path
  
  cat("master_path_\n")
  cat(sprintf(as.character(master_path)))
  cat("\n")
  
  #### Read Prefix table ----
  
 # filename<-paste(type,"_PREFIXES_Table_EXPANDED.txt",sep="")
  
  Prefix.Table<-read.table(file=opt$PREFIXES_Table_EXPANDED, sep="\t", header=T, stringsAsFactors = F)
  
  cat("Prefix.Table_\n")
  str(Prefix.Table)
  cat("\n")
  
  
  
  #### LOOP TO ASSIGN FASTQ FILES ----
  
 
  MASTER_ROUTE<-as.character(master_path)
  
  SAMPLES<-unique(Prefix.Table$sample)
  
  list_gather_DEF<-list()
  
  for(i in 1:length(SAMPLES))
  {
    
    sample<-SAMPLES[i]
    
    cat("sample_\n")
    cat(sprintf(as.character(sample)))
    cat("\n")
    
    Prefix.Table_sel<-Prefix.Table[which(Prefix.Table$sample == sample),]
    
    cat("Prefix.Table_sel_\n")
    str(Prefix.Table_sel)
    cat("\n")
    
  # quit(status = 1)
    list_gather<-list()
    
    for(k in 1:length(Prefix.Table_sel$Sample_Name))
    {
      PREFIX<-Prefix.Table_sel$Sample_Name[k]
      
      
      PREFIX<-gsub("_","-",PREFIX)
      
      cat("PREFIX_\n")
      cat(sprintf(as.character(PREFIX)))
      cat("\n")
      
      Prefix.Table_sel_PREFIX<-Prefix.Table_sel[k,]
      
      cat("Prefix.Table_sel_PREFIX_\n")
      str(Prefix.Table_sel_PREFIX)
      cat("\n")
      
      ### It should start with the PREFIX
      
      FILES_sel<-FILES$File[which(FILES$Prefix == PREFIX)]
      
      cat("FILES_sel_\n")
      str(FILES_sel)
      cat("\n")
      
      pivotal_point<-length(FILES_sel)
      
      cat("pivotal_point_\n")
      cat(sprintf(as.character(pivotal_point)))
      cat("\n")
      
      # quit(status = 1)
      
      if(pivotal_point == 2)
      {
        R1_file<-FILES_sel[grep("R1_[0-9]+\\.fastq.gz", FILES_sel,perl=T)]
        R2_file<-FILES_sel[grep("R2_[0-9]+\\.fastq.gz", FILES_sel,perl=T)]
        
        R1_file<-paste(MASTER_ROUTE,R1_file,sep='')
        R2_file<-paste(MASTER_ROUTE,R2_file,sep='')
        
        Partial<-as.data.frame(cbind(Prefix.Table_sel_PREFIX,R1_file,R2_file))
        
        # cat("Partial_1\n")
        # str(Partial)
        # cat("\n")
        
        list_gather[[k]]<-Partial
        # 
        # quit(status = 1)
        
       
        
        
      }else{
        
        if(pivotal_point == 1)
        {
          R1_file<-FILES_sel[grep("R1_[0-9]+\\.fastq.gz", FILES_sel,perl=T)]
          R2_file<-FILES_sel[grep("R2_[0-9]+\\.fastq.gz", FILES_sel,perl=T)]
          
          if(length(R1_file) >0)
          {
            R1_file<-paste(MASTER_ROUTE,R1_file,sep='')
            R2_file<-"NA"
            
          }else{
            
            R2_file<-paste(MASTER_ROUTE,R2_file,sep='')
            R1_file<-"NA"
          }
          
          Partial<-as.data.frame(cbind(Prefix.Table_sel_PREFIX,R1_file,R2_file))
          
          # cat("Partial_2\n")
          # str(Partial)
          # cat("\n")
          
          list_gather[[k]]<-Partial

          # quit(status = 1)
          
         
          
        }else{
          
          print("WARNING_FILES_PREFIX")
          # quit(status = 1)
          
          PREFIX<-paste(gsub("-","_",PREFIX),'_',sep='')
          
          cat("PREFIX_\n")
          cat(sprintf(as.character(PREFIX)))
          cat("\n")
          # 
          # FILES$Prefix<-gsub("_[^_]+$","",FILES$File)
          # 
          # cat("FILES_\n")
          # str(FILES)
          # cat("\n")
          
          ### It should start with the PREFIX
          
          FILES_sel<-FILES$File[grep(paste('^',PREFIX,sep=''),FILES$File)]
          
          cat("FILES_sel_\n")
          str(FILES_sel)
          cat("\n")
          
          pivotal_point<-length(FILES_sel)
          
          cat("pivotal_point_\n")
          cat(sprintf(as.character(pivotal_point)))
          cat("\n")
          
          if(pivotal_point == 2)
          {
            R1_file<-FILES_sel[grep("R1_[0-9]+\\.fastq.gz", FILES_sel,perl=T)]
            R2_file<-FILES_sel[grep("R2_[0-9]+\\.fastq.gz", FILES_sel,perl=T)]
            
            R1_file<-paste(MASTER_ROUTE,R1_file,sep='')
            R2_file<-paste(MASTER_ROUTE,R2_file,sep='')
            
            Partial<-as.data.frame(cbind(Prefix.Table_sel_PREFIX,R1_file,R2_file))
            
            # cat("Partial_1\n")
            # str(Partial)
            # cat("\n")
            
            list_gather[[k]]<-Partial
            # 
            # quit(status = 1)
            
            
            
            
          }else{
            
            if(pivotal_point == 1)
            {
              R1_file<-FILES_sel[grep("R1_[0-9]+\\.fastq.gz", FILES_sel,perl=T)]
              R2_file<-FILES_sel[grep("R2_[0-9]+\\.fastq.gz", FILES_sel,perl=T)]
              
              if(length(R1_file) >0)
              {
                R1_file<-paste(MASTER_ROUTE,R1_file,sep='')
                R2_file<-"NA"
                
              }else{
                
                R2_file<-paste(MASTER_ROUTE,R2_file,sep='')
                R1_file<-"NA"
              }
              
              Partial<-as.data.frame(cbind(Prefix.Table_sel_PREFIX,R1_file,R2_file))
              
              # cat("Partial_2\n")
              # str(Partial)
              # cat("\n")
              
              list_gather[[k]]<-Partial
              
              # quit(status = 1)
              
              
              
            }else{
              
              
              
          
              PREFIX<-gsub("^[^_]+_","",PREFIX)
              
              cat("PREFIX_3\n")
              cat(sprintf(as.character(PREFIX)))
              cat("\n")
              # 
              # FILES$Prefix<-gsub("_[^_]+$","",FILES$File)
              # 
              # cat("FILES_\n")
              # str(FILES)
              # cat("\n")
              
              ### It should start with the PREFIX
              
              FILES_sel<-FILES$File[grep(paste('^',PREFIX,sep=''),FILES$File)]
              
              cat("FILES_sel_\n")
              str(FILES_sel)
              cat("\n")
              
              pivotal_point<-length(FILES_sel)
              
              cat("pivotal_point_\n")
              cat(sprintf(as.character(pivotal_point)))
              cat("\n")
              
              if(pivotal_point == 2)
              {
                R1_file<-FILES_sel[grep("R1_[0-9]+\\.fastq.gz", FILES_sel,perl=T)]
                R2_file<-FILES_sel[grep("R2_[0-9]+\\.fastq.gz", FILES_sel,perl=T)]
                
                R1_file<-paste(MASTER_ROUTE,R1_file,sep='')
                R2_file<-paste(MASTER_ROUTE,R2_file,sep='')
                
                Partial<-as.data.frame(cbind(Prefix.Table_sel_PREFIX,R1_file,R2_file))
                
                # cat("Partial_1\n")
                # str(Partial)
                # cat("\n")
                
                list_gather[[k]]<-Partial
                # 
                # quit(status = 1)
                
                
                
                
              }else{
                
                if(pivotal_point == 1)
                {
                  R1_file<-FILES_sel[grep("R1_[0-9]+\\.fastq.gz", FILES_sel,perl=T)]
                  R2_file<-FILES_sel[grep("R2_[0-9]+\\.fastq.gz", FILES_sel,perl=T)]
                  
                  if(length(R1_file) >0)
                  {
                    R1_file<-paste(MASTER_ROUTE,R1_file,sep='')
                    R2_file<-"NA"
                    
                  }else{
                    
                    R2_file<-paste(MASTER_ROUTE,R2_file,sep='')
                    R1_file<-"NA"
                  }
                  
                  Partial<-as.data.frame(cbind(Prefix.Table_sel_PREFIX,R1_file,R2_file))
                  
                  # cat("Partial_2\n")
                  # str(Partial)
                  # cat("\n")
                  
                  list_gather[[k]]<-Partial
                  
                  # quit(status = 1)
                  
                  
                  
                }else{
                  
                  print("UNABLE_TO_RESCUE_PREFIX\n")
                  print("FILE might not be part of the table\n")
                  
                  quit(status = 1)
                }
                
             
                
              }
                
          
              
              
            }
          
          
          
       
          }# else pivotal_point == 2
          
        }# else pivotal_point == 1 ORIGINAL
      }# else pivotal_point == 2 ORIGINAL
    } # k
    
    TABLE_files = unique(as.data.frame(data.table::rbindlist(list_gather, fill = T)))
    
    cat("TABLE_files\n")
    str(TABLE_files)
    cat("\n")
    
    list_gather_DEF[[i]]<-TABLE_files
    
    # if(sample == "H3")
    # {
    #   quit(status = 1)
    #   
    # }
    # 
    
  } #i
  
  
  TABLE_files = unique(as.data.frame(data.table::rbindlist(list_gather_DEF, fill = T)))
  
  cat("TABLE_files\n")
  str(TABLE_files)
  cat("\n")
  
  
  
  # quit(status = 1)
  
  
  DEF<-merge(Prefix.Table, TABLE_files,
             by=colnames(Prefix.Table),
             all=T)
  
  cat("DEF\n")
  str(DEF)
  cat("\n")
  
  
  DEF_NO_NA<-DEF[!is.na(DEF$R1_file),]
  
  cat("DEF_NO_NA\n")
  str(DEF_NO_NA)
  cat("\n")
  
  
  
  
  #### SAVE ----
  
  setwd
  
  filename<-paste(type,"_Table_of_files.txt",sep="")
  
  write.table(DEF_NO_NA, file=filename,
              sep="\t", quote=F,
              row.names = F)
  
  # quit(status = 1)
  
  cat("THE END 2\n")
  
  # quit(status=1)
}

genie.meta_and_flash_input_generator = function(option_list)
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
  
  #### Read Prefix table ----
  
  filename<-paste(type,"_PREFIXES_Table_EXPANDED.txt",sep="")
  
  Prefix.Table<-read.table(file=filename, sep="\t", header=T, stringsAsFactors = F)
  
  cat("Prefix.Table_\n")
  str(Prefix.Table)
  cat("\n")
  
  #### Read input_regions ----
  
  #filename<-paste(type,"_input_regions.tsv",sep="")
  
  input_regions<-read.table(file=opt$input_regions, sep="\t", header=T, stringsAsFactors = F)
  
  colnames(input_regions)[which(colnames(input_regions)  == "HGNC")]<-"hgnc"
  
  
  cat("input_regions_\n")
  str(input_regions)
  cat("\n")
  
  #### Read Table_of_files ----
  
  filename<-paste(type,"_Table_of_files.txt",sep="")
  
  Table_of_files<-read.table(file=filename, sep="\t", header=T, stringsAsFactors = F)
  
  cat("Table_of_files_\n")
  str(Table_of_files)
  cat("\n")
  
  colnames(Table_of_files)[which(colnames(Table_of_files)  == "HGNC")]<-"hgnc"
  
  # quit(status=1)
  
  #### HK Table_of_files ----
  
  Table_of_files$name<-Table_of_files$hgnc
  
  #Table_of_files$name[which(Table_of_files$name == "EROS")]<-"eros"
  #Table_of_files_subset<-Table_of_files[-which(Table_of_files$name == "SAV1"),]
  
  Table_of_files$replicate<-gsub("DNA","",Table_of_files$Replicate)
  Table_of_files$replicate<-gsub("_","",Table_of_files$replicate)
  Table_of_files$replicate_full<-paste(Table_of_files$name,
                                              Table_of_files$Sample_Name,
                                              sep="_")
  
  Table_of_files$type<-gsub("_[0-9]+","",Table_of_files$Replicate)
  
  colnames(Table_of_files)[which(colnames(Table_of_files) == "Sample_Name")]<-"prefix"
  colnames(Table_of_files)[which(colnames(Table_of_files) == "Sample_ID")]<-"plate"
  colnames(Table_of_files)[which(colnames(Table_of_files) == "index")]<-"Barcode_FWD"
  colnames(Table_of_files)[which(colnames(Table_of_files) == "index2")]<-"Barcode_REV"
  colnames(Table_of_files)[which(colnames(Table_of_files) == "R1_file")]<-"Read_1_file"
  colnames(Table_of_files)[which(colnames(Table_of_files) == "R2_file")]<-"Read_2_file"
  colnames(Table_of_files)[which(colnames(Table_of_files) == "id")]<-"rsid"
  
  
  Table_of_files$PCR<-"NA"
  Table_of_files$bam<-"NA"
  
  # Table_of_files$prefix<-gsub("_","-",Table_of_files$prefix)
  
  cat("Table_of_files_\n")
  str(Table_of_files)
  cat("\n")
  
  # quit(status=1)
  
  Header_Jeremy_annotation_file<-c("name","replicate","replicate_full","rsid","hgnc","sample","type","PCR","prefix","plate","Barcode_FWD","Barcode_REV","Read_1_file","Read_2_file","bam","gDNA_1","gDNA_2","gDNA_3","cDNA_1","cDNA_2","cDNA_3","cDNA_4","cDNA_5","cDNA_6")
  
  Header_Jeremy_annotation_file_subset<-Header_Jeremy_annotation_file[which(Header_Jeremy_annotation_file%in%colnames(Table_of_files))]
  
  cat("Header_Jeremy_annotation_file_subset\n")
  cat(sprintf(as.character(Header_Jeremy_annotation_file_subset)))
  cat("\n")
  
  #### regionsFile & fastqMetaFile ----
  
  regionsFile = input_regions#[-which(input_regions$name == "SAV1"),]
  
  fastqMetaFile = Table_of_files[,which(colnames(Table_of_files)%in%Header_Jeremy_annotation_file_subset)]
  
  cat("fastqMetaFile_1\n")
  str(fastqMetaFile)
  cat("\n")
  
  fastqMetaFile<-fastqMetaFile[,Header_Jeremy_annotation_file_subset]
  
  cat("fastqMetaFile_2\n")
  str(fastqMetaFile)
  cat("\n")
  
  #### Merge & Parameter calculation ----
  
  readLenParam = 150
  sdParam = NULL
  
  if (is.null(readLenParam) || is.na(readLenParam)) {
    readLenParam = 150
  } else {
    readLenParam = strtoi(readLenParam)
  }
  if (is.null(sdParam) || is.na(sdParam)) {
    sdParam = 20
  } else {
    sdParam = strtoi(sdParam)
  }
  
  regions.df = regionsFile
  
  meta.df = fastqMetaFile
  
  df = meta.df %>% dplyr::left_join(regions.df, by="name") %>%
    dplyr::mutate(amplicon_size = end - start + 1) %>%
    group_by(replicate_full) %>%
    dplyr::mutate(read_len = readLenParam, amplicon_sd = sdParam,
                  expected_overlap = max(5, 2*read_len - amplicon_size)) %>%
    dplyr::mutate(output_file = paste0(replicate_full, ".extendedFrags.fastq.gz"),
                  min_overlap = getMinOverlap(amplicon_size, read_len),
                  max_mismatch_dens = min(0.15, 0.1 + 1 / expected_overlap)) %>%
    ungroup()
  
  df$max_mismatch_dens = sprintf("%.3g", df$max_mismatch_dens)
  
  df = df %>% dplyr::select(replicate_full, Read_1_file, Read_2_file, read_len, amplicon_size, amplicon_sd, min_overlap, max_mismatch_dens, output_file)
  
  
  #### SAVE ----
  
  filename<-paste(type,"_genie.meta.tsv",sep="")
  
  write.table(meta.df, file=filename, sep="\t", quote=F, 
              row.names=F, col.names=T, eol="\n")
  
  filename<-paste(type,"_flash_input.tsv",sep="")
  
  write.table(df, file=filename, sep="\t", quote=F, 
              row.names=F, col.names=T, eol="\n")
  
  
}


################################  REBOOT CREATE REPLICATES FILE ---------------


### 244 FLASH YES/NO + ALIGNMENT
### 250 All the scripts and the genIE functions (Jeremy's)

printList = function(l, prefix = "    ") {
  list.df = data.frame(val_name = names(l), value = as.character(l))
  list_strs = apply(list.df, MARGIN = 1, FUN = function(x) { paste(x, collapse = " = ")})
  cat(paste(paste(paste0(prefix, list_strs), collapse = "\n"), "\n"))
}

getMinOverlap = function(amplicon_size, read_len) {
  minOverlap = 10
  if (2*read_len - amplicon_size < 10) {
    minOverlap = 2*read_len - amplicon_size - 1
    if (minOverlap < 4) {
      minOverlap = 4
    }
  }
  minOverlap
}


#### main script ----

main = function() {
  cmd_line = commandArgs()
  cat("Command line:\n")
  cat(paste(gsub("--file=", "", cmd_line[4], fixed=T),
            paste(cmd_line[6:length(cmd_line)], collapse = " "),
            "\n\n"))
  option_list <- list(
    make_option(c("--master_file"), type="character", default=NULL, 
                metavar="FILE.txt", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--manifest"), type="character", default=NULL, 
                metavar="FILE.txt", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--equivalence"), type="character", default=NULL, 
                metavar="FILE.txt", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--fasta"), type="character", default=NULL, 
                metavar="FILE.txt", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--PREFIXES_Table_EXPANDED"), type="character", default=NULL, 
                metavar="FILE.txt", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--FILES"), type="character", default=NULL, 
                metavar="FILE.txt", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--input_regions"), type="character", default=NULL, 
                metavar="FILE.txt", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--master_path"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--parameters"), type="character", default=NULL, 
                metavar="FILE.txt", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="filename", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
    
  )
  parser = OptionParser(usage = "136_Rscript_edition_Add_Expression_List.R
                        --master_file FILE.txt
                        --equivalence_genIE2 FILE.txt
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
 
   PREFIXES_Table_EXPANDED_generator(opt)
   Table_of_files_generator(opt)
  genie.meta_and_flash_input_generator(opt)
  
}


###########################################################################

system.time( main() )
