#Time:2022.10.29
#Author:Dai Yulong
#组蛋白输入文件格式转换

ARGV<-commandArgs(TRUE)
if (length(ARGV)<2) {
  cat("Rscript preAnalysisHistone.R  input outpath\n")
  q()
}

input<-as.character(ARGV[1])
pn<-as.character(ARGV[2])

if(FALSE){
  input<-"E:/projectResearch/07histonePipline/inputfile.txt"
  pn<-"E:/projectResearch/07histonePipline"
}

suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(stringr))
suppressMessages(library(openxlsx))

if (!dir.exists(pn)) {
  dir.create(pn)
}

if(str_ends(input,'csv')){
  dat <- read.csv(input,check.names = FALSE)
}else if (str_ends(input,'xlsx')) {
  dat <- read.xlsx(input,1)
}else{
  dat <- read.delim(input,check.names = FALSE)
}

index <- which(grepl(")",dat[[1]]))

all_dat <- data.frame()
for (i in 1:length(index)) {
  
  if (i+1 <= length(index)) {
   
    block <- dat[index[i]:(index[i+1]-1),]
    block <- block %>% add_column(class=block[1,1],.after = 1)
  
    all_dat <- rbind(all_dat,block)
  }else if(i+1 > length(index)){
    block <- dat[index[i]:nrow(dat),]
    block <- block %>% add_column(class=block[1,1],.after = 1)

    all_dat <- rbind(all_dat,block)
  }
  all_dat <- all_dat[complete.cases(all_dat),]
}

all_dat <- separate(all_dat,class,c("peptide","type"),"\\(",convert = T)
all_dat <- separate(all_dat,type,c("type",NA,NA),"\\)",convert = T)


write.csv(all_dat,paste0(pn, "/infile.csv"),row.names = F)
