#!/usr/bin/env Rscript
#Time:2022.11.04
#Author: daiyulong
#function: 读取参数文件

ARGV<-commandArgs(TRUE)
if (length(ARGV)<3) {
  cat("Rscript readParaHistone.R inputPara input outpath\n")
  q()
}

inputPara<-as.character(ARGV[1])
input<-as.character(ARGV[2])
pn<-as.character(ARGV[3])

if(FALSE){
  inputPara<-'E:/projectResearch/07histonePipline/parameter.xlsx'
  input<-'E:/projectResearch/07histonePipline/infile.csv'
  pn<-'E:/projectResearch/07histonePipline'
}

suppressMessages(library(openxlsx))
suppressMessages(library(readxl))
suppressMessages(library(ggsci))



if (!dir.exists(pn)) {
  dir.create(pn)
}

mycols<- c(pal_lancet("lanonc")(9)[1:8], pal_jama("default")(7), pal_nejm("default")(8), pal_jco("default")(10))
shape.vec<- rep(21:25, 10)

# sg<-read_excel(inputPara,sheet=1,na="NA")
sg<-read.xlsx(inputPara,sheet=1,na.strings = "NA")
sg<-as.data.frame(sg)
colnames(sg)[1:2]<-c("Sample","Group")
rownames(sg)<-sg[[1]]
sg$Shape<-shape.vec[as.integer(as.factor(sg$Group))]
sg$Color<-mycols[as.integer(as.factor(sg$Group))]
write.table(sg,paste0(pn,"/groupFile.txt"),row.names = FALSE,quote=FALSE,sep="\t")

# comp<-read_excel(inputPara,sheet=2,na ="NA")
comp<-read.xlsx(inputPara,sheet=2,na.strings ="NA")
write.table(comp,paste0(pn,"/comparison.txt"),row.names= FALSE, quote=FALSE,sep="\t")

if(stringr::str_ends(input,'xlsx')){
  dd<-read.xlsx(input,sheet=1)
}else if (stringr::str_ends(input,'csv')) {
  dd <- read.csv(input,check.names = FALSE)
}else{
  dd<-read.delim(input,check.names = FALSE)
}


for(i in 1:NROW(comp)){
  #i=1
  method = tolower(as.character(comp[i,1]))
  if(method=='ttest' | method=='pairedttest'){
    g1=as.character(comp[i,2])
    g2=as.character(comp[i,3])
    stringr::str_trim(g1)
    stringr::str_trim(g2)
    sample<-c(sg$Sample[sg$Group==g1],sg$Sample[sg$Group==g2])
    #创建目录
    if(!dir.exists(paste0(pn,"/",g1,".vs.",g2))){
      dir.create(paste0(pn,"/",g1,".vs.",g2), recursive = TRUE)
    }
    #生成input文件
    compFile<-dd[,c(colnames(dd)[1],sample)] ## modified 2022-09-19
    write.table(compFile,paste0(pn,"/",g1,".vs.",g2,"/input.txt"), row.names = FALSE, quote=FALSE, sep="\t")
    comgsg<-sg[sample,]
    comgsg$Type<-1
    comgsg$Type[comgsg$Group==g2]<-0
    comgsg<-comgsg[,c(1,2,5,3,4)]
    write.table(comgsg, paste0(pn,'/',g1,".vs.",g2,"/group.txt"),row.names=FALSE,quote=FALSE, sep="\t")
  }else if(method=='anova'){
    anovaGroup=as.character(comp[i,2])
    stringr::str_trim(anovaGroup)
    groups<-stringr::str_split(anovaGroup,';')[[1]]
    sample<-sg$Sample[sg$Group %in% groups]
    output <- stringr::str_replace_all(anovaGroup,';','_')
    compFile <- dd[, c(colnames(dd)[1],sample)]
    if(!dir.exists(paste0(pn,"/Anova/",output))){
      dir.create(paste0(pn,"/Anova/",output), recursive =TRUE)
    }
    write.table(compFile, paste0(pn,"/Anova/", output, '/input.txt'), row.names=FALSE, quote=FALSE, sep='\t')
    comgsg<-sg[sample,]
    write.table(comgsg,paste0(pn,'/Anova/',output,'/group.txt'), row.names=FALSE, quote=FALSE, sep='\t')
  }

}

