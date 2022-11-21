#!/usr/bin/Rscript

#多种检验方法的整合

ARGS<-commandArgs(TRUE)
options(stringsAsFactors = FALSE)

if(length(ARGS)!=6){
  cat("Rscript multitest.R input outfold groupFile method FCcutoff pcutoff\n
method: stt: student t test\n
        wtt: Welch t test\n
        ptt: paired t test\n
        mwut: Mann-Whitney U test\n
        ft: F test\n
		fc: fold change
		anova: one-way anova test
		kwt: Kruskal-Wallis test
      ")
  q()
}
suppressMessages(library(stringr))

methods=c('stt','wtt','ptt','mwut','ft', 'fc', 'anova', 'kwt')

input<-as.character(ARGS[1])
ginput<-as.character(ARGS[2])
pn   <-as.character(ARGS[3])
method<-as.character(ARGS[4])
FCcutoff<-as.numeric(ARGS[5])
pcutoff<-as.numeric(ARGS[6])

if(FALSE){
  input<-'E:/projectResearch/07histonePipline/Mode.vs.Control/input.txt'
  ginput<-'E:/projectResearch/07histonePipline/Mode.vs.Control/group.txt'
  pn<-'E:/projectResearch/07histonePipline/Mode.vs.Control'
  method<-'stt'
  FCcutoff=2
  pcutoff=0.05
}


if(method %in% methods){
  cat("Method is OK")
}else{
  cat(paste0(method,":Error: ***********Method parameter is error"))
  q()
}

if(!dir.exists(pn)){
  dir.create(pn, recursive = TRUE)
}

loadData<-function(input){
  if(str_ends(input,".csv")){
    df<-read.csv(input,  stringsAsFactors = FALSE, check.names =FALSE)
    rownames(df)<-df[[1]]
  }else if(str_ends(input,".xlsx")){
    df<-xlsx::read.xlsx(input, sheetIndex = 1)
    rownames(df)<-df[[1]]
  }else{
    df<-read.delim(input, sep="\t", stringsAsFactors = FALSE, check.names = FALSE)
    rownames(df)<-df[[1]]
  }
  return(df)
}

df<-loadData(input)

minValue<- min(df[,2:ncol(df)])

for(i in 2:ncol(df)){
  df[[i]][df[[i]]==minValue]<-NA
}

# df <- df[complete.cases(df),]

gd <- loadData(ginput)
samples<-rownames(gd)
cl <- gd$Type
groups<- as.character(gd$Group)

for(s in samples){
  df[,s]<-as.double(df[,s])
}

m<-as.matrix(df[, samples])
mbak<-as.matrix(df[, samples])

cv<-function(x){
  100*sd(x)/mean(x)
}

se<-function(x){
  sd(x)/sqrt(length(x))
}

# student t test
multi.ttest<-function(x){
  errorflag<-try(t.test(x[which(cl==0)],x[which(cl==1)],var.equal=TRUE),silent=TRUE)
  if('try-error' %in% class(errorflag)){
    NA
  } else{
    t.test(x[which(cl==0)],x[which(cl==1)],var.equal=TRUE)$p.value
  }
}

# paried student t test
multi.pairedttest<-function(x){
  errorflag<-try(t.test(x[which(cl==0)],x[which(cl==1)],var.equal=TRUE, paired=TRUE),silent=TRUE)
  if('try-error' %in% class(errorflag)){
    NA
  }else{
    t.test(x[which(cl==0)],x[which(cl==1)],var.equal=TRUE, paired=TRUE)$p.value
  }
}

# Welch t test
multi.welchttest<-function(x){
  errorflag<-try(t.test(x[which(cl==0)],x[which(cl==1)],var.equal=TRUE),silent=TRUE)
  if('try-error' %in% class(errorflag)){
    NA
  }else{
    t.test(x[which(cl==0)],x[which(cl==1)],var.equal=FALSE)$p.value
  }
}

# Mann-whitney U test
multi.wilcox.test<-function(x){
  errorflag<-try(wilcox.test(x[which(cl==0)],x[which(cl==1)],exact=FALSE,correct=FALSE),silent=TRUE)
  if('try-error' %in% class(errorflag)){
    NA
  }else{
    wilcox.test(x[which(cl==0)],x[which(cl==1)],exact=FALSE,correct=FALSE)$p.value
  }
}

# F test
multi.var.test<-function(x){
  errorflag<-try(var.test(x[which(cl==0)],x[which(cl==1)]), silent=TRUE)
  if('try-error' %in% class(errorflag)){
    NA
  }else{
    var.test(x[which(cl==0)],x[which(cl==1)])$p.value
  }
}

# Fold Change
multi.fc<-function(x){
  mean(x[which(cl==1)]) / mean(x[which(cl==0)])
}

# one-way anova
multi.anova<-function(x){
  A<-groups
  d <- data.frame(x, A)
  aov.mis <- aov(x~A, data = d)
  summary(aov.mis)[[1]][,5][1]
}

# Kruskal-Wallis test
multi.kwt<-function(x){
  A <- groups
  d <- data.frame(x, A)
  y1 <- kruskal.test(x~A, data=d)
  y1$p.value
}


if(method=='stt'){
  df$P.value <- apply(m,1,multi.ttest)
  df$FC <- apply(mbak,1,multi.fc)
  df$Log2FC<-log2(df$FC)
}else if(method == 'wtt'){
  df$P.value <- apply(m,1, multi.welchttest)
  df$FC <- apply(mbak,1,multi.fc)
  df$Log2FC<-log2(df$FC)
}else if(method == 'ptt'){
  df$P.value <- apply(m,1, multi.pairedttest)
  df$FC <- apply(mbak,1,multi.fc)
  df$Log2FC<-log2(df$FC)
}else if(method == 'mwut'){
  df$P.value <- apply(m,1, multi.wilcox.test)
  df$FC <- apply(mbak,1,multi.fc)
  df$Log2FC<-log2(df$FC)
}else if(method == 'ft'){
  df$P.value <- apply(m,1, multi.var.test)
  df$FC <- apply(mbak,1,multi.fc)
  df$Log2FC<-log2(df$FC)
}else if(method == 'fc'){
  df$FC <- apply(mbak,1,multi.fc)
  df$Log2FC<-log2(df$FC)
}else if(method == 'anova'){
  df$P.value <- apply(m,1, multi.anova)
}else if(method == 'kwt'){
  df$P.value <- apply(m, 1, multi.kwt)
}else{
  cat("Method Error\n")
  q()
}

for(groupname in levels(as.factor(groups))){
  meanname <- paste0("Mean.", groupname)
  df[,meanname] <- apply(m[,groups==groupname], 1, mean)
  sdname <- paste0("SD.", groupname)
  df[,sdname] <- apply(m[,groups==groupname], 1, sd)
  sename <- paste0("SE.", groupname)
  df[,sename] <- apply(m[,groups==groupname], 1, se)
  cvname <- paste0("CV.", groupname)
  df[,cvname] <- apply(m[,groups==groupname], 1, cv)
}

if("P.value" %in% colnames(df)){
  df$FDR <- p.adjust(df$P.value, method='fdr')
  #library(qvalue)
  #df$qvalue<-qvalue(df$P.value)
}


df$threshold<-"NoSig"
df$threshold[which(df$Log2FC>=log2(FCcutoff) & df$P.value <=pcutoff)]<-"Up"
df$threshold[which(df$Log2FC<=log2(1/FCcutoff) & df$P.value <=pcutoff)]<-"Down"


write.table(df,file=paste0(pn,"/test.txt"), sep='\t', row.names=FALSE, quote=FALSE)
write.table(df[df$threshold=="Up" | df$threshold=="Down",],file=paste0(pn,"/DEP.txt"), sep='\t', row.names=FALSE, quote=FALSE)
write.table(df[df$threshold=="Up",],file=paste0(pn,"/Up.txt"), sep='\t', row.names=FALSE, quote=FALSE)
write.table(df[df$threshold=="Down",],file=paste0(pn,"/Down.txt"), sep='\t', row.names=FALSE, quote=FALSE)


