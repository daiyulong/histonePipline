#!/usr/bin/env Rscript
# Objective : draw corr,heatmap,boxplot..
# Author: daiyulong
# Time: 2022/11/07
#######################################
ARGS <- commandArgs(TRUE)

if(length(ARGS) < 5){
  print(ARGS)
  cat("Usage: matrixDrawHistone.R input_file group_file outputfold ncol4facet ylab scaleMethod [TRUE|FALSE delete min value]")
  q()
}

suppressMessages(library(tidyverse))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))

input <- as.character(ARGS[1])
inputgroup <- as.character(ARGS[2])
pn    <- as.character(ARGS[3])
ylab <-as.character(ARGS[4])
scaleMethod<-as.character(ARGS[5])

if(FALSE){
  input="E:/projectResearch/07histonePipline/infile.csv"
  inputgroup="E:/projectResearch/07histonePipline/groupFile.txt"
  pn="E:/projectResearch/07histonePipline"
  ylab='Intensity'
  scaleMethod='log2'
  delMin='TRUE'
}

if(length(ARGS)==6){
	delMin <- as.character(ARGS[6])
	if(delMin == 'TRUE'){
		delMin=TRUE
	}else{
		delMin=FALSE
	}
}else{
	delMin<-FALSE
}


if(scaleMethod!='none'){
	ylab=Hmisc::capitalize(paste0(scaleMethod," Scaled ",ylab))
}
###################################
ggsci="lancet"
xlab = "Sample"

##################################

# if(!dir.exists(pn)){
#   dir.create(pn,recursive = T)
# }


if(str_ends(input,'csv')){
  df <- read.csv(input,sep = ",",check.names = FALSE)
}else{
  df <- read.delim(input,check.names = FALSE)
}

if(str_ends(inputgroup,'csv')){
  dg<-read.csv(inputgroup,check.names = FALSE)
}else{
  dg <- read.delim(inputgroup,check.names = FALSE)
}
colnames(dg)[1:2]<- c("Sample","Group")
dg<-dg[,1:2]

# df<-df[,c(1,which(colnames(df) %in% dg[[1]]))]
# bakdf<-df

#####绘制饼图########
ac <- nrow( df[grepl("ac",df[[1]]),])
me1 <- nrow( df[grepl("me1",df[[1]]),])
me2 <- nrow( df[grepl("me2",df[[1]]),])
me3 <- nrow( df[grepl("me3",df[[1]]),])

datt <- data.frame("Acetyl"=ac,"Monomethyl"=me1,"Dimethyl"=me2,"Trimethyl"=me3)
datt <- as.data.frame(t(datt))
datt$type <- rownames(datt)

datt <- datt[,c(2,1)]
colnames(datt) <- c("type","n")

#分组统计并计算百分比
datpie <- datt %>% 
  mutate(percent = 100 * n/sum(n))


p <- ggplot(datpie,aes(x="",y=n,fill=type,group=type))+
  geom_bar(width = 1,stat = "identity")+
  coord_polar("y",start = 0)+
  
  theme_void()+theme(legend.position = "right")+
  geom_text(aes(1.3,label =paste0(round(percent, 1),"%")),position =  position_stack(vjust=0.5), size=2)
p

ggsave(paste0(pn,"/pie.png"),units = "in",dpi = 300)
ggsave(paste0(pn,"/pie.pdf"))


################
df<-df[,c(1,which(colnames(df) %in% dg[[1]]))]
bakdf<-df

minValue<- min(df[,2:ncol(df)])
if(delMin){
  for(i in 2:ncol(df)){
    df[[i]][df[[i]]==minValue]<-NA
  }
}
cn <- colnames(df)

width = (NCOL(df)-1)*0.8+2
height = 6

ldf<-melt(df,id.vars=cn[1],measure.vars=cn[2:length(cn)],variable.name="Sample",value.name="Intensity")
names(ldf)<-c("ID","Sample","Intensity")


d<-merge(ldf, dg, by='Sample', all=TRUE)
colnames(d)
d<-d[!is.na(d$Intensity) & d$Intensity!=Inf & d$Intensity!=-Inf,]


if(delMin){
	d<-d[!is.na(d$Intensity),]
}

if(scaleMethod == "log2"){
	if(any(d$Intensity==0,na.rm=TRUE)){
		cat("Use log2(x+1)\n")
		scaleMethod="log2plus";
		d$Intensity<- -log2(d$Intensity+1)
	}else{
		d$Intensity<- -log2(d$Intensity)
	}
}else if(scaleMethod == "log2plus"){
	d$Intensity<- -log2(d$Intensity+1)
}else if(scaleMethod == "log10"){
	if(any(d$Intensity==0)){
		cat("Use log10(x+1)\n")
		scaleMethod="log10plus";
		d$Intensity<- log10(d$Intensity+1)
	}else{
		d$Intensity<- log10(d$Intensity)
	}
}

dataConvert<-function(d,method){
  multi.zscore<-function(x){
    (x-mean(x))/sd(x)
  }
  
  if(method == 'log2'){
    d<-cbind(d[[1]],-log2(d[,2:NCOL(d)]))
  }else if(method == "log2plus"){
    d<-cbind(d[[1]],-log2(d[,2:NCOL(d)]+1))
  }else if(method == 'log10'){
    d<-cbind(d[[1]],log10(d[,2:NCOL(d)]))
  }else if(method == 'zscore'){
    d<-cbind(d[[1]], apply(d[,2:NCOL(d)], 1, multi.zscore))
  }else if(method == 'log10plus'){
    d<-cbind(d[[1]], log10(d[,2:NCOL(d)]+1))
  }
  return(d)
}

scaledf<-dataConvert(df,scaleMethod)

#plotBoxplotGroup
plotBoxplotGroup2 <- function(df, plotPoint=TRUE, plotMeanPoint=TRUE, plotBox=TRUE, outimage="", title="", xlab="", ylab="", filllabel="",
                             colorlabel="", palette="Set2", width=8, height=6, units='in', dpi=300, ggsci=NA){
  library(ggplot2)
  library(cowplot)
  library(ggsci)
  mycols<- c(ggsci::pal_lancet("lanonc")(9)[1:8], ggsci::pal_jama("default")(7), ggsci::pal_nejm("default")(8), ggsci::pal_jco("default")(10))
  
  colnames(df)<-c("Sample","Group","Value")
  df$Group<-as.factor(df$Group)
  df$Sample<-as.factor(df$Sample)
  GroupNumber<-length(levels(factor(df$Group)))
  SimpleNumber<-length(unique(df$Sample))
  
  if(SimpleNumber>20){
    width = SimpleNumber*0.5+2
    barwidth=0.7
    pointsize=1
  }else{
    barwidth=0.5
    pointsize=2
  }
  if(width>40){
    width=40
  }
  
  df$Value[is.na(df$Value)]<-0
  
  
  if(plotPoint){
    p<-ggplot(df,aes(x=Sample, y=Value))+
      geom_boxplot(fill=NA,width=barwidth,outlier.colour=NA)+
      ggtitle(title)+xlab(xlab)+ylab(ylab)+labs(color=colorlabel)+
      #scale_fill_brewer(palette=palette)+
      theme_cowplot()
    p<-p+geom_quasirandom(aes(colour=Group),width=0.2)
    
    if(is.na(ggsci)){
      p<-p+scale_color_brewer(palette=palette)
    }else{
      p<-addGGsci(p,'color', ggsci)
    }
    
  }else{
    p<-ggplot(df,aes(x=Sample,y=Value,fill=Group))+
      geom_boxplot(width=barwidth,outlier.colour="darkred")+
      ggtitle(title)+xlab(xlab)+ylab(ylab)+labs(fill=filllabel)+
      theme_cowplot()
    if(is.na(ggsci)){
      p<-p+scale_fill_manual(values=mycols[1:GroupNumber])
    }else{
      p<-addGGsci(p,'fill',ggsci)
    }
  }
  if(plotMeanPoint){
    p<-p+ stat_summary(fun=mean, colour="#ECECFF", geom="point", shape=18, size=pointsize, show.legend = FALSE)
  }
  
  p<-p+theme(axis.text.x=element_text(angle=45,hjust=1))
  print(p)
  
  if(outimage!=""){
    ggsave(outimage,dpi=dpi,width=width,height = height, units = units)
  }
  return(p)
}

#plotViolinGroup
plotViolinGroup2 <- function(df, plotPoint=FALSE, plotMeanPoint=FALSE, plotBox=TRUE, outimage="", title="", xlab="", ylab="", filllabel="",
                            colorlabel="", palette="Set2", width=8, height=6, units='in', dpi=300, ggsci=NA){
  colnames(df)<-c("Sample","Group","Value")
  df$Value[is.na(df$Value)]<-0
  
  library(ggplot2)
  library(cowplot)
  library(ggsci)
  mycols<- c(ggsci::pal_lancet("lanonc")(9), ggsci::pal_jama("default")(7), ggsci::pal_nejm("default")(8), ggsci::pal_jco("default")(10))
  GroupNumber<-length(levels(factor(df$Group)))
  SimpleNumber<-length(unique(df$Sample))
  
  if(SimpleNumber>20){
    width = SimpleNumber*0.5+2
    barwidth=0.7
    pointsize=1
  }else{
    barwidth=0.5
    pointsize=2
  }
  print(SimpleNumber)
  print(barwidth)
  
  if(width>40){
    width=40
  }
  
  if(plotPoint){
    p<-ggplot(df,aes(x=Sample,y=Value))+
      #geom_violin(width=.8,trim=FALSE)+
      #geom_boxplot(fill=NA,width=0.5,outlier.colour=NA)+
      ggtitle(title)+xlab(xlab)+ylab(ylab)+labs(color=colorlabel)+
      #scale_fill_brewer(palette=palette)+
      theme_cowplot()
    
    p<-p+geom_quasirandom(aes(colour=Group),width=0.2)
    if(is.na(ggsci)){
      p<-p+scale_color_brewer(palette=palette)
    }else{
      p<-addGGsci(p,'color', ggsci)
    }
  }else{
    p<-ggplot(df,aes(x=Sample,y=Value,fill=Group))+
      geom_violin(width=.8,trim=FALSE)+
      #geom_boxplot(width=0.5,outlier.colour="darkred")+
      ggtitle(title)+xlab(xlab)+ylab(ylab)+labs(fill=filllabel)+
      theme_cowplot()
    if(is.na(ggsci)){
      p<-p+scale_fill_manual(values=mycols[1:GroupNumber])
    }else{
      p<-addGGsci(p,'fill',ggsci)
    }
  }
  if(plotBox){
    if(plotPoint){
      p<-p+geom_boxplot(fill=NA,width=0.2,outlier.colour=NA)
    }else{
      p<-p+geom_boxplot(fill='#FCFCFC',width=0.2,outlier.colour=NA)
    }
    
  }
  if(plotMeanPoint){
    p<-p+ stat_summary(fun=mean, colour="#ECECFF", geom="point", shape=18, size=pointsize, show.legend = FALSE)
  }
  if(length(df$Group[!duplicated(df$Group)]) == 1){
    p<-p + theme(legend.position = "none",axis.text.x=element_text(angle=45,hjust=1))
  }else{
    p<-p+theme(axis.text.x=element_text(angle=45,hjust=1))
  }
  
  print(p)
  
  if(outimage!=""){
    ggsave(outimage,dpi=dpi,width=width,height = height, units = units)
  }
  return(p)
}

#boxplot分析图
plotBoxplotGroup2(d[,c("Sample","Group","Intensity")], outimage = paste0(pn,"/boxplot.png"),
                 plotPoint=FALSE, width=width, height=4, ylab=ylab)
plotBoxplotGroup2(d[,c("Sample","Group","Intensity")], outimage = paste0(pn,"/boxplot.pdf"),
                 plotPoint=FALSE, width=width, height=4, ylab=ylab)

#violin分析图
plotViolinGroup2(d[,c("Sample","Group","Intensity")], outimage = paste0(pn,"/violin.png"),plotMeanPoint=TRUE,
                plotPoint=FALSE, width=width, height=4, ylab=ylab)
plotViolinGroup2(d[,c("Sample","Group","Intensity")], outimage = paste0(pn,"/violin.pdf"),plotMeanPoint=TRUE,
                plotPoint=FALSE, width=width, height=4, ylab=ylab)


#相关性分析图
mdf <- scaledf[2:ncol(scaledf)]
rownames(mdf)<-as.character(scaledf[[1]])
m<-as.matrix(mdf)
#rownames(mdf)<-as.vector(df[,1])
mcor <- cor(m,use='complete.obs')

SampleNumber<-length(unique(d$Sample))
if(SampleNumber>30){
  tl.cex=0.5
}else if(SampleNumber>20){
  tl.cex=0.6
}else if(SampleNumber>15){
  tl.cex=0.8
}else{
  tl.cex=1
}

library(corrplot)
if(SampleNumber<15){
  palette_1 <- RColorBrewer::brewer.pal(n=11, name = "RdYlGn") 
  palette_2 <- rev(palette_1)
  png(paste0(pn,"/corr.png"), width=7, height=7, units='in', res=300)
  corrplot(mcor,
           method = "circle",
           number.cex = 0.7,
           tl.cex = tl.cex,
           tl.srt = 45,
           #method = "ellipse",
           addCoef.col = "blue",
           number.digits=2,
           type = "lower", col = palette_2, diag = FALSE, mar = c(1,1,1,1))
  dev.off()
  cairo_pdf(paste0(pn,"/corr.pdf"), width=7, height=7)
  corrplot(mcor,
           method = "circle",
           number.cex = 0.7,
           tl.cex = tl.cex,
           tl.srt = 45,
           #method = "ellipse",
           addCoef.col = "blue",
           number.digits=2,
           type = "lower", col = palette_2, diag = FALSE, mar = c(1,1,1,1))
  dev.off()
  
  
}else{
  palette_1 <- RColorBrewer::brewer.pal(n=11, name = "RdYlGn") 
  palette_2 <- rev(palette_1)
  png(paste0(pn,"/corr.png"), width=7, height=7, units='in', res=300)
  corrplot(mcor,
           method = "pie",
           #number.cex = 0.7,
           tl.cex = tl.cex,
           tl.srt = 45,
           #method = "ellipse",
           #addCoef.col = "blue",
           type = "lower", col = palette_2, diag = FALSE, mar = c(1,1,1,1))
  dev.off()
  cairo_pdf(paste0(pn,"/corr.pdf"), width=7, height=7)
  corrplot(mcor,
           method = "pie",
           #number.cex = 0.7,
           tl.cex = tl.cex,
           tl.srt = 45,
           #method = "ellipse",
           #addCoef.col = "blue",
           
           type = "lower", col = palette_2, diag = FALSE, mar = c(1,1,1,1))
  dev.off()
  
}


# 无scale的相关性图
mdf <- df[2:ncol(df)]
rownames(mdf)<-as.character(df[[1]])
m<-as.matrix(mdf)
outfold<-paste0(pn, '/../images/NonScale')
if(!dir.exists(outfold)){
  dir.create(outfold, recursive = TRUE)
}

write.csv(m, paste0(outfold,"/cor_input.csv"))
mcor <- cor(m, use='complete.obs')
write.csv(mcor,paste0(outfold,'/cor.csv'))


SampleNumber<-length(unique(d$Sample))
if(SampleNumber>30){
  tl.cex=0.5
}else if(SampleNumber>20){
  tl.cex=0.6
}else if(SampleNumber>15){
  tl.cex=0.8
}else{
  tl.cex=1
}

library(corrplot)
if(SampleNumber<15){
  palette_1 <- RColorBrewer::brewer.pal(n=11, name = "RdYlGn") 
  palette_2 <- rev(palette_1)
  png(paste0(outfold,"/corr.png"), width=7, height=7, units='in', res=300)
  corrplot(mcor,
           method = "circle",
           number.cex = 0.7,
           tl.cex = tl.cex,
           tl.srt = 45,
           #method = "ellipse",
           addCoef.col = "blue",
           number.digits=2,
           type = "lower", col = palette_2, diag = FALSE, mar = c(1,1,1,1))
  dev.off()
  cairo_pdf(paste0(outfold,"/corr.pdf"), width=7, height=7)
  corrplot(mcor,
           method = "circle",
           number.cex = 0.7,
           tl.cex = tl.cex,
           tl.srt = 45,
           #method = "ellipse",
           addCoef.col = "blue",
           number.digits=2,
           type = "lower", col = palette_2, diag = FALSE, mar = c(1,1,1,1))
  dev.off()
  
  
}else{
  palette_1 <- RColorBrewer::brewer.pal(n=11, name = "RdYlGn") 
  palette_2 <- rev(palette_1)
  png(paste0(outfold,"/corr.png"), width=7, height=7, units='in', res=300)
  corrplot(mcor,
           method = "pie",
           #number.cex = 0.7,
           tl.cex = tl.cex,
           tl.srt = 45,
           #method = "ellipse",
           #addCoef.col = "blue",
           type = "lower", col = palette_2, diag = FALSE, mar = c(1,1,1,1))
  dev.off()
  cairo_pdf(paste0(outfold,"/corr.pdf"), width=7, height=7)
  corrplot(mcor,
           method = "pie",
           #number.cex = 0.7,
           tl.cex = tl.cex,
           tl.srt = 45,
           #method = "ellipse",
           #addCoef.col = "blue",
           
           type = "lower", col = palette_2, diag = FALSE, mar = c(1,1,1,1))
  dev.off()
  
}


checkValue<-function(x){
  if(sum(x==minValue) > length(x)*0.8){
    return(FALSE)
  }else{
    return(TRUE)
  }
}

#heatmap绘图
library(pheatmap)
annotation_col <- dg
rownames(annotation_col)<-dg[[1]]
annotation_col$Sample<-NULL

mmdf<-bakdf[,2:ncol(bakdf)]
rownames(mmdf)<-as.character(bakdf[[1]])

if(delMin){
  flag <- apply(mmdf, 1, checkValue)
  mdf<-mmdf[flag,]
}else{
  mdf<-mmdf
}
border=NA


if(NROW(mdf)>50){
  try(
    pheatmap(mdf,
             filename=paste0(pn,"/heatmap.png"),
             cellwidth=12,
             cellheight=3,
             scale = 'row',
             cluster_cols = TRUE,
             cluster_rows = TRUE,
             show_colnames = TRUE,
             annotation_col=annotation_col,
             #color =colorRampPalette(c("blue","white","red"))(51),
             color =colorRampPalette(c("blue","snow2","red"))(51),

             border=border,
             clustering_distance_rows = "euclidean",#correlation,"euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
             clustering_distance_cols  = "correlation",#"euclidean",#"correlation",
             clustering_method = "average",#"ward.D", "ward.D2", "single", "complete", "average" (= UPGMA),
             #                            "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)
             show_rownames = F)
  )
  
  try(
    pheatmap(mdf,
             filename=paste0(pn,"/heatmap_noClusterCol.png"),
             cellwidth=12,
             cellheight=3,
             scale = 'row',
             cluster_cols = FALSE,
             cluster_rows = TRUE,
             show_colnames = TRUE,
             annotation_col=annotation_col,
             #color =colorRampPalette(c("blue","white","red"))(51),
             color =colorRampPalette(c("blue","snow2","red"))(51),
             
             border=border,
             clustering_distance_rows = "euclidean",#correlation,"euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
             clustering_distance_cols  = "correlation",#"euclidean",#"correlation",
             clustering_method = "average",#"ward.D", "ward.D2", "single", "complete", "average" (= UPGMA),
             #                            "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)
             show_rownames = F)
  )
  
  try(
  pheatmap(mdf,
           filename=paste0(pn,"/heatmap.pdf"),
           cellwidth=12,
           cellheight=3,
           scale = 'row',
           cluster_cols = TRUE,
           cluster_rows = TRUE,
           show_colnames = TRUE,
           annotation_col=annotation_col,
           #color =colorRampPalette(c("blue","white","red"))(51),
           color =colorRampPalette(c("blue","snow2","red"))(51),

           border=border,
           clustering_distance_rows = "euclidean",#correlation,"euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
           clustering_distance_cols  = "euclidean",#"correlation",#"correlation",
           clustering_method = "average",#"complete",#"ward.D", "ward.D2", "single", "complete", "average" (= UPGMA),
           #                            "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)
           show_rownames = FALSE)
  )
  try(
  pheatmap(mdf,
           filename=paste0(pn,"/heatmap_noClusterCol.pdf"),
           cellwidth=12,
           cellheight=3,
           scale = 'row',
           cluster_cols = FALSE,
           cluster_rows = TRUE,
           show_colnames = TRUE,
           annotation_col=annotation_col,
           #color =colorRampPalette(c("blue","white","red"))(51),
           color =colorRampPalette(c("blue","snow2","red"))(51),
           
           border=border,
           clustering_distance_rows = "euclidean",#correlation,"euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
           clustering_distance_cols  = "euclidean",#"correlation",#"correlation",
           clustering_method = "average",#"complete",#"ward.D", "ward.D2", "single", "complete", "average" (= UPGMA),
           #                            "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)
           show_rownames = FALSE)
  
  )
}else{
  try(
    pheatmap(mdf,
             filename=paste0(pn,"/heatmap.png"),
             cellwidth=12,
             height=0.5*NROW(mdf),
             scale = 'row',
             cluster_cols = TRUE,
             cluster_rows = TRUE,
             show_colnames = TRUE,
             annotation_col=annotation_col,
             #color =colorRampPalette(c("blue","white","red"))(51),
             color =colorRampPalette(c("blue","snow2","red"))(51),
             fontsize_row = 6,
             border=border,
             clustering_distance_rows = "euclidean",#correlation,"euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
             clustering_distance_cols  = "correlation",#"euclidean",#"correlation",
             clustering_method = "average",#"ward.D", "ward.D2", "single", "complete", "average" (= UPGMA),
             #                            "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)
             show_rownames = T)
  )
  
  try(
    pheatmap(mdf,
             filename=paste0(pn,"/heatmap_noClusterCol.png"),
             cellwidth=12,
             height=0.5*NROW(mdf),
             scale = 'row',
             cluster_cols = FALSE,
             cluster_rows = TRUE,
             show_colnames = TRUE,
             annotation_col=annotation_col,
             #color =colorRampPalette(c("blue","white","red"))(51),
             color =colorRampPalette(c("blue","snow2","red"))(51),
             fontsize_row = 6,
             border=border,
             clustering_distance_rows = "euclidean",#correlation,"euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
             clustering_distance_cols  = "correlation",#"euclidean",#"correlation",
             clustering_method = "average",#"ward.D", "ward.D2", "single", "complete", "average" (= UPGMA),
             #                            "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)
             show_rownames = T)
  )
  try(
  pheatmap(mdf,
           filename=paste0(pn,"/heatmap.pdf"),
           cellwidth=12,
           height=8,
           scale = 'row',
           cluster_cols = TRUE,
           cluster_rows = TRUE,
           show_colnames = TRUE,
           annotation_col=annotation_col,
           #color =colorRampPalette(c("blue","white","red"))(51),
           color =colorRampPalette(c("blue","snow2","red"))(51),
           fontsize_row = 6,
           border=border,
           clustering_distance_rows = "euclidean",#correlation,"euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
           clustering_distance_cols  = "euclidean",#"correlation",#"correlation",
           clustering_method = "average",#"ward.D", "ward.D2", "single", "complete", "average" (= UPGMA),
           #                            "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)
           show_rownames = TRUE)
  )
  try(
  pheatmap(mdf,
           filename=paste0(pn,"/heatmap_noClusterCol.pdf"),
           cellwidth=12,
           height=8,
           scale = 'row',
           cluster_cols = FALSE,
           cluster_rows = TRUE,
           show_colnames = TRUE,
           annotation_col=annotation_col,
           #color =colorRampPalette(c("blue","white","red"))(51),
           color =colorRampPalette(c("blue","snow2","red"))(51),
           fontsize_row = 6,
           border=border,
           clustering_distance_rows = "euclidean",#correlation,"euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
           clustering_distance_cols  = "euclidean",#"correlation",#"correlation",
           clustering_method = "average",#"ward.D", "ward.D2", "single", "complete", "average" (= UPGMA),
           #                            "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)
           show_rownames = TRUE)
  )
}


