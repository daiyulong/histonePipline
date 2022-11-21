
ARGS<-commandArgs(TRUE)
options(stringsAsFactors = FALSE)

if(length(ARGS)!=4){
  cat("Rscript drawStackedBar.R input  infile  groupFile outfold\n")
  q()
}
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(ggsci))

input<-as.character(ARGS[1])
infile<-as.character(ARGS[2])
ginput   <-as.character(ARGS[3])
pn<-as.character(ARGS[4])


if(FALSE){
  input<-'E:/projectResearch/07histonePipline/Mode.vs.Control/test.txt'
  infile <- "E:/projectResearch/07histonePipline/infile.csv"
  ginput<-'E:/projectResearch/07histonePipline/Mode.vs.Control/group.txt'
  pn<-'E:/projectResearch/07histonePipline/Mode.vs.Control'
}


###输入文件形如test.txt，但是不能有NA ###
prodata<-read.delim(input,check.names = F)
alldat <- read.csv(infile,check.names = F)
dg <- read.delim(ginput,check.names = F)

pn <- paste0(pn,"/StackedBar")
if (!dir.exists(pn)) {
  dir.create(pn,recursive = T)
}
prodata <- merge(alldat[,1:2],prodata,by = "Group")

p_list <- unique(prodata[[2]])

for (i in 1:length(p_list)) {
  #i=5
  pname <- p_list[i]
  print(pname)
  prodata2 <- prodata[which(prodata[[2]]==pname),]
  prodata2 <- prodata2[,c(1,which(colnames(prodata2) %in% dg[[1]]))]
  
  #删除数据含有非修饰（unmod）的行数据
  prodata2 <- prodata2[!grepl("unmod",prodata2[[1]]),]
  
  #删除数据中含有NA的行数据
  pro_data=prodata2
  pro_data <- pro_data[complete.cases(pro_data),]
  
  #判断绘图数据是否为空
  if (NROW(pro_data)==0) {
    print("绘图数据为空，无法绘图\n")
    next
  }
  
  mycom<-matrix(nrow = nrow(pro_data) ,ncol=3)
  mycom[,1]<-pro_data[,1]
  
  
  i=1
  while(i<nrow(pro_data)+1){
    x<-pro_data[i,2:4]
    y<-pro_data[i,5:7]
    ttest<-t.test(x,y)
    mycom[i,2]<-ttest$p.value
    if(mycom[i,2]<0.001){
      mycom[i,3]=3
    }else if(mycom[i,2]<0.01){
      mycom[i,3]=2
    }else if(mycom[i,2]<0.05){
      mycom[i,3]=1
    }else{
      mycom[i,3]=0
    }
    i=i+1
  }
  

  data1 <- matrix(nrow = nrow(pro_data) ,ncol= 3  )
  colnames(data1) = c("Sample","Mode","Control")
  data1[,1]=pro_data[,1]
  data1[,2]=rowMeans(pro_data[,2:4])
  data1[,3]=rowMeans(pro_data[,5:7])
  
  
  data2<-data.frame(data1) 
  data3 <- melt(data2,id="Sample")
  
  
  i=1
  while (i <= nrow(pro_data) ) {
    data3$sd[i]=sd(pro_data[i,2:4])
    data3$sd[i+nrow(pro_data)]=sd(pro_data[i,5:7])
    i=i+1
  }
  
  colnames(data3)[2] = "Group"
  fin_data <- data3 %>% 
      group_by(Group) %>% 
    arrange(Group,Sample) %>% 
    mutate(cum_num = cumsum(value)) %>% 
    ungroup()
  
  
  fin_data$mycom1<-NA
  fin_data$mycom2<-NA
  fin_data$mycom3<-NA
  
  i=2
  while(i<nrow(mycom)+1){
    if(mycom[i,3]==3){
      fin_data$mycom3[i]=(fin_data$cum_num[i]+fin_data$cum_num[i-1])/2
    }else if(mycom[i,3]==2){
      fin_data$mycom2[i]=(fin_data$cum_num[i]+fin_data$cum_num[i-1])/2
    }else if(mycom[i,3]==1){
      fin_data$mycom1[i]=(fin_data$cum_num[i]+fin_data$cum_num[i-1])/2
    }
    i=i+1
  }
  
  if(mycom[1,3]==3){
    fin_data$mycom3[1]=(fin_data$cum_num[1]+0)/2
  }else if(mycom[1,3]==2){
    fin_data$mycom2[1]=(fin_data$cum_num[1]+0)/2
  }else if(mycom[1,3]==1){
    fin_data$mycom1[1]=(fin_data$cum_num[1]+0)/2
  }
  
  mydata <- fin_data
  mydata$value <- as.numeric(mydata$value)
  mydata$Group <- as.character(mydata$Group)
  # sd列为标准差
  # cum_num列为value累加值（用于误差棒绘制）
  # com列为对应数量星星所在位置
  datcheck <- mydata %>% 
    summarise_all(~all(is.na(.)))
  
  
  p<-ggplot(data=mydata,aes(x=Group,y=value,fill=Sample))+
    labs(x='',y='Proportion of Total Area')+
    scale_y_continuous(expand=c(0,0))+
    geom_bar(stat="identity",colour = "black",width = 0.5,position = position_stack(reverse = TRUE))+
    geom_errorbar(aes(ymin=cum_num,ymax=cum_num+sd),width=0.1)+
    coord_flip()+
    theme_classic()+
    theme(axis.title.x = element_text(size = 10),axis.text =element_text(size = 10),legend.title = element_blank())+
    scale_fill_igv()
  
  
  if (datcheck[["mycom1"]]==FALSE){
    p <- p+annotate("text", x = mydata$Group, y =mydata$mycom1, label ="*",vjust=0.8)
  }else{
    p <- p
  }
  
  if (datcheck[["mycom2"]]==FALSE){
    p <- p+annotate("text", x = mydata$Group, y =mydata$mycom2, label ="**",vjust=0.8)
  }else{
    p <- p
  }
  
  if (datcheck[["mycom3"]]==FALSE){
    p <- p+annotate("text", x = mydata$Group, y =mydata$mycom3, label ="***",vjust=0.8)
  }else{
    p <- p
  }
  
  p
  
  ggsave(paste0(pn,"/",pname, "_StackedBar.png"),units = "in",dpi = 300)
  
}
