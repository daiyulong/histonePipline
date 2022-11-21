#drawHistone
ARGS<-commandArgs(TRUE)
options(stringsAsFactors = FALSE)

if(length(ARGS)!=4){
  cat("Rscript drawerrorBar.R infile testinput  groupFile outfold\n")
  q()
}
suppressMessages(library(tidyr))
suppressMessages(library(reshape2))

infile<-as.character(ARGS[1])
input<-as.character(ARGS[2])
ginput   <-as.character(ARGS[3])
pn<-as.character(ARGS[4])


if (FALSE) {
  infile<-"E:/projectResearch/07histonePipline/infile.csv"
  input<-"E:/projectResearch/07histonePipline/Mode.vs.Control/test.txt"
  ginput<-"E:/projectResearch/07histonePipline/Mode.vs.Control/group.txt"
  pn <- "E:/projectResearch/07histonePipline"
}

d <- read.csv(infile,sep = ",",check.names = F)
tdat <- read.delim(input,header = T,check.names = F)
dg <- read.delim(ginput,header = T,check.names = F)
#pn <- "E:/projectResearch/07histonePipline"

pn <- paste0(pn,"/ErrordBar")
if (!dir.exists(pn)) {
  dir.create(pn,recursive = T)
}

tdat <- merge(d[,1:2],tdat,by = "Group")
pep_list <- unique(tdat[[2]])

for (i in 1:length(pep_list)) {
  name <- pep_list[i]
  d1 <- tdat[which(tdat[[2]]==name),c(1,2,which(grepl("Mean.",colnames(tdat))))]
  dd1 <- melt(d1,id = c("Group","peptide"),variable.name="class",value.name = "meanValue")
  
  d2 <- tdat[which(tdat[[2]]==name),c(1,2,which(grepl("SD.",colnames(tdat))))]
  dd2 <- melt(d2,id = c("Group","peptide"),variable.name="class",value.name = "sdValue")
  
  dd <- cbind(dd1,dd2["sdValue"])
  dd <- separate(dd,class,c(NA,"class"),sep = "\\.")
  dd <- dd[complete.cases(dd),]
  
  library(ggplot2)
  library(ggpubr)
  library(ggthemes)
  theme_set(theme_base())
  dodge <- position_dodge(width=.9)
  
  p <- ggplot(data=dd) +
    geom_bar(aes(x=Group, y=meanValue, fill=class), 
             stat="identity", position=dodge) +
    geom_errorbar(aes(x=Group, ymin=meanValue-sdValue, ymax=meanValue+sdValue, color=class), 
                  stat="identity", position=dodge, width=.3)+
    labs(y = "Relative abundance", x = "") +
    theme_bw()+
    theme(axis.text.x = element_text(size = 8,angle = 40,vjust = 0.7,hjust = 0.7),
          plot.title = element_text(hjust=0.5),
          legend.title = element_blank())+
    labs(title = name)
    
  p
  
  ggsave(paste0(pn,"/",name,"_errorBar.png"),width = 7,height = 7,units = "in",dpi = 300)
  
}


