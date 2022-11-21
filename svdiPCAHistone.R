#!/usr/bin/Rscript
#Author: Yu Hong
#PCA分析，使用SVDimpute算法，可以做Cross-Valiation, 可以容忍大量的缺失数据 , 缺失数据>10%

args<-commandArgs(TRUE)

if(length(args)!=4){
  cat("\nUsage: Rscript svdiPCA.R input groupfile outprefix scalemethod[none|pareto|vector|uv]\n")
  q()
}

suppressMessages(library(pcaMethods))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(ggsci))
suppressMessages(library(stringr))

#data(metaboliteData)
scalemethods<-c("none", "pareto", "vector", "uv")

mycols<- c(pal_lancet("lanonc")(9)[1:8], pal_jama("default")(7), pal_nejm("default")(8), pal_jco("default")(10))

input<-as.character(args[1])
groupfile<-as.character(args[2])
pn<-as.character(args[3])
sm <- as.character(args[4])


if(FALSE){
  input<-"E:/projectResearch/07histonePipline/infile.csv"
  groupfile<-"E:/projectResearch/07histonePipline/groupFile.txt"
  pn<-"E:/projectResearch/07histonePipline"
  sm<-'uv'
}

if(sm %in% scalemethods){
  print(sm)
}else{
  print("scale methods error, only support: none, pareto, vector, uv")
  q()
}

# pn <- paste0(pn,"/SampleAnalysis")
# 
# if(!dir.exists(pn)){
#   dir.create(pn,recursive = TRUE)
# }


###################define##################
width=8
height=6
dpi=300
units='in'
scree.num=5

###########################################
if (str_ends(input,"csv")) {
  dd <- read.csv(input,check.names = FALSE)
}else if (str_ends(input,"xlsx")) {
  dd <- read.xlsx(input,1)
}else{
  dd <- read.delim(input,row.names = 1,check.names = FALSE,quote = "")
}
# dd<-read.delim(input, row.names = 1, check.names = FALSE, quote="")
gdf<-read.delim(groupfile,row.names=1)

if(NROW(gdf)<3){
  cat("Sample Number < 3, Cannot do PCA\n\n")
  q()
}

dd<-dd[,rownames(gdf)]
tdd<-t(dd)

#对数据进行scale
md<-prep(tdd, scale=sm, center=TRUE)
rownames(md)<-rownames(tdd)
colnumb<-NCOL(md)

nPcs<- NROW(dd)

if(nPcs > 5){
  nPcs=5
}
print(nPcs)

resSVDI <- pca(md, method="svdImpute", center=FALSE, nPcs=nPcs)

#q2SVDI <- Q2(resSVDI, md, fold=7,  type="impute")
if(NCOL(md)< 7){
  cvfold <- NCOL(md)
}else{
  cvfold <-7
}
q2SVDI <- Q2(resSVDI, md, fold=cvfold)

#d <- merge(md, scores(resSVDI), by=0)
d<-scores(resSVDI)
d<-as.data.frame(d)
d$Group<-gdf[rownames(d),1]
d$Sample<-rownames(d)
rownames(d)<-d$Sample
z<-sprintf("%.2f",resSVDI@R2*100)

row_numb=length(q2SVDI)
outdd<-data.frame(Component=dimnames(resSVDI@scores)[[2]][1:row_numb], R2=resSVDI@R2[1:row_numb], R2cum=resSVDI@R2cum[1:row_numb], Q2=q2SVDI)
outdd$Q2cum<-cumsum(outdd$Q2)
write.table(outdd,paste0(pn,"/pca_Variance.xls"), sep='\t',row.names = FALSE, quote=FALSE)


Cairo::Cairo(file = paste0(pn,"/pca_scree.png"), unit="in", dpi=dpi, width=width, height=height, type="png", bg="white");

pcvars<-resSVDI@R2[1:row_numb]
names(pcvars)<-dimnames(resSVDI@scores)[[2]][1:row_numb]
cumvars<-resSVDI@R2cum[1:row_numb]
names(cumvars)<-dimnames(resSVDI@scores)[[2]][1:row_numb]
ylims <- range(c(pcvars,cumvars));
extd<-(ylims[2]-ylims[1])/10
miny<- ifelse(ylims[1]-extd>0, ylims[1]-extd, 0);
maxy<- ifelse(ylims[2]+extd>1, 1.0, ylims[2]+extd);
par(mar=c(5,5,6,3));
plot(pcvars, type='l', col='blue', main='Scree plot', xlab='PC index', ylab='Variance explained', ylim=c(miny, maxy), axes=F)
text(pcvars, labels =paste(100*round(pcvars,3),'%'), adj=c(-0.3, -0.5), srt=45, xpd=T)
points(pcvars, col='red');

lines(cumvars, type='l', col='green')
text(cumvars, labels =paste(100*round(cumvars,3),'%'), adj=c(-0.3, -0.5), srt=45, xpd=T)
points(cumvars, col='red');

abline(v=1:scree.num, lty=3);
axis(2);
axis(1, 1:length(pcvars), 1:length(pcvars));
dev.off();

write.table(d,paste0(pn,"/pca_scores.xls"), row.names = FALSE, quote = FALSE, sep="\t")

png(paste(pn,"/biplot.png",sep=""),res=dpi,width=width,height = height, units = units)
biplot(resSVDI,cex=.1)
dev.off()



group_numb = length(unique(d$Group))
shape.vec<- (1:group_numb + 14)  %% 26

ddtype<-gdf[rownames(d),c("Group", "Shape","Color")]
ddtype<-ddtype[!duplicated(ddtype$Group),]
shape.man<-ddtype$Shape
names(shape.man)<-ddtype$Group

color.man<-ddtype$Color
names(color.man)<-ddtype$Group

################第1，2主成分散点图####################
##label
p<-ggplot(d,aes(x=PC1,y=PC2))+
  stat_ellipse(type="norm",geom="polygon",alpha=0.5,color="black", fill=NA, size=1)+
  geom_vline(xintercept = 0, color="black")+
  geom_hline(yintercept = 0, color="black")+
  geom_point(aes(fill=Group, shape=Group),size=4, alpha=0.9)+
  geom_text_repel(label=d$Sample,size=3, colour='black', force=1, max.overlaps = Inf,
                  #min.segment.length = 0,
                  nudge_y=(range(d$PC2)[2] - range(d$PC2)[1])/20 * (ifelse(d$PC2>0,1,-1)))+
  xlab(paste("PC1(",z[1],"%)",sep=""))+
  ylab(paste("PC2(",z[2],"%)",sep=""))+
  scale_shape_manual(values=shape.man)+
  scale_color_manual(values=color.man)+
  scale_fill_manual(values=color.man)+
  stat_ellipse(aes(fill=Group),type="norm",geom="polygon",alpha=0.1,color=NA)+
  #guides(fill="none")+
  theme_bw()+
  theme(panel.grid = element_blank())

ggsave(paste0(pn,"/PCA1.2.png"),plot=p,dpi=dpi,width=width,height = height, units =units)
ggsave(paste0(pn,"/PCA1.2.pdf"),plot=p,dpi=dpi,width=width,height = height, units =units)

##no labels
p<-ggplot(d,aes(x=PC1,y=PC2))+
  stat_ellipse(type="norm",geom="polygon",alpha=0.5,color="black", fill=NA, size=1)+
  geom_vline(xintercept = 0, color="black")+
  geom_hline(yintercept = 0, color="black")+
  geom_point(aes(fill=Group, shape=Group),size=4, alpha=0.9)+
  #geom_text_repel(label=d$Sample,size=3, colour='black', force=1,
  #                nudge_y=(range(d$PC2)[2] - range(d$PC2)[1])/20 * (ifelse(d$PC2>0,1,-1)))+
  xlab(paste("PC1(",z[1],"%)",sep=""))+
  ylab(paste("PC2(",z[2],"%)",sep=""))+
  scale_shape_manual(values=shape.man)+
  scale_color_manual(values=color.man)+
  scale_fill_manual(values=color.man)+
  stat_ellipse(aes(fill=Group),type="norm",geom="polygon",alpha=0.1,color=NA)+
  #guides(fill="none")+
  theme_bw()+
  theme(panel.grid = element_blank())


ggsave(paste0(pn,"/PCA1.2_nolabel.png"),plot=p,dpi=dpi,width=width,height = height, units =units)
ggsave(paste0(pn,"/PCA1.2_nolabel.pdf"),plot=p,dpi=dpi,width=width,height = height, units =units)
#png(paste("PCA/",pn,".PCA1.2.3.png",sep=""),res=dpi,width=height,height = height, units = units)


################第1，2主成分散点图####################
##label
p<-ggplot(d,aes(x=PC1,y=PC2))+
  #stat_ellipse(type="norm",geom="polygon",alpha=0.5,color="black", fill=NA, size=1)+
  #geom_vline(xintercept = 0, color="black")+
  #geom_hline(yintercept = 0, color="black")+
  geom_point(aes(fill=Group, shape=Group),size=4, alpha=0.9)+
  geom_text_repel(label=d$Sample,size=3, colour='black', force=1, max.overlaps = Inf,#min.segment.length = 0,
                  nudge_y=(range(d$PC2)[2] - range(d$PC2)[1])/20 * (ifelse(d$PC2>0,1,-1)))+
  xlab(paste("PC1(",z[1],"%)",sep=""))+
  ylab(paste("PC2(",z[2],"%)",sep=""))+
  scale_shape_manual(values=shape.man)+
  scale_color_manual(values=color.man)+
  scale_fill_manual(values=color.man)+
  #stat_ellipse(aes(fill=Group),type="norm",geom="polygon",alpha=0.1,color=NA)+
  #guides(fill="none")+
  theme_bw()#+
  #theme(panel.grid = element_blank())

ggsave(paste0(pn,"/PCA1.2.s2.png"),plot=p,dpi=dpi,width=width,height = height, units =units)
ggsave(paste0(pn,"/PCA1.2.s2.pdf"),plot=p,dpi=dpi,width=width,height = height, units =units)

##no labels
p<-ggplot(d,aes(x=PC1,y=PC2))+
  #stat_ellipse(type="norm",geom="polygon",alpha=0.5,color="black", fill=NA, size=1)+
  #geom_vline(xintercept = 0, color="black")+
  #geom_hline(yintercept = 0, color="black")+
  geom_point(aes(fill=Group, shape=Group),size=4, alpha=0.9)+
  #geom_text_repel(label=d$Sample,size=3, colour='black', force=1,
  #                nudge_y=(range(d$PC2)[2] - range(d$PC2)[1])/20 * (ifelse(d$PC2>0,1,-1)))+
  xlab(paste("PC1(",z[1],"%)",sep=""))+
  ylab(paste("PC2(",z[2],"%)",sep=""))+
  scale_shape_manual(values=shape.man)+
  scale_color_manual(values=color.man)+
  scale_fill_manual(values=color.man)+
  #stat_ellipse(aes(fill=Group),type="norm",geom="polygon",alpha=0.1,color=NA)+
  #guides(fill="none")+
  theme_bw()#+
  #theme(panel.grid = element_blank())


ggsave(paste0(pn,"/PCA1.2_nolabel.s2.png"),plot=p,dpi=dpi,width=width,height = height, units =units)
ggsave(paste0(pn,"/PCA1.2_nolabel.s2.pdf"),plot=p,dpi=dpi,width=width,height = height, units =units)
#png(paste("PCA/",pn,".PCA1.2.3.png",sep=""),res=dpi,width=height,height = height, units = units)


################第1，2主成分散点图####################
##label
p<-ggplot(d,aes(x=PC1,y=PC2))+
  stat_ellipse(type="norm",geom="polygon",alpha=0.5,color="black", fill=NA, size=1)+
  geom_vline(xintercept = 0, color="black")+
  geom_hline(yintercept = 0, color="black")+
  geom_point(aes(fill=Group, shape=Group),size=4, alpha=0.9)+
  geom_text_repel(label=d$Sample,size=3, colour='black', force=1, max.overlaps = Inf,#min.segment.length = 0,
                  nudge_y=(range(d$PC2)[2] - range(d$PC2)[1])/20 * (ifelse(d$PC2>0,1,-1)))+
  xlab(paste("PC1(",z[1],"%)",sep=""))+
  ylab(paste("PC2(",z[2],"%)",sep=""))+
  scale_shape_manual(values=shape.man)+
  scale_color_manual(values=color.man)+
  scale_fill_manual(values=color.man)+
  #stat_ellipse(aes(fill=Group),type="norm",geom="polygon",alpha=0.1,color=NA)+
  #guides(fill="none")+
  theme_bw()#+
#theme(panel.grid = element_blank())

ggsave(paste0(pn,"/PCA1.2.s3.png"),plot=p,dpi=dpi,width=width,height = height, units =units)
ggsave(paste0(pn,"/PCA1.2.s3.pdf"),plot=p,dpi=dpi,width=width,height = height, units =units)

##no labels
p<-ggplot(d,aes(x=PC1,y=PC2))+
  stat_ellipse(type="norm",geom="polygon",alpha=0.5,color="black", fill=NA, size=1)+
  geom_vline(xintercept = 0, color="black")+
  geom_hline(yintercept = 0, color="black")+
  geom_point(aes(fill=Group, shape=Group),size=4, alpha=0.9)+
  #geom_text_repel(label=d$Sample,size=3, colour='black', force=1,
  #                nudge_y=(range(d$PC2)[2] - range(d$PC2)[1])/20 * (ifelse(d$PC2>0,1,-1)))+
  xlab(paste("PC1(",z[1],"%)",sep=""))+
  ylab(paste("PC2(",z[2],"%)",sep=""))+
  scale_shape_manual(values=shape.man)+
  scale_color_manual(values=color.man)+
  scale_fill_manual(values=color.man)+
  #stat_ellipse(aes(fill=Group),type="norm",geom="polygon",alpha=0.1,color=NA)+
  #guides(fill="none")+
  theme_bw()#+
#theme(panel.grid = element_blank())


ggsave(paste0(pn,"/PCA1.2_nolabel.s3.png"),plot=p,dpi=dpi,width=width,height = height, units =units)
ggsave(paste0(pn,"/PCA1.2_nolabel.s3.pdf"),plot=p,dpi=dpi,width=width,height = height, units =units)

##label

#三维散点图
library(scatterplot3d)
png(paste(pn,"/PCA123.png",sep=""),res=dpi,width=width,height = height, units = units)

colorDefine<- function(v, mycols){
  v<-as.factor(v)
  levels<-levels(v)
  N <- length(levels)
  vv<-as.character(v)
  for(i in 1:N){
    vv[vv==levels[i]] <- mycols[i]
  }
  return(vv)
}

colors<-colorDefine(d$Group, mycols)
s3d<-scatterplot3d(d[,c("PC1","PC2","PC3")],xlab=paste0("PC1(",z[1],"%)"),ylab=paste0("PC2(",z[2],"%)"),zlab=paste0("PC3(",z[3],"%)"),
                   #color=colors,
                   color=gdf[rownames(d),"Color"],
                   pch=gdf[rownames(d),"Shape"],
                   bg=gdf[rownames(d),"Color"],
                   type="h")
text(s3d$xyz.convert(d[,c("PC1","PC2","PC3")]),labels=d$Sample,cex= 0.6)
dev.off()


Plot3D <- function(x, y = NULL, z = NULL, color = par("col"), pch = NULL,
                   main = NULL, sub = NULL, xlim = NULL, ylim = NULL, zlim = NULL,
                   xlab = NULL, ylab = NULL, zlab = NULL, scale.y = 1, angle = 40,
                   axis = TRUE, tick.marks = TRUE, label.tick.marks = TRUE,
                   x.ticklabs = NULL, y.ticklabs = NULL, z.ticklabs = NULL,
                   y.margin.add = 0, grid = TRUE, box = TRUE, lab = par("lab"),
                   lab.z = mean(lab[1:2]), type = "p", highlight.3d = FALSE,
                   mar = c(5, 3, 4, 3) + 0.1, col.axis = par("col.axis"),
                   col.grid = "grey", col.lab = par("col.lab"), cex.symbols = par("cex"),
                   cex.axis = 0.8 * par("cex.axis"), cex.lab = par("cex.lab"),
                   font.axis = par("font.axis"), font.lab = par("font.lab"),
                   lty.axis = par("lty"), lty.grid = 2, lty.hide = 1,
                   lty.hplot = par("lty"), log = "", ...)
# log not yet implemented
{
  ## Uwe Ligges <ligges@statistik.tu-dortmund.de>,
  ## http://www.statistik.tu-dortmund.de/~ligges
  ##
  ## For MANY ideas and improvements thanks to Martin Maechler!!!
  ## Parts of the help files are stolen from the standard plotting functions in R.
  
  mem.par <- par(mar = mar)
  x.scal <- y.scal <- z.scal <- 1
  xlabel <- if (!missing(x)) deparse(substitute(x))
  ylabel <- if (!missing(y)) deparse(substitute(y))
  zlabel <- if (!missing(z)) deparse(substitute(z))
  
  ## color as part of `x' (data.frame or list):
  if(!is.null(d <- dim(x)) && (length(d) == 2) && (d[2] >= 4))
    color <- x[,4]
  else if(is.list(x) && !is.null(x$color))
    color <- x$color
  
  ## convert 'anything' -> vector
  xyz <- xyz.coords(x=x, y=y, z=z, xlab=xlabel, ylab=ylabel, zlab=zlabel,
                    log=log)
  if(is.null(xlab)) { xlab <- xyz$xlab; if(is.null(xlab)) xlab <- "" }
  if(is.null(ylab)) { ylab <- xyz$ylab; if(is.null(ylab)) ylab <- "" }
  if(is.null(zlab)) { zlab <- xyz$zlab; if(is.null(zlab)) zlab <- "" }
  
  if(length(color) == 1)
    color <- rep(color, length(xyz$x))
  else if(length(color) != length(xyz$x))
    stop("length(color) ", "must be equal length(x) or 1")
  
  angle <- (angle %% 360) / 90
  yz.f <- scale.y * abs(if(angle < 1) angle else if(angle > 3) angle - 4 else 2 - angle)
  yx.f <- scale.y * (if(angle < 2) 1 - angle else angle - 3)
  if(angle > 2) { ## switch y and x axis to ensure righthand oriented coord.
    temp <- xyz$x; xyz$x <- xyz$y; xyz$y <- temp
    temp <- xlab;  xlab <- ylab;   ylab <- temp
    temp <- xlim;  xlim <- ylim;   ylim <- temp
  }
  angle.1 <- (1 < angle && angle < 2) || angle > 3
  angle.2 <- 1 <= angle && angle <= 3
  dat <- cbind(as.data.frame(xyz[c("x","y","z")]), col = color)
  
  n <- nrow(dat);
  y.range <- range(dat$y[is.finite(dat$y)])
  
  ### 3D-highlighting / colors / sort by y
  if(type == "p" || type == "h") {
    y.ord <- rev(order(dat$y))
    dat <- dat[y.ord, ]
    if(length(pch) > 1)
      if(length(pch) != length(y.ord))
        stop("length(pch) ", "must be equal length(x) or 1")
    else pch <- pch[y.ord]
    daty <- dat$y
    daty[!is.finite(daty)] <- mean(daty[is.finite(daty)])
    if(highlight.3d && !(all(diff(daty) == 0)))
      dat$col <- rgb(seq(0, 1, length = n) * (y.range[2] - daty) / diff(y.range), g=0, b=0)
  }
  
  ### optim. axis scaling
  p.lab <- par("lab")
  ## Y
  y.range <- range(dat$y[is.finite(dat$y)], ylim)
  y.prty <- pretty(y.range, n = lab[2],
                   min.n = max(1, min(.5 * lab[2], p.lab[2])))
  y.scal <- round(diff(y.prty[1:2]), digits = 12)
  y.add <- min(y.prty)
  dat$y <- (dat$y - y.add) / y.scal
  y.max <- (max(y.prty) - y.add) / y.scal
  
  x.range <- range(dat$x[is.finite(dat$x)], xlim)
  x.prty <- pretty(x.range, n = lab[1],
                   min.n = max(1, min(.5 * lab[1], p.lab[1])))
  x.scal <- round(diff(x.prty[1:2]), digits = 12)
  dat$x <- dat$x / x.scal
  x.range <- range(x.prty) / x.scal
  x.max <- ceiling(x.range[2])
  x.min <-   floor(x.range[1])
  if(!is.null(xlim)) {
    x.max <- max(x.max, ceiling(xlim[2] / x.scal))
    x.min <- min(x.min,   floor(xlim[1] / x.scal))
  }
  x.range <- range(x.min, x.max)
  ## Z
  z.range <- range(dat$z[is.finite(dat$z)], zlim)
  z.prty <- pretty(z.range, n = lab.z,
                   min.n = max(1, min(.5 * lab.z, p.lab[2])))
  z.scal <- round(diff(z.prty[1:2]), digits = 12)
  dat$z <- dat$z / z.scal
  z.range <- range(z.prty) / z.scal
  z.max <- ceiling(z.range[2])
  z.min <-   floor(z.range[1])
  if(!is.null(zlim)) {
    z.max <- max(z.max, ceiling(zlim[2] / z.scal))
    z.min <- min(z.min,   floor(zlim[1] / z.scal))
  }
  z.range <- range(z.min, z.max)
  
  ### init graphics
  plot.new()
  if(angle.2) {x1 <- x.min + yx.f * y.max; x2 <- x.max}
  else        {x1 <- x.min; x2 <- x.max + yx.f * y.max}
  plot.window(c(x1, x2), c(z.min, z.max + yz.f * y.max))
  temp <- strwidth(format(rev(y.prty))[1], cex = cex.axis/par("cex"))
  if(angle.2) x1 <- x1 - temp - y.margin.add
  else        x2 <- x2 + temp + y.margin.add
  plot.window(c(x1, x2), c(z.min, z.max + yz.f * y.max))
  if(angle > 2) par("usr" = par("usr")[c(2, 1, 3:4)])
  usr <- par("usr") # we have to remind it for use in closures
  title(main, sub, ...)
  
  ### draw axis, tick marks, labels, grid, ...
  xx <- if(angle.2) c(x.min, x.max) else c(x.max, x.min)
  if(grid) {
    ## grids
    ###################
    # XY wall
    i <- x.min:x.max;
    segments(i, z.min, i + (yx.f * y.max), yz.f * y.max + z.min,
             col = col.grid, lty = lty.grid);
    
    i <- 0:y.max;
    segments(x.min + (i * yx.f), i * yz.f + z.min,
             x.max + (i * yx.f), i * yz.f + z.min,
             col = col.grid, lty = lty.grid);
    
    ######################
    # XZ wall
    # verticle lines
    temp <- yx.f * y.max;
    temp1 <- yz.f * y.max;
    i <- (x.min + temp):(x.max + temp);
    segments(i, z.min + temp1, i, z.max + temp1,
             col = col.grid, lty = lty.grid);
    
    # horizontal lines
    i <- (z.min + temp1):(z.max + temp1);
    segments(x.min + temp, i, x.max + temp, i,
             col = col.grid, lty = lty.grid)
    
    
    ##################
    # YZ wall
    # horizontal lines
    i <- xx[2]:x.min;
    mm <- z.min:z.max;
    segments(i, mm, i + temp, mm + temp1,
             col = col.grid, lty = lty.grid);
    # verticle lines
    i <- 0:y.max;
    segments(x.min + (i * yx.f), i * yz.f + z.min,
             xx[2] + (i * yx.f), i * yz.f + z.max,
             col = col.grid, lty = lty.grid)
    
    
    # make the axis into solid line
    segments(x.min, z.min, x.min + (yx.f * y.max), yz.f * y.max + z.min,
             col = col.grid, lty = lty.hide);
    segments(x.max, z.min, x.max + (yx.f * y.max), yz.f * y.max + z.min,
             col = col.axis, lty = lty.hide);
    segments(x.min + (y.max * yx.f), y.max * yz.f + z.min,
             x.max + (y.max* yx.f), y.max * yz.f + z.min,
             col = col.grid, lty = lty.hide);
    segments(x.min + temp, z.min + temp1, x.min + temp, z.max + temp1,
             col = col.grid, lty = lty.hide);
    segments(x.max + temp, z.min + temp1, x.max + temp, z.max + temp1,
             col = col.axis, lty = lty.hide);
    segments(x.min + temp, z.max + temp1, x.max + temp, z.max + temp1,
             col = col.axis, lty = lty.hide);
    segments(xx[2], z.max, xx[2] + temp, z.max + temp1,
             col = col.axis, lty = lty.hide);
  }
  if(axis) {
    if(tick.marks) { ## tick marks
      xtl <- (z.max - z.min) * (tcl <- -par("tcl")) / 50
      ztl <- (x.max - x.min) * tcl / 50
      mysegs <- function(x0,y0, x1,y1)
        segments(x0,y0, x1,y1, col=col.axis, lty=lty.axis)
      ## Y
      i.y <- 0:y.max
      mysegs(yx.f * i.y - ztl + xx[1], yz.f * i.y + z.min,
             yx.f * i.y + ztl + xx[1], yz.f * i.y + z.min)
      ## X
      i.x <- x.min:x.max
      mysegs(i.x, -xtl + z.min, i.x, xtl + z.min)
      ## Z
      i.z <- z.min:z.max
      mysegs(-ztl + xx[2], i.z, ztl + xx[2], i.z)
      
      if(label.tick.marks) { ## label tick marks
        las <- par("las")
        mytext <- function(labels, side, at, ...)
          mtext(text = labels, side = side, at = at, line = -.5,
                col=col.lab, cex=cex.axis, font=font.lab, ...)
        ## X
        if(is.null(x.ticklabs))
          x.ticklabs <- format(i.x * x.scal)
        mytext(x.ticklabs, side = 1, at = i.x)
        ## Z
        if(is.null(z.ticklabs))
          z.ticklabs <- format(i.z * z.scal)
        mytext(z.ticklabs, side = if(angle.1) 4 else 2, at = i.z,
               adj = if(0 < las && las < 3) 1 else NA)
        ## Y
        temp <- if(angle > 2) rev(i.y) else i.y ## turn y-labels around
        if(is.null(y.ticklabs))
          y.ticklabs <- format(y.prty)
        else if (angle > 2)
          y.ticklabs <- rev(y.ticklabs)
        text(i.y * yx.f + xx[1],
             i.y * yz.f + z.min, y.ticklabs,
             pos=if(angle.1) 2 else 4, offset=1,
             col=col.lab, cex=cex.axis/par("cex"), font=font.lab)
      }
    }
    
    ## axis and labels
    
    mytext2 <- function(lab, side, line, at)
      mtext(lab, side = side, line = line, at = at, col = col.lab,
            cex = cex.lab, font = font.axis, las = 0)
    ## X
    lines(c(x.min, x.max), c(z.min, z.min), col = col.axis, lty = lty.axis)
    mytext2(xlab, 1, line = 1.5, at = mean(x.range))
    ## Y
    lines(xx[1] + c(0, y.max * yx.f), c(z.min, y.max * yz.f + z.min),
          col = col.axis, lty = lty.axis)
    mytext2(ylab, if(angle.1) 2 else 4, line= 0.5, at = z.min + y.max * yz.f)
    
    ## Z
    lines(xx[c(2,2)], c(z.min, z.max), col = col.axis, lty = lty.axis)
    mytext2(zlab, if(angle.1) 4 else 2, line= 1.5, at = mean(z.range))
    
  }
  
  ### plot points
  x <- dat$x + (dat$y * yx.f)
  z <- dat$z + (dat$y * yz.f)
  col <- as.character(dat$col)
  if(type == "h") {
    z2 <- dat$y * yz.f + z.min
    segments(x, z, x, z2, col = col, cex = cex.symbols, lty = lty.hplot, ...)
    points(x, z, type = "p", col = col, pch = pch, cex = cex.symbols, ...)
  }
  else points(x, z, type = type, col = col, pch = pch, cex = cex.symbols, ...)
  
  ### box-lines in front of points (overlay)
  if(axis && box) {
    lines(c(x.min, x.max), c(z.max, z.max),
          col = col.axis, lty = lty.axis)
    lines(c(0, y.max * yx.f) + x.max, c(0, y.max * yz.f) + z.max,
          col = col.axis, lty = lty.axis)
    lines(xx[c(1,1)], c(z.min, z.max), col = col.axis, lty = lty.axis)
  }
  
  
  # par(mem.par) # we MUST NOT set the margins back
  ### Return Function Object
  ob <- ls() ## remove all unused objects from the result's enviroment:
  rm(list = ob[!ob %in% c("angle", "mar", "usr", "x.scal", "y.scal", "z.scal", "yx.f",
                          "yz.f", "y.add", "z.min", "z.max", "x.min", "x.max", "y.max",
                          "x.prty", "y.prty", "z.prty")])
  rm(ob)
  invisible(list(
    xyz.convert = function(x, y=NULL, z=NULL) {
      xyz <- xyz.coords(x, y, z)
      if(angle > 2) { ## switch y and x axis to ensure righthand oriented coord.
        temp <- xyz$x; xyz$x <- xyz$y; xyz$y <- temp
      }
      y <- (xyz$y - y.add) / y.scal
      return(list(x = xyz$x / x.scal + yx.f * y,
                  y = xyz$z / z.scal + yz.f * y))
    },
    points3d = function(x, y = NULL, z = NULL, type = "p", ...) {
      xyz <- xyz.coords(x, y, z)
      if(angle > 2) { ## switch y and x axis to ensure righthand oriented coord.
        temp <- xyz$x; xyz$x <- xyz$y; xyz$y <- temp
      }
      y2 <- (xyz$y - y.add) / y.scal
      x <- xyz$x / x.scal + yx.f * y2
      y <- xyz$z / z.scal + yz.f * y2
      mem.par <- par(mar = mar, usr = usr)
      on.exit(par(mem.par))
      if(type == "h") {
        y2 <- z.min + yz.f * y2
        segments(x, y, x, y2, ...)
        points(x, y, type = "p", ...)
      }
      else points(x, y, type = type, ...)
    },
    plane3d = function(Intercept, x.coef = NULL, y.coef = NULL,
                       lty = "dashed", lty.box = NULL, ...){
      if(!is.atomic(Intercept) && !is.null(coef(Intercept))) Intercept <- coef(Intercept)
      if(is.null(lty.box)) lty.box <- lty
      if(is.null(x.coef) && length(Intercept) == 3){
        x.coef <- Intercept[if(angle > 2) 3 else 2]
        y.coef <- Intercept[if(angle > 2) 2 else 3]
        Intercept <- Intercept[1]
      }
      mem.par <- par(mar = mar, usr = usr)
      on.exit(par(mem.par))
      x <- x.min:x.max
      ltya <- c(lty.box, rep(lty, length(x)-2), lty.box)
      x.coef <- x.coef * x.scal
      z1 <- (Intercept + x * x.coef + y.add * y.coef) / z.scal
      z2 <- (Intercept + x * x.coef +
               (y.max * y.scal + y.add) * y.coef) / z.scal
      segments(x, z1, x + y.max * yx.f, z2 + yz.f * y.max, lty = ltya, ...)
      y <- 0:y.max
      ltya <- c(lty.box, rep(lty, length(y)-2), lty.box)
      y.coef <- (y * y.scal + y.add) * y.coef
      z1 <- (Intercept + x.min * x.coef + y.coef) / z.scal
      z2 <- (Intercept + x.max * x.coef + y.coef) / z.scal
      segments(x.min + y * yx.f, z1 + y * yz.f,
               x.max + y * yx.f, z2 + y * yz.f, lty = ltya, ...)
    },
    
    wall3d = function(Intercept, x.coef = NULL, y.coef = NULL,
                      lty = "dashed", lty.box = NULL, ...){
      if(!is.atomic(Intercept) && !is.null(coef(Intercept))) Intercept <- coef(Intercept)
      if(is.null(lty.box)) lty.box <- lty
      if(is.null(x.coef) && length(Intercept) == 3){
        x.coef <- Intercept[if(angle > 2) 3 else 2]
        y.coef <- Intercept[if(angle > 2) 2 else 3]
        Intercept <- Intercept[1]
      }
      mem.par <- par(mar = mar, usr = usr)
      on.exit(par(mem.par))
      x <- x.min:x.max
      ltya <- c(lty.box, rep(lty, length(x)-2), lty.box)
      x.coef <- x.coef * x.scal
      z1 <- (Intercept + x * x.coef + y.add * y.coef) / z.scal
      z2 <- (Intercept + x * x.coef +
               (y.max * y.scal + y.add) * y.coef) / z.scal
      segments(x, z1, x + y.max * yx.f, z2 + yz.f * y.max, lty = ltya, ...)
      y <- 0:y.max
      ltya <- c(lty.box, rep(lty, length(y)-2), lty.box)
      y.coef <- (y * y.scal + y.add) * y.coef
      z1 <- (Intercept + x.min * x.coef + y.coef) / z.scal
      z2 <- (Intercept + x.max * x.coef + y.coef) / z.scal
      segments(x.min + y * yx.f, z1 + y * yz.f,
               x.max + y * yx.f, z2 + y * yz.f, lty = ltya, ...)
    },
    box3d = function(...){
      mem.par <- par(mar = mar, usr = usr)
      on.exit(par(mem.par))
      lines(c(x.min, x.max), c(z.max, z.max), ...)
      lines(c(0, y.max * yx.f) + x.max, c(0, y.max * yz.f) + z.max, ...)
      lines(c(0, y.max * yx.f) + x.min, c(0, y.max * yz.f) + z.max, ...)
      lines(c(x.max, x.max), c(z.min, z.max), ...)
      lines(c(x.min, x.min), c(z.min, z.max), ...)
      lines(c(x.min, x.max), c(z.min, z.min), ...)
    }
  ))
}

cls<-factor(d$Group, levels=unique(d$Group))
cols <- hyplot::GetColorSchema(cls)
legend.nm <- unique(as.character(cls))
uniq.cols <- unique(cols)
pchs <- as.numeric(cls)+1 %% 26
uniq.pchs <- unique(pchs)
xlabel <- paste("PC1", "(", z[1], "%)")
ylabel <- paste("PC2", "(", z[2], "%)")
zlabel <- paste("PC3", "(", z[3], "%)")
##################主成分分析3D图
png(paste(pn,"/PCA123_s2.png",sep=""),res=dpi,width=8,height = 8, units = units)

Plot3D(d$PC1, d$PC2, d$PC3, xlab= xlabel, ylab=ylabel,
       zlab=zlabel, angle =60, color=gdf[rownames(d),"Color"], pch=gdf[rownames(d),"Shape"], box=F)
#legend("topleft", legend = legend.nm, pch=uniq.pchs %% 26, col=uniq.cols)

if(group_numb < 6){
  legend("topleft", legend = legend.nm, pch=shape.man,col=color.man);
}else if (group_numb < 10){
  legend("topleft", legend = legend.nm, pch=shape.man, col=color.man, cex=0.75);
}else if (group_numb < 20){
  legend("topleft",legend = legend.nm, pch=shape.man, col=color.man, cex=0.5);
}else{
  legend("topleft",legend = legend.nm, pch=shape.man, col=color.man, cex=0.5,
         xjust=1, yjust=0, xpd=T, ncol=5)
}

dev.off()

