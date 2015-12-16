# R script to place the raw data for figure 2

require(ggplot2)
require(grid)
require(gridExtra)
require(dplyr)

fnPlotAllCOs <- function(allCOData,dcoData,chromosome) {
  stepSize <- 1000000
  
  plot <- ggplot(allCOData,aes(x=AvePosition,y=yPosition,xmin=0,ymax=4,color=factor(Class),shape=factor(Class))) +
    geom_point(position="dodge",size=3) +
    scale_color_manual(values = c("tco"="red", "dco"="cadetblue2", "sco"="blue1", "nco"="darkgoldenrod3","fake"="gray94")) +
    scale_shape_manual(values = c("tco"=124, "dco"=124,"sco"=124,"nco"=124,"fake"=19)) +
    theme_bw() +
    expand_limits(x = 0, y = -10) +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(.01, .01)) +
    theme(panel.background = element_blank(), 
          panel.border = element_blank(),
          legend.position="none", 
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_blank(), 
          plot.margin=unit(c(0,0,0,0),"mm"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.background = element_blank(),
          panel.border = element_blank(),
          axis.text.x=element_text(size=6,vjust=1)
    ) +
    scale_x_discrete (
      breaks=(round(seq(0, max(allCOData$AvePosition), by = stepSize),1)), 
      labels=round((seq(0, max(allCOData$AvePosition), by = stepSize) / stepSize),0)
    ) 
  
  currDcoData <- subset(dcoData,Chr == chromosome)
  
  yPos = -1
  for(i in 1:nrow(currDcoData)) {
    plot <- plot + annotate("segment", x=currDcoData[i,]$fiveprime, xend=currDcoData[i,]$threeprime, y=yPos, yend=yPos, colour="blue", size=.75)
    yPos <- yPos - 1
  }  
  return(plot)
}

## Release 6
y<-annotation_custom(grobTree(textGrob("y",x=.015,y=.05,gp=gpar(fontsize=10))))  #356,509..361,245
w<-annotation_custom(grobTree(textGrob("w",x=.118,y=.05,gp=gpar(fontsize=10)))) #2,790,599..2,796,466 
cv<-annotation_custom(grobTree(textGrob("cv",x=.241,y=.05,gp=gpar(fontsize=10))))#5,690,002..5,693,084
v<-annotation_custom(grobTree(textGrob("v",x=.464,y=.05,gp=gpar(fontsize=10)))) #10,923,972..10,925,631
g<-annotation_custom(grobTree(textGrob("g",x=.583,y=.05,gp=gpar(fontsize=10))))  #13,727,204..13,736,278
f<-annotation_custom(grobTree(textGrob("f",x=.733,y=.05,gp=gpar(fontsize=10))))  #17,232,942..17,280,964
B<-annotation_custom(grobTree(textGrob("B",x=.750,y=.05,gp=gpar(fontsize=10)))) #17334,529-17,537,993
car<-annotation_custom(grobTree(textGrob("car",x=.831,y=.05,gp=gpar(fontsize=10)))) #

al<-annotation_custom(grobTree(textGrob("al",x=0.016,y=.05,gp=gpar(fontsize=10))))
dp<-annotation_custom(grobTree(textGrob("dp",x=0.193,y=.05,gp=gpar(fontsize=10))))
cl<-annotation_custom(grobTree(textGrob("cl",x=0.235,y=.05,gp=gpar(fontsize=10))))
b<-annotation_custom(grobTree(textGrob("b",x=0.588,y=.05,gp=gpar(fontsize=10))))
pr<-annotation_custom(grobTree(textGrob("pr",x=0.854,y=.05,gp=gpar(fontsize=10))))
lt<-annotation_custom(grobTree(textGrob("lt",x=0.975,y=.05,gp=gpar(fontsize=10))))

cn<-annotation_custom(grobTree(textGrob("cn",x=0.308,y=.05,gp=gpar(fontsize=10))))
sca<-annotation_custom(grobTree(textGrob("sca",x=0.506,y=.05,gp=gpar(fontsize=10))))
px<-annotation_custom(grobTree(textGrob("px",x=0.891,y=.05,gp=gpar(fontsize=10))))
bw<-annotation_custom(grobTree(textGrob("bw",x=0.931,y=.05,gp=gpar(fontsize=10))))

ru<-annotation_custom(grobTree(textGrob("ru",x=0.049,y=.05,gp=gpar(fontsize=10))))
h<-annotation_custom(grobTree(textGrob("h",x=0.309,y=.05,gp=gpar(fontsize=10))))
D<-annotation_custom(grobTree(textGrob("D",x=0.504,y=.05,gp=gpar(fontsize=10))))
th<-annotation_custom(grobTree(textGrob("th",x=0.571,y=.05,gp=gpar(fontsize=10))))
st<-annotation_custom(grobTree(textGrob("st",x=0.587,y=.05,gp=gpar(fontsize=10))))

p<-annotation_custom(grobTree(textGrob("p",x=0.27,y=.05,gp=gpar(fontsize=10))))
cu<-annotation_custom(grobTree(textGrob("cu",x=0.349,y=.05,gp=gpar(fontsize=10))))
ss<-annotation_custom(grobTree(textGrob("ss",x=0.511,y=.05,gp=gpar(fontsize=10))))
sr<-annotation_custom(grobTree(textGrob("sr",x=0.565,y=.05,gp=gpar(fontsize=10))))
e<-annotation_custom(grobTree(textGrob("e",x=0.662,y=.05,gp=gpar(fontsize=10))))
Pr<-annotation_custom(grobTree(textGrob("Pr",x=0.813,y=.05,gp=gpar(fontsize=10))))
ca<-annotation_custom(grobTree(textGrob("ca",x=0.929,y=.05,gp=gpar(fontsize=10))))

coData <- read.csv("~/simr/2016_CO_NCO_Paper/CO_NCO_data/co-ncoDetail.r6.tsv",sep = "\t", header = T)
dcoData <- read.csv("~/simr/2016_CO_NCO_Paper/CO_NCO_data//dcoData.r6.tsv",sep = "\t", header = T)
chr <- c("chrX","chr2L","chr2R","chr3L","chr3R")
size <- c(23542271,23513712,25286936,28110227,32079331)
chrSizes <- data.frame(chr,size)

# Subset by batch or dad
coData <- subset(coData,Batch == 1)
dcoData <- subset(dcoData,Batch == 1)

coData <- subset(coData,Batch == 2)
dcoData <- subset(dcoData,Batch == 2)

# Assign values to yPosition based on Class column
coData$yPosition <- 0;
coData$yPosition[coData$Class=="sco"]<-1
coData$yPosition[coData$Class=="dco"]<-2
coData$yPosition[coData$Class=="nco"]<-3

newRow1 <- data.frame(Stock="Z1",Chr="chrX",AvePosition=23542271,Class="fake",yPosition=0,Dad="na",Batch=0)
newRow2 <- data.frame(Stock="Z1",Chr="chr2L",AvePosition=23513712,Class="fake",yPosition=0,Dad="na",Batch=0)
newRow3 <- data.frame(Stock="Z1",Chr="chr2R",AvePosition=25286936,Class="fake",yPosition=0,Dad="na",Batch=0)
newRow4 <- data.frame(Stock="Z1",Chr="chr3L",AvePosition=28110227,Class="fake",yPosition=0,Dad="na",Batch=0)
newRow5 <- data.frame(Stock="Z1",Chr="chr3R",AvePosition=32079331,Class="fake",yPosition=0,Dad="na",Batch=0)

coData <- rbind(coData,newRow1,newRow2,newRow3,newRow4,newRow5)

chrX <- fnPlotAllCOs(filter(coData,Chr=="chrX"),dcoData,"chrX")
chrX <- chrX + y + w + cv + v + g + f + B + car
chr2L <- fnPlotAllCOs(filter(coData,Chr=="chr2L"),dcoData,"chr2L")
chr2L <- chr2L + al + dp + cl + b + pr + lt
chr2R <- fnPlotAllCOs(filter(coData,Chr=="chr2R"),dcoData,"chr2R")
chr2R <- chr2R + cn + sca + px + bw
chr3L <- fnPlotAllCOs(filter(coData,Chr=="chr3L"),dcoData,"chr3L")
chr3L <- chr3L + ru + h + D + th + st
chr3R <- fnPlotAllCOs(filter(coData,Chr=="chr3R"),dcoData,"chr3R")
chr3R <- chr3R + p + cu + ss + sr + e + Pr + ca

# This raw plot is then made prettier in Illustrator
grid.arrange(chrX,chr2L,chr2R,chr3L,chr3R,ncol=1)