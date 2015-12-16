## This script generates a per-stock distribution of SCO, DCO, TCO, and NCO events.

## Author: Danny Miller danrdanny <shift-2> gmail <dot> com

## Note that one TCO is depicted as a DCO but is red
## The output of this script is made much prettier in Illustrator by hand

## Release 6
coData <- read.csv("~/simr/2016_CO_NCO_Paper/CO_NCO_data/co-ncoDetail.r6.tsv",sep = "\t", header = T)
# Format of file:
#Stock  Chr	Pos	Gap	Class	Count
#w3.1	chrX	2432603	58	co	0
#cs2.1	chrX	2495362	694	co	1
#cs13.1	chrX	3470244	4177	dco	0
#cs15.1	chrX	3905685	841	nco	0

# Run the graphs once for wStocks and once for csStocks
stocks <- c("w3.1","w3.4","w3.5","w3.6","w3.8","w3.9","w3.11","w3.12","w3.13","w3.14","w3.15","w3.16","w3.17","w3.18","w3.21","w3.24","w3.25","w3.26","w4.1","w4.2","w4.3","w4.4","w4.5","w4.6","w4.7","w4.8","w4.9","w4.11","w4.12","w4.13","w4.15","w4.16","w4.17","w5.1","w5.2","w5.3","w5.4","w5.5","w5.6","w5.7","w5.8","w5.9","w5.11","w5.12","w5.13","w5.14","w5.15","w5.16","w5.17","w11.1","w11.2","w11.3","w11.4","w11.5","w11.6","w11.7","w11.8","w11.10","w11.11","w11.12","w11.13","w12.1","w12.2","w12.3","w12.4","w12.5","w12.6","w12.7","w12.9","w12.10","w12.11","w12.12","w12.13","w13.1","w13.2","w13.4","w13.5","w13.6","w13.7","w13.8","w13.9","w13.10","w13.11","w13.12","w13.13","w15.1","w15.2","w15.3","w15.4","w15.5","w15.6","w15.7","w15.8","w15.9","w15.10","w15.11","w15.12","w15.13")
stocks <- c("cs1.1","cs1.2","cs1.3","cs1.4","cs1.5","cs1.6","cs1.7","cs1.8","cs2.1","cs2.2","cs2.3","cs2.4","cs2.5","cs2.6","cs2.7","cs2.8","cs5.1","cs5.2","cs5.3","cs5.4","cs5.5","cs5.6","cs5.7","cs5.8","cs6.1","cs6.2","cs6.3","cs6.4","cs6.5","cs6.6","cs6.7","cs6.8","cs7.1","cs7.2","cs7.3","cs7.4","cs7.5","cs7.6","cs7.7","cs7.8","cs8.1","cs8.2","cs8.3","cs8.4","cs8.5","cs8.6","cs8.7","cs8.8","cs9.1","cs12.1","cs12.2","cs12.3","cs12.4","cs12.5","cs12.6","cs12.7","cs12.8","cs12.9","cs12.10","cs12.11","cs12.12","cs12.13","cs12.14","cs12.15","cs12.16","cs12.17","cs12.18","cs13.1","cs13.2","cs13.3","cs13.4","cs13.5","cs13.6","cs13.7","cs13.8","cs13.9","cs13.10","cs13.11","cs13.12","cs13.13","cs13.14","cs13.15","cs13.16","cs13.17","cs13.18","cs14.1","cs14.2","cs14.3","cs14.4","cs14.5","cs14.6","cs14.7","cs14.8","cs14.9","cs14.10","cs14.11","cs14.12","cs14.13")

chr <- c("chrX","chr2L","chr2R","chr3L","chr3R")
size <- c(23542271,23513712,25286936,28110227,32079331)
chrSizes <- data.frame(chr,size)

# Assign values to yPosition based on Class column
coData$yPosition <- 0;
coData$yPosition[coData$Class=="sco"]<-1
coData$yPosition[coData$Class=="dco"]<-2
coData$yPosition[coData$Class=="nco"]<-3

stockData<-data.frame(Stock=character(),Chr=character(),AvePosition=integer(),Class=character(),Dad=character(),Batch=integer(),yPosition=integer(),stringsAsFactors=FALSE)

for (currStock in stocks) {
  print(currStock)
  
  stockData <- rbind(stockData,filter(coData,Stock==currStock))
  #stockData <- filter(coData,Stock==currStock)
  
  newRow1 <- data.frame(Stock=currStock,Chr="chrX",AvePosition=23542271,Class="fake",Dad="na",Batch=0,yPosition=0)
  newRow2 <- data.frame(Stock=currStock,Chr="chr2L",AvePosition=23513712,Class="fake",Dad="na",Batch=0,yPosition=0)
  newRow3 <- data.frame(Stock=currStock,Chr="chr2R",AvePosition=25286936,Class="fake",Dad="na",Batch=0,yPosition=0)
  newRow4 <- data.frame(Stock=currStock,Chr="chr3L",AvePosition=28110227,Class="fake",Dad="na",Batch=0,yPosition=0)
  newRow5 <- data.frame(Stock=currStock,Chr="chr3R",AvePosition=32079331,Class="fake",Dad="na",Batch=0,yPosition=0)
  
  stockData <- rbind(stockData,newRow1,newRow2,newRow3,newRow4,newRow5)
  stockData$yPosition<-1
}
stockData$Chr <- factor(stockData$Chr, levels = c("chrX", "chr2L", "chr2R", "chr3L", "chr3R"))

stepSize <- 1000000
a <- ggplot(data = stockData, aes(x=AvePosition, y=yPosition,color=factor(Class),shape=factor(Class)))
a <- a + geom_point(size = 4)
a <- a + scale_color_manual(values = c("tco"="red1", "dco"="cadetblue2", "sco"="blue1", "nco"="darkgoldenrod3","fake"="gray94")) 
a <- a + scale_shape_manual(values = c("tco"=15,"dco"=15,"sco"=19,"nco"=124,"fake"=20))  
a <- a + theme_bw()
a <- a + scale_x_continuous(expand = c(0, 0)) 
a <- a + scale_y_continuous(expand = c(0, 0))
a <- a + theme(panel.background = element_rect(fill='gray94'),
               panel.border = element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.grid.major.x = element_blank(),
               panel.grid.minor.x = element_blank(),
               legend.position="none", 
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank(),
               axis.title.y = element_blank(),
               axis.title.x = element_blank(), 
               
               axis.text.x=element_text(size=6,vjust=1)
)
a <- a + facet_grid(Stock ~ Chr)
a