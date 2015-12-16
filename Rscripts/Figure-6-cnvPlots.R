require(cowplot) # for doing multi-graph plots

fnPlotlog2Depth <- function(data) {
  stepSize <- 1000000

  plot <- ggplot() + 
    geom_point(data=data, map=aes(x=Pos, y=log2depth), size=1) +
    theme_bw() +
    scale_x_discrete(breaks=(round(seq(min(data$Pos), max(data$Pos), by = stepSize),1)), labels=round((seq(min(data$Pos), max(data$Pos), by = stepSize) / stepSize),0)) +
    scale_y_continuous( limits = c(-1.5,1.5), expand = c(0,0) ) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  
  return(plot)
}

chrList <- c("chrX","chr2L","chr2R","chr3L","chr3R","chr4")

# Shared CNVs
stock <- "cs14.10"
stock <- "cs13.10"
stock <- "cs12.2"
stock <- "w3.16"

# De novo stocks
stock <- "cs13.12"
stock <- "cs14.5"
stock <- "cs7.5"
stock <- "cs12.3"

# To generate plots run each stock by hand
fullFile <- paste("~/simr/2016_CO_NCO_Paper/CO_NCO_data/cnvData/",stock,".cnv", sep="")
data <- read.csv(fullFile,sep = "\t", header = T)

chrXplot <- fnPlotlog2Depth(subset(data,Chr=="chrX"))
chr2Lplot <- fnPlotlog2Depth(subset(data,Chr=="chr2L"))
chr2Rplot <- fnPlotlog2Depth(subset(data,Chr=="chr2R"))
chr3Lplot <- fnPlotlog2Depth(subset(data,Chr=="chr3L"))
chr3Rplot <- fnPlotlog2Depth(subset(data,Chr=="chr3R"))

print(plot_grid(chrXplot, chr2Lplot, chr2Rplot, chr3Lplot, chr3Rplot, rel_heights = c(1,1,1,1,1), labels = c("chrX","chr2L","chr2R","chr3L","chr3R"), nrow=5, ncol=1))
