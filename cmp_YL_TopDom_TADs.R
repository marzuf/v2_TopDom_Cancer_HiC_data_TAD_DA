require(foreach)
require(doMC)
registerDoMC(4)
chromo <- "chr1"

all_chromos <- paste0("chr", c(1:22))

binSizeKb <- 40

cell_line <- "ENCSR489OCU_NCI-H460"


all_TopDom_dt <- foreach(chromo = all_chromos, .combine='rbind') %dopar% {
  inFile <- file.path(paste0(cell_line, "_", binSizeKb, "kb"), "FINAL_DOMAINS", paste0(cell_line, "_", chromo, "_TopDom_", binSizeKb, "kb_final_domains.txt"))
  stopifnot(file.exists(inFile))
  in_dt <- read.delim(inFile, stringsAsFactors = FALSE, header=FALSE, col.names=c("chromo", "start", "end"))
  in_dt
}
all_TopDom_dt$caller <- "TopDom"


all_YL_dt <- foreach(chromo = all_chromos, .combine='rbind') %dopar% {
  inFile <- file.path("../v2_Yuanlong_Cancer_HiC_data_TAD_DA", paste0(cell_line, "_", binSizeKb, "kb"), "FINAL_DOMAINS", paste0(cell_line, "_", chromo, "_YL_", binSizeKb, "kb_final_domains.txt"))
  stopifnot(file.exists(inFile))
  in_dt <- read.delim(inFile, stringsAsFactors = FALSE, header=FALSE, col.names=c("chromo", "start", "end"))
  in_dt
}
all_YL_dt$caller <- "YL"

aggByChr_TopDom <- aggregate(start~chromo+caller, data=all_TopDom_dt, FUN=length)
colnames(aggByChr_TopDom)[3] <- "nbrTADs"

aggByChr_YL <- aggregate(start~chromo+caller, data=all_YL_dt, FUN=length)
colnames(aggByChr_YL)[3] <- "nbrTADs"

all_aggByChr_dt <- rbind(aggByChr_TopDom, aggByChr_YL)



all_dt <- rbind(all_TopDom_dt, all_YL_dt)
all_dt$domainSize <- all_dt$end - all_dt$start + 1
all_dt$domainSize_log10 <- log10(all_dt$domainSize)


my_cols <- c(pal_jama()(5)[c(3, 2)])

plotTit <- paste0(cell_line, " - all chromos")
mySub <- paste0(binSizeKb, " kb HiC data")
fontFamily <- "Hershey"
legTitle <- ""

custom_p <- function(x) {
  x2 <- x +   
    guides(color=FALSE)+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    theme(
      text = element_text(family=fontFamily),
      panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
      panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
      panel.background = element_rect(fill = "transparent"),
      panel.grid.major.x =  element_blank(),
      panel.grid.minor.x =  element_blank(),
      axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
      axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
      axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
      axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5),
      plot.title = element_text(hjust=0.5, size = 16, face="bold"),
      plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
      legend.title = element_text(face="bold")
    ) 
  
}

pNbr <- ggdensity(all_aggByChr_dt,
                  title="# TADs/chromosome",
                 x = "nbrTADs",
                 y = "..density..",
                 # combine = TRUE,                  # Combine the 3 plots
                 xlab = "# TADs by chromo", 
                 # add = "median",                  # Add median line. 
                 rug = FALSE,                      # Add marginal rug
                 color = "caller", 
                 fill = "caller",
                 palette = "jco"
) +
  scale_color_manual(values=my_cols)+
  ggtitle(plotTit, subtitle = mySub)+
  scale_fill_manual(values=my_cols)  +
  labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density")  

pNbr <- custom_p(pNbr)
pNbr
  

pSize <- ggdensity(all_dt,
                   title="Domain size [bp, log10]",
                  x = "domainSize_log10",
                  y = "..density..",
                  # combine = TRUE,                  # Combine the 3 plots
                  xlab = "# TADs by chromo", 
                  # add = "median",                  # Add median line. 
                  rug = FALSE,                      # Add marginal rug
                  color = "caller", 
                  fill = "caller",
                  palette = "jco"
) +
  scale_color_manual(values=my_cols)+
  ggtitle(plotTit, subtitle = mySub)+
  scale_fill_manual(values=my_cols)  +
  labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density")  

pSize <- custom_p(pSize)
pSize





outFile <- file.path(outFolder, paste0("signif_all_meanCorr_dist_density.", plotType))
ggsave(p2b, file=outFile, height=myHeight, width=myWidth)
cat(paste0("... written: ", outFile, "\n"))

