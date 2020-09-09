require(foreach)
require(doMC)
registerDoMC(40)
require(ggpubr)
require(ggplot2)
require(ggsci)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

# Rscript cmp_YL_TopDom_TADs.R

outFolder <- "CMP_YL_TOPDOM_TADS"
dir.create(outFolder, recursive = TRUE)
myHeight <- 5
myWidth <- 7
plotType <- "svg"


plotTypeB <- "png"
myHeightB <- myWidthB <- 400

chromo <- "chr1"

all_chromos <- paste0("chr", c(1:22))

binSizeKb <- 40

cell_line <- "ENCSR489OCU_NCI-H460"

buildTable1 <- F
buildTable2 <- F


all_cell_lines <- c(
  "ENCSR489OCU_NCI-H460",
  "ENCSR444WCZ_A549",
  "LG1",
  "Barutcu_MCF-10A",
  "Barutcu_MCF-7",
  "GM12878",
  "K562",
  "Rao_HCT-116_2017",
  "LG2",
  "LI",
  "PA2",
  "PA3",
  "GSE109229_BT474",
  "GSE118514_22Rv1",
  "GSE118588_Panc_beta",
  "ENCSR079VIJ_G401",
  "ENCSR312KHQ_SK-MEL-5",
  "ENCSR346DCU_LNCaP",
  "ENCSR401TBQ_Caki2",
  "ENCSR444WCZ_A549",
  "ENCSR489OCU_NCI-H460",
  "ENCSR549MGQ_T47D",
  "ENCSR862OGI_RPMI-7951",
  "GSE105194_cerebellum",
  "GSE105194_spinal_cord",
  "GSE105318_DLD1",
  "GSE105381_HepG2",
  "GSE118514_RWPE1",
  "Panc1_rep12",
  "ENCSR504OTV_transverse_colon",
  "GSE109229_SKBR3",
  "GSE99051_786_O",
  "HMEC"
)

cell_line = all_cell_lines[1]

if(buildTable1) {
  
  all_TopDom_dt <- foreach(cell_line = all_cell_lines, .combine='rbind') %dopar% {
    chr_dt <- read.delim(file.path(paste0(cell_line, "_", binSizeKb, "kb"), "genes2tad", paste0("all_assigned_regions.txt")), 
                         col.names=c("chromo", "region", "start", "end"),
                         stringsAsFactors = FALSE, header=F)
    stopifnot(chr_dt$end > chr_dt$start)
    
    
    chr_dt <- chr_dt[grepl("_TAD", chr_dt$region),]
    # chr_dt <- foreach(chromo = all_chromos, .combine='rbind') %dopar% {
    #   inFile <- file.path(paste0(cell_line, "_", binSizeKb, "kb"), "FINAL_DOMAINS", paste0(cell_line, "_", chromo, "_TopDom_", binSizeKb, "kb_final_domains.txt"))
    #   if(!file.exists(inFile))cat(paste0(inFile, "\n"))
    #   stopifnot(file.exists(inFile))
    #   in_dt <- read.delim(inFile, stringsAsFactors = FALSE, header=FALSE, col.names=c("chromo", "start", "end"))
    #   in_dt
    # }  
    chr_dt$cell_line <- cell_line
    chr_dt
  }
  all_TopDom_dt$caller <- "TopDom"
  
  all_YL_dt <- foreach(cell_line = all_cell_lines, .combine='rbind') %dopar% {
  # all_YL_dt <- foreach(chromo = all_chromos, .combine='rbind') %dopar% {
  #   inFile <- file.path("../v2_Yuanlong_Cancer_HiC_data_TAD_DA", paste0(cell_line, "_", binSizeKb, "kb"), "FINAL_DOMAINS", paste0(cell_line, "_", chromo, "_YL_", binSizeKb, "kb_final_domains.txt"))
  #   stopifnot(file.exists(inFile))
  #   
    
    if(cell_line == "GM12878") {
      in_dt <- read.delim(file.path("../v2_Yuanlong_Cancer_HiC_data_TAD_DA_GM12878/", paste0(cell_line, "_", binSizeKb, "kb"), "genes2tad", "all_assigned_regions.txt"),
                          col.names=c("chromo", "region", "start", "end"),
                          stringsAsFactors = FALSE, header=F)
      
    } else {
      in_dt <- read.delim(file.path("../v2_Yuanlong_Cancer_HiC_data_TAD_DA", paste0(cell_line, "_", binSizeKb, "kb"), "genes2tad", "all_assigned_regions.txt"),
                          col.names=c("chromo", "region", "start", "end"),
                          stringsAsFactors = FALSE, header=F)
      
    }
    stopifnot(in_dt$end > in_dt$start)
    
    
    in_dt <- in_dt[grepl("_TAD", in_dt$region),]
    # in_dt <- read.delim(inFile, stringsAsFactors = FALSE, header=FALSE, col.names=c("chromo", "start", "end"))
    in_dt$cell_line <- cell_line
    in_dt
  }
  all_YL_dt$caller <- "YL"
  
  
  outFile <- file.path(outFolder, "all_TopDom_dt.Rdata")
  save(all_TopDom_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
  outFile <- file.path(outFolder, "all_YL_dt.Rdata")
  save(all_YL_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
} else {
  all_TopDom_dt <- get(load(file.path(outFolder, "all_TopDom_dt.Rdata")))
  all_YL_dt <- get(load(file.path(outFolder, "all_YL_dt.Rdata")))
}

aggByChr_TopDom <- aggregate(start~chromo+caller+cell_line, data=all_TopDom_dt, FUN=length)
colnames(aggByChr_TopDom)[4] <- "nbrTADs"

aggByChr_YL <- aggregate(start~chromo+caller+cell_line, data=all_YL_dt, FUN=length)
colnames(aggByChr_YL)[4] <- "nbrTADs"

all_aggByChr_dt <- rbind(aggByChr_TopDom, aggByChr_YL)


all_dt <- rbind(all_TopDom_dt, all_YL_dt)
all_dt$domainSize <- all_dt$end - all_dt$start + 1
all_dt$domainSize_log10 <- log10(all_dt$domainSize)


my_cols <- c(pal_jama()(5)[c(3, 2)])

plotTit <- paste0(cell_line, " - all chromos")
mySub <- paste0(binSizeKb, " kb HiC data - all datasets (n=", length(unique(all_dt$cell_line)),")")
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
  return(x2)
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
# pNbr
outFile <- file.path(outFolder, paste0("signif_cmp_YL_TopDom_domainNbr_density.", plotType))
ggsave(pNbr, file=outFile, height=myHeight, width=myWidth)
cat(paste0("... written: ", outFile, "\n"))



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
# pSize


outFile <- file.path(outFolder, paste0("cmp_YL_TopDom_domainSize_density.", plotType))
ggsave(pSize, file=outFile, height=myHeight, width=myWidth)
cat(paste0("... written: ", outFile, "\n"))

all_TopDom_dt$tad_size_log10 <- log10(all_TopDom_dt$end - all_TopDom_dt$start + 1)
aggMeanSize_TopDom_dt <- aggregate(tad_size_log10~cell_line + caller + chromo, FUN=mean, data = all_TopDom_dt)
all_YL_dt$tad_size_log10 <- log10(all_YL_dt$end - all_YL_dt$start + 1)
aggMeanSize_YL_dt <- aggregate(tad_size_log10~cell_line + caller + chromo, FUN=mean, data = all_YL_dt)
aggMeanSize_dt <- merge(aggMeanSize_TopDom_dt, aggMeanSize_YL_dt, by=c("cell_line", "chromo"), suffixes = c("_TopDom", "_YL"))

my_x <- aggMeanSize_dt[,"tad_size_log10_YL"]
my_y <- aggMeanSize_dt[,"tad_size_log10_TopDom"]

outFile <- file.path(outFolder, paste0("tad_size_mean_dotplot.", plotTypeB))
do.call(plotTypeB, list(outFile, height=myHeightB, width=myWidthB))
plot(
  x = my_x,
  y = my_y,
  xlab = "tad_size_log10_YL",
  ylab = "tad_size_log10_TopDom",
  pch=16,
  main="mean TAD size [log10] by chromo by cell line"
)
mtext(side=3, text = paste0("all DS  - n = ", length(unique(aggMeanSize_dt$cell_line)), " - (tot = ", nrow(aggMeanSize_dt), ")"))
addCorr(x=my_x,y=my_y, bty="n", legPos = "topright")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

        
aggNbr_TopDom_dt <- aggregate(tad_size_log10~cell_line + caller + chromo, FUN=length, data = all_TopDom_dt)
aggNbr_YL_dt <- aggregate(tad_size_log10~cell_line + caller + chromo, FUN=length, data = all_YL_dt)
colnames(aggNbr_TopDom_dt)[4] <- "tad_nbr"
colnames(aggNbr_YL_dt)[4] <- "tad_nbr"
aggNbr_TopDom_dt$tad_nbr_log10 <- log10(aggNbr_TopDom_dt$tad_nbr)
aggNbr_YL_dt$tad_nbr_log10 <- log10(aggNbr_YL_dt$tad_nbr)

aggNbr_dt <- merge(aggNbr_TopDom_dt, aggNbr_YL_dt, by=c("cell_line", "chromo"), suffixes = c("_TopDom", "_YL"))

my_x <- aggNbr_dt[,"tad_nbr_log10_YL"]
my_y <- aggNbr_dt[,"tad_nbr_log10_TopDom"]

outFile <- file.path(outFolder, paste0("tad_nbr_dotplot.", plotTypeB))
do.call(plotTypeB, list(outFile, height=myHeightB, width=myWidthB))

plot(
  x = my_x,
  y = my_y,
  xlab = "tad_nbr_log10_YL",
  ylab = "tad_nbr_log10_TopDom",
  pch=16,
  main="# TADs [log10] by chromo by cell line"
)
mtext(side=3, text = paste0("all DS  - n = ", length(unique(aggNbr_dt$cell_line)), " - (tot = ", nrow(aggNbr_dt), ")"))
addCorr(x=my_x,y=my_y, bty="n", legPos = "topleft")

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



if(buildTable2) {
  source("get_moc.R")
  
  all_MoC_dt <- foreach(cell_line = all_cell_lines, .combine='rbind') %dopar% {
    
    
    
    if(cell_line == "GM12878") {
      yl_all_dt <- read.delim(file.path("../v2_Yuanlong_Cancer_HiC_data_TAD_DA_GM12878/", paste0(cell_line, "_", binSizeKb, "kb"), "genes2tad", "all_assigned_regions.txt"),
                          col.names=c("chromo", "region", "start", "end"),
                          stringsAsFactors = FALSE, header=F)
      
    } else {
      yl_all_dt <- read.delim(file.path("../v2_Yuanlong_Cancer_HiC_data_TAD_DA", paste0(cell_line, "_", binSizeKb, "kb"), "genes2tad", "all_assigned_regions.txt"),
                          col.names=c("chromo", "region", "start", "end"),
                          stringsAsFactors = FALSE, header=F)
      
    }
    yl_all_dt <- yl_all_dt[grepl("_TAD", yl_all_dt$region),]
    
    
    td_all_dt <- read.delim(file.path(paste0(cell_line, "_", binSizeKb, "kb"), "genes2tad", paste0("all_assigned_regions.txt")), 
                         col.names=c("chromo", "region", "start", "end"),
                         stringsAsFactors = FALSE, header=F)
    td_all_dt <- td_all_dt[grepl("_TAD", td_all_dt$region),]
    
    
    
    chr_dt <- foreach(chromo = all_chromos, .combine='rbind') %do% {
      # tdFile <- file.path(paste0(cell_line, "_", binSizeKb, "kb"), "FINAL_DOMAINS", paste0(cell_line, "_", chromo, "_TopDom_", binSizeKb, "kb_final_domains.txt"))
      # stopifnot(file.exists(tdFile))
      # ylFile <- file.path("../v2_Yuanlong_Cancer_HiC_data_TAD_DA", paste0(cell_line, "_", binSizeKb, "kb"), "FINAL_DOMAINS", paste0(cell_line, "_", chromo, "_YL_", binSizeKb, "kb_final_domains.txt"))
      # stopifnot(file.exists(ylFile))
      
      td_dt <- td_all_dt[td_all_dt$chromo == chromo, c("chromo", "start", "end")]
      stopifnot(nrow(td_dt) > 0)
      stopifnot(length(unique(td_dt$chromo)) == 1)
      
      yl_dt <- yl_all_dt[yl_all_dt$chromo == chromo, c("chromo", "start", "end")]
      stopifnot(nrow(yl_dt) > 0)
      stopifnot(length(unique(yl_dt$chromo)) == 1)
      
      chr_moc <- get_MoC_dt(dt1=td_dt, dt2=yl_dt, chrSize=NULL)
      out_dt <- data.frame(
        chromo = chromo,
        moc = chr_moc,
        stringsAsFactors = FALSE
      )
      rownames(out_dt)<-NULL
      out_dt
    }
    chr_dt$cell_line = cell_line
    chr_dt
  }
  outFile <- file.path(outFolder, "all_MoC_dt.Rdata")
  save(all_MoC_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  all_MoC_dt <- get(load(file.path(outFolder, "all_MoC_dt.Rdata")))
}

plotType <- "svg"
myHeight <- 5
myWidth <- 5
plotCex <- 1.2
jitterCol <-  "blue"

agg_MoC <- aggregate(moc~cell_line, FUN=mean, data=all_MoC_dt)
agg_MoC <- agg_MoC[order(agg_MoC$moc, decreasing = TRUE),]

all_MoC_dt$chromo <- factor(all_MoC_dt$chromo, levels=all_chromos)

all_MoC_dt$cell_line <- factor(all_MoC_dt$cell_line, levels=agg_MoC$cell_line)


outFile <- file.path(outFolder, paste0("YL_TopDom_MoC_byChromo_boxplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
# par(mar = par()$mar + c(10,3,0,0))
boxplot(moc ~ chromo, outline=FALSE,
        data = all_MoC_dt, main = "YL - TopDom MoC by chromo", 
        xlab="", ylab="MoC", cex.main=plotCex, cex.lab=plotCex, cex.axis=plotCex, las=2)
stripchart(moc~chromo, vertical = TRUE, data = all_MoC_dt,
           method = "jitter", add = TRUE, pch = 20, col = jitterCol)
# legend("topright", legend=legText, bty="n", cex=0.9)
mtext(side=2, text="MoC", cex=plotCex, line=5)
mtext(side=3, text = paste0("all DS (# = ", length(unique(all_MoC_dt$cell_line)), ")"))
foo <- dev.off()

cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, paste0("YL_TopDom_MoC_byDS_boxplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
# par(mar = par()$mar + c(10,3,0,0))
boxplot(moc ~ cell_line, outline=FALSE,
        data = all_MoC_dt, main = "YL - TopDom MoC by DS", 
        xlab="", ylab="MoC", cex.main=plotCex, cex.lab=plotCex, cex.axis=plotCex, las=2)
stripchart(moc~cell_line, vertical = TRUE, data = all_MoC_dt,
           method = "jitter", add = TRUE, pch = 20, col = jitterCol)
# legend("topright", legend=legText, bty="n", cex=0.9)
mtext(side=2, text="MoC", cex=plotCex, line=5)
mtext(side=3, text = paste0("all chromo (# = ", length(unique(all_MoC_dt$chromo)), ")"))
foo <- dev.off()

cat(paste0("... written: ", outFile, "\n"))