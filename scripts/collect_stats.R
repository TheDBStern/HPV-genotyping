#!/usr/bin/env Rscript

library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

files <- list.files(args[1])
refs <- unique(unlist(lapply(files, function (x) paste(strsplit(x,'\\.')[[1]][1:2],collapse='.'))))

dfout <- data.frame()
for (ref in refs){
    meandepth <- read.table(paste(args[1],'/',ref,'.bam.mean_depth',sep=''))[[1]]
    breadth <- read.table(paste(args[1],'/',ref,'.bam.0.breadth',sep=''))[[1]]
    breadth10 <- read.table(paste(args[1],'/',ref,'.bam.10.breadth',sep=''))[[1]]
    species <- str_remove(read_file(paste(args[1],'/',ref,'.species',sep='')), '\n')
    pwid <- read.table(paste(args[1],'/',ref,'.bam.pwid',sep=''))[[1]]

    df <- data.frame("Ref"=ref,"Species"=species,"MeanDepth"=meandepth,"PercentCovered"=breadth,"PercentCovered10"=breadth10,"PairwiseID"=pwid)
    dfout <- rbind(dfout,df)

    depthdat <- read.table(paste(args[1],'/',ref,'.perbase_depth',sep=''))
    colnames(depthdat) <- c("chr","pos","cov")

    #calculate mean coverage in 100bp windows
    win_size <- 50
    nwindows <- ceiling(nrow(depthdat) / win_size)
    starts <- seq(1,nrow(depthdat),by=win_size)
    ends <- c(seq(win_size,nrow(depthdat),by=win_size),nrow(depthdat))

    winres <- data.frame()
    for (i in 1:length(starts)){
        start <- starts[i]
        end <- ends[i]
        meancov <- depthdat %>%
            slice(start:end) %>%
            summarize(meancov=mean(cov))
        df <- data.frame(start,end,meancov)
        winres <- rbind(winres,df)
    }

    write.table(dfout,args[2],sep='\t',quote=F,row.names=F)

    ggplot(winres, aes(start, meancov)) +
	   #geom_area() +
       geom_rect(aes(xmin=start,xmax=end,ymin=0,ymax=meancov),
                colour="#7FB3D5", fill="#7FB3D5", size=0.5, alpha=1) +
       labs(x="Position",y="Coverage",
            title = paste(ref,":",species,"\nMean Depth:",meandepth,"\nPercent Covered (min 10 reads):",breadth10,"\nPairwise ID:",pwid,sep='')) +
       theme_minimal(base_size = 8)
    ggsave(paste(args[1],'/',ref,'.coverage_stats.pdf',sep=''),
            width = 6, height = 3, units = "in")
    #dev.off()
}
