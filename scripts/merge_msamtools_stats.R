library(dplyr)

files <- read.table('../files.txt')

combine_res <- function(filename){
    abund <- read.table(paste(filename,'.pave.profile.txt',sep=''),h=T)
    colnames(abund) <- c('RelAb')
    abund$Ref <- rownames(abund)
    stats <- read.delim(paste(filename,'.pave.stats.txt',sep=''),h=T)
    res <- merge(stats, abund,by='Ref') %>%
            arrange(desc(RelAb))
    write.table(res,paste(filename,'.results.txt',sep=''),sep='\t',quote=F,row.names=F)
}

sapply(files[,1],combine_res)
