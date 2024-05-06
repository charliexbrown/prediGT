args = commandArgs(trailingOnly=TRUE)

library(qqman)
library(dplyr)

if (length(args)==0) {
  stop("usage is args[0], input file and output_path.n", call.=FALSE)
}

## define output files
manhattan_file<- file.path(args[2], "manhattan.png")
qqfile <- file.path(args[2], "qqplot.png")

## read plink results df
plink.result <- read.delim(args[1], header = TRUE, comment.char = '$')
## rename cols
colnames(plink.result) <- c("CHR", "POS", "ID", "REF", "ALT", "A1", "TEST", "OBS_CT", "OR", "SE", "Z_STAT", "P")

## remove snp where p is na
plink.result<- plink.result[!is.na(plink.result$P),]

## function to retain top 5000 snps based on pvalue
sample_gwas_result <- function(gwas_result, p_cutoff = 1e-4, n_snp = 5000, by_p = TRUE) {
  gwas_result <- gwas_result %>% filter(TEST == 'ADD')
  gwas_result_list <- split(gwas_result, gwas_result$P)
  
   if (by_p){
    snps_selected <- sample(1:nrow(gwas_result_list[[2]]), n_snp,p=1/gwas_result_list[[2]]$P, replace=TRUE)
   } else{ 
     snps_selected <- sample(1:nrow(gwas_result_list[[2]]), n_snp, replace=TRUE)
   }
  
  gwas_result_list[[2]] <- gwas_result_list[[2]][snps_selected, ]
  bind_rows(gwas_result_list)
}

small_plink_result <- sample_gwas_result(plink.result)
small_plink_result %>% arrange(P) %>% head 


names(small_plink_result)[1] <- 'CHR'

## transform chrX to number
small_plink_result$chrno <- with(small_plink_result, as.numeric(ifelse(CHR=='X', '23',CHR)))

## remove duplicated snps and with Na pval
small_plink_result<- subset(small_plink_result, !is.na(P) & !duplicated(ID))

#manhattanplot
png(manhattan_file, type="cairo")
manhattan(small_plink_result, bp="POS", chr="chrno", snp="ID", annotatePval = 0.01, annotateTop= FALSE,  genomewideline = -log10(1e-08), col = c("blue4", "orange3"))
dev.off()

## qqplot
png(qqfile, type="cairo")
qq(plink.result$P)
dev.off()
