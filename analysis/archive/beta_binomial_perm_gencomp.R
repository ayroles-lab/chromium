library(data.table)
library(aod)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
# file_number <- as.integer(args[1])
file_number = 1

alt_count_files <- list.files(path="/Genomics/ayroleslab2/scott/git/chromium/data/counts_matrices/", pattern = "alternate_counts_*", full.names=TRUE)
ref_count_files <- list.files(path="/Genomics/ayroleslab2/scott/git/chromium/data/counts_matrices/", pattern = "reference_counts_*", full.names=TRUE)

# fread is faster than read.table
alt_count_data <- read.table(alt_count_files[file_number], header=TRUE, sep="\t", row.names=1)
ref_count_data <- read.table(ref_count_files[file_number], header=TRUE, sep="\t", row.names=1)

both_count_data <- alt_count_data + ref_count_data

# count data to use for beta binomial modeling
# alt=fread('/Genomics/ayroleslab2/shared/longevity/final_files/29Jun20_merged_alt_counts_allCHR.txt')
# both=fread('/Genomics/ayroleslab2/shared/longevity/final_files/29Jun20_merged_both_counts_allCHR.txt')
info=read_csv('/Genomics/ayroleslab2/scott/git/chromium/individual_metadata.csv')
info$name = gsub('2428__','',info$indiv)

info = info[order(info$name),]
both_count_data = both_count_data[,order(colnames(both_count_data))]
alt_count_data = alt_count_data[,order(colnames(alt_count_data))]
pvals <- matrix(nrow=dim(alt_count_data)[1],ncol=2)

library(tictoc)
number_of_loci = dim(alt_count_data)[1]

# pb = txtProgressBar(min=0,max=dim(alt_count_data)[1],style=3)
output_mat = matrix(nrow=number_of_loci,ncol=5)

# which(output_mat[,"SNP"] == "2R:19643940")
# 
colnames(output_mat) <- c("SNP","obs_ctrl","pval_ctrl","obs_crvi","pval_crvi")
tic()
# i = which(output_mat[,"SNP"] == "2R:19643940")

for (i in 1:number_of_loci){

	dt = data.table(depth = unlist(both_count_data[i,],use.names=FALSE), alt_depth = unlist(alt_count_data[i,],use.names=FALSE), sex = info$sex, treatment = info$treatment)
	dt = dt[dt$depth > 0,]

    dt$treatment <- factor(dt$treatment,levels=c("T0","Control","3mM"))

    mod1 <- tryCatch(betabin(cbind(alt_depth, depth - alt_depth) ~ treatment + sex, ~1,data=dt,method="L-BFGS-B"), error=function(x){})
    
    # Add try-catch blocks to handle errors in accessing the summary attributes
    obs_ctrl <- tryCatch(attributes(summary(mod1))$Coef[2,1], error=function(x) NA)
    pval_ctrl <- tryCatch(attributes(summary(mod1))$Coef[2,4], error=function(x) NA)

    obs_crvi <- tryCatch(attributes(summary(mod1))$Coef[3,1], error=function(x) NA)
    pval_crvi <- tryCatch(attributes(summary(mod1))$Coef[3,4], error=function(x) NA)
    
    # output_mat[i,] = c(rownames(alt_count_data)[i],obs_ctrl,pval_ctrl,obs_crvi,pval_crvi)

    print(paste0("Number " ,i, " of ", dim(alt_count_data)[1], " (", 100*round(i/number_of_loci,5), "%)"))
}
toc()

# output_mat[output_mat[,"SNP"] == "2R:19643940",]
# both_count_data[output_mat[,"SNP"] == "2R:19643940",]

write.table(output_mat, file=paste0("/Genomics/ayroleslab2/scott/git/chromium/data/tmp_betabin_gc_",file_number,".csv"), sep=",", row.names=FALSE, col.names=TRUE)


# # Get number of nas in ref and alt for snp
# sum(is.na(alt_count_data[output_mat[,"SNP"] == "2R:19643940",]))
# sum(is.na(ref_count_data[output_mat[,"SNP"] == "2R:19643940",]))

# dim(ref_count_data[output_mat[,"SNP"] == "2R:19643940",])
# summary(t(ref_count_data[output_mat[,"SNP"] == "2R:19643940",]))

# ref_count_data[output_mat[,"SNP"] == "2R:19643940",]

# table(info[which(!is.na(ref_count_data[output_mat[,"SNP"] == "2R:19643940",])),c("treatment","sex")])
