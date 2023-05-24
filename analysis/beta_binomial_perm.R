library(data.table)
library(aod)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
file_number <- as.integer(args[1])
file_number = 1

alt_count_files <- list.files(path="/Genomics/ayroleslab2/scott/git/chromium/data/counts_matrices/", pattern = "alternate_counts*", full.names=TRUE)
ref_count_files <- list.files(path="/Genomics/ayroleslab2/scott/git/chromium/data/counts_matrices/", pattern = "reference_counts_*", full.names=TRUE)

# fread is faster than read.table
alt_count_data <- read.table(alt_count_files[file_number], header=TRUE, sep="\t", row.names=1)
ref_count_data <- read.table(ref_count_files[file_number], header=TRUE, sep="\t", row.names=1)

both_count_data <- alt_count_data + ref_count_data


info=read_csv('/Genomics/ayroleslab2/scott/git/chromium/individual_metadata.csv')
info$name = gsub('2428__','',info$indiv)

info = info[order(info$name),]
both_count_data = both_count_data[,order(colnames(both_count_data))]
alt_count_data = alt_count_data[,order(colnames(alt_count_data))]

# Filter snp data to only include SNPs with non-zero alts
filter <- rowSums(alt_count_data,na.rm=TRUE) == 0 | apply(alt_count_data,1,function(x) all(is.na(x)))
alt_count_data <- alt_count_data[!filter,]
both_count_data <- both_count_data[!filter,]

coef_betabin <- matrix(nrow=dim(alt_count_data)[1],ncol=2)


B <- 10
plouf <- CJ(it = 1:B)
library(tictoc)
number_of_loci = dim(alt_count_data)[1]
# number_of_loci = 10
# pb = txtProgressBar(min=0,max=dim(alt_count_data)[1],style=3)
output_mat = matrix(nrow=number_of_loci,ncol=8)
colnames(output_mat) <- c("obs_ctrl","pval_ctrl","obs_crvi","pval_crvi","obs_contrast","pval_contrast","emp_fdr_pvalthresh_05","significant")
tic()
i= 347
for (i in 1:number_of_loci){
	# setTxtProgressBar(pb,i)
    # i=3538
	dt = data.table(depth = unlist(both_count_data[i,],use.names=FALSE), alt_depth = unlist(alt_count_data[i,],use.names=FALSE), sex = info$sex, treatment = info$treatment)
	dt = dt[dt$depth > 0,]

    dt$treatment <- factor(dt$treatment,levels=c("T0","Control","3mM"))

    mod1 <- tryCatch(betabin(cbind(alt_depth, depth - alt_depth) ~ treatment + sex, ~1,data=dt,method="L-BFGS-B",hessian=TRUE), error=function(x){})
    
    # Add try-catch blocks to handle errors in accessing the summary attributes
    obs_ctrl <- tryCatch(attributes(summary(mod1))$Coef[2,1], error=function(x) NA)
    pval_ctrl <- tryCatch(attributes(summary(mod1))$Coef[2,4], error=function(x) NA)

    obs_crvi <- tryCatch(attributes(summary(mod1))$Coef[3,1], error=function(x) NA)
    pval_crvi <- tryCatch(attributes(summary(mod1))$Coef[3,4], error=function(x) NA)


    coef_betabin[i,] <- tryCatch(attributes(summary(mod1))$Coef[2:3,1], error=function(x) c(NA, NA))

    dt_no_t0 <- dt[dt$treatment %in% c("Control","3mM")]

    mod2 <- tryCatch(betabin(cbind(alt_depth, depth - alt_depth) ~ treatment + sex, ~1,data=dt_no_t0,method="L-BFGS-B",hessian=TRUE), error=function(x){})
    obs_contrast <- tryCatch(attributes(summary(mod2))$Coef[2,1], error=function(x) NA)
    pval_contrast <- tryCatch(attributes(summary(mod2))$Coef[2,4], error=function(x) NA)
    
    f <- function(x) tryCatch(attributes(summary(x))$Coef[2,4], error=function(x) NA)
    dt_no_t0[, .(fitted = f(betabin(cbind(alt_depth, depth - alt_depth) ~ treatment + sex , random = ~1, data = .SD, method="L-BFGS-B",hessian=TRUE)))]

    Ntrial <- plouf[,.(depth = dt_no_t0$depth, alt_depth =dt_no_t0$alt_depth, sex = dt_no_t0$sex, treatment = sample(dt_no_t0$treatment)), by = .(it)]

    # # betabin
    res <- Ntrial[, .(tryCatch(f(betabin(cbind(alt_depth, depth - alt_depth) ~ treatment + sex, ~1,data=.SD,method="L-BFGS-B",hessian=TRUE)), error=function(x){ NA })),
        by = .(it)]

    fdr_pval <- tryCatch(quantile(res$V1, 0.05), error=function(x) NA)

    sig <- tryCatch(as.numeric(pval_contrast < fdr_pval), error = function(x) NA)
    
    output_mat[i,] = c(obs_ctrl,pval_ctrl,obs_crvi,pval_crvi,obs_contrast,pval_contrast,fdr_pval,sig)

    print(paste0("Number " ,i, " of ", dim(alt_count_data)[1], " (", 100*round(i/dim(alt_count_data)[1],5), "%)"))
}
toc()

write.table(output_mat, file=paste0("/Genomics/ayroleslab2/scott/git/chromium/data/betabin_contrast_w_perm/betabin_perm_output_10perm_fdr05_",file_number,".csv"), sep=",", row.names=FALSE, col.names=TRUE)
