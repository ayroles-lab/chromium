library(data.table)
library(aod)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
file_number <- as.integer(args[1])
# file_number = 1

alt_count_files <- list.files(path="/Genomics/ayroleslab2/scott/git/chromium/data/counts_matrices/", pattern = "alternate_counts*", full.names=TRUE)
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

coef_betabin <- matrix(nrow=dim(alt_count_data)[1],ncol=3)
coef_glm <- matrix(nrow=dim(alt_count_data)[1],ncol=3)


B <- 10000
plouf <- CJ(it = 1:B)
library(tictoc)
# pb = txtProgressBar(min=0,max=dim(alt_count_data)[1],style=3)
for (i in 1:dim(alt_count_data)[1]){
	# setTxtProgressBar(pb,i)
	dt = data.table(depth = unlist(both_count_data[i,],use.names=FALSE), alt_depth = unlist(alt_count_data[i,],use.names=FALSE), sex = info$sex, treatment = info$treatment)
	dt = dt[dt$depth > 0,]

    dt$treatment <- factor(dt$treatment,levels=c("T0","Control","3mM"))

    mod1 <- tryCatch(betabin(cbind(alt_depth, depth - alt_depth) ~ treatment + sex, ~1,data=dt), error=function(x){})
    obs_ctrl = attributes(summary(mod1))$Coef[2,1]
    obs_crvi = attributes(summary(mod1))$Coef[3,1]

    coef_betabin[i,] = attributes(summary(mod1))$Coef[2:4,1]

    mod_glm <- tryCatch(glm(cbind(alt_depth, depth - alt_depth) ~ treatment + sex, data=dt, family = binomial), error=function(x){})
    coef_glm[i,] = coef(mod_glm)[2:4]
    
    f = function(x) x@fixed.param[2:4]
    dt[,.(fitted = f(betabin(cbind(alt_depth, depth - alt_depth) ~ treatment + sex , random = ~1, data = .SD)))]
    
    Ntrial <- plouf[,.(depth = dt$depth, alt_depth =dt$alt_depth, sex = dt$sex, treatment = sample(dt$treatment)), by = .(it)]

    # # betabin
    # res <- Ntrial[, .(tryCatch(f(betabin(cbind(alt_depth, depth - alt_depth) ~ treatment + sex, ~1,data=.SD)), error=function(x){ NA })),
    #     by = .(it)] 

    f_glm = function(x) coef(x)[3:4]
    res = Ntrial[,.(glm = f_glm(glm(cbind(alt_depth, depth - alt_depth) ~ sex + treatment, family = binomial, data = .SD))), by = .(it)]


    p_ctrl = sum((abs(res$glm[seq(1,nrow(res),by=2)]) > abs(coef(mod_glm)[2]))/B)
    p_crvi = sum((abs(res$glm[seq(2,nrow(res),by=2)]) > abs(coef(mod_glm)[3]))/B)
    pvals[i,] = c(p_ctrl,p_crvi)
    
    print(paste0("Number " ,i, " of ", dim(alt_count_data)[1], " (", round(i/dim(alt_count_data)[1],5), "%)"))
}

st=format(Sys.time(), "%Y%m%d_%H%M")

output_table = data.table(snp=rownames(alt_count_data)[1:length(pvals)], pvalues_ctrl = as.vector(pvals[,1]),
    pvalues_crvi = as.vector(pvals[,2]),coef_betabin_control=coef_betabin[,1],coef_betabin_crvi = coef_betabin[,2],coef_betabin_sex_male=coef_betabin[,3],
    coef_glm_control=coef_glm[,1],coef_glm_crvi = coef_glm[,2],coef_glm_sex_male=coef_betabin[,3])

write.table(output_table,paste0(toString(st),'_glm_permutation_pvals_10000iter_',toString(file_number),'.csv'),sep=',',row.names=F)
