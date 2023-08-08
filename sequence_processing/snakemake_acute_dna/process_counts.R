library(data.table)

args <- commandArgs(trailingOnly = TRUE)
ad_file <- args[1]
dp_file <- args[2]
alt_counts_file <- args[3]
ref_counts_file <- args[4]

ad <- fread(ad_file)
dp <- fread(dp_file)
ad <- unique(ad)
dp <- unique(dp)

ad$loc_id <- paste(ad$V1, ad$V2, sep='_')
dp$loc_id <- paste(dp$V1, dp$V2, sep='_')

both <- merge(ad, dp, by='loc_id')
both2 <- subset(both, (V3.x + V4) == V3.y)

ad2 <- ad[which(ad$loc_id %in% both2$loc_id), ]

setkey(ad2, V1, V2)
names <- unique(ad2[, V2])
sites <- unique(ad2[, V1])

ref_counts <- as.data.table(unique(ad2[, V1]))
setnames(ref_counts, "V1", "site")
for(i in 1:length(names)) {
  ref_counts[, names[i] := ad2[.(as.factor(ref_counts[, site]), names[i]), V3]]
}

alt_counts <- as.data.table(unique(ad2[, V1]))
setnames(alt_counts, "V1", "site")
for(i in 1:length(names)) {
  alt_counts[, names[i] := ad2[.(as.factor(alt_counts[, site]), names[i]), V4]]
}

write.table(alt_counts, alt_counts_file, row.names=F, sep='\t', quote=F)
write.table(ref_counts, ref_counts_file, row.names=F, sep='\t', quote=F)
