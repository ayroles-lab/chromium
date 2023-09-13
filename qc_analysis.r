install.packages("tidyverse")
library(tidyverse)
library(ggplot2)
df <- read_tsv("barcode_splitter_from_htseq.out")
df$Percent = as.numeric(sub("%", "", df$Percent))
ggplot(df, aes(x=Percent)) +
    geom_histogram(
binwidth=.01,
                   colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") +
    ggtitle("Percent of Reads Mapped to Barcode Distribution")
ggsave("percent_reads_by_barcode_all.png", width=16, height=4)

ggplot(df, aes(x=Percent)) +
    geom_histogram(aes(y=..density..),
                   binwidth=.001,
                   colour="black", fill="white") +
    # geom_density(alpha=.2, fill="#FF6666") +
    xlim(0, 0.05) +
    ggtitle("Percent of Reads Mapped to Barcode Distribution")
ggsave("percent_reads_by_barcode_0_to_1pct.png", width=16, height=4)


