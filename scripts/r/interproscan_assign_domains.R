library(tidyverse)
library(plyranges)
library(parallel)

args <- commandArgs(trailingOnly = T)
in_file <- args[1]
#eval_cutoff <- as.numeric(args[2])
n_cores <- as.numeric(args[2])
out_file <- args[3]
#in_file <- "data/tsv/interproscan_PF00474_interproID.tsv.gz"
#n_cores <- 8


# Read tsv
df_interproscan <- read.csv(
  in_file, colClasses = c("character", "character", "integer", "integer", "numeric", "character", "character"), 
  header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL", "-")) %>% as_tibble()
# Read InterPro short name list
df_interpro_short <- read.csv("db/interpro/interpro_short_name_list.tsv.gz", header = TRUE, sep = "\t", quote = "", na.strings=c("", "NULL", "-")) %>% as_tibble()

# Filter out interpro_acc == NA
# e_val is no longer available because it also contains "scores" (https://interproscan-docs.readthedocs.io/en/latest/UserDocs.html#output-formats)
# Nest by prot_acc & interpro_acc
df_interproscan_nest <- 
  df_interproscan %>% filter(!is.na(interpro_acc)) %>%
  group_by(prot_acc, seq_len, interpro_acc, interpro_desc) %>% nest() %>% 
  left_join(
    df_interpro_short %>% dplyr::rename(interpro_acc = "interpro_id", interpro_protein_n = "protein_count"), 
    by = "interpro_acc") #%>% 
  #select(prot_acc,seq_len,interpro_acc,data,protein_count,short_name,type,interpro_desc)

# Assign range_num
df_interproscan_range_num <- do.call(
  rbind,
  mclapply(
    1:dim(df_interproscan_nest)[1], function (i) {
      #i <- 8
      #i <- 1363
      rng <- df_interproscan_nest[i,]$data[[1]] %>% as_iranges()
      rng_reduced <- reduce_ranges(rng) %>% mutate(range_num = seq(1, length(start)))
      bind_cols(
        df_interproscan_nest[i,c("prot_acc","seq_len","interpro_acc")],
        rng_reduced %>% as.data.frame() %>% as_tibble(),
        df_interproscan_nest[i,c("short_name","type","interpro_desc","interpro_protein_n")],
      ) %>% ungroup()
    }, mc.cores = n_cores
  )
)

# Write tsv
write.table(df_interproscan_range_num, out_file, append=FALSE, quote=FALSE, sep = "\t", row.names=FALSE)
