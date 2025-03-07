library(httr2)
library(jsonlite)
library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)

# Set vriavles
in_file <- args[1]
# in_file <- "db/interpro_api/PF00864.tsv.gz"
out_file <- args[2]

# Load uniprot protein tsv
df_protein <- read.csv(in_file, header = TRUE, sep = "\t", quote = "\"", na.strings=c("", "NULL")) %>% as_tibble()

#### Extract protein info by API ####
# Set empty df
df_result_unify <- tibble()
# API loop
for (i in 1:dim(df_protein)[1]) {
  acc <- df_protein$accession[i]
  source_db <- df_protein$source_database[i]
  print(paste(i,acc,source_db))
  
  #### InterPro ####
  # API endpoint
  url_interpro <- paste0("https://www.ebi.ac.uk/interpro/api/entry/InterPro/protein/", source_db, "/", acc, 
                "?page_size=200&extra_fields=short_name")
  # create and send a request
  interpro_response <- request(url_interpro) %>% req_perform()
  if (!resp_status_desc(interpro_response) == "No Content") {
    # get content in text format
    interpro_content <- resp_body_string(interpro_response)
    # convert JSON data to list format
    interpro_data <- fromJSON(interpro_content)
    # interpro_data %>% str()
    # Extract results
    df_interpro_metadata <- interpro_data$results$metadata[c("accession","name","source_database","type","integrated")]
    df_interpro_extra_fields <- interpro_data$results$extra_fields
  } else {
    df_interpro_metadata <- data.frame(matrix(ncol = 5, nrow = 0))
    colnames(df_interpro_metadata) <- c("accession","name","source_database","type","integrated")
    df_interpro_extra_fields <- data.frame(matrix(ncol = 1, nrow = 0))
    colnames(df_interpro_extra_fields) <- c("short_name")
  }
  
  #### Pfam ####
  # API endpoint
  url_pfam <- paste0("https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/", source_db, "/", acc, 
                         "?page_size=200&extra_fields=short_name")
  # create and send a request
  pfam_response <- request(url_pfam) %>% req_perform()
  if (!resp_status_desc(pfam_response) == "No Content") {
    # get content in text format
    pfam_content <- resp_body_string(pfam_response)
    # convert JSON data to list format
    pfam_data <- fromJSON(pfam_content)
    # pfam_data %>% str()
    # Extract results
    df_pfam_metadata <- pfam_data$results$metadata[c("accession","name","source_database","type","integrated")]
    df_pfam_extra_fields <- pfam_data$results$extra_fields
  } else {
    df_pfam_metadata <- data.frame(matrix(ncol = 5, nrow = 0))
    colnames(df_pfam_metadata) <- c("accession","name","source_database","type","integrated")
    df_pfam_extra_fields <- data.frame(matrix(ncol = 1, nrow = 0))
    colnames(df_pfam_extra_fields) <- c("short_name")
  }
  
  #### Combine result tsv ####
  df_results <- rbind(
    cbind(df_interpro_metadata, df_interpro_extra_fields),
    cbind(df_pfam_metadata, df_pfam_extra_fields)
  ) %>% as_tibble()
  df_result_unify <- rbind(df_result_unify, cbind(acc, df_results))
  
  # Sleep 1 sec
  Sys.sleep(1)
}

#### Export tsv ####
# Rename col names
colnames(df_result_unify) <- 
  c("accession", "dbEntry", "dbEntry_name", "dbEntry_source", "dbEntry_type", "dbEntry_integrated", "dbEntry_short_name")
# Export
write.table(df_result_unify, out_file, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)
