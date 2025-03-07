library(httr2)
library(jsonlite)
library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)

# Set vriavles
interpro_id <- args[1]
#interpro_id <- "PF00864"
entry_id <- args[2]
#entry_id <- "pfam"
out_file <- args[3]

# API endpoint url (page_size limit = 200)
url <- paste0("https://www.ebi.ac.uk/interpro/api/protein/uniprot/entry/", entry_id, "/", 
              interpro_id, "?page_size=200&extra_fields=sequence,is_fragment")

df_result_unify <- tibble()
while (!is.null(url)) {
  # create and send a request
  response <- request(url) %>% req_perform()
  # get content in text format
  content <- resp_body_string(response)
  # convert JSON data to list format
  data <- fromJSON(content)
  #data %>% str()
  # Extract results
  metadata_cols <- colnames(data$results$metadata)
  df_metadata <- data$results$metadata[,metadata_cols[!metadata_cols %in% c("source_organism")]]
  df_source_organism <- data$results$metadata[,c("source_organism")]
  #df_entries <- do.call(rbind, data$results$entries)
  df_extra_fields <- data$results$extra_fields
  # Combine result tsv
  df_result_unify <- rbind(df_result_unify, cbind(df_metadata, df_source_organism, df_extra_fields)) %>% as_tibble()
  #df_result_unify <- rbind(df_result_unify, cbind(df_metadata, df_source_organism, df_entries, df_extra_fields)) %>% as_tibble()
  # Set next URL
  url <- data$`next`
  # Sleep 3 sec
  Sys.sleep(3)
}

# Export tsv
write.table(df_result_unify, out_file, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)
