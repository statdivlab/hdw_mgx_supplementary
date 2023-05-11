##  Summarize contigs assembled 

contigs_summary_files <- list.files(path = "/Users/paulinetrinh/Documents/HDW/metagenomic/HMP_comparison/hdw_contigs_summary", pattern = "*.txt", full.names = T)

formatting_fn <- function(x) {
  read_tsv(x) %>% 
    column_to_rownames("contigs_db")
}

summary_info <- sapply(contigs_summary_files, formatting_fn, simplify=FALSE)  
summary_tbl <- do.call(cbind,summary_info)

hdw_contigs_summary_files <- list.files(path = "/Users/paulinetrinh/Documents/HDW/metagenomic/HMP_comparison/hdw_contigs_summary/HDW_files", pattern = "*.txt", full.names = T)
hdw_summary_info <- sapply(hdw_contigs_summary_files, formatting_fn, simplify=FALSE)  
hdw_summary_tbl <- do.call(cbind,hdw_summary_info)


all_summary_tbl <- merge(hdw_summary_tbl, summary_tbl, by = 'row.names', all = TRUE)
