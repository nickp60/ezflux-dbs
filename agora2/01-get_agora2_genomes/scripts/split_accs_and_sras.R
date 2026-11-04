library(tidyverse)


if (!file.exists("41587_2022_1628_MOESM3_ESM.xlsx")) download.file("https://static-content.springer.com/esm/art%3A10.1038%2Fs41587-022-01628-0/MediaObjects/41587_2022_1628_MOESM3_ESM.xlsx", "41587_2022_1628_MOESM3_ESM.xlsx")
metaraw <- readxl::read_excel("./41587_2022_1628_MOESM3_ESM.xlsx", sheet = 2, skip = 3) %>%
  rename(
      genome_link = `Genome link`,
      ref_seq_id = `RefSeq ID`
  )
  #janitor::clean_names()

meta <- metaraw %>%
  mutate(x = gsub(".*(GC._.*?)_.*", "\\1", genome_link),
         x= ifelse(grepl("https://www.ncbi.nlm.nih.gov/assembly/",x),
                   gsub("/", "", gsub("https://www.ncbi.nlm.nih.gov/assembly/(.*?)", "\\1", x)), x),
         x = ifelse(grepl("ftp", x), ref_seq_id, x),
         sra = ifelse(grepl("fastq.gz", genome_link), basename(gsub(";.*", "", genome_link)), ""),
         sra = gsub("_1.fastq.gz", "", sra)
  )

table(is.na(meta$x))
meta$x
meta$sra
table(!is.na(meta$x) | meta$sra != "" )
meta %>% filter(is.na(x) & sra == "")

cols_of_interest = c("MicrobeID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# we filter for sra == "" for accs to cover cases like GCF_000369125.1 where we have both an assembly and an sra
accs_out <- meta %>% filter(sra == "") %>% select(all_of(cols_of_interest), x) %>% filter(!is.na(x)) %>%
group_by(x) %>% slice(1)
write.table(accs_out, "agora2_accs", row.names = FALSE, quote=F, sep="\t")
meta %>% filter(is.na(x))  %>% select(all_of(cols_of_interest), sra) %>% rename(x=sra) %>% filter(x != "" ) %>%
    write.table("agora2_sras", row.names = FALSE, quote = F, sep="\t")
