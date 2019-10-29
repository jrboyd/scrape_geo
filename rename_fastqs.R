library(data.table)
library(magrittr)
library(openxlsx)
geo_dt = fread("GSE63525.csv", header = F)
colnames(geo_dt) = c("GSM", "cell_line", "hic_id", "type", "cell_desc", "srr")

supp_dt = fread("GSE63525_supp_info.csv")
colnames(supp_dt) = c("hic_id", "cell_line", "bio_rep", "hic_type", "digest", "fixation", "biotin", "a", "b", "c", "d")

geo_dt = merge(geo_dt, supp_dt, by = "hic_id")
# pub_df = openxlsx::read.xlsx("pnas.1518552112.st01.xlsx", sheet = 6, startRow = 2)
# pub_df = pub_df[!is.na(df$Genotype),]
# pub_dt = as.data.table(pub_df)
# pub_dt[, hic_id := tstrsplit(Library, "-", keep = 3)]


# dt = merge(geo_dt, pub_dt)
setwd("~/HiC-Pro/inputs/raohuntley")
fastq_f = dir(pattern = ".fastq$")
file_dt = data.table(file_path = fastq_f, srr = basename(fastq_f) %>% sub(".fastq", "", .))
dt = merge(geo_dt, file_dt, by = "srr")
dim(dt)
# dt$cell_line = 
# sub(" (CCL186)", "", dt$cell_line)
dt$cell_line = sapply(strsplit(dt$cell_line.x, " "), function(x)x[1])
for(i in seq_len(nrow(dt))){
  old_f = dt[i, file_path]
  srr = sapply(strsplit(old_f, "\\."), function(x)x[1])
  new_f = paste0(dt[i, paste(cell_line, hic_id, GSM, sep = "_")], "_", srr, "_", dt$digest, ".fastq")
  cmd = paste("mv", old_f, new_f)
  print(cmd)
  # system(cmd)
}
