

setwd("~/R/scrape_geo")
source("functions_scrape_geo.R")
library(data.table)
# DUMP_SCRIPT = "scripts/fastq_dump_wrapper.sh"
DUMP_SCRIPT = "~/R/scrape_geo/scripts/fastq_dump_wrapper.PE.sh"
DUMP_SCRIPT = normalizePath(DUMP_SCRIPT)
stopifnot(file.exists(DUMP_SCRIPT))
SPLIT_SCRIPT = "scripts/split_PE_fastq.sh"
SPLIT_SCRIPT = normalizePath(SPLIT_SCRIPT)
stopifnot(file.exists(SPLIT_SCRIPT))

gse_hits = "GSE161501"


chara = list("cancer type", "cell type")
# Characteristics	

gse_mats = lapply(gse_hits, function(gse){
  mat = scrape_geo(GSE_id = gse,
                   do_srr = TRUE, debug = F, 
                   key_str = c("Title", rep("Chara", length(chara))), 
                   key_idx = append(list(5), chara))#, "cell line", "chip_antibody", "donor_id", "donor_sex"))
  mat
})
sapply(gse_mats, nrow)
# Specify GSE here

odir = paste0("fastqs_GSE161501_Diedrich_ALL_ATAC_and_RNA")
dir.create(odir)
odir = normalizePath(odir)
setwd(odir)

mat = as.matrix(rbindlist(lapply(gse_mats, as.data.frame)))

# recommend saving full matrix to avoid having to scrape again.
write.table(mat, file = "scrape_matrix.txt", sep = "\t")

# transform into cleaner format with better colnames
dt = as.data.table(mat)
setnames(dt, paste0("V", seq(3+length(chara))), c("gse", "gsm", "title", unlist(chara)))
for(ch in chara){
  dt[[ch]] = sub(".+: ", "", dt[[ch]])
}
# recommend saving this processed matrix too.
fwrite(dt, file = "scrape_matrix.processed.txt", sep = "\t")

dt[, file_name := paste(paste0(sub("_", ".", title)), gsm, srr, "fastq", sep = ".")]

# These are the critical inputs for the actual fastq-dump call
fastq_names = dt$file_name
srrs = dt$srr

for(i in seq_along(srrs)){
  cmd = paste("qsub -N fastq-dump -v PATH=$PATH -o", odir, "-e", odir, DUMP_SCRIPT, srrs[i], fastq_names[i])
  jid = system(cmd, intern = TRUE)
  jid = strsplit(jid, " ")[[1]][3]
  message(jid)  
  # cmd2 = paste("qsub -N fastq-split -v PATH=$PATH -o", odir, "-e", odir, 
  #              "-hold_jid", jid,
  #              SPLIT_SCRIPT, fastq_names[i], 100)
  # system(cmd2, intern = TRUE)
}
 q
