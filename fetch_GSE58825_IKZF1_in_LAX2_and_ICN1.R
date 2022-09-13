

setwd("~/R/scrape_geo")
source("functions_scrape_geo.R")
library(data.table)
DUMP_SCRIPT = "scripts/fastq_dump_wrapper.sh"
DUMP_SCRIPT = normalizePath(DUMP_SCRIPT)
stopifnot(file.exists(DUMP_SCRIPT))
SPLIT_SCRIPT = "scripts/split_PE_fastq.sh"
SPLIT_SCRIPT = normalizePath(SPLIT_SCRIPT)
stopifnot(file.exists(SPLIT_SCRIPT))

gse_hits = "GSE58825"


chara = list("cell type", "chip-ab", "replicate")
# Characteristics	
# cell type: Xenograft-derived BCR-ABL1+ (Ph+) human pre-B ALL
# chip-ab: anti-Ikaros, clone H-100 (Santa Cruz, cat # sc-13039 X)
#  replicate: ICN1-replicate1 (H-100 Ab)

gse_mats = lapply(gse_hits, function(gse){
  mat = scrape_geo(GSE_id = gse,
                   do_srr = TRUE, debug = F, 
                   key_str = c("Title", rep("Chara", length(chara))), 
                   key_idx = append(list(5), chara))#, "cell line", "chip_antibody", "donor_id", "donor_sex"))
  mat
})




sapply(gse_mats, nrow)
# Specify GSE here

odir = paste0("fastqs_GSE58825_IKZF1_in_LAX2_and_ICN1")
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

# Now we need to create parallel vectors of final fastq names (fastq_names) and SRR ids (srrs)
# Do any selection/filtering here.
# All useful elements should be carried forward to the fastq name including GSM/SRR

dt[, cell := ifelse(grepl("ICN1", title), "ICN1", "LAX2")]
dt$cell %>% table
dt[, mark := ifelse(grepl("input", title), "input", "IKZF1")]
dt$rep = c("Hrep1", "Hrep1", "rep1", "Crep1", "Hrep2", "Crep1", "Crep2", "rep1")
stopifnot(dt[, .N, .(cell, mark, rep)]$N < 2)

dt[, file_name := paste(paste(cell, mark, rep, sep = "_"), gsm, srr, "fastq", sep = ".")]

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

