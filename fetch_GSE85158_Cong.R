setwd("~/R/scrape_geo")
source("functions_scrape_geo.R")
library(data.table)
DUMP_SCRIPT = "scripts/fastq_dump_wrapper.sh"
DUMP_SCRIPT = normalizePath(DUMP_SCRIPT)
stopifnot(file.exists(DUMP_SCRIPT))
SPLIT_SCRIPT = "scripts/split_PE_fastq.sh"
SPLIT_SCRIPT = normalizePath(SPLIT_SCRIPT)
stopifnot(file.exists(SPLIT_SCRIPT))

# Specify GSE here
gse = "GSE85158"

# navigate to a GSM and determine what Characteristics are available
# cell line: AU565
# cell type: Breast cancer cell line
# treatment: Control
# chip antibody: H3K4me3 (Abcam, ab8580, lot GR1902371-1 and GR224370-1)
chara = list("cell line", "cell type", 'treatment', "chip antibody")
mat = scrape_geo(GSE_id = gse,
                 do_srr = TRUE, debug = F, 
                 key_str = c("Title", rep("Chara", length(chara))), 
                 key_idx = append(list(5), chara))#, "cell line", "chip_antibody", "donor_id", "donor_sex"))
head(mat)



odir = paste0("fastqs_", gse, "_Cong")
dir.create(odir)
odir = normalizePath(odir)
setwd(odir)

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

dt[, c("mark") := tstrsplit(`chip antibody`, " ", keep = 1)]
dt$mark %>% table

dt = dt[mark %in% c("H3K4me3", 'H3K27me3', "none")]
dt$`cell line` %>% unique %>% paste(collapse = "\n") %>% message
dt = dt[`cell line` %in% c("MCF10A", "MCF7", "MB231")]
dt

dt[, file_name := paste0(sub("Input", "input", gsub("\\.", "_", title)), ".", gsm, ".fastq")]

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
