

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

gse = "GSE69149"


chara = list("cancer type", "cell type")
# Characteristics	

mat = scrape_geo(GSE_id = gse,
                 do_srr = TRUE, debug = F, 
                 key_str = c("Title", rep("Chara", length(chara))), 
                 key_idx = append(list(5), chara))#, "cell line", "antibody", "replicate"))
mat

# Specify GSE here

odir = paste0("fastqs_GSE69149_Sokolova16_U2OS_and_RPE_NPAT_and_HINF_chipseq")
dir.create(odir)
odir = normalizePath(odir)
setwd(odir)

write.table(mat, file = "scrape_matrix.txt", sep = "\t")

# transform into cleaner format with better colnames
dt = as.data.table(rbind(mat))
setnames(dt, paste0("V", seq(3+length(chara))), c("gse", "gsm", "title", unlist(chara)))
for(ch in chara){
  dt[[ch]] = sub(".+: ", "", dt[[ch]])
}
# recommend saving this processed matrix too.
fwrite(dt, file = "scrape_matrix.processed.txt", sep = "\t")

dt = dt[grepl("ChIP", title)]
dt = dt[!(grepl("E2F1", title) | grepl("CASP", title))]
dt[, c("method", "target", "cell") := tstrsplit(title, "_")]
dt[, cell := sub("hTERT-", "", cell)]
dt[, rep := "rep1"]
dt[grepl("IgG", target), rep := target]
dt[grepl("IgG", target), target := "input"]

dt[, title := paste(cell, target, rep, sep = "_")]

dt[, file_name := paste(title, gsm, srr, "fastq", sep = ".")]

# These are the critical inputs for the actual fastq-dump call
fastq_names = dt$file_name
srrs = dt$srr

for(i in seq_along(srrs)){
  cmd = paste("qsub -N fastq-dump -v PATH=$PATH -o", odir, "-e", odir, DUMP_SCRIPT, srrs[i], fastq_names[i])
  jid = system(cmd, intern = TRUE)
  jid = strsplit(jid, " ")[[1]][3]
  message(jid)  
}
