source("functions_scrape_geo.R")
library(data.table)

DUMP_SCRIPT = "scripts/fastq_dump_wrapper.sh"
DUMP_SCRIPT = normalizePath(DUMP_SCRIPT)
stopifnot(file.exists(DUMP_SCRIPT))
SPLIT_SCRIPT = "scripts/split_PE_fastq.sh"
SPLIT_SCRIPT = normalizePath(SPLIT_SCRIPT)
stopifnot(file.exists(SPLIT_SCRIPT))


gse = "GSE40129" #will be combined with GSE19465
chara = list("cell type", 'chip antibody', "stimulated with")
mat = scrape_geo(GSE_id = gse, 
                 do_srr = TRUE, debug = F, key_str = c("Title", rep("Chara", length(chara))), 
                 key_idx = append(list(5), chara))#, "cell line", "chip_antibody", "donor_id", "donor_sex"))
head(mat)
# mat = cbind(mat, rownames(mat))
dt = as.data.table(mat)
setnames(dt, c("V1", "V2", "V3", "V4", "V5", "V6"), c("gse", "gsm", "title", unlist(chara)))
for(ch in chara){
  dt[[ch]] = sub(".+: ", "", dt[[ch]])
}

odir = paste0("fastqs_", gse, "_Carroll_ESR1_enhancers")
dir.create(odir)
odir = normalizePath(odir)
setwd(odir)

dt[, file_name := file.path(odir, paste0(title, ".fastq"))]
fastq_names = dt$file_name
srrs = dt$srr

for(i in seq_along(srrs)){
  cmd = paste("qsub -N fastq-dump -v PATH=$PATH -o", odir, "-e", odir, DUMP_SCRIPT, srrs[i], fastq_names[i])
  jid = system(cmd, intern = TRUE)
  jid = strsplit(jid, " ")[[1]][3]
  
  # cmd2 = paste("qsub -N fastq-split -v PATH=$PATH -o", odir, "-e", odir, 
  #              "-hold_jid", jid,
  #              SPLIT_SCRIPT, fastq_names[i], 100)
  # system(cmd2, intern = TRUE)
}
