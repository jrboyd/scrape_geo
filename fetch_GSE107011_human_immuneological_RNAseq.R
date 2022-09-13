setwd("~/R/scrape_geo")
source("functions_scrape_geo.R")
library(data.table)
DUMP_SCRIPT = "scripts/fastq_dump_wrapper.sh"
DUMP_SCRIPT = normalizePath(DUMP_SCRIPT)
stopifnot(file.exists(DUMP_SCRIPT))
SPLIT_SCRIPT = "scripts/split_PE_fastq.sh"
SPLIT_SCRIPT = normalizePath(SPLIT_SCRIPT)
stopifnot(file.exists(SPLIT_SCRIPT))


gse = "GSE107011"
chara = list("cell type", 'disease status', 'gender')
mat = scrape_geo(GSE_id = gse, 
                 do_srr = TRUE, debug = F, key_str = c("Title", rep("Chara", length(chara))), 
                 key_idx = append(list(5), chara))#, "cell line", "chip_antibody", "donor_id", "donor_sex"))
head(mat)
mat
# mat = cbind(mat, rownames(mat))
dt = as.data.table(mat)
setnames(dt, paste0("V", seq(3+length(chara))), c("gse", "gsm", "title", unlist(chara)))
for(ch in chara){
  dt[[ch]] = sub(".+: ", "", dt[[ch]])
}

dt[is.na(`cell type`), `cell type` := "human"]
dt[, title := gsub("\\(", "", title)]
dt[, title := gsub("\\)", "", title)]
dt[, title := gsub('/', "-", title)]
dt[grepl("-", title)]
# dt = dt[`cell type` == "human"]

table(dt$`disease status` )

dt[, file_name := file.path(odir, paste0(paste("human", title, sep = "_"), ".fastq"))]
dt[, root := paste("human", title, sep = "_")]

fwrite(dt[, .(srr, root)], "../../GSE107011_monaco_human_immune.csv")

odir = paste0("fastqs_", gse, "_RNAseq")
dir.create(odir)
odir = normalizePath(odir)
setwd(odir)




fastq_names = dt$file_name
srrs = dt$srr

for(i in seq_along(srrs)){
  cmd = paste("qsub -N fastq-dump -v PATH=$PATH -o", odir, "-e", odir, DUMP_SCRIPT, srrs[i], fastq_names[i])
  jid = system(cmd, intern = TRUE)
  jid = strsplit(jid, " ")[[1]][3]
  message(jid)
  # cmd2 = paste("qsub -N fastq-split -v PATH=$PATH -o", odir, "-e", odir,
  #              # "-hold_jid", jid,
  #              SPLIT_SCRIPT, fastq_names[i], 50)
  # # message(cmd2)
  # system(cmd2, intern = TRUE)
}
