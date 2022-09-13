setwd("~/R/scrape_geo")
source("functions_scrape_geo.R")
library(data.table)
DUMP_SCRIPT = "scripts/fastq_dump_wrapper.sh"
DUMP_SCRIPT = normalizePath(DUMP_SCRIPT)
stopifnot(file.exists(DUMP_SCRIPT))
SPLIT_SCRIPT = "scripts/split_PE_fastq.sh"
SPLIT_SCRIPT = normalizePath(SPLIT_SCRIPT)
stopifnot(file.exists(SPLIT_SCRIPT))


gse = "GSE53595"
chara = list("cell type", 'strain')
mat = scrape_geo(GSE_id = gse, gsm_override = c("GSM1296535", "GSM1296537", "GSM1296572", "GSM1296573", "GSM1296574", "GSM1296575", "GSM1296576", "GSM1296577", "GSM1296578", "GSM1296579"),
                 do_srr = TRUE, debug = F, key_str = c("Title", rep("Chara", length(chara))), 
                 key_idx = append(list(5), chara))#, "cell line", "chip_antibody", "donor_id", "donor_sex"))
head(mat)
mat

mat[grepl("ChIP", mat[,3]),]

mat = mat[grepl("ChIP", mat[,3]),]



# mat = cbind(mat, rownames(mat))
dt = as.data.table(mat)
setnames(dt, paste0("V", seq(3+length(chara))), c("gse", "gsm", "title", unlist(chara)))
for(ch in chara){
  dt[[ch]] = sub(".+: ", "", dt[[ch]])
}


odir = paste0("fastqs_", gse, "_Ikaros_proB")
dir.create(odir)
odir = normalizePath(odir)
setwd(odir)


dt[, c("mark", "status") := tstrsplit(title, " ", keep = c(1, 4))]
dt$mouse = sub(":.+", "", dt$title)
dt[, stage := "matB"]

dt[, status := c("Ikzf1(∆/+)" = "haplo", "Ikzf1(∆/–)" = "del", "9102" = "WT", "10195" = "WT")[status]]
dt$tissue = c(rep("lymphnode", 3), "spleen")
dt[, rep := paste0("rep", seq(.N)), .(mark, status)]
dt[,  rep := paste0("rep", seq(.N)), .(stage, tissue)]
dt$mark = "Ikaros"

# dt[, file_name := file.path(odir, paste0("Busslinger_", "proB_", ifelse(grepl("Ikaros", title), "IKAROS", "input"), "_R1.fastq"))]
dt[, file_name := paste(paste("Busslinger", stage, mark, rep, sep = "_"), mouse, "fastq", sep = ".")]


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
