setwd("~/R/scrape_geo")
source("functions_scrape_geo.R")
library(data.table)
DUMP_SCRIPT = "scripts/fastq_dump_wrapper.sh"
DUMP_SCRIPT = normalizePath(DUMP_SCRIPT)
stopifnot(file.exists(DUMP_SCRIPT))
SPLIT_SCRIPT = "scripts/split_PE_fastq.sh"
SPLIT_SCRIPT = normalizePath(SPLIT_SCRIPT)
stopifnot(file.exists(SPLIT_SCRIPT))


gse = "GSE38046"
chara = list('strain', "genotype", "tissue", "antibody", "extraction method")
mat = scrape_geo(GSE_id = gse, 
                 do_srr = TRUE, debug = F, key_str = c("Title", rep("Chara", length(chara))), 
                 key_idx = append(list(5), chara))#, "cell line", "chip_antibody", "donor_id", "donor_sex"))
head(mat)
mat

dim(mat)

mat = mat[grepl("ChIP", mat[,3]), ]
mat = mat[!grepl("Pax5", mat[,3]), ]

# mat = cbind(mat, rownames(mat))
dt = as.data.table(mat)
setnames(dt, paste0("V", seq(3+length(chara))), c("gse", "gsm", "title", unlist(chara)))
for(ch in chara){
  dt[[ch]] = sub(".+: ", "", dt[[ch]])
}


odir = paste0("fastqs_", gse, "_Trupke_Ikaros_proB")
dir.create(odir)
odir = normalizePath(odir)
setwd(odir)

dt[, c("mark", "stage", "mouse") := tstrsplit(title, " ", keep = c(1, 3, 4))]
dt[mouse == "B", "mouse" := tstrsplit(title, " ", keep = c(5))]
dt[mouse == "input", "mouse" := tstrsplit(title, " ", keep = c(6))]
dt[, mouse := sub("\\[", "", mouse)]
dt[, mouse := sub("\\]", "", mouse)]
dt[, mouse := sub(",", "", mouse)]
dt[, stage := sub('-', '', stage)]
dt[, rep := seq_len(.N), .(mouse)]
dt[, rep := paste0("rep", rep)]
dt[, name := paste("Trupke", stage, mark, mouse, rep, sep = "_")]

dt[, file_name := file.path(odir, paste0(name, ".fastq"))]

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



if(FALSE){
  pool_dir = file.path(odir, "combined")
  dir.create(pool_dir)
  
  files = dir(odir, pattern = "fastq$")
  pool_dt = data.table(file = files)
  pool_dt = pool_dt[, tstrsplit(file, "_"), .(file)]
  pool_dt[, pool_group := paste(V1, V2, V3, V4, sep = "_")]
  pool_dt[, final_group := paste(V1, V2, V3, sep = "_")]
  
  pdt = unique(pool_dt[, .(final_group, pool_group)])
  pdt[, rep := paste0("rep", seq_len(.N)), .(final_group)]
  
  pool_dt = merge(pool_dt, pdt[, .(pool_group, rep)], by = "pool_group")
  pool_dt[, final_file := paste0(pool_group, "_", rep, ".fastq")]
  todo = split(pool_dt$file, pool_dt$final_file)
  for(final_f in names(todo)){
    cmd = paste("cat", paste(todo[[final_f]], collapse = " "), ">", file.path(pool_dir, final_f))
    system(cmd)
  }
}
