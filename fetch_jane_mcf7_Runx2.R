source("functions_scrape_geo.R")

gse = "GSE86538" #will be combined with GSE19465
chara = list("cell type", "cell subtype", "condition")
mat = scrape_geo(GSE_id = gse, 
                 # gsm_todo = "GSM537640", 
                 do_srr = TRUE, debug = F, key_str = c("Title", rep("Chara", length(chara))), 
                 key_idx = append(list(5), chara))#, "cell line", "chip_antibody", "donor_id", "donor_sex"))
mat = cbind(mat, rownames(mat))

save(mat, file = paste0(gse, ".save"))

odir = paste0(gse, "_fastqs")
dir.create(odir)
setwd(odir)
colnames(mat) = c("description", chara, "srr", "gsm")
library(data.table)
dt = as.data.table(mat)
for(ch in chara){
  dt[[ch]] = sub(".+: ", "", dt[[ch]])
}

dt = dt[grepl("RNA", description)]
dt[is.na(condition), condition := "control"]
dt[`cell subtype` == "Tamoxifen Resistance", `cell subtype` := "tamR"]
dt[`cell subtype` == "RUNX2 DOX inducible model", `cell subtype` := "Runx2Dox"]
dt[condition == "Full Medium + Doxycycline", condition := "wDox"]
dt[condition == "Full Medium", condition := "noDox"]
dt[, rep := paste0("R", seq(.N)), by = .(`cell subtype`, condition)]

dt[, file_name := paste0(paste(`cell type`, `cell subtype`, condition, rep, sep = '_'), ".fastq")]

fastq_names = dt$file_name
srrs = dt$srr
for(i in seq_along(srrs)){
  jid = system(paste("qsub -v PATH=$PATH ~/scripts/fastq_dump_wrapper.sh", srrs[i], fastq_names[i]), intern = TRUE)
  jid = strsplit(jid, " ")[[1]][3]
  # system(paste("qsub -hold_jid", jid, "-v PATH=$PATH ~/scripts/split_PE_fastq.sh", fastq_names[i], 100), intern = TRUE)
}

# for(i in which(grepl("HiC", fastq_names))){
#   system(paste("qsub -v PATH=$PATH ~/scripts/split_PE_fastq.sh", fastq_names[i], 102), intern = TRUE)
# }

system("qstat")

read.table(fastq_names[1], nrows = 8, sep = "\n")

# for(i in seq_along(srrs)){
#   jid = system(paste("qsub -v PATH=$PATH ~/scripts/split_PE_fastq.sh", file.path("cut_fastq", fastq_names[i]), 100), intern = TRUE)
#   jid = strsplit(jid, " ")[[1]][3]
#   print(jid)
# }

dt_cd[, rep := paste0("R", seq(.N)), .(donor_id, chip_antibody)]
out_dt = dt_cd[, .(file_name, sub("RO ", "CD34-", donor_id), chip_antibody, rep)]
top_dt = data.table(lines = c(paste0("in=", getwd()), paste0("out=", "/slipstream/galaxy/uploads/working/qc_framework/output_", gse)))
qc_f = paste0("qc_config_", gse, ".csv")
fwrite(top_dt, qc_f, col.names = FALSE)
fwrite(out_dt, qc_f, col.names = FALSE, append = TRUE)

setwd('..')
