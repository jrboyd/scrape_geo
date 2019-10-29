source("functions_scrape_geo.R")
library(data.table)
gse = "GSE55502" #will be combined with GSE19465
chara = list("cell type")
mat = scrape_geo(GSE_id = gse, 
                 do_srr = TRUE, debug = F, key_str = c("Title", rep("Chara", length(chara))), 
                 key_idx = append(list(5), chara))#, "cell line", "chip_antibody", "donor_id", "donor_sex"))
mat = cbind(mat, rownames(mat))
dt = as.data.table(mat)
colnames(dt) = c("description", unlist(chara), "srr", "gsm")

odir = paste0(gse, "_rodrigo_H9_fastqs2")
dir.create(odir)
setwd(odir)


library(data.table)
dt = as.data.table(mat)
for(ch in chara){
  dt[[ch]] = sub(".+: ", "", dt[[ch]])
}
dt[, cell := "H9norm"]
dt[grepl("PAX6", `cell type`), cell := "H9pax6"]
dt[, experiment := tstrsplit(description, "_", keep = 1)]
dt[, cycle := tstrsplit(description, "_", keep = 3)]
dt[is.na(cycle), cycle := "async"]
dt[, rep := paste0("rep", seq(.N)), .(cell, experiment, cycle)]
dt = dt[rep != "rep2"]
dt[, file_name := paste0(paste(sep = "_", cell, experiment, cycle, rep), ".", srr, ".", gsm, ".fastq")]
# dt = dt[grepl("RNA", description)]
# dt[is.na(condition), condition := "control"]
# dt[`cell subtype` == "Tamoxifen Resistance", `cell subtype` := "tamR"]
# dt[`cell subtype` == "RUNX2 DOX inducible model", `cell subtype` := "Runx2Dox"]
# dt[condition == "Full Medium + Doxycycline", condition := "wDox"]
# dt[condition == "Full Medium", condition := "noDox"]
# dt[, rep := paste0("R", seq(.N)), by = .(`cell subtype`, condition)]

# dt[, cell := sub("_.+", "", id)]
# dt[, rep := seq(.N), by = .(`cell`)]
# dt[, file_name := paste0(paste(cell, mark, rep, sep = '_'), ".",gsm, ".fastq")]
# dt[, file_name := sub("replicate", "rep", file_name)]

fastq_names = dt$file_name
srrs = dt$srr
for(i in seq_along(srrs)){
  jid = system(paste("qsub -v PATH=$PATH ~/scripts/fastq_dump_wrapper.sh", srrs[i], fastq_names[i]), intern = TRUE)
  jid = strsplit(jid, " ")[[1]][3]
  # system(paste("qsub -hold_jid", jid, "-v PATH=$PATH ~/scripts/split_PE_fastq.sh", fastq_names[i], 100), intern = TRUE)
}
qc_dt = dt[, .(file_name, cell, cycle, experiment, rep)]
qc_dt = qc_dt[cycle != "async"]
qc_dt[experiment == "Input", experiment := "input"]
# fwrite(qc_dt, "qc_config_hESC.csv")
# # for(i in which(grepl("HiC", fastq_names))){
# #   system(paste("qsub -v PATH=$PATH ~/scripts/split_PE_fastq.sh", fastq_names[i], 102), intern = TRUE)
# # }
# 
system("qstat")
# 
# read.table(fastq_names[1], nrows = 8, sep = "\n")
# 
# # for(i in seq_along(srrs)){
# #   jid = system(paste("qsub -v PATH=$PATH ~/scripts/split_PE_fastq.sh", file.path("cut_fastq", fastq_names[i]), 100), intern = TRUE)
# #   jid = strsplit(jid, " ")[[1]][3]
# #   print(jid)
# # }
# 
# dt_cd[, rep := paste0("R", seq(.N)), .(donor_id, chip_antibody)]
# out_dt = dt_cd[, .(file_name, sub("RO ", "CD34-", donor_id), chip_antibody, rep)]
top_dt = data.table(lines = c(paste0("in=", getwd()), 
                              paste0("out=", "/slipstream/galaxy/uploads/working/qc_framework/output_", gse, "_rodrigo_H9")))
qc_f = paste0("qc_config_", gse, ".csv")
fwrite(top_dt, qc_f, col.names = FALSE)
fwrite(qc_dt, qc_f, col.names = FALSE, append = TRUE)
# 
# setwd('..')
