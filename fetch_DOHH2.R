source("functions_scrape_geo.R")
library(data.table)
gse = "DOHH2" #will be combined with GSE19465
chara = list("antibody", "line")
mat1 = scrape_geo(GSE_id = "GSE86744", 
                 do_srr = TRUE, debug = FALSE, 
                 key_str = c("Title", rep("Chara", length(chara))), 
                 key_idx = append(list(5), chara))

mat2 = scrape_geo(GSE_id = "GSE86701", 
                 do_srr = TRUE, debug = FALSE, 
                 key_str = c("Title", rep("Chara", length(chara))), 
                 key_idx = append(list(5), chara))

mat3 = scrape_geo(GSE_id = "GSE86678", 
                 do_srr = TRUE, debug = FALSE, 
                 key_str = c("Title", rep("Chara", length(chara))), 
                 key_idx = append(list(5), chara))

mat = rbind(mat1, mat2, mat3)
dt = as.data.table(mat)
colnames(dt) = c("gse", "gsm", "title", unlist(chara), "srr")

odir = paste0(gse, "_DOHH2_ChIPseq")
dir.create(odir)


save(mat, dt, file = file.path(odir, "scrape_result.save"))

setwd(odir)
library(data.table)

dt[, mark := sub(".+: ", "", antibody)]
dt[mark == "Control", mark := "input"]
dt[grepl("me3", mark), mark := toupper(mark)]
dt[, cell := sub(".+: ", "", line)]

dt[, rep := paste0("R", as.numeric(factor(gsm))), .(mark)]

dt[, file_name := paste0(cell, "_", mark, "_", rep, ".", seq_len(.N),  ".fastq")]

fwrite(dt, file = file.path(odir, "scrape_result.csv"))

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
