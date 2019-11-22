source("functions_scrape_geo.R")
chara = list("line", "cell_type", "donor_sex", "donor_age", "donor_id", "chip_antibody")

gse = "CD34" #will be combined with GSE19465

# skip = "GSM537640"
mat1 = scrape_geo(GSE_id = "GSE17312", 
                 # gsm_todo = "GSM537640", 
                 do_srr = TRUE, debug = F, key_str = c("Title", rep("Chara", length(chara))), 
                 key_idx = append(list(5), chara))#, "cell line", "chip_antibody", "donor_id", "donor_sex"))
# mat = cbind(mat, rownames(mat))

mat2 = scrape_geo(GSE_id = "GSE19465", 
                 # gsm_todo = "GSM537640", 
                 do_srr = TRUE, debug = F, key_str = c("Title", rep("Chara", length(chara))), 
                 key_idx = append(list(5), chara))#
# mat2 = cbind(mat2, rownames(mat2))
mat = rbind(mat1, mat2)

save(mat, file = "GSE17312_and_GSE19465
     .save")

patient_mat = mat[grepl("(01517)|(01536)|(01549)|(01562)", mat[,8]),]

mat[grepl("H3K4me3.+Primary", mat[,3]),]


mat_k27 = mat[grepl("H3K27me3.+Mobilized CD34", mat[,3]),]
mat_k27[grepl("01517", mat_k27[,8]),]
mat_k27[grepl("01536", mat_k27[,8]),]
mat_k27[grepl("01549", mat_k27[,8]),]
mat[grepl("WCE.+Mobilized CD34", mat[,3]),]
mat[grepl("WCE", mat[,3]),]
#STOPPED HERE
# write.table(mat, file = "GSE63018_srrs.csv", col.names = FALSE, quote = FALSE, sep = ",", row.names = FALSE)
# mat[,2] = sub("Biological_Replicate", "rep", mat[,1])

odir = paste0(gse, "_fastqs")
dir.create(odir)
setwd(odir)
colnames(mat) = c("description", chara, "srr", "gsm")
library(data.table)
dt = as.data.table(mat)
for(ch in chara){
  dt[[ch]] = sub(".+: ", "", dt[[ch]])
}
dt[, chip_antibody := sub("Me", "me", chip_antibody)]
dt[, chip_antibody := sub("Ac", "ac", chip_antibody)]
dt[grepl("WCE", description), chip_antibody := "input"]

chk = dcast(dt[grepl("CD34", description) & srr != "na"][order(donor_id)], "donor_id~chip_antibody")
complete = chk[H3K27me3 > 0 & H3K4me3 > 0 & input > 0]$donor_id

dt_cd = dt[donor_id %in% complete & srr != "na"][chip_antibody %in% c("H3K27me3", "H3K4me3", "input")][order(chip_antibody)][order(donor_id)]
dt_cd[, file_name := paste0(paste("CD34", sub("RO ", "", donor_id), chip_antibody, gsm, srr, sep = '_'), ".fastq")]

gse = scrape_gse(dt_cd$gsm)
gse = sapply(gse, function(x)x[1])
gse = cbind(gsm = names(gse), gse)
tmp = dt_cd[, .(donor_id, chip_antibody, gsm, srr)]
tmp = merge(tmp, gse)
tmp = tmp[, .(donor_id, chip_antibody, gse, gsm, srr)]
fwrite(tmp, "CD34_supplemental.txt")

fastq_names = dt_cd$file_name
srrs = dt_cd$srr
for(i in seq_along(srrs)){
  jid = system(paste("qsub -v PATH=$PATH ~/scripts/fastq_dump_wrapper.sh", srrs[i], fastq_names[i]), intern = TRUE)
  jid = strsplit(jid, " ")[[1]][3]
  # system(paste("qsub -hold_jid", jid, "-v PATH=$PATH ~/scripts/split_PE_fastq.sh", fastq_names[i], 100), intern = TRUE)
}

# for(i in which(grepl("HiC", fastq_names))){
#   system(paste("qsub -v PATH=$PATH ~/scripts/split_PE_fastq.sh", fastq_names[i], 102), intern = TRUE)
# }

system("qstat")

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
