source("functions_scrape_geo.R")

mat = scrape_geo(GSE_id = "GSE63018", do_srr = TRUE)
mat = cbind(mat, rownames(mat))
# write.table(mat, file = "GSE63018_srrs.csv", col.names = FALSE, quote = FALSE, sep = ",", row.names = FALSE)

dir.create("GSE63018_fastqs")
setwd("GSE63018_fastqs")
fastq_names = paste0("GSE63018.", gsub(" ", "_", mat[,1]), ".patient_", sub(".+: ", "", mat[,2]), ".fastq")
srrs = mat[,3]
# for(i in seq_along(srrs)){
#   jid = system(paste("qsub -v PATH=$PATH ~/scripts/fastq_dump_wrapper.sh", srrs[i], fastq_names[i]), intern = TRUE)
#   jid = strsplit(jid, " ")[[1]][3]
#   system(paste("qsub -hold_jid", jid, "-v PATH=$PATH ~/scripts/split_PE_fastq.sh", fastq_names[i], 100), intern = TRUE)
# }

for(i in seq_along(srrs)){
  jid = system(paste("qsub -v PATH=$PATH ~/scripts/split_PE_fastq.sh", file.path("cut_fastq", fastq_names[i]), 100), intern = TRUE)
  jid = strsplit(jid, " ")[[1]][3]
  print(jid)
}



setwd('..')
