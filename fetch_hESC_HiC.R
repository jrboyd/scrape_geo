source("functions_scrape_geo.R")

mat = scrape_geo(GSE_id = "GSE86821", do_srr = TRUE, debug = F, key_str = c("Title", "cell line"))
mat = cbind(mat, rownames(mat))
# write.table(mat, file = "GSE63018_srrs.csv", col.names = FALSE, quote = FALSE, sep = ",", row.names = FALSE)

mat[,2] = sub("Biological_Replicate", "rep", mat[,1])

dir.create("GSE86821_fastqs")
setwd("GSE86821_fastqs")
fastq_names = paste0("GSE86821.", gsub(" ", "_", mat[,2]), ".fastq")
srrs = mat[,3]
for(i in seq_along(srrs)){
  jid = system(paste("qsub -v PATH=$PATH ~/scripts/fastq_dump_wrapper.sh", srrs[i], fastq_names[i]), intern = TRUE)
  jid = strsplit(jid, " ")[[1]][3]
  system(paste("qsub -hold_jid", jid, "-v PATH=$PATH ~/scripts/split_PE_fastq.sh", fastq_names[i], 100), intern = TRUE)
}

for(i in which(grepl("HiC", fastq_names))){
  system(paste("qsub -v PATH=$PATH ~/scripts/split_PE_fastq.sh", fastq_names[i], 102), intern = TRUE)
}

system("qstat")

# for(i in seq_along(srrs)){
#   jid = system(paste("qsub -v PATH=$PATH ~/scripts/split_PE_fastq.sh", file.path("cut_fastq", fastq_names[i]), 100), intern = TRUE)
#   jid = strsplit(jid, " ")[[1]][3]
#   print(jid)
# }



setwd('..')
