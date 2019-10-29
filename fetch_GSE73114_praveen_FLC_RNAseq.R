source("functions_scrape_geo.R")

GSE = "GSE73114"

mat = scrape_geo(GSE_id = GSE, do_srr = TRUE)
mat = cbind(mat, rownames(mat))
# write.table(mat, file = paste0(GSE, "_srrs.csv"), col.names = FALSE, quote = FALSE, sep = ",", row.names = FALSE)

odir = paste0(GSE, "_fastqs")
dir.create(odir)
setwd(odir)

fastq_names = paste0(GSE, ".", gsub(" ", "_", mat[,3]), ".", mat[,5], ".fastq")
srrs = mat[,4]
# for(i in seq_along(srrs)){
#   jid = system(paste("qsub -v PATH=$PATH ~/scripts/fastq_dump_wrapper.sh", srrs[i], fastq_names[i]), intern = TRUE)
#   jid = strsplit(jid, " ")[[1]][3]
#   jid2 = system(paste("qsub -hold_jid", jid, "-v PATH=$PATH ~/scripts/split_PE_fastq.sh", fastq_names[i], 96), intern = TRUE)
# }
vi 


setwd('..')