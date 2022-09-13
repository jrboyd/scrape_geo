

setwd("~/R/scrape_geo")
source("functions_scrape_geo.R")
library(data.table)
# DUMP_SCRIPT = "scripts/fastq_dump_wrapper.sh"
DUMP_SCRIPT = "~/R/scrape_geo/scripts/fastq_dump_wrapper.PE.sh"
DUMP_SCRIPT = normalizePath(DUMP_SCRIPT)
stopifnot(file.exists(DUMP_SCRIPT))
SPLIT_SCRIPT = "scripts/split_PE_fastq.sh"
SPLIT_SCRIPT = normalizePath(SPLIT_SCRIPT)
stopifnot(file.exists(SPLIT_SCRIPT))

gse = "GSE141139"


chara = list("cancer type", "cell type")
# Characteristics	

mat = scrape_geo(GSE_id = gse,
                 do_srr = TRUE, debug = F, 
                 key_str = c("Title", rep("Chara", length(chara))), 
                 key_idx = append(list(5), chara))#, "cell line", "antibody", "replicate"))
mat

is_hic = grepl("Hi-C", mat[,3])
is_chip = !grepl("Hi-C", mat[,3]) & !grepl("RNA", mat[,3])

mat.hic = mat[is_hic,]
mat.chip = mat[is_chip,]
# Specify GSE here

odir = paste0("fastqs_GSE141139_Kang20_U2OS_and_RPE_ChIP_and_HiC")
dir.create(odir)
odir = normalizePath(odir)
setwd(odir)

write.table(mat, file = "scrape_matrix.txt", sep = "\t")

# transform into cleaner format with better colnames
dt = as.data.table(rbind(mat.hic, mat.chip))
setnames(dt, paste0("V", seq(3+length(chara))), c("gse", "gsm", "title", unlist(chara)))
for(ch in chara){
  dt[[ch]] = sub(".+: ", "", dt[[ch]])
}
# recommend saving this processed matrix too.
fwrite(dt, file = "scrape_matrix.processed.txt", sep = "\t")

dt[grepl("Hi-C", title), title := paste0("U2OS_", title)]
dt$title
dt$title = sub("-Hi-C-", "_HIC_", dt$title)
dt$title = gsub(" ", "_", dt$title)
dt$title = sub("-ChIP-seq-", "_chipseq_", dt$title)
dt$title = gsub("-", "_", dt$title)
dt$title = sub("Input", "input", dt$title)
dt$title = sub("/", "-", dt$title)

dt[, file_name := paste(title, gsm, srr, "fastq", sep = ".")]

# These are the critical inputs for the actual fastq-dump call
fastq_names = dt$file_name
srrs = dt$srr

for(i in seq_along(srrs)){
  cmd = paste("qsub -N fastq-dump -v PATH=$PATH -o", odir, "-e", odir, DUMP_SCRIPT, srrs[i], fastq_names[i])
  jid = system(cmd, intern = TRUE)
  jid = strsplit(jid, " ")[[1]][3]
  message(jid)  
}
