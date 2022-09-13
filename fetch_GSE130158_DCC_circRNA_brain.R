

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

gse_hits = "GSE130158"


chara = list("tissue", "cell markers")
# Characteristics	

gse_mats = lapply(gse_hits, function(gse){
  mat = scrape_geo(GSE_id = gse,
                   do_srr = TRUE, debug = F, 
                   key_str = c("Title", rep("Chara", length(chara))), 
                   key_idx = append(list(5), chara))#, "cell line", "chip_antibody", "donor_id", "donor_sex"))
  mat
})




sapply(gse_mats, nrow)
# Specify GSE here

odir = paste0("fastqs_GSE130158_DCC_circRNA")
dir.create(odir)
odir = normalizePath(odir)
setwd(odir)

mat = as.matrix(rbindlist(lapply(gse_mats, as.data.frame)))

# recommend saving full matrix to avoid having to scrape again.
write.table(mat, file = "scrape_matrix.txt", sep = "\t")

# transform into cleaner format with better colnames
dt = as.data.table(mat)
setnames(dt, paste0("V", seq(3+length(chara))), c("gse", "gsm", "title", unlist(chara)))
for(ch in chara){
  dt[[ch]] = sub(".+: ", "", dt[[ch]])
}
# recommend saving this processed matrix too.
fwrite(dt, file = "scrape_matrix.processed.txt", sep = "\t")

# 1) bone marrow CD34+ cells: CD34+CD38neglinneg (HSC), CD34+CD45RA+CD38+CD10neg CD62Lhi linneg (LMPP), (CD34+ CD38+CD10+CD45RA+ linneg (CLP) and CD34+ CD38+CD19+linneg (BCP); 
# 2) thymic CD34+ cells: CD34+CD7negCD1anegCD4negCD8neg (Thy1), CD34+CD7+CD1anegCD4negCD8neg (Thy2), CD34+CD7+CD1a+CD4negCD8neg (Thy3); and 
# 3) thymic CD34neg cells: CD4+CD8+ (Thy4), CD3+CD4+CD8neg (Thy5), and CD3+CD4neg CD8+ (Thy6) (Supplemental Fig.1b)
# CD34+CD38neglinneg (HSC)
# CD34+CD45RA+CD38+CD10neg CD62Lhi linneg (LMPP)
# CD34+CD38+CD10+CD45RA+ linneg (CLP)
# CD34+CD38+CD19+linneg (BCP)
# CD34+CD7negCD1anegCD4negCD8neg (Thy1)
# CD34+CD7+CD1anegCD4negCD8neg (Thy2)
# CD34+CD7+CD1a+CD4negCD8neg (Thy3)
# CD4+CD8+ (Thy4)
# CD3+CD4+CD8neg (Thy5)

# Now we need to create parallel vectors of final fastq names (fastq_names) and SRR ids (srrs)
# Do any selection/filtering here.
# All useful elements should be carried forward to the fastq name including GSM/SRR


dt[, file_name := paste(paste0(title), gsm, srr, "fastq", sep = ".")]
dt = dt[(tissue %in% c("brain", "gut") & ! grepl("30", title))]

# These are the critical inputs for the actual fastq-dump call
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

