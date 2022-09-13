

setwd("~/R/scrape_geo")
source("functions_scrape_geo.R")
library(data.table)
DUMP_SCRIPT = "scripts/fastq_dump_wrapper.sh"
DUMP_SCRIPT = normalizePath(DUMP_SCRIPT)
stopifnot(file.exists(DUMP_SCRIPT))
SPLIT_SCRIPT = "scripts/split_PE_fastq.sh"
SPLIT_SCRIPT = normalizePath(SPLIT_SCRIPT)
stopifnot(file.exists(SPLIT_SCRIPT))

json_content = jsonlite::read_json("ENCODE_MCF7_chipseq.json")
json_content
json_content

rec_search  = function(x, str){
  # message(str)
  lapply(seq_along(x), function(i){
    if(is.list(x[[i]])){
      lapply(x[[i]], rec_search, str = paste(str, ifelse(is.null(names(x)[i]), i, names(x)[i]), sep = ","))
    }else{
      if(grepl("GSE", x[[i]])){
        message("\n----start")
        message(str)
        message("\n")
        paste(str, x[[i]], sep = "\n")  
        browser()
      }else{
        list()
      }
      
    }
  })
}

gse_hits = lapply(json_content[["@graph"]], function(x){
  x$dbxrefs[grepl("GSE", x$dbxrefs)]
})
gse_hits = unique(unlist(gse_hits))
gse_hits = sub("GEO:", "", gse_hits)

chara = list("line", "antibody", "passage number", "lab")
gse_mats = lapply(gse_hits, function(gse){
  mat = scrape_geo(GSE_id = gse,
                   do_srr = TRUE, debug = F, 
                   key_str = c("Title", rep("Chara", length(chara))), 
                   key_idx = append(list(5), chara))#, "cell line", "chip_antibody", "donor_id", "donor_sex"))
  mat
})




sapply(gse_mats, nrow)
# Specify GSE here

# antibody: H3K9me2
# line: MCF-7
# biomaterial_type: immortalized cell line
# description: mammary epithelium, doubling time ~ 29hoursATCC HTB-22
# biosample encode accession: ENCBS789UPK (SAMN05733837)
# age: 69 year
# dev stage: adult
# Sex: female
# health state: breast cancer (adenocarcinoma)
# passage number: 146
# lab: Bradley Bernstein, Broad


odir = paste0("fastqs_ENCODE")
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

# Now we need to create parallel vectors of final fastq names (fastq_names) and SRR ids (srrs)
# Do any selection/filtering here.
# All useful elements should be carried forward to the fastq name including GSM/SRR

dt = rbind(
  dt[grepl("^H[0-9]K", antibody)],
  dt[grepl("Control", antibody)],
  dt[antibody %in% c("ESR1", "BRCA2", "CTCF", "GATA3", "FOXA1", "JUN")]
)

dt[, mark := antibody]
dt[mark == "Control", mark := "input"]
dt$cell = "MCF7"
# dt[, rep := paste0("rep", seq_len(.N)), .(cell, mark)]

tmp = unique(dt[, .(cell, mark, gsm)])
tmp[, .N, .(cell, mark)]
tmp[, rep := paste0("rep", seq_len(.N)), .(cell, mark)]
# dt$rep = NULL
dt = merge(dt, tmp, by = c("cell", "mark", "gsm"))

dt[, file_name := paste(paste(cell, mark, rep, sep = "_"), gsm, srr, "fastq", sep = ".")]

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

### AFTER FINISHED
files = dir(pattern = "fastq$")
dt_pool = data.table(file = files)
dt_pool[, c("name", 'gsm', "srr", "suff") := tstrsplit(file, "\\.") ]
dt_pool = dt_pool[, .(srr = paste(srr, collapse = "_"), to_pool = paste(file, collapse = " "), len = .N), .(name, gsm)]
dt_pool

pool_dir = "pooled"
dir.create(pool_dir)
dt_pool[, new_file := paste(name, gsm, srr, "fastq", sep = ".")]
i = 1
for(i in seq_len(nrow(dt_pool))){
  len = dt_pool$len[i]
  to_pool = dt_pool$to_pool[i]  
  new_file = file.path(pool_dir, dt_pool$new_file[i])
  if(len == 1){
    cmd = paste("mv", to_pool, new_file)  
  }else{
    cmd = paste("cat", to_pool, ">", new_file, "&& rm", to_pool)  
  }
  message(cmd)
  system(cmd)
}

to_del = dt[!lab %in% c("Michael Snyder, Stanford", "Bradley Bernstein, Broad")]$gsm %>% unique
sapply(to_del, function(del){
  f = dir(pool_dir, pattern = del, full.names = TRUE)
  file.remove(f)
})

files = dir(pool_dir, full.names = TRUE)
dt_qc = data.table(file = files)
dt_qc[, c("cell", "mark", "rep") := tstrsplit(basename(file), "[_\\.]", keep = 1:3)]
dt_qc
dt_qc[, file := basename(file)]
fwrite(dt_qc, file.path(pool_dir, "qc_ENCODE_MCF7.body.csv"), sep = ",", col.names = FALSE)
