library(magrittr)
get_right = function(key, txt){
  split = sapply(strsplit(txt, key), function(x)paste(x[-1], collapse = key))
  return(split)
}
get_left = function(key, txt){
  split = sapply(strsplit(txt, key), function(x)x[1])
  return(split)
}
get_mid = function(key_left, key_right, txt){
  right_split = get_right(key_left, txt)
  mid_split = get_left(key_right, right_split)
  return(mid_split)
}

get_attrib = function(gsm_lines, key, w = NULL, debug = FALSE){
  k = which(grepl(key, gsm_lines))
  attrib_line = gsm_lines[k+1]
  split = strsplit(attrib_line, split = "[\\\"><]")[[1]]
  if(debug) print(split)
  if(is.null(w)){
    split[which(split == "/td") - 1]  
  }else if(is.character(w)){
    split[which(grepl(w, split))[1]] 
  }else{
    split[w]  
  }
}

scrape_gse = function(
  gsms
){
  base_url = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="
  gse = (parallel::mclapply(gsms, function(g){
    url = paste0(base_url, g)
    con = curl::curl(url = url)
    gsm_lines = readLines(con)
    close(con)
    gse_line = gsm_lines[grepl("GSE", gsm_lines)]
    gse = regmatches(gse_line, regexpr("GSE[0-9]+", gse_line))
    gse
  }))
  names(gse) = gsms
  gse
}

readLines.recursive = function(url, wait = .5){
  success = FALSE
  tryCatch({
    con = curl::curl(url = url)
    lines = readLines(con)  
    close(con)
    success = TRUE
  }, error = function(e){
    message(e)
    
  })  
  if(!success){
    lines = readLines.recursive(url, wait)
  }
  lines
}

#' Title
#'
#' @return
#' @export
#' @import curl
#' @examples
scrape_geo = function(GSE_id, 
                      gsm_override = NULL,
                      gsm_todo = NULL,
                      gsm_skip = NULL,
                      do_srr = TRUE, 
                      debug = FALSE, 
                      key_str = c("Source name", "Chara", "Description"), 
                      key_idx = as.list(rep(5, length(key_str)))){
  stopifnot(length(key_str) == length(key_idx))
  
  
  
  if(!is.null(gsm_override)){
    gsms = gsm_override
  }else{
    gse_url = paste0("http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", GSE_id)
    gse_lines = readLines.recursive(gse_url)
    
    keep = which(grepl("GSM", gse_lines))
    gsm_lines = gse_lines[keep]
    # gsms = get_mid(key_left = "geoaxema_recenter)\">", key_right = "</a></t", gsm_lines)
    if(!is.null(gsm_todo)){
      gsms = gsm_todo
    }else{
      gsms = sapply(gsm_lines, function(x){
        m = regexpr(pattern = "GSM[0-9]+", text = x)
        regmatches(x = x, m = m)
      })
    }
    gsms = setdiff(gsms, gsm_skip)
  }
  names(gsms) = NULL
  
  mat = matrix("", nrow = 0, ncol = 2+length(key_str) + sum(do_srr))
  # rownames(mat) = gsms
  base_url = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="
  for(g in gsms){
    print(g)
    url = paste0(base_url, g)
    message(url)
    
    
    

    gsm_lines = readLines.recursive(url)
    
    if(debug) browser()
    
    vals = sapply(seq_along(key_str), function(ks){
      get_attrib(gsm_lines, key_str[ks], key_idx[[ks]], debug = debug)
    })
    if(do_srr){
      srx_line = gsm_lines[grepl("href", gsm_lines) & grepl("SRX", gsm_lines)][1]
      if(is.na(srx_line)){
        srr = "na"
      }else{
        m = regexpr(pattern = "http.+term=SRX[0-9]+", text = srx_line)
        srx_url = regmatches(x = srx_line, m = m)
        
        srx_lines = readLines.recursive(srx_url)
        m = gregexpr(pattern = "SRR[0-9]+", text = srx_lines)
        srr = regmatches(x = srx_lines, m = m) %>% unlist %>% unique
      }
      newL = cbind(matrix(rep(c(GSE_id, g, vals), length(srr)), nrow = length(srr), byrow = TRUE), srr)
      print(newL)
      mat = rbind(mat, newL)
      #mat[g, ] = newL
    }else{
      mat = rbind(mat, vals)
      # mat[g, ] = vals
    }
    #   close(srx_con)
    #   
    #   m = regexpr(pattern = "SRX[0-9].+", text = ftp_line)
    #   ftp = sub("\"", "", regmatches(x = ftp_line, m = m))
    #   ftp_con = curl::curl(url = ftp)
    #   readLines(ftp_con)
  }
  return(mat)
}


# mat = scrape_geo("GSE74716")
# mat = mat[9:12, 2:1]
# mat[,2] = c("MCF10A_input_R1", 
#             "MCF10A_BRG1_R1",
#             "MCF10A_input_R2",
#             "MCF10A_BRG1_R2")
# write.csv(mat, file = "GSE74716_srrs.csv")
### then: 
### qsub -cwd ~/../home/joeboyd/scripts/fetch_srr.sh ~/../home/joeboyd/R/scrape_geo/GSE74716_srrs.csv

