readMSP = function(file,
                             threads = 4) {
  msp.data <- readr::read_lines(file)
  if (length(grep("BEGIN IONS", msp.data)) > 0) {
    msp.data <- msp.data[msp.data != ""]
    temp.idx1 <- grep("BEGIN IONS", msp.data)
    temp.idx2 <- grep("END IONS", msp.data)
    if (length(temp.idx2) < length(temp.idx1)) {
      temp.idx2 <- c(temp.idx2, length(msp.data))
    }
    
    temp.idx <-
      purrr::map2(.x = temp.idx1, temp.idx2, function(x, y) {
        c(x + 1, y - 1)
      })
    # future::plan(strategy = future::multisession, workers = threads)
    ms2_spec <- purrr::map(
      .x = temp.idx,
      .f = function(x) {
        temp_spec <- msp.data[x[1]:x[2]]
        temp_spec <- temp_spec
        spec_info <-
          temp_spec[stringr::str_detect(temp_spec, "[A-Za-z]")]
        spec <-
          temp_spec[!stringr::str_detect(temp_spec, "[A-Za-z]")]
        spec_info <- stringr::str_split(spec_info, "\\=") %>%
          do.call(rbind, .)
        mz <-
          as.numeric(spec_info[grep("MASS|MZ", spec_info[, 1]), 2])
        rt <-
          as.numeric(spec_info[grep("RT|RETETION", spec_info[, 1]), 2])
        spec_info <- c(mz = mz, rt = rt)
        
        spec <- purrr::map(
          .x = spec,
          .f = function(x) {
            stringr::str_split(x, " ")[[1]] %>% as.numeric()
          }
        ) %>%
          do.call(rbind, .)
        
        spec <-
          spec %>% as.matrix()
        
        colnames(spec) <- c("mz", "intensity")
        
        spec <- list(info = spec_info, spec = spec)
        spec
      }
    )
    return(ms2_spec)
  } else{
    # n.tot <- length(msp.data)
    n.null <- which(msp.data == '')
    
    temp.idx1 <- c(1, n.null[-length(n.null)])
    temp.idx2 <- n.null - 1
    
    temp.idx <- data.frame(temp.idx1, temp.idx2,
                           stringsAsFactors = FALSE)
    temp.idx <- apply(temp.idx, 1, list)
    
    temp.idx <- lapply(temp.idx, unlist)
    
    # n.spec <- which(grepl('^\\d', msp.data))
    # n.info <- seq(n.tot)[-c(n.spec, n.null)]
    
    # pbapply::pboptions(style = 1)
    info.spec =
      purrr::map(
        temp.idx,
        .f =
          function(idx) {
            temp.msp.data <- msp.data[idx[1]:idx[2]]
            
            temp.msp.data <-
              temp.msp.data[temp.msp.data != ""]
            info.idx <-
              grep("[A-Za-z]", temp.msp.data)
            temp.info <- temp.msp.data[info.idx]
            temp.info <-
              strsplit(temp.info, split = ":")
            temp.info <- do.call(rbind, temp.info)
            temp.info <- data.frame(temp.info,
                                    stringsAsFactors = FALSE)
            temp.info[, 2] <-
              stringr::str_trim(temp.info[, 2])
            colnames(temp.info) <-
              rownames(temp.info) <- NULL
            rownames(temp.info) <- temp.info[, 1]
            temp.info <- temp.info[, -1, drop = FALSE]
            
            temp.spec <- temp.msp.data[-info.idx]
            
            if (length(temp.spec) != 0) {
              if (length(grep(" ", temp.spec[1])) == 1) {
                temp.spec <- strsplit(temp.spec, split = ' ')
              }
              
              if (length(grep("\t", temp.spec[1])) == 1) {
                temp.spec <- strsplit(x = temp.spec, split = "\t")
              }
              
              temp.spec <- do.call(rbind, temp.spec)
              temp.spec <- data.frame(temp.spec,
                                      stringsAsFactors = FALSE)
              colnames(temp.spec) <-
                c('mz', 'intensity')
              rownames(temp.spec) <- NULL
              temp.spec$mz <-
                as.numeric(as.character(temp.spec$mz))
              temp.spec$intensity <-
                as.numeric(temp.spec$intensity)
              temp.spec <-
                temp.spec[temp.spec$intensity != 0,]
            } else{
              temp.spec <- NULL
            }
            
            list('info' = temp.info,
                 'spec' = temp.spec)
          }
      )
    
    mz.idx <- grep("[Mm][Zz]", rownames(info.spec[[1]][[1]]))
    rt.idx <-
      grep("Time$|TIME$|time$|^RT|^rt|^Rt", rownames(info.spec[[1]][[1]]))
    
    ##fix bug in msp data from metAnalyzer
    if (length(rt.idx) == 0) {
      cat(crayon::yellow("The msp data are from MetAnalyzer software.\n"))
      rt.idx <-
        grep("NAME|Name|name", rownames(info.spec[[1]][[1]]))
      ##rt.idx is the name of peak
      info.spec <- lapply(info.spec, function(x) {
        info <- x[[1]]
        mz <- as.numeric(info[mz.idx, 1])
        rt <- as.character(info[rt.idx, 1])
        info <- c(mz, rt)
        names(info) <- c("mz", "rt")
        x[[1]] <- info
        x
      })
    } else{
      info.spec <- purrr::map(info.spec, function(x) {
        info <- x[[1]]
        mz <- as.numeric(info[mz.idx, 1])
        rt <- as.numeric(info[rt.idx, 1])
        info2 = info[-c(mz.idx, rt.idx), , drop = FALSE]
        info1 <- c(mz, rt)
        info = c(info1, info2[, 1])
        names(info) <- c("mz", "rt", rownames(info2))
        x[[1]] <- info
        x
      })
    }
    
  }
  
  remove.idx <-
    which(unlist(lapply(info.spec, function(x)
      is.null(x[[2]]))))
  if (length(remove.idx) > 0) {
    info.spec <- info.spec[-remove.idx]
  }
  
  info.spec <- info.spec
}


WriteMSP = function(compID, content, filename) {
  sink(filename)
  info <- content$info
  spec <- content$spec
  head <- names(info)
  cat('compID: ', compID, '\n', sep = '')
  for (i in seq_along(head)){
    cat(head[i], ': ', info[i], '\n', sep = '')
  }
  cat('Num Peaks: ', nrow(spec), '\n', sep = '')
  for (nr in seq(nrow(spec))) {
    cat(paste(spec[nr, ], collapse = ' '), '\n', sep = '')
  }
  cat('\n')
  sink()
}

