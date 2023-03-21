# read in indexes for all samples -----
allindex <- function(dir = sgms, pattern = '*.raw$', named = TRUE){
  fn <- list.files(dir, pattern = pattern)
  listIndex <- lapply(file.path(dir, fn), function(p){readIndex(p)})
  if (named == TRUE) {
    names(listIndex) <- stringr::str_replace(fn, '\\.raw$', '')
  } else if (length(named) == legnth(fn)){
    names(listIndex) <- named
  } else {
    stop('named list if not the same length as ')
  }
  return(listIndex)
  # allindex <- listIndex
  # saveRDS(listIndex, file.path(local, 'rds/allindex.rds'))
}

# modify the tidymod table for peptide centric table that is compatible with downstream.
peptideCentric <- function(tbl){
  hs <- UniProt.ws::UniProt.ws(9606)
  geneSymbol <- UniProt.ws::select(x = hs, keys = tbl %>% tidyr::separate_rows(Proteins) %>% dplyr::pull(Proteins), columns = c('gene_primary'), keytype = 'UniProtKB')
  d <- dplyr::left_join(tbl, geneSymbol %>% dplyr::select(Entry, `Gene.Names..primary.`), by = c('Proteins' = 'Entry')) %>% 
    dplyr::mutate(uniqueID = paste(`Gene.Names..primary.`, type, `Positions within proteins`, sep = '_')) %>% 
    # dplyr::select(uniqueID, Proteins, `Modified sequence`, `MS/MS scan number`, `Best localization raw file`, `Charge.evidence`, `m/z`, `MS/MS m/z`, `Calibrated retention time`) %>% 
    tidyr::nest(data = -c(`Modified sequence`, `MS/MS scan number`)) %>% # `Modified sequence` is the unique identifier. 
    dplyr::mutate(uniqueID = purrr::map_chr(data, ~{.x %>% dplyr::pull(uniqueID) %>% paste0(collapse = '; ' )}), 
                  Proteins = purrr::map_chr(data, ~{.x %>% dplyr::pull(Proteins) %>% paste0(collapse = '; ' )}),
                  `Best localization raw file` = purrr::map_chr(data, ~{.x %>% dplyr::pull(`Best localization raw file`) %>% unique(.) %>% paste0(collapse = '; ' )}),
                  charge = purrr::map_dbl(data, ~{.x %>% dplyr::pull(`Charge.evidence`) %>% unique(.) %>% as.numeric(.)}),
                  mz = purrr::map_dbl(data, ~{.x %>% dplyr::pull(`m/z`) %>% unique(.) %>% as.numeric(.)}),
                  isolatemz = purrr::map_dbl(data, ~{.x %>% dplyr::pull(`MS/MS m/z`) %>% unique(.) %>% as.numeric(.)}),
                  rt = purrr::map(data, ~{.x %>% dplyr::pull(`Calibrated retention time`) %>% unique(.) %>% unlist(.)}),
                  rn = dplyr::row_number()) %>% 
    dplyr::select(-data) %>% 
    dplyr::mutate(pngLab = uniqueID, 
                  pngLab = ifelse(stringr::str_detect(pngLab, pattern = ";"), stringr::str_extract(pngLab, '^[^;]+;'), pngLab),
                  geneName = stringr::str_extract(pngLab, '^[^_]+_') %>% trimws(., whitespace = "_"),
                  rest = stringr::str_remove(pngLab, '^[^_]+_'),
                  site = stringr::str_extract(rest, '\\d+') %>% trimws(., whitespace = "_"),
                  mod = stringr::str_remove(rest, '\\d+') %>% trimws(., whitespace = ";") %>% trimws(., whitespace = "_"),
                  uniprot = ifelse(stringr::str_detect(Proteins, pattern = ";"), stringr::str_extract(Proteins, '^[^;]+;'), Proteins) %>% trimws(., whitespace = ";"))
  proteinName <- UniProt.ws::select(x = hs, keys = d  %>% dplyr::pull(uniprot), columns = c('protein_name'), keytype = 'UniProtKB')
  dplyr::left_join(d, proteinName, by = c('uniprot' = 'Entry'))
}

# finding corresponding MS1 for MS2 52162 and 96907 ------
scanInfo <- function(index, ms2scanNum, colname = NULL){
  ms2row <- index[index$scan == ms2scanNum,]
  if (is.null(colname)){
    return(ms2row)
  } else {
    stopifnot(colname %in% c('masterScan', 'rtinseconds', 'precursorMass', 'charge', 'monoisotopicMz'))
    # stopinnot(index$MSOrder == 'Ms') # enable query RT info from MS1
    return(ms2row[[colname]])
  }
}


# extract unmodified sequence -----
unmodSeq <- function(modseq_ori, maxlength = 5){
  subseq <- lapply(c("\\|[^\\|]+_", "\\|[^\\|]+\\|", "_[^\\|]+\\|"), function(pattern){
    s <- stringr::str_replace_all(modseq_ori, modname_pattern(bracket = TRUE), '|')
    stringr::str_extract_all(s, pattern = pattern)}) |> unlist() |> trimws(whitespace = '(_|\\|)') 
  unmodSeqs <- subseq[nchar(subseq) >= maxlength]
  l <- length(unmodSeqs)
  if (l > 0){
    unmodSeqs <- paste0('_', unmodSeqs, "_")
    names(unmodSeqs) <- paste0('unmod', seq_len(l))
    return(unmodSeqs)
  } 
}


# find scan number from the split search ----
scanNumSearching <- function(pc_spectbl_split, ms2snlist = FALSE, indextbl, raw = FALSE, mz = isolatemz, ppm = 5, rttol = 3.5){
  # ppm and rttol cutoff are chosen by comparing against maxquant scan number from Raw search 
  t <- pc_spectbl_split %>% 
    dplyr::mutate(scanNumList = purrr::pmap(list(indextbl, {{ mz }}, rt), ~{
      ..1[..1[['precursorMass']] > (..2-..2*ppm/1000000) & 
            ..1[['precursorMass']] < (..2+..2*ppm/1000000) &
            (..1[['rtinseconds']]/60) > (..3[1] - rttol) & 
            (..1[['rtinseconds']]/60) < (..3[2] + rttol), 'scan'] %>% .[!is.na(.)]})) 
  if (raw) {
    # if raw search results are used, you can also obtain a list of spectra for identification 
    t <- t %>% 
      dplyr::mutate(scanNumList = purrr::map2(`MS/MS scan number`, scanNumList, ~{c(.x, unlist(.y)) %>% unique(.) %>% unlist(.)}))
  }
  t <- t %>% 
    dplyr::mutate(l = lengths(scanNumList),
                  i = ceiling(l/2),
                  scanNum = purrr::map2_int(scanNumList, i, ~{.x[.y] %>% ifelse(length(.) == 0, NA, .) %>% as.integer(.)}))
  switch(as.character(ms2snlist), 
         'TRUE' = t %>% dplyr::pull(scanNumList), 
         'T' = t %>% dplyr::pull(scanNumList), 
         'FALSE' = t %>% dplyr::pull(scanNum),
         'F' = t %>% dplyr::pull(scanNum), 
         'table' = t)
}



# constand and modification table ----
unmodified = cbind(AA="-", mono=0.0, modname = 'unmodified',  desc="unmodified") # 0
mgo_MGH <- cbind(AA="R",mono=54.01057, modname = 'mgo_MGH', desc="mgo wo H2O, R only, H2OC2") # 1
mgo_CEA <- cbind(AA='R',mono=72.02113, modname = 'mgo_CEA', desc="mgo w H2O, R only, H4O2C3") # 2
mgo_CEL <- cbind(AA='K',mono=72.02113, modname = 'mgo_CEL', desc="mgo w H2O, K only, H4O2C3") # 3
`N-glyceroyl(K)` <- cbind(AA='K',mono=88.01604, modname = 'N-glyceroyl(K)', desc="mgo w H2O, K only, H4O2C3") # 4
`N-phosphoglyceroyl(K)` <- cbind(AA='K',mono=167.98237, modname = 'N-phosphoglyceroyl(K)', desc="mgo w H2O, K only, H4O2C3") # 5
`bp(Ywhc)` <- cbind(AA= c('Y,W,H,C'), mono=361.14601, modname = 'bp(Ywhc)', desc="bp(Ywhc)") # 6
`Acetyl (Protein N-term)` <- cbind(AA= 'N-term', mono=42.01056, modname = 'Acetyl (Protein N-term)', desc="Acetylation of the protein N-terminus") # 7
`Oxidation (M)` <- cbind(AA= 'M', mono=15.99491, modname = 'Oxidation (M)', desc="Oxidation") # 8
modtbl <- as.data.frame(rbind(unmodified, mgo_MGH, mgo_CEA, mgo_CEL, `N-glyceroyl(K)`, `N-phosphoglyceroyl(K)`, 
                              `bp(Ywhc)`, `Acetyl (Protein N-term)`, `Oxidation (M)`), stringsAsFactors = FALSE) # unmodified as 0,
modname_pattern <- function(modname = modtbl$modname, bracket = FALSE){
  if (bracket){
    modname <- paste0('(', modname, ')')
  }
  paste0(modname, collapse = '|') |> 
    stringr::str_replace_all(pattern = '\\(', replacement = '\\\\(') |> 
    stringr::str_replace_all(pattern = '\\)', replacement = '\\\\)') # modtype pattern make it a grep-able with ()bracket
}
H <- 1.007825
C <- 12.000000
O <- 15.994915
N <- 14.003074
P <- 30.973763

# function for fragmentation ---- 
HCD_Ion <- function(b, y){
  y2 = (y + H) / 2
  b2 = (b + H) / 2
  y3 = (y + 2 * H) / 3
  b3 = (b + 2 * H) / 3
  return(cbind(b = b, `b2+` = b2, `b3+` = b3,  y = y, `y2+` = y2, `y3+` = y3))
}


# extract sequence w\o modification-----
modoutput <- function(modseq_ori){
  trimws(modseq_ori, whitespace = '_') |> stringr::str_replace_all(modname_pattern(bracket = TRUE), '')
}

# obtain fragmentation pattern, only works for sequence with one modifications -----
fragmentIon <- function(modseq_ori){
  # taking care of N-acetylation bacause it is at the start before any sequence. 
  Nterm_Acetylation <- FALSE
  if (grepl('\\(Acetyl \\(Protein N-term\\)\\)', modseq_ori)){
    modseq_ori <- stringr::str_replace(modseq_ori, '\\(Acetyl \\(Protein N-term\\)\\)', '') 
    Nterm_Acetylation <- TRUE
  }
  modseq <- trimws(modseq_ori, whitespace = '_')
  allmod <- stringr::str_extract_all(modseq, modname_pattern()) |> unlist()
  modpos <- c()
  for (m in allmod){
    modpos <- append(modpos, stringr::str_locate(modseq, modname_pattern(bracket = TRUE))[1, 'start']-1)
    modseq <- stringr::str_replace(modseq, modname_pattern(bracket = TRUE), '')
  }
  modstring <- paste0(rep(0, stringr::str_length(modseq)), collapse = "")
  for (i in seq_along(modpos)){
    stringr::str_sub(modstring, start = modpos[i], end = modpos[i]) <- as.character(which(modtbl$modname %in% allmod[i])-1)
  }
  if (Nterm_Acetylation) {
    stringr::str_sub(modstring, start = 1, end = 1) <- as.character(which(modtbl$modname == 'Acetyl (Protein N-term)')-1)
  } 
  protViz::fragmentIon(modseq, FUN = HCD_Ion, modified = modstring, modification = modtbl$mono)[[1]]
}


# precursor charge number (side legend of MS1 plot ) or isotopicPattern ----
# sequence (w mod), isotope number calculated until, fragIon used to find Max peak around (fragmented) precursor 
# cha: the current precursor charge number which is going to fragment; maxcha: the max change number display on the side legend of 
precursorMz <- function(modseq_ori, mode = c("chargeNum", "isoPattern"), fragIon = NULL, cha, maxcha = cha, isoNumber = 4, named = TRUE){
  seq <- modoutput(modseq_ori) # sequence wo modification
  # matching existing pattern from modtbl what you see the original sequence. 
  modtype <- stringr::str_extract_all(modseq_ori, pattern = modname_pattern()) |> unlist()
  modsum <- modtbl$mono[modtbl$modname %in% modtype] |> as.numeric() |> sum()
  fw <- protViz:: parentIonMass(seq) + modsum
  if (mode == "chargeNum"){
    # charge number 
    clen <- seq_len(maxcha)
    chargeNum <- (fw + (clen-1) * H)/clen # fw is the parent ion mass!! 
    chargeLab <- c('[M+H]+', '[M+2H]+2', '[M+3H]+3', '[M+4H]+4', '[M+5H]+5','[M+6H]+6','[M+7H]+7')
    names(chargeNum) <- chargeLab[clen]
  }
  if (mode == "isoPattern"){
    # isotopic pattern of a particular precursor 
    stopifnot(!is.null(fragIon))
    iso <- seq_len(isoNumber)
    isoPattern <- (fw + (cha-1)*H)/cha + c(0, iso)/cha # fw is the parent ion mass, not formular weight!! 
    isoLab <- c('M', 'M+1', 'M+2', 'M+3', 'M+4', 'M+5', 'M+6', 'M+7')
    mzerror <- isoPattern - fragIon # indicate which precursor got fragmented
    isoLab[mzerror == min(abs(mzerror))] <- paste0(isoLab[mzerror == min(abs(mzerror))], '(F)')
    names(isoPattern) <- isoLab[c(0, iso) + 1]
  }
  # output
  output <- switch(mode, chargeNum = chargeNum, isoPattern = isoPattern)
  return(output)
}


# measured/actual isotope pattern 1) with strongest local signal intensity within the theoretical ppm window ----
measuredIsotopes <- function(locx, locy, theoretical, otype = c('x', 'y'),  ppm = 20){
  lapply(theoretical, function(n){
    xinRegion <- locx > n - n * ppm * 10^(-6) & locx < n + n * ppm * 10^(-6)
    if (any(xinRegion) && any(locy[xinRegion] != 0)){ 
      # exclude 2 cases 1) xinRegion all FALSE 2) locy[xinRegion] all zero (max will return the same length if all elements are 0) 
      switch(otype, x = locx[xinRegion & locy == max(locy[xinRegion])], y = locy[xinRegion & locy == max(locy[xinRegion])])
    }
  }) |> unlist() 
}


# replace rawrr:::plot.rawrrSpectrum to show precursor and its isotope -----
plot.rawrrSpectrum <- function (x, modseq_ori, cha, ppm = 20, relative = TRUE, centroid = FALSE, SN = FALSE, legend = TRUE, 
                                diagnostic = FALSE, mmz, isolatemz = NULL, ...){
  if (is.null(mmz)){
    xlim <- x$massRange
    ylim <- c(0, 1.2 * max(x$centroid.intensity,na.rm = TRUE))
  } else {
    upperlim <- switch(as.character(cha), '2' = 3.5, '3' = 2.5, '4' = 2, '5' = 2, '6' = 2) # calculate the upper lim of ms1 based on charge number
    xlim <- c(floor(mmz), floor(mmz) + upperlim) # expand x around monoisotopic peak
    # length(x$centroid.mZ); length(x$centroid.intensity)
    reducedCintensity <- x$centroid.intensity[which(x$centroid.mZ > xlim[1] & x$centroid.mZ < xlim[2])]
    ylim <- c(0, 1.2 * max(reducedCintensity, na.rm = TRUE)) # y is based on current window.
    # length(x$mZ); length(x$intensity)
    mZidx <- which(x$mZ > xlim[1] & x$mZ < xlim[2])
    ylim_notCentroid <- c(0, 1.19 * max(x$intensity[mZidx]))
  }
  stopifnot(is.rawrrSpectrum(x))
  if (centroid) {
    stopifnot(x$centroidStream)
    if (SN) {
      plot(x = x$centroid.mZ, y = x$centroid.intensity/x$noise, 
           type = "h", xlim = xlim, xlab = "Centroid m/z", 
           ylab = "Centroid Signal/Noise", frame.plot = FALSE, 
           ...)
    } else {
      plot(x = x$centroid.mZ, y = x$centroid.intensity, 
           type = "h", xlim = xlim, ylim = ylim, xlab = "Centroid m/z", 
           ylab = "Centroid Intensity", frame.plot = FALSE, 
           ...)
      if (all(c("charges", "resolutions") %in% names(x))) {
        n <- length(x$centroid.intensity)
        if (n > 10) 
          n <- 10
        i <- order(x$centroid.intensity, decreasing = TRUE)[seq_len(n)]
        text(x = x$centroid.mZ[i], y = x$centroid.intensity[i], 
             pos = 3, labels = paste(format(x$centroid.mZ[i], 
                                            nsmall = 4), "\nz = ", x$charges[i], "\nR = ", 
                                     x$resolutions[i]), cex = 0.5)
      }}
    } else { # centroid = FALSE
    if (relative) { # relative = TRUE
      plot(x = x$mZ[mZidx], y = x$intensity[mZidx]/max(x$intensity[mZidx]), type = "h", 
           xlim = xlim, xlab = "m/z", ylim = c(0, 1.15), ylab = "Relative Intensity", 
           frame.plot = FALSE)#, ...)
    } else {
      plot(x = x$mZ[mZidx], y = x$intensity[mZidx], type = "h", xlim = xlim, 
           xlab = "m/z", ylim = ylim_notCentroid, ylab = "Intensity", frame.plot = FALSE, 
           ...)
    }
    }
  
  # annotating precursor plot 
  theoretical <- precursorMz(modseq_ori, mode = 'isoPattern', fragIon = mmz, cha = cha) 
  measured_locx <- measuredIsotopes(x$mZ[mZidx], x$intensity[mZidx]/max(x$intensity[mZidx]), theoretical = theoretical, otype = 'x', ppm = ppm)
  measured_locy <- measuredIsotopes(x$mZ[mZidx], x$intensity[mZidx]/max(x$intensity[mZidx]), theoretical = theoretical, otype = 'y', ppm = ppm)
  axis(3, theoretical, names(theoretical), las = 2) # upper box: theoretical isotopic number
  axis(1, measured_locx, sprintf("%.4f", measured_locx), las = 2, cex.axis = 0.8) # lower box: max m/z within 20ppm of theoretical m/z
  points(measured_locx, measured_locy, col = "blue", pch = 22)
  chargeNum <- precursorMz(modseq_ori, mode = 'chargeNum', cha = cha)
  legend("right", sprintf("% 10.4f   %s", unname(chargeNum), names(chargeNum)), title = "Precursor Ions", bty = "n", cex = 0.65) # precursor charge number
  box()
  if (legend) {
    legend("topleft", paste(c("Scan#: ", "Scan Type: ", "RT [min]: ", 
                              #"Base peak mass [m/z]: ", "Base peak intensity: ", 
                              "TIC: "), c(x$title, x$scanType, round(x$rtinseconds/60, digits = 2), 
                                          # format(x$basePeak[1], nnsmall = 4), format(x$basePeak[2], scientific = TRUE), 
                                          format(x$TIC, scientific = TRUE))), bty = "n", cex = 0.5)
  }
  if (diagnostic) {
    legend("left", legend = paste(c("Injection time [ms]: ", 
                                    "Max. Injection time [ms]: ", "AGC target: ", "Resolution: "), 
                                  c(x$`Ion Injection Time (ms)`, x$`Max. Ion Time (ms)`, 
                                    format(x$`AGC Target`, scientific = TRUE), format(x$`FT Resolution`, 
                                                                                      scientific = TRUE))), bty = "n", cex = 0.5, 
           text.col = "grey")
  }
  invisible(x)
}


# modified protViz::psm to have 1) correctly formated fragmentation ion pattern -----
PSM <- function (sequence, spec, FUN = defaultIon, plot = TRUE, 
                 fi = fragmentIon(sequence, FUN = FUN)[[1]], fragmentIonError = 0.6){
  n <- nchar(sequence)
  pim <- fi$y[nrow(fi)]
  by.mZ <- numeric()
  by.label <- character()
  fi.names <- names(fi)
  for (i in 1:ncol(fi)) {
    by.mZ <- c(by.mZ, fi[, i])
    # by.label <- c(by.label, paste(fi.names[i], 1:n, sep = ""))
    prodIons <- rep(fi.names[i], n)
    stringr::str_sub(prodIons, 2,1) <- paste0(1:n,']')
    by.label <- c(by.label, trimws(prodIons, whitespace = ']'))
  }
  NN <- findNN_(q = by.mZ, vec = spec$mZ)
  mZ.error <- spec$mZ[NN] - by.mZ
  if (plot == TRUE) {
    plot(mZ.error ~ spec$mZ[NN], pch = 22, ylim = c(-5 * fragmentIonError, 5 * fragmentIonError))
    abline(h = fragmentIonError, col = "grey")
    abline(h = -fragmentIonError, col = "grey")
    abline(h = 0, col = "grey", lwd = 2)
    plot(mZ.error[mZ.error.idx <- order(mZ.error)], main = paste("Error of", 
                                                                 sequence, "(parent ion mass =", round(pim, 2), "Da)"), 
         ylim = c(-5 * fragmentIonError, 5 * fragmentIonError), 
         pch = 22, sub = paste("The error cut-off is", fragmentIonError, 
                               "Da (grey line)."))
    abline(h = fragmentIonError, col = "grey")
    abline(h = -fragmentIonError, col = "grey")
    abline(h = 0, col = "grey", lwd = 2)
    text(1:length(by.label), mZ.error[mZ.error.idx], by.label[mZ.error.idx], 
         cex = 0.75, pos = 3)
    hits = (abs(mZ.error) < fragmentIonError)
    nHits <- sum(hits)
    sumMZerror = round(sum(abs(mZ.error[hits])), 2)
    avgMZerror = round(sumMZerror/nHits, 2)
    cover = round(nHits/(nrow(fi) * ncol(fi)), 2)
    legend("topleft", paste(c("nHits", "sumMZerror", "avgMZerror", 
                              "cover"), as.character(c(nHits, sumMZerror, avgMZerror, 
                                                       cover)), sep = "="))
  }
  return(list(mZ.Da.error = mZ.error, mZ.ppm.error = 1e+06 * 
                mZ.error/by.mZ, idx = NN, label = by.label, score = -1, 
              sequence = sequence, fragmentIon = fi))
}


# modified peakplot function 1)adding ppm option -----
# oriSeq = modseq_ori 
# spec = prodIon
# itol = 5
# unit = 'ppm'
# fi = fragmentIon(oriSeq)
# peptideSequence = modoutput(oriSeq)
# sub = paste(trimws(oriSeq, whitespace = "_"), spec$title,sep = " / ")
# pattern.abc = "[abc].*"; pattern.xyz = "[xyz].*";ion.axes = TRUE
peakPlot <- function (oriSeq, spec, peptideSequence = modoutput(oriSeq), fi = fragmentIon(oriSeq),
                      FUN = defaultIon, sub = paste(trimws(oriSeq, whitespace = "_"), spec$title,sep = " / "), 
                      itol = 20, unit = c('ppm', 'da'), pattern.abc = "[abc].*", pattern.xyz = "[xyz].*", 
                      ion.axes = TRUE, ...) {
  m <- PSM(peptideSequence, spec, FUN, fi = fi, plot = FALSE)
  max.intensity <- max(spec$intensity, na.rm = TRUE)
  plot(spec$mZ, spec$intensity, xlab = "m/z", ylab = "intensity", 
       type = "h", sub = sub, axes = "F", ...)
  if (unit == 'ppm'){
    LABEL.abc <- (abs(m$mZ.ppm.error) < itol) & (regexpr(pattern.abc, m$label) > 0)
    LABEL.xyz <- (abs(m$mZ.ppm.error) < itol) & (regexpr(pattern.xyz, m$label) > 0)
  } else if (unit == 'da'){
    LABEL.abc <- (abs(m$mZ.Da.error) < itol) & (regexpr(pattern.abc, m$label) > 0)
    LABEL.xyz <- (abs(m$mZ.Da.error) < itol) & (regexpr(pattern.xyz, m$label) > 0)
  } else{
    stop('unkown unit')
  }
  if (ion.axes) {
    if (length(m$idx[LABEL.abc]) > 0) {
      axis(1, spec$mZ[m$idx[LABEL.abc]], m$label[LABEL.abc], 
           las = 2)
    }
    axis(2)
    if (length(m$idx[LABEL.xyz]) > 0) {
      axis(3, spec$mZ[m$idx[LABEL.xyz]], m$label[LABEL.xyz], 
           col.axis = "blue", las = 2)
    }
  } else {
    axis(1)
    axis(2)
    a.at <- spec$mZ[m$idx[LABEL.abc | LABEL.xyz]]
    a.label <- m$label[LABEL.abc | LABEL.xyz]
    if (length(a.at) > 0) {
      axis(3, a.at, a.label, col.axis = "black", las = 2)
    } else {
      print("WARNING")
      print(a.at)
      print(a.label)
    }
  }
  box()
  axis(4, seq(0, max.intensity, length = 6), seq(0, 100, length = 6))
  # protViz:::.peakplot.label(spec = spec, match = m, itol = itol, pch = 22)
  # peakplot.label(spec = spec, m =m, itol = itol, unit = unit, pch = 22)
  points(spec$mZ[m$idx[LABEL.abc]], spec$intensity[m$idx[LABEL.abc]], col = "black", pch = 22, ...)
  points(spec$mZ[m$idx[LABEL.abc]], spec$intensity[m$idx[LABEL.abc]], col = "black", type = "h")
  points(spec$mZ[m$idx[LABEL.xyz]], spec$intensity[m$idx[LABEL.xyz]], col = "blue", pch = 22, ...)
  points(spec$mZ[m$idx[LABEL.xyz]], spec$intensity[m$idx[LABEL.xyz]], col = "blue", type = "h")
  sortedFragmentIonsTable <- data.frame(label = c(m$label[LABEL.abc],m$label[LABEL.xyz]), 
                                        mass = c(spec$mZ[m$idx[LABEL.abc]],spec$mZ[m$idx[LABEL.xyz]]))
  if (nrow(sortedFragmentIonsTable) > 0) {
    sortedFragmentIonsTable <- sortedFragmentIonsTable[order(sortedFragmentIonsTable$mass), ]
    legend("right", sprintf("% 10.4f   %s", sortedFragmentIonsTable$mass, 
                            sortedFragmentIonsTable$label), title = "Fragment Ions", 
           bty = "n", cex = 0.65)
  }
  # return(m)
}

# zoomed-in the chromatogram bound by rtmin-offset and rtmin+offset. exclude zero intensities -----
# (when plot using `elimiate0 = TRUE` for aesthetic purpose; when calculating auc, using FALSE. TODO check with xcalibur integration)
chromatoZoomed <- function(x, rtmin, offset = 10, eliminate0 = TRUE){
  for (i in seq_len(length(x))){
    if(eliminate0){
      idx <- which(x[[i]]$times > rtmin - offset & x[[i]]$times < rtmin + offset & x[[i]]$intensities > 0)
    } else {
      idx <- which(x[[i]]$times > rtmin - offset & x[[i]]$times < rtmin + offset)
    }
    x[[i]]$times <- x[[i]]$times[idx]
    x[[i]]$intensities <- x[[i]]$intensities[idx]
  }
  return(x)
}

# integrate AUC using trapezoid method between a manually chosen "integalab" -----
auc <- function(x, integalab){
  AUC <- vector(mode = 'double', length = length(x))
  for (i in seq_len(length(x))){
    inten_InArea <- x[[i]]$intensities[x[[i]]$times >= integalab[1] & x[[i]]$times <= integalab[2]]
    time_InArea <- x[[i]]$times[x[[i]]$times >= integalab[1] & x[[i]]$times <= integalab[2]]
    AUC[i] <- lapply(seq_len(length(inten_InArea) - 1), function(n){
      (inten_InArea[n] + inten_InArea[n+1]) * (time_InArea[n+1] - time_InArea[n]) / 2 
    })|> unlist() |> sum()
  }
  return(AUC)
}

# replace rawrr:::plot.rawrrChromatogram  -----
plot.rawrrChromatogramSet <- function (x_ori, rtmin, offset = 10, integalab, wAUC = FALSE, diagnostic = FALSE, ...){
  if(!is.null(rtmin)){
    x <- chromatoZoomed(x = x_ori, rtmin = rtmin, offset = offset, eliminate0 = TRUE)
  }
  stopifnot(attr(x, "class") == "rawrrChromatogramSet")
  if (attr(x, "type") == "xic") {
    plot(0, 0, type = "n", 
         xlim = range(unlist(lapply(x, function(o) {o$times}))), 
         ylim = range(unlist(lapply(x, function(o) {o$intensities}))), 
         frame.plot = FALSE, xlab = "Retention Time [min]",
         ylab = "Intensities", ...)
    cm <- hcl.colors(length(x), "Set 2")
    mapply(function(o, co) {
      lines(o$times, o$intensities, col = co)
    }, x, cm)
    if (wAUC){
      monoisotopicMass <- lapply(x, function(o) {o$mass}) |> unlist()
      AUC <- chromatoZoomed(x = x_ori, rtmin = rtmin, offset = offset, eliminate0 = FALSE) |> auc(integalab = integalab)
      legend("topleft", bg = '#FFFFFF80', legend = sprintf("% 8.4f   %.0f", monoisotopicMass, AUC), 
             col = cm, pch = 16, title = "AUC bounded by Pink", cex = 0.75)
      ymax <- max(lapply(x, function(n) n$intensities) |> unlist(), na.rm = TRUE)
      polygonx <- c(integalab, integalab[2], integalab[1])
      polygony <- c(0,0,ymax,ymax)
      polygon(polygonx, polygony, col = '#FFC0CB40', border = '#FFFFFF00')
    } else {
      legend("topleft", bg = '#FFFFFF80', legend = as.character(sapply(x, function(o) {o$mass})), 
             col = cm, pch = 16, title = "target mass [m/z]", cex = 0.75)
    }
    if (diagnostic) {
      legend("topright", 
             legend = paste(c("File: ", "Filter: ","Type: ", "Tolerance: "), 
                            c(basename(attr(x,"file")), attr(x, "filter"), attr(x, "type"), attr(x, "tol"))), 
             bty = "n", cex = 0.75, text.col = "black")
    }
  }
  invisible(x)
}


# not work for  modification peptides because mod means all of that kind will be modified -----
# find precursor mass, p----
# mod = 'mgo_MGH';mtbl = modtbl
# modTo26letter <- function(mod, mtbl = modtbl){
#   
#   modaa <- mtbl$AA[mtbl$modname == mod]
#   modmass <- mtbl$mono[mtbl$modname == mod] |> as.numeric()
#   data(AA)
#   AA <- dplyr::left_join(data.frame(LETTERS), AA, by = c('LETTERS' = 'letter1'))
#   mmm <- AA$Monoisotopic+(AA$LETTERS == modaa)*modmass
#   replace(mmm, is.na(mmm), 0L)
# }
# modTo26letter('mgo_MGH')


# (not used) isotope pattern extracted from index (find precursor and charge from index, not the actual peak)----
# extractedIsotopes <- function(index, ms2scanNum, isoNumber = 4){
#   mmz <- scanInfo(index, ms2scanNum, 'monoisotopicMz')
#   cha <- scanInfo(index, ms2scanNum, 'charge')
#   isotopicIons <- mmz + c(0, seq_len(4))/cha
#   return(isotopicIons)
# }
# measured <- extractedIsotopes(high2Index, ms2_scanNum) # 535.9305 536.2638 536.5972 536.9305 537.2638
