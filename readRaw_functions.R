# finding corresponding MS1 for MS2 52162 and 96907 ------
scanInfo <- function(index, ms2scanNum, colname = c('masterScan', 'rtinseconds', 'precursorMass', 'charge', 'monoisotopicMz')){
  ms2row <- index[index$scan == ms2scanNum,]
  stopifnot(colname %in% c('masterScan', 'rtinseconds', 'precursorMass', 'charge', 'monoisotopicMz'))
  # stopinnot(index$MSOrder == 'Ms') # enable query RT info from MS1
  ms2row[[colname]]
}

# modification table ----
unmodified = cbind(AA="-", mono=0.0, modname = 'unmodified',  desc="unmodified") # 0
mgo_MGH <- cbind(AA="R",mono=54.01057, modname = 'mgo_MGH', desc="mgo wo H2O, R only, H2OC2") # 1
mgo_CEA <- cbind(AA='R',mono=72.02113, modname = 'mgo_CEA', desc="mgo w H2O, R only, H4O2C3") # 2
mgo_CEL <- cbind(AA='K',mono=72.02113, modname = 'mgo_CEL', desc="mgo w H2O, K only, H4O2C3") # 3
`N-glyceroyl(K)` <- cbind(AA='K',mono=88.01604, modname = 'N-glyceroyl(K)', desc="mgo w H2O, K only, H4O2C3") # 4
`N-phosphoglyceroyl(K)` <- cbind(AA='K',mono=167.98237, modname = 'N-phosphoglyceroyl(K)', desc="mgo w H2O, K only, H4O2C3") # 5
`bp(Ywhc)` <- cbind(AA= c('Y,W,H,C'), mono=361.14601, modname = 'bp(Ywhc)', desc="bp(Ywhc)") # 6
`Acetyl (Protein N-term)` <- cbind(AA= 'N-term', mono=42.01056, modname = 'Acetyl (Protein N-term)', desc="Acetylation of the protein N-terminus") # 7
`Oxidation (M)` <- cbind(AA= 'M', mono=15.99491, modname = 'Oxidation (M)', desc="Oxidation") # 8
modtbl <- as.data.frame(rbind(unmodified, mgo_MGH, mgo_CEA, mgo_CEL, `N-glyceroyl(K)`, `N-phosphoglyceroyl(K)`, `bp(Ywhc)`, `Acetyl (Protein N-term)`, `Oxidation (M)`), stringsAsFactors = FALSE) # unmodified as 0,

# extract mod sequence -----
modoutput <- function(modseq_ori){
  # the 2nd trimws will be useful in 
  trimws(modseq_ori, whitespace = '_') |> stringr::str_replace_all('\\([^\\)]+\\)', '') |> trimws(whitespace = '\\)') 
}

# function for fragmentation ---- 
HCD_Ion <- function(b, y){
  C <- 12.000000
  H <- 1.007825
  O <- 15.994915
  N <- 14.003074
  P <- 30.973763
  y2 = (y + H) / 2
  b2 = (b + H) / 2
  return(cbind(b = b, `b2+` = b2, y = y, `y2+` = y2))
}

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

# obtain fragmentation pattern, only works for sequence with one modifications -----
fragmentIon <- function(modseq_ori){
  Nterm_Acetylation <- FALSE
  if (grepl('\\(Acetyl \\(Protein N-term\\)\\)', modseq_ori)){
    modseq_ori <- stringr::str_replace(modseq_ori, '\\(Acetyl \\(Protein N-term\\)\\)', '') 
    Nterm_Acetylation <- TRUE
  }
  modseq <- trimws(modseq_ori, whitespace = '_')
  modrepl <- stringr::str_replace_all(modseq, '\\([^\\)]+\\)', '')
  modpos <- stringr::str_locate_all(modseq, '\\([^\\)]+\\)') |> lapply(function(n) n[1]) |> unlist()-1
  modtype <- stringr::str_extract_all(modseq, pattern = "(?<=\\()(\\w+)(?=\\))")|> lapply(function(n) n[1]) |> unlist()
  modstring <- paste0(rep(0, nchar(modrepl)), collapse = "")
  substring(modstring, first = modpos, last = modpos+1) <- as.character(which(modtbl$modname == modtype)-1)
  if (Nterm_Acetylation) {
    substring(modstring, first = 1, last = 2) <- as.character(which(modtbl$modname == 'Acetyl (Protein N-term)')-1)
    protViz::fragmentIon(modrepl, FUN = HCD_Ion, modified = modstring, modification = modtbl$mono)[[1]] # N_term argument is not working
  } else {
    protViz::fragmentIon(modrepl, FUN = HCD_Ion, modified = modstring, modification = modtbl$mono)[[1]]
  }
}


# (not used) isotope pattern extracted from index (find precusor and charge from index, not the actual peak)
extractedIsotopes <- function(index, ms2scanNum, isoNumber = 4){
  mmz <- scanInfo(index, ms2scanNum, 'monoisotopicMz')
  cha <- scanInfo(index, ms2scanNum, 'charge')
  isotopicIons <- mmz + c(0, seq_len(4))/cha
  return(isotopicIons)
}
# measured <- extractedIsotopes(high2Index, ms2_scanNum) # 535.9305 536.2638 536.5972 536.9305 537.2638

# precursor charge number
precursorChargeNum <- function(modseq_ori, cha, named = TRUE){
  seq <- modoutput(modseq_ori) # sequence wo modification
  modseq <- trimws(modseq_ori, whitespace = '_')
  modtype <- stringr::str_extract(modseq, pattern = "(?<=\\()(\\w+)(?=\\))")
  fw <- protViz:: parentIonMass(seq) + as.numeric(modtbl$mono[modtype == modtbl$modname])
  H <- 1.007825
  precursor <- vector(mode = 'double', length = cha)
  for (c in seq_len(cha)){
    precursor[c] <- (fw + (c-1) * H)/c
  }
  if(named){
    precursor <- round(precursor, digits = 4)
    precLab <- c('[M+H]+', '[M+2H]+2', '[M+3H]+3', '[M+4H]+4', '[M+5H]+5')
    paste0(precLab[seq_len(cha)], '=', precursor)
  } else {
    return(precursor)
  }
}

# theoretic isotope pattern ----
# (sequence (w mod), charge from Index$charge, number of isotope calculated, fragIon measure Max peak closest to Index'$monoisotopicMz)
theoreticalIsotopes <- function(modseq_ori, cha, isoNumber = 4, fragIon){
  seq <- modoutput(modseq_ori) # sequence wo modification
  modseq <- trimws(modseq_ori, whitespace = '_')
  modtype <- stringr::str_extract(modseq, pattern = "(?<=\\()(\\w+)(?=\\))")
  fw <- protViz:: parentIonMass(seq) + as.numeric(modtbl$mono[modtype == modtbl$modname])
  H <- 1.007825
  isotopicIons <- (fw + (cha - 1)*H)/cha + c(0, seq_len(isoNumber))/cha 
  isoLab <- c('M', 'M+1', 'M+2', 'M+3', 'M+4', 'M+5', 'M+6', 'M+7')
  mzerror <- isotopicIons - fragIon
  isoLab[mzerror == min(abs(mzerror))] <- paste0(isoLab[mzerror == min(abs(mzerror))], '(F)')
  names(isotopicIons) <- isoLab[c(0, seq_len(isoNumber)) + 1]
  return(isotopicIons) # named isotopic pattern, c('M', 'M+1', 'M+2', 'M+3', 'M+4')
}

# measured/actual isotope pattern 1) with strongest local signal intensity within the theoretical ppm window
measuredIsotopes <- function(locx, locy, theoretical, otype = c('x', 'y'),  ppm = 20){
  lapply(theoretical, function(n){
    xinRegion <- locx > n - n * ppm * 10^(-6) & locx < n + n * ppm * 10^(-6)
    switch(otype, x = locx[locy == max(locy[xinRegion])], y = locy[locy == max(locy[xinRegion])])
  }) |> unlist() |>  round(digits = 4)
}


# replace rawrr:::plot.rawrrSpectrum to show precursor and its isotope -----
plot.rawrrSpectrum <- function (x, modseq_ori, cha, ppm = 20, relative = TRUE, centroid = FALSE, SN = FALSE, legend = TRUE, 
                                diagnostic = FALSE, mmz, ...){
  if (is.null(mmz)){
    xlim <- x$massRange
    ylim <- c(0, 1.2 * max(x$centroid.intensity))
  } else {
    xlim <- c(floor(mmz), floor(mmz)+2.5) # expand x around monoisotopic peak
    # length(x$centroid.mZ); length(x$centroid.intensity)
    reducedCintensity <- x$centroid.intensity[which(x$centroid.mZ > xlim[1] & x$centroid.mZ < xlim[2])]
    ylim <- c(0, 1.2 * max(reducedCintensity)) # y is based on current window.
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
           frame.plot = FALSE, ...)
    } else {
      plot(x = x$mZ[mZidx], y = x$intensity[mZidx], type = "h", xlim = xlim, 
           xlab = "m/z", ylim = ylim_notCentroid, ylab = "Intensity", frame.plot = FALSE, 
           ...)
    }
    }
  
  # annotating precursor plot 
  theoretical <- theoreticalIsotopes(modseq_ori = modseq_ori, cha = cha, fragIon = mmz)
  measured_locx <- measuredIsotopes(x$mZ[mZidx], x$intensity[mZidx]/max(x$intensity[mZidx]), theoretical = theoretical, otype = 'x')
  measured_locy <- measuredIsotopes(x$mZ[mZidx], x$intensity[mZidx]/max(x$intensity[mZidx]), theoretical = theoretical, otype = 'y')
  axis(3, theoretical, names(theoretical), las = 2) # upper box: theoretical isotopic pattern
  axis(1, measured_locx, measured_locx, las = 2, cex.axis = 0.8) # lower box: max m/z within 20ppm of theoretical m/z
  points(measured_locx, measured_locy, col = "blue", pch = 22)
  legend("right", precursorChargeNum(modseq_ori, cha = cha), title = "Precursor Ions", bty = "n", cex = 0.65) # precursor charge number
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


# modified protViz::psm to have 1) correctly formated fragmentation ion pattern 
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
    legend("right", sprintf("% 10.3f   %s", sortedFragmentIonsTable$mass, 
                            sortedFragmentIonsTable$label), title = "Fragment Ions", 
           bty = "n", cex = 0.65)
  }
  # return(m)
}









theoretical   # 535.9319 536.2653 536.5986 536.9319 537.2653 
measured_locx # 535.9290 536.2656 536.5991 536.9330 537.2671 

high2Chr <- readChromatogram(file.path(sgms, 'high2.raw'), mass = theoretical, tol = 20)

function (rawfile, mass = NULL, tol = 10, filter = "ms", type = "xic") 
{
  .isAssemblyWorking()
  .checkRawFile(rawfile)
  stopifnot(type %in% c("xic", "bpc", "tic"))
  if (type == "xic") {
    if (is.null(mass)) {
      stop("No mass vector is provided.")
    }
    e <- .rawrrSystem2Source(rawfile, input = mass, rawrrArgs = sprintf("xic %f %s", 
                                                                        tol, shQuote(filter)))
    rv <- lapply(e$chromatogram, function(x) {
      attr(x, "filename") <- rawfile
      attr(x, "type") <- "xic"
      class(x) <- "rawrrChromatogram"
      x
    })
    if (is.null(rv[[1]]$times)) {
      errmsg <- c("The extraction of xic(s) failed for an unknown reason.", 
                  "\nPlease check the System Requirements.")
      stop(errmsg)
    }
  }
  else {
    rv <- .readChromatogramTicBpc(rawfile, filter, type)
  }
  attr(rv, "filter") <- filter
  attr(rv, "filename") <- rawfile
  if (type == "xic") {
    attr(rv, "type") <- "xic"
    attr(rv, "tol") <- tol
    class(rv) <- "rawrrChromatogramSet"
  }
  else if (type == "tic") {
    attr(rv, "type") <- "tic"
    class(rv) <- "rawrrChromatogram"
  }
  else {
    attr(rv, "type") <- "bpc"
    class(rv) <- "rawrrChromatogram"
  }
  rv
}



