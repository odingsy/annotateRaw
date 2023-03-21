# DONE Non of N-acetyl is correct. correcting this 
# DONE two or more modification implement delete/extract using non-greedy regex match
# DONE exploring scan number independent mode, guessing scan number from maxquant 'MS/MS m/z'. 
# TODO come up a function enable semi-automating `integalab` selection. 
# TODO testing the current auc vs auc from xcalibur
# TODO integration multiple sample (w multiple CV) 
# TODO identify and quantify constituting unmodifed peptides. 
# TODO for precursor matches, instead of of find maximum around precursor (mmz), should do maxquant output (theoretical precursor mass from sequence and charge )vs raw index (actual precursor mass vs charge)
# TODO testing mass tolerance in by protVis_findNN_ is described. 


# log
# 230316: outpng_raw_5ppmAllSamples, outpng_split_5ppmAllSamples
# scanNumSearching: ppm = 5, rttol = 3.5
# plot: ppm = 5
# peakplot: itol = 5
# readChromatogram: tol = 10

# 230320
# run until: "077_ctrllong1_unmod2_c2_55708"


library(rawrr)
library(protViz)
library(magrittr)
library(parallel)


# directory -----
local <- '' # local path, code, output 
sgms <- '' # .raw file dir
source(file.path(local,'readRaw_functions.R'))
# source('https://raw.githubusercontent.com/odingsy/annotateRaw/main/readRaw_functions.R')


# read in index ----
# allindex <- allindex(dir = sgms, pattern = '*.raw$', named = TRUE)
allindex <- readRDS(file.path(local, 'rds/allindex.rds'))

tbl <- readRDS(file.path(local, 'rds/tbl_modpng.rds'))
combiPlot <- function(i){
  modseq_ori <- tbl[[i, 'allseq']]
  cha <- tbl[[i, 'cha']]
  indexName <- tbl[[i, 'indexnames']]
  index <- allindex[[indexName]]
  scanNum <- tbl[[i, 'scanNumList']] %>% unlist()
  peptideNum <- tbl[[i, 'rn']]
  modtype <- tbl[[i, 'allseq_name']]
  isolatemz <- tbl[[i, 'isolatemz']]
  if (length(scanNum) == 0) {next}
    for (j in seq_along(scanNum)){
    tryCatch({
      # all parameters 
      ms2_scanNum <- scanNum[j] # tbl[[i, 'MS/MS scan number']]
      rawFile <- file.path(sgms, paste0(indexName, '.raw'))
      # integalab = c(83.37490, 84.47190)
      wAUC = FALSE # whether chromatogram will calculate AUC
      fn <- paste(sprintf("%03.0f", peptideNum), stringr::str_remove(indexName, '_'), modtype, paste0('c', cha), ms2_scanNum, sep = '_')
      print(fn)
      
      # graphing parameters
      graphics.off()
      png(file = file.path(local, "outpng_raw", paste0(fn, ".png")), width = 12, height = 7, units = 'in', res = 300)
      # layout setting
      layout.matrix <- matrix(c(3, 3, 1, 2), byrow = TRUE, nrow = 2)
      layout(mat = layout.matrix,
             heights = c(2, 3), # Heights of the two rows
             widths = c(2, 2)) # Widths of the two columns
      # layout.show(3)
      
      # precursor isotopic profile. 
      x <- readSpectrum(rawFile, scan = scanInfo(index, ms2_scanNum, 'masterScan'))[[1]]
      mmz <- scanInfo(index, ms2_scanNum, 'precursorMass')
      plot(x, modseq_ori, cha = cha, relative = TRUE, SN = FALSE, mmz = mmz, ppm = 5) # xlim take monoisotopic and derived xaxis window around this value.
      
      # PSM using peakplot 
      prodIon <- readSpectrum(rawFile, scan = ms2_scanNum)[[1]]
      peakPlot(oriSeq = modseq_ori, spec = prodIon, itol = 5, unit = 'ppm')
      
      # chromatograms of the isotopic pattern 
      theoretical <- precursorMz(modseq_ori, mode = 'isoPattern', fragIon = mmz, cha = cha)
      high2Chr <- readChromatogram(rawFile, mass = theoretical, tol = 10)
      rtms1 <- scanInfo(index, ms2_scanNum, 'rtinseconds')/60
      plot.rawrrChromatogramSet(high2Chr, rtmin = rtms1, offset = 10, integalab = integalab, wAUC = wAUC)
      dev.off()
    }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
}}

 
mclapply(1:nrow(tbl), combiPlot, mc.cores = detectCores() - 1)



