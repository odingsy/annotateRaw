library(rawrr)
library(protViz)
library(magrittr)



# directory -----
basePath <- ''
local <- file.path(basePath, '') # local where output png exists
sgms <- '' # .raw file dir
maxquantTXTs <- file.path(local, '') # containing maxquant txt output
source(file.path(basePath,''))

# read in index ----
# allindex <- allindex(dir = sgms, pattern = '.raw$', named = TRUE) # allindex(dir = sgms, pattern = '*.raw$', named = TRUE)
# saveRDS(allindex, file.path(local, 'rds/allindex.rds'))
allindex <- readRDS(file.path(local, 'rds/allindex.rds'))
indexnames <- stringr::str_replace(list.files(sgms, pattern = '.raw$'), '.raw$', '')

# Run this for manual input of peptide sequence -----
tbl <- tblContructor(peptideSeq = c('_TALQEVY(bp(Ywhc))TLAEHR_', '_AAEAAASAY(bp(Ywhc))YNPGNPHNVYMPTSQPPPPPYYPPEDK_'), 
                     indexnames = indexnames, Modification = 'mgo_MGH',ppm  = 20, rt = 2, cha = 3)

# Run this to extract information from maxquant output ----
tbl <- peptideCentric(maxquantTXTs = maxquantTXTs, modtbl = 'xxx.txt')
tbl <- tbl %>% 
  tidyr::expand_grid(., tibble::tibble(indexnames = indexnames)) %>% 
  dplyr::mutate(indextbl = purrr::map(indexnames, ~{allindex[[.]]})) %>% 
  scanNumSearching(., ms2snlist = 'table', indextbl = indextbl, raw = FALSE, mz = `MS/MS m/z`, ppm = ppm,rt = `Retention time`, rttol = 7) %>%
  dplyr::select(-indextbl)


for (i in seq_len(nrow(tbl))){
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
      pngPath <- file.path(local, "outpng", paste0(fn, ".png"))
      if (file.exists(pngPath)) {next}
      graphics.off()
      png(file = pngPath, width = 12, height = 7, units = 'in', res = 300)
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
      plot.rawrrChromatogramSet(high2Chr, rtmin = rtms1, offset = 10, isotopeName = names(theoretical), integalab = integalab, wAUC = wAUC)
      dev.off()
    }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}
