library(rawrr)
library(protViz)
library(magrittr)



# directory -----
rawdir <- '' # .raw file dir
maxqdir <- '' # maxquant dir
outpngdir <- '' # output png files
source(file.path(outpngdir,'readRaw_functions.R'))

# read in index ----
# allindex <- extractIndex(dir = rawdir, pattern = '.raw$', named = TRUE) # allindex(dir = rawdir, pattern = '*.raw$', named = TRUE)
# saveRDS(allindex, file.path(outpngdir, 'rds/savedIndices.rds'))
allindex <- readRDS(file.path(outpngdir, 'rds/savedIndices.rds'))
indexnames <- stringr::str_replace(list.files(rawdir, pattern = '.raw$'), '.raw$', '')

# Method1: Manually input of a peptide sequence -----
tbl <- modsManual(peptideSeq = c('_VASLR(mgo_MGH)ETYGDMADCCEK_', '_FLGDR(mgo_MGH)DFNQLSSR_'), 
                     indexnames = indexnames, Modification = 'mgo_MGH',ppm  = 20, rt = 2, cha = 3)


# Method2: extract information from maxquant output ----
tbl <- modsFromMaxquant(maxquantTXTs = maxqdir, modtbl = 'mgo_MGHSites.txt')


# Main: extracting info from the raw file ----
future::plan(future::multisession, workers = 50)
t <- tbl %>% 
  head(50) %>% 
  tidyr::unnest(scanNumList) %>% 
  dplyr::mutate(index = furrr::future_map(indexnames, ~allindex[[.]]),
                rawFile = file.path(rawdir, paste0(indexnames, '.raw')),
                wAUC = FALSE) %>% 
  dplyr::rename(peptideNum = rn, modtype = Modification, ms2_scanNum = scanNumList) %>% 
  dplyr::mutate(fn = paste(sprintf("%03.0f", peptideNum), stringr::str_remove(indexName, '_'), modtype, paste0('c', cha), ms2_scanNum, sep = '_'),
                pngPath = file.path(outpngdir, "outpng", paste0(fn, ".png")),
                x = furrr::future_pmap(list(rawFile, index, ms2_scanNum), ~{readSpectrum(..1, scan = scanInfo(..2, ..3, 'masterScan'))[[1]]}),
                mmz = furrr::future_pmap(list(index, ms2_scanNum), ~scanInfo(..1, ..2, 'precursorMass')),
                prodIon = furrr::future_pmap(list(rawFile, ms2_scanNum), ~{readSpectrum(..1, scan = ..2)[[1]]}),
                theoretical = furrr::future_pmap(list(allseq, mmz, cha), ~{precursorMz(..1, mode = 'isoPattern', fragIon = ..2, cha = ..3)}), 
                high2Chr = furrr::future_pmap(list(rawFile, theoretical), ~{readChromatogram(..1, mass = ..2, tol = 10)}),
                rtms1 = furrr::future_pmap(list(index, ms2_scanNum), ~{scanInfo(..1, ..2, 'rtinseconds')/60}))


purrr::pwalk(t, plotMS_arg)



