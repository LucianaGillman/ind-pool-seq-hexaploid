setwd("/Users/mundomac/Desktop/2024/Doctorado/primer_paper/script/dartr")
library(adegenet)
library (dartR)
library('poppr')
library(PopGenReport)
library(networkD3)
library('ggplot2')
##1 tengo que meter la matriz al bior entera 

data_br <- dartR::gl.read.dart(filename = "S4 Table.csv",
                       lastmetric = 15, ind.metafile = 'metadata.csv')

pools <- read.csv('pools_ind.csv')
pools$Pooles60 <- rep(60)

k=0
staditics_list <- list()
fst_list <- list()


for (indi in c(20,30,40,50,60)){
  p <- paste('Pooles',indi,sep = "")
  poolesf <- pools %>% filter(get(p) == indi)
  data_pool <- dartR::gl.keep.ind(data_br, poolesf$muestras, recalc = T, mono.rm = T, verbose = NULL)
  gl_data <- dartR::gl.filter.callrate(data_pool, method='loc', threshold = 0.9, verbose = 3)
  k=k+1
  gl_data_maf<- dartR::gl.filter.maf(
    gl_data,
    threshold = 0.05,
    by.pop = FALSE,
    recalc = FALSE,
    save2tmp = FALSE,
    verbose = NULL)
  
  #calculo de Ge y Fis  
  estadisticas <-gl.report.heterozygosity(gl_data_maf)
  
  #Calucluo de Ae
  geni_60 <- dartR::gl2gi(gl_data_maf)
  al_rich <- allel.rich(geni_60)
  estadisticas$Ae <-al_rich$mean.richness 
  
  ## reportar alelos privados una población contra el resto
  al_priv_porpop <- gl.report.pa(gl_data_maf, method = 'one2rest') #calculates the number of private alleles and fixed alleles between pairs of populations
  estadisticas$priv_one2rest <- al_priv_porpop$priv1
  
  estadisticas_name<- paste("estadisticas_",indi,"ind.csv", sep="")
  write.csv(estadisticas, estadisticas_name)
  staditics_list[[k]] <- estadisticas
  
  #Fst
  fst <- gl.fst.pop(gl_data_maf)
  Fst_name<- paste("Fst_",indi,"ind.csv", sep="")
  write.csv (fst$Fsts, Fst_name)
  fst_list <- fst$Fsts
  
}

saveRDS(staditics_list,"staditics_list.rds")
saveRDS(fst_list,"fst_list.rds")



