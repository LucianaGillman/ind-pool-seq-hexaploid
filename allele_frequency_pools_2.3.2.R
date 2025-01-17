# Script to calculate allele frequencies for different pool combinations

# Load required libraries
library(tidyverse)  # For data manipulation and visualization
library(openxlsx)   # For reading and writing Excel files

# Set working directory
# Update this path as needed to point to the location of your data
setwd("/Users/mundomac/Desktop/2024/Doctorado/primer_paper/respuesta/scripts/correr/2.3")

# Load the dataset
countfile = read.csv("/Users/mundomac/Desktop/2024/Doctorado/primer_paper/respuesta/Material_suplementario/S3_Table.csv")

# Assign column names from the third row of the dataset
colnames(countfile) <- c(countfile[3,])
# Select only columns ending with "_1" or "_2", which correspond to two depths
countfileS <- countfile%>% select(ends_with('_1'),ends_with('_2'))

#add names to the pools
nombre_pool <- c ("24.20.1_1M", "24.20.2_1M","24.30.1_1M", "24.30.2_1M","24.40.1_1M",
                  "24.40.2_1M","24.50.1_1M", "24.50.2_1M","24.60.1_1M", "24.60.2_1M",
                  "28.20.1_1M", "28.20.2_1M","28.30.1_1M", "28.30.2_1M","28.40.1_1M",
                  "28.40.2_1M","28.50.1_1M", "28.50.2_1M","28.60.1_1M", "28.60.2_1M",
                  "50.20.1_1M", "50.20.2_1M","50.30.1_1M", "50.30.2_1M","50.40.1_1M",
                  "50.40.2_1M","50.50.1_1M", "50.50.2_1M","50.60.1_1M", "50.60.2_1M",
                  "87.20.1_1M", "87.20.2_1M","87.30.1_1M", "87.30.2_1M","87.40.1_1M",
                  "87.40.2_1M","87.50.1_1M", "87.50.2_1M","87.60.1_1M", "87.60.2_1M",
                  "88.20.1_1M", "88.20.2_1M","88.30.1_1M", "88.30.2_1M","88.40.1_1M",
                  "88.40.2_1M","88.50.1_1M", "88.50.2_1M","88.60.1_1M", "88.60.2_1M",
                  "24.20.1_2M", "24.20.2_2M","24.30.1_2M", "24.30.2_2M","24.40.1_2M",
                  "24.40.2_2M","24.50.1_2M", "24.50.2_2M","24.60.1_2M", "24.60.2_2M",
                  "28.20.1_2M", "28.20.2_2M","28.30.1_2M", "28.30.2_2M","28.40.1_2M",
                  "28.40.2_2M","28.50.1_2M", "28.50.2_2M","28.60.1_2M", "28.60.2_2M",
                  "50.20.1_2M", "50.20.2_2M","50.30.1_2M", "50.30.2_2M","50.40.1_2M",
                  "50.40.2_2M","50.50.1_2M", "50.50.2_2M","50.60.1_2M", "50.60.2_2M",
                  "87.20.1_2M", "87.20.2_2M","87.30.1_2M", "87.30.2_2M","87.40.1_2M",
                  "87.40.2_2M","87.50.1_2M", "87.50.2_2M","87.60.1_2M", "87.60.2_2M",
                  "88.20.1_2M", "88.20.2_2M","88.30.1_2M", "88.30.2_2M","88.40.1_2M",
                  "88.40.2_2M","88.50.1_2M", "88.50.2_2M","88.60.1_2M", "88.60.2_2M")

colnames(countfileS) <- nombre_pool

# Remove rows and convert all values to numeric for analysis
countfilef <- countfileS [-c(1:5),]
countfilef <- mutate_all(countfilef, function(x) as.numeric(as.character(x)))
rownames(countfilef) <- countfile[-c(1:5),1]

# calculo de frecuencias alelicas, promedio los duplicados y elimino los SNP con frecuencias menores de 0.05.
#armo csv con las tres profundidades

# Main loop for calculating allele frequencies across different accessions, pool sizes, MAF, and missing data 
for (ac in c('24','28','50','87','88')){
  for (i in c(0.1,0.2,0.3,0.4,0.5,0.6)) {
    pooles_ac <- cbind(AlleleID= rownames(countfilef),countfilef %>%select(starts_with(ac)))
    # Calculate missing data and filter rows with missing data proportion above threshold
    df_new <- apply(pooles_ac[, -1],2, function(x) head(x,-1)/(head(x,-1)+tail(x,-1)))
    toDelete <- seq(1, nrow(df_new), 2)
    df_=df_new[ toDelete ,]
    df<- as.data.frame(df_[rep(1:nrow(as.data.frame(df_)),1,each=2),])
    rownames(df) <- pooles_ac[,1]
    cuentamiss = function(datos) sum(is.na(datos))
    datos_perdidos = apply(df, 1, cuentamiss )/length(df[1,])
    pooles_dp = as.data.frame(pooles_ac[ - which(datos_perdidos>i), ])
    
    # Combine data from depths 1M and 2M
    pooles_1y2 <- pooles_dp[,-1]
    colnames(pooles_1y2)<-gsub(".1_","_",colnames(pooles_1y2))
    colnames(pooles_1y2)<-gsub(".2_","_",colnames(pooles_1y2))
    dat1y2 <- as.data.frame(t(rowsum(t(pooles_1y2), group = colnames(pooles_1y2), na.rm = T)))
    data1y2 <-dat1y2
    colnames(dat1y2)<-gsub("_1M","",colnames(dat1y2))
    colnames(dat1y2)<-gsub("_2M","",colnames(dat1y2))
    dat3 <- as.data.frame(t(rowsum(t(dat1y2), group = colnames(dat1y2), na.rm = T)))
    colnames(dat3) <- paste(colnames(dat3),"3M",sep="_")
    dtmp=as.data.frame(cbind(AlleleID=rownames(data1y2),data1y2, dat3)) 
    
    snp_data_sep<- separate(dtmp,col = 1 ,into = c("Clon","Allele" ,"ID"), sep="-", remove=F)
    snp_data_sep$Allele <- sub("", "Allele_1", snp_data_sep$Allele)
    snp_data_sep$Allele[snp_data_sep$Allele != 'Allele_1'] <- 'Allele_2'
    
    # Combine Allele1 and Allele2 into a single column
    snp_data_long <- snp_data_sep %>%
     pivot_longer(cols=5:dim(snp_data_sep)[2],
                   names_to = "Sample",
                   values_to = "Conteo")
    
    snp_data_wider<- snp_data_long[,-1] %>%
    pivot_wider(names_from = "Allele",
                values_from = "Conteo")
   
    # Calculate allele frequencies
    snp_data_wider<- snp_data_wider%>%
      mutate(Allele1_frec= Allele_1/(Allele_1+Allele_2))%>%
      mutate(Allele2_frec= Allele_2/(Allele_1+Allele_2))

    snp_data_long_wid <- snp_data_wider%>%
      pivot_longer(cols= 6:7,
                   names_to = "Allele",
                   values_to = "Frequence")
   
   # Filter SNPs by minor allele frequency
    for (mf in c(0.05, 0.01)){
      mfs <-  1-mf
      Datafreqs_longd <- snp_data_long_wid %>% filter(Frequence >= mf)%>% filter(Frequence <= mfs)
      Datafreqs_longd$Allele<-gsub("Allele","", Datafreqs_longd$Allele)
      Datafreqs_longd$Allele<-gsub("_frec","", Datafreqs_longd$Allele)
      
      #identificar cada columna por separado
      freqs <- separate(
        Datafreqs_longd,
        col= Sample,
        into = c('pool','depth'),
        remove = T,
        sep = "_",
        convert = FALSE)
      freqs <- separate(
        freqs,
        col= pool,
        into = c('Accesion','Tamano'),
        remove = T,
        convert = FALSE)
      t <- sub("0.","",i)
      # Save filtered data to CSV files for further analysis
      n <- paste('/Users/mundomac/Desktop/2024/Doctorado/primer_paper/respuesta/scripts/correr/2.3/pool_matrix/Datafreqs_1_2_3M_SV', ac,'ac_', t,'0dp_',mf,"maf",".csv", sep = "") 
      write_csv(freqs, n)
    }
  }
}

# End of script

