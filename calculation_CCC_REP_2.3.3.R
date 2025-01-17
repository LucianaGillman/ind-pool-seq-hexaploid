# This script calculates the Representativity (%) and concordance correlation coefficient (CCC) between allele frequencies
# of pooled and individuals samples.

# Set the working directory to the location of input and output files
# Load required libraries
library(tidyverse)  # For data manipulation and visualization
library(epiR)       # For calculating CCC

setwd("/Users/mundomac/Desktop/2024/Doctorado/primer_paper/respuesta/scripts/correr/2.3")

# Initialize empty lists to store CCC results and SNP counts
ccc_list_est <- list()  # Estimated CCC values
ccc_list_low <- list()  # Lower confidence interval of CCC
ccc_list_upper <- list()  # Upper confidence interval of CCC
SNP_comun <- numeric()  # Number of Common SNPs between pooled and individuals data sets
SNP_ind <- numeric()  # Number of SNPs from individuals
nameslist <- list()  # List to store descriptive names for outputs

# Counter for tracking iterations
k = 0

# Iterate through accessions, sizes, depths, and minor allele frequencies (MAFs)

for (ac in c(24, 28, 50, 87, 88)) {  # Loop through accessions
  for (tm in c(20, 30, 40, 50, 60)) {  # Loop through size categories
    for (dp in c(10, 20, 30, 40, 50, 60)) {  # Loop through depth thresholds
      for (maf in c(0.05, 0.01)) {  # Loop through MAF thresholds
        # Load individual allele frequency data
        dpii <- dp/100  # Convert depth to decimal format
        ni <- paste('/Users/mundomac/Desktop/2024/Doctorado/primer_paper/respuesta/scripts/correr/2.3/ind_matrix/freqs_SVi_',tm,'ind_',ac,"ac_",dpii,'dp.csv', sep = "")
        frec_ind <- read.csv(ni)
        
        # Adjust Allele column to match the pool data format
        frec_ind$Allele <- c(2,frec_ind$Allele[-dim(frec_ind)[1]])
        
        # Load pool allele frequency data
        np <- paste('/Users/mundomac/Desktop/2024/Doctorado/primer_paper/respuesta/scripts/correr/2.3/pool_matrix/Datafreqs_1_2_3M_SV',ac,'ac_',dp,'dp_',maf,'maf','.csv', sep = "")
        frec_pool <- read.csv(np)
        
        # Filter pool data by size and depth
        for (de in c ("1M","2M","3M")){ # Loop through depths (1M, 2M, 3M)
          frec_pool_filt2 <- frec_pool %>% filter(Tamano==tm,depth ==de)
          
          # Increment counter
          k=k+1
          
          # Merge individual and pool data on common columns
          join_ip <- inner_join(frec_ind, frec_pool_filt2, by = c('Clon', 'Allele', "ID"))
          
          # Calculate CCC using epiR
          ccc <- epi.ccc(join_ip$freq_ind,join_ip$Frequence, 
                         ci = "z-transform", conf.level = 0.95, rep.measure = FALSE)
          
          # Store CCC estimates and SNP counts
          ccc_list_est[[k]] <- ccc$rho.c[1]  # CCC estimate
          ccc_list_low[[k]] <- ccc$rho.c[2]  # Lower CI
          ccc_list_upper[[k]] <- ccc$rho.c[3]  # Upper CI
          SNP_comun[[k]] <- dim(join_ip)[1] / 2  # Common SNPs (divide by 2 as counts are paired)
          SNP_ind[[k]] <- dim(frec_ind)[1] / 2  # SNPs from individuals
          
          # Generate descriptive name for the current iteration
          nombre <- paste(tm, '_',de,'_',ac, "_", dp,'_',maf,sep = "") 
          nameslist[[k]] <- nombre
            
        }
      }
    }
  }
}

# Combine results into a single data frame for export
CCC_Rep <-as.data.frame(cbind(
  name = unlist(nameslist), 
  Number_SNPs_common= unlist(SNP_comun),
  Representativity = SNP_comun*100/SNP_ind,
  CCC=unlist(ccc_list_est),
  CCClow=unlist(ccc_list_low),
  CCCupper=unlist(ccc_list_upper)
))

# Split descriptive names into separate columns for easier interpretation
CCC_sep <- separate(
  CCC_Rep,
  col= name,
  into = c('Sample Size','Sequence Depth','Accession',"Missing Data", "MAF pool"),
  remove = T,
  sep = "_",
  convert = FALSE)

write.csv(CCC_sep, 'CCC_Rep.csv')
