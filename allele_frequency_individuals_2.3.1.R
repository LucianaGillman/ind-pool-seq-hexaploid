#This present script is for calculating the allele frequencies across accessions 
#form individual samples with score data set, varying sample size and missing data
#With this script you going to save the frequencies of each accession, sample size 
#and missing data in a file. 

# Set the working directory where your scripts and output files are located
setwd('/Users/mundomac/Desktop/2024/Doctorado/primer_paper/respuesta/scripts/correr/2.3')

# Load necessary libraries for genetic data analysis and data manipulation
library(adegenet)  # For genetic data analysis
library(tidyverse) # For data manipulation and visualization

# Read the main SNP data file. Ensure the file path is correct, and handle missing values with `na.strings`
datos <- read.csv("/Users/mundomac/Desktop/2024/Doctorado/primer_paper/respuesta/Material_suplementario/S2_Table.csv", 
                  stringsAsFactors = FALSE, na.strings = "na")

# Adjust column names for better readability and use the correct rows to name the columns
colnames(datos) <- c(datos[5, 1:30], datos[3, 31:dim(datos)[2]])
datos_tra <- datos[-c(1:5), -c(3:30)] # Remove irrelevant rows and columns for cleaner data

# Separate the AlleleID into individual components: Clon, Allele, and ID
snp_data_sep <- separate(datos_tra, col = 1, into = c("Clon", "Allele", "ID"), sep = "-", remove = FALSE)
snp_data_sep$Allele <- sub("", "Allele_1", snp_data_sep$Allele) # Label first allele as Allele_1
snp_data_sep$Allele[snp_data_sep$Allele != 'Allele_1'] <- 'Allele_2' # Label the rest as Allele_2

# Convert the data to long format for easier manipulation
snp_data_long <- snp_data_sep %>%
  pivot_longer(cols = 6:dim(snp_data_sep)[2], 
               names_to = "Sample", 
               values_to = "SNPs")

# Group data by Clon, ID, and Sample, and combine SNP values for each combination
snp_data_j <- snp_data_long %>%
  group_by(Clon, ID, Sample) %>%
  mutate(SNPuno = paste(SNPs, collapse = " "))

# Function to standardize SNP encoding: Convert "0 1" to "1", "1 0" to "0", "1 1" to "2", and missing data to NA
substitute_01 <- function(df, column) {
  df[[column]] <- sub("0 1", "1", df[[column]])
  df[[column]] <- sub("1 0", "0", df[[column]]) 
  df[[column]] <- sub("1 1", "2", df[[column]])
  df[[column]] <- sub("- -", NA, df[[column]])
  return(df)
}
df <- substitute_01(snp_data_j, "SNPuno")

# Filter for Allele_1 to simplify subsequent analyses
filter <- df %>%
  filter(Allele == "Allele_1")

# Convert the filtered data to a wide format, with each sample as a column
data_one_row <- pivot_wider(filter[,-c(2,3,4,5,7)], 
                            names_from = "Sample", 
                            values_from = SNPuno)

# Load pool data and add an additional column for Pooles60
pools <- read.csv('/Users/mundomac/Desktop/2024/Doctorado/primer_paper/respuesta/scripts/github/pools_ind.csv')
pools <- cbind(pools, Pooles60 = rep(60))

# Main loop for calculating allele frequencies across different accessions, pool sizes, and missing data 
for (indi in c(20, 30, 40, 50, 60)) { # Iterate over pool sizes
  p <- paste('Pooles', indi, sep = "")
  pooles_mezclaf <- pools %>% filter(get(p) == indi) # Filter for current pool size
  names_mezclaf <- names(data_one_row)[(names(data_one_row) %in% pooles_mezclaf$muestras)]
  
  for (ac in c(24, 28, 50, 87, 88)) { # Iterate over accessions
    pooles_ac <- pooles_mezclaf %>% filter(Accesion == ac) # Filter by accession
    names_ac <- names_mezclaf[(names_mezclaf %in% pooles_ac$muestras)]
    ind_ac <- cbind(AlleleID = data_one_row[,1], data_one_row[, names_ac])
    
    # Calculate missing data percentage
    cuentamiss <- function(datos) sum(is.na(datos))
    FreqMiss <- apply(ind_ac[,-1], 1, cuentamiss) / length(ind_ac[1,])
    
    for (i in seq(0.1, 0.6, by = 0.1)) { # Filter by missing data threshold
      datos <- ind_ac %>% filter(FreqMiss < i)
      
      # Calculate allele frequencies
      frecuence <- apply(datos[,-1], 1, function(snp) { 
        fAA <- length(which(snp == 0))
        faa <- length(which(snp == 1))
        fAa <- length(which(snp == 2))
        q <- (faa + fAa / 2) / (faa + fAa + fAA)
        return(q)
      })
      
      # Prepare frequency data for export
      freq <- as.data.frame(cbind(SNP = datos[,1], Freq1 = frecuence, Freq2 = 1 - frecuence))
      Datafreqs <- freq %>% filter(Freq1 >= 0.05) %>% filter(Freq1 <= 0.95) # Filter by MAF >= 0.05
      
      freqs <- separate(Datafreqs, col = SNP, sep = "-", into = c('Clon', "2", 'ID'), 
                        remove = FALSE, convert = FALSE)
      freqlong <- pivot_longer(freqs[,-3], cols = 4:5, names_to = "Allele", values_to = "freq_ind")
      freqlong$Allele <- gsub("Freq", "", freqlong$Allele)
      
      # Save frequency data to a CSV file
      alname <- paste("/Users/mundomac/Desktop/2024/Doctorado/primer_paper/respuesta/scripts/correr/2.3/ind_matrix/freqs_SVi_", 
                      indi, "ind_", ac, "ac_", i, "dp.csv", sep = "")
      write.csv(freqlong, alname)
    }
  }
}
