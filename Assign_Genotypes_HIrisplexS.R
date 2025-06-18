########Things you need to set#################
# Set working directory
setwd("N:\projects\wgs_pipeline_testing_6180\users\wmj412\testing\analysis\8A-D2-SS-root")

# Set directory where you want to save the text files
output_dir <- "N:\projects\wgs_pipeline_testing_6180\users\wmj412\testing"

#Set filtering criteria

DP=4
GQ=10

#############################################

# Read VCF files from the folder
files <- list.files(path = getwd(), pattern = "\\.vcf$", full.names = TRUE)

# Load necessary package
# install.packages("tidyverse")
library(tidyverse)


# Function to read a VCF file
read_vcf <- function(file) {
  # Read the file
  lines <- readLines(file)
  
  # Filter out header lines (starting with ##)
  data_lines <- lines[!grepl("^##", lines)]
  
  # Separate header (starting with #) and data
  header <- data_lines[grepl("^#", data_lines)]
  data <- data_lines[!grepl("^#", data_lines)]
  
  # Read the data into a data frame
  if (length(data) > 0) {
    data <- read.table(text = data, header = FALSE, stringsAsFactors = FALSE)
    colnames(data) <- strsplit(header, "\t")[[1]]
  } else {
    data <- data.frame()
  }
  
  return(data)
}

# Read all VCF files into a list
vcf_list <- lapply(files, read_vcf)

# Assign names of files to the list
names(vcf_list) <- basename(files)

# Function to rename column #CHROM to CHROM
rename_column <- function(df) {
  colnames(df)[colnames(df) == "#CHROM"] <- "CHROM"
  return(df)
}

# Apply the function to each element (data frame) in the list
vcf_list <- lapply(vcf_list, rename_column)
# Extract GT, DP, and GQ for each SNP
extract_field <- function(df, index) {
  sapply(strsplit(df[, 10], ":"), function(x) x[index])
}

for(i in seq_along(vcf_list)) {
  vcf_list[[i]]$GT <- extract_field(vcf_list[[i]], 1)
  vcf_list[[i]]$DP <- as.numeric(extract_field(vcf_list[[i]], 3))
  vcf_list[[i]]$GQ <- as.numeric(extract_field(vcf_list[[i]], 4))
}

# Define the function to create the new column based on GT
create_genotype_column <- function(df) {
  df$GT_Allel <- apply(df, 1, function(row) {
    gt <- row["GT"]
    ref <- row["REF"]
    alt <- strsplit(row["ALT"], ",")[[1]][1]
    
    if (gt == "0/0") {
      return(paste0(ref, ref))
    } else if (gt == "1/1") {
      return(paste0(alt, alt))
    } else if (gt == "0/1") {
      return(paste0(ref, alt))
    } else {
      return(NA)
    }
  })
  
  # Change GT_Allel to NA if DP < DP-value selected or GQ < GQ-value selected
  df$GT_Allel[df$DP < DP | df$GQ < GQ] <- NA
  
  return(df)
}

# Process each element in the list to create the new column
vcf_list <- lapply(vcf_list, create_genotype_column)

# Assign the genotypes to the HIrisPlex-S rsIDs

H_POS <- data.frame(CHROM = paste0("chr", c(
  11, 9, 16, 11, 15, 16, 6, 15, 12, 14, 15, 11, 15, 15, 15, 15, 5,
  14, 15, 15, 16, 16, 16, 16, 16, 16, 16, 15, 20, 14, 5, 16, 16, 6,
  20, 20, 15, 9, 16, 16, 16
)),
POS = c(
  89178528, 16858086, 89919722, 89284793, 28111713, 89919683,
  396321, 28026629, 88934558, 92307319, 28120472, 89277878, 48134287,
  28042975, 27942626, 28285036, 33951588, 92416482, 27985172, 27951891,
  89919436, 89919510, 89919709, 89919736, 89920138, 89919714, 89919532,
  28208069, 34630286, 92334859, 33958854, 89317317, 89917970, 457748,
  34077942, 34197406, 28251049, 12709305, 89957798, 89919746, 89919344),
  
  HIrisPlexS_ID = c(
    "rs1042602", "rs10756819", "rs1110400", "rs1126809", "rs1129038",
    "rs11547464", "rs12203592", "rs12441727", "rs12821256", "rs12896399",
    "rs12913832", "rs1393350", "rs1426654", "rs1470608", "rs1545397",
    "rs1667394", "rs16891982", "rs17128291", "rs1800407", "rs1800414",
    "rs1805005", "rs1805006", "rs1805007", "rs1805008", "rs1805009",
    "rs201326893", "rs2228479", "rs2238289", "rs2378249", "rs2402130",
    "rs28777", "rs3114908", "rs3212355", "rs4959270", "rs6059655",
    "rs6119471", "rs6497292", "rs683", "rs8051733", "rs885479",
    "rs312262906")
)

# Create a list to store the updated data frames
HIrisPlexS_list <- list()

# Iterate through vcf_files and perform merge
for (i in seq_along(vcf_list)) {
  
  # Merge H_POS with vcf_list[[i]] using CHROM and POS as keys
  merged_df <- merge(H_POS, vcf_list[[i]][c("CHROM", "POS", "GT", "DP", "GQ", "GT_Allel")], by = c("CHROM", "POS"))
  
  # Assign name from vcf_files to merged_df
  name <- names(vcf_list)[i]
  names(merged_df) <- make.names(c("CHROM", "POS", "HIrisPlexS_ID", "GT", "DP", "GQ", "GT_Allel"), unique = TRUE)
  
  # Reorder columns so HIrisPlexS_ID and GT_Allel are last
  merged_df <- merged_df[, c(setdiff(names(merged_df), c("HIrisPlexS_ID", "GT_Allel")), "HIrisPlexS_ID", "GT_Allel")]
  
  # Store merged_df in merged_list
  HIrisPlexS_list[[name]] <- merged_df
}

#order the list according to HIrisPlex-S upload order
order_ids <- c(
  "rs312262906", "rs11547464", "rs885479", "rs1805008", "rs1805005",
  "rs1805006", "rs1805007", "rs1805009", "rs201326893", "rs2228479",
  "rs1110400", "rs28777", "rs16891982", "rs12821256", "rs4959270",
  "rs12203592", "rs1042602", "rs1800407", "rs2402130", "rs12913832",
  "rs2378249", "rs12896399", "rs1393350", "rs683", "rs3114908",
  "rs1800414", "rs10756819", "rs2238289", "rs17128291", "rs6497292",
  "rs1129038", "rs1667394", "rs1126809", "rs1470608", "rs1426654",
  "rs6119471", "rs1545397", "rs6059655", "rs12441727", "rs3212355",
  "rs8051733"
)

# Order merged_list according to order_ids

for (i in seq_along(HIrisPlexS_list)) {
  # Reorder rows based on order_ids
  HIrisPlexS_list[[i]] <- HIrisPlexS_list[[i]][match(order_ids, HIrisPlexS_list[[i]]$HIrisPlexS_ID), ]
}


# Iterate through HIrisPlexS_list and write each element to a text file
#for (name in names(HIrisPlexS_list)) {
  # File name based on the name of the element in HIrisPlexS_list
#  filename <- paste0(output_dir, name, ".txt")
  
  # Write the data frame to a text file with sep="\" and no row names
#  write.table(HIrisPlexS_list[[name]], file = filename, sep = "\t", row.names = FALSE, quote=FALSE)
  
#  cat("File", filename, "written.\n")
#}

# Initialize an empty list to store GT_Allel columns
GT_Allel_list <- list()

# Extract GT_Allel columns from each data frame in HIrisPlexS_list
for (i in seq_along(HIrisPlexS_list)) {
  GT_Allel_list[[i]] <- HIrisPlexS_list[[i]]$GT_Allel
}

# Transpose GT_Allel_list into a row
HIrisPlexS_upload <- as.data.frame(t(sapply(GT_Allel_list, as.character)))

# Insert names(HIrisPlexS_list) as the first column
HIrisPlexS_upload <- cbind(Name = names(HIrisPlexS_list), HIrisPlexS_upload)

colnames(HIrisPlexS_upload) <- c("sampleid", "rs312262906_A", "rs11547464_A", "rs885479_T", "rs1805008_T", "rs1805005_T",
                  "rs1805006_A", "rs1805007_T", "rs1805009_C", "rs201326893_A", "rs2228479_A", "rs1110400_C",
                  "rs28777_C", "rs16891982_C", "rs12821256_G", "rs4959270_A", "rs12203592_T", "rs1042602_T",
                  "rs1800407_A", "rs2402130_G", "rs12913832_T", "rs2378249_C", "rs12896399_T", "rs1393350_T",
                  "rs683_G", "rs3114908_T", "rs1800414_C", "rs10756819_G", "rs2238289_C", "rs17128291_C",
                  "rs6497292_C", "rs1129038_G", "rs1667394_C", "rs1126809_A", "rs1470608_A", "rs1426654_G",
                  "rs6119471_C", "rs1545397_T", "rs6059655_T", "rs12441727_A", "rs3212355_A", "rs8051733_C")

# Change alleles to minor allele count based on the HIrisPlex-S annotation (https://walshlab.sitehost.iu.edu/hpsconvertsonline.R)

HIrisPlexS_upload$rs312262906_A<-gsub('C/C', '0', HIrisPlexS_upload$rs312262906_A)
HIrisPlexS_upload$rs312262906_A<-gsub('CC', '0', HIrisPlexS_upload$rs312262906_A)
HIrisPlexS_upload$rs312262906_A<-gsub('A/A', '2', HIrisPlexS_upload$rs312262906_A)
HIrisPlexS_upload$rs312262906_A<-gsub('AA', '2', HIrisPlexS_upload$rs312262906_A)
HIrisPlexS_upload$rs312262906_A<-gsub('C/A', '1', HIrisPlexS_upload$rs312262906_A)
HIrisPlexS_upload$rs312262906_A<-gsub('CA', '1', HIrisPlexS_upload$rs312262906_A)
HIrisPlexS_upload$rs312262906_A<-gsub('A/C', '1', HIrisPlexS_upload$rs312262906_A)
HIrisPlexS_upload$rs312262906_A<-gsub('AC', '1', HIrisPlexS_upload$rs312262906_A)
HIrisPlexS_upload$rs312262906_A<-gsub('C', 'NA', HIrisPlexS_upload$rs312262906_A)
HIrisPlexS_upload$rs312262906_A<-gsub('A', 'NA', HIrisPlexS_upload$rs312262906_A)
HIrisPlexS_upload$rs312262906_A<-gsub('G', 'NA', HIrisPlexS_upload$rs312262906_A)
HIrisPlexS_upload$rs312262906_A<-gsub('T', 'NA', HIrisPlexS_upload$rs312262906_A)
HIrisPlexS_upload$rs312262906_A<-gsub('./.', 'NA', HIrisPlexS_upload$rs312262906_A)
HIrisPlexS_upload$rs11547464_A<-gsub('G/G', '0', HIrisPlexS_upload$rs11547464_A)
HIrisPlexS_upload$rs11547464_A<-gsub('GG', '0', HIrisPlexS_upload$rs11547464_A)
HIrisPlexS_upload$rs11547464_A<-gsub('A/A', '2', HIrisPlexS_upload$rs11547464_A)
HIrisPlexS_upload$rs11547464_A<-gsub('AA', '2', HIrisPlexS_upload$rs11547464_A)
HIrisPlexS_upload$rs11547464_A<-gsub('G/A', '1', HIrisPlexS_upload$rs11547464_A)
HIrisPlexS_upload$rs11547464_A<-gsub('GA', '1', HIrisPlexS_upload$rs11547464_A)
HIrisPlexS_upload$rs11547464_A<-gsub('A/G', '1', HIrisPlexS_upload$rs11547464_A)
HIrisPlexS_upload$rs11547464_A<-gsub('AG', '1', HIrisPlexS_upload$rs11547464_A)
HIrisPlexS_upload$rs11547464_A<-gsub('G', 'NA', HIrisPlexS_upload$rs11547464_A)
HIrisPlexS_upload$rs11547464_A<-gsub('A', 'NA', HIrisPlexS_upload$rs11547464_A)
HIrisPlexS_upload$rs11547464_A<-gsub('C', 'NA', HIrisPlexS_upload$rs11547464_A)
HIrisPlexS_upload$rs11547464_A<-gsub('T', 'NA', HIrisPlexS_upload$rs11547464_A)
HIrisPlexS_upload$rs11547464_A<-gsub('./.', 'NA', HIrisPlexS_upload$rs11547464_A)
HIrisPlexS_upload$rs885479_T<-gsub('G/G', '0', HIrisPlexS_upload$rs885479_T)
HIrisPlexS_upload$rs885479_T<-gsub('GG', '0', HIrisPlexS_upload$rs885479_T)
HIrisPlexS_upload$rs885479_T<-gsub('A/A', '2', HIrisPlexS_upload$rs885479_T)
HIrisPlexS_upload$rs885479_T<-gsub('AA', '2', HIrisPlexS_upload$rs885479_T)
HIrisPlexS_upload$rs885479_T<-gsub('G/A', '1', HIrisPlexS_upload$rs885479_T)
HIrisPlexS_upload$rs885479_T<-gsub('GA', '1', HIrisPlexS_upload$rs885479_T)
HIrisPlexS_upload$rs885479_T<-gsub('A/G', '1', HIrisPlexS_upload$rs885479_T)
HIrisPlexS_upload$rs885479_T<-gsub('AG', '1', HIrisPlexS_upload$rs885479_T)
HIrisPlexS_upload$rs885479_T<-gsub('G', 'NA', HIrisPlexS_upload$rs885479_T)
HIrisPlexS_upload$rs885479_T<-gsub('A', 'NA', HIrisPlexS_upload$rs885479_T)
HIrisPlexS_upload$rs885479_T<-gsub('C', 'NA', HIrisPlexS_upload$rs885479_T)
HIrisPlexS_upload$rs885479_T<-gsub('T', 'NA', HIrisPlexS_upload$rs885479_T)
HIrisPlexS_upload$rs885479_T<-gsub('./.', 'NA', HIrisPlexS_upload$rs885479_T)
HIrisPlexS_upload$rs1805008_T<-gsub('C/C', '0', HIrisPlexS_upload$rs1805008_T)
HIrisPlexS_upload$rs1805008_T<-gsub('CC', '0', HIrisPlexS_upload$rs1805008_T)
HIrisPlexS_upload$rs1805008_T<-gsub('T/T', '2', HIrisPlexS_upload$rs1805008_T)
HIrisPlexS_upload$rs1805008_T<-gsub('TT', '2', HIrisPlexS_upload$rs1805008_T)
HIrisPlexS_upload$rs1805008_T<-gsub('C/T', '1', HIrisPlexS_upload$rs1805008_T)
HIrisPlexS_upload$rs1805008_T<-gsub('CT', '1', HIrisPlexS_upload$rs1805008_T)
HIrisPlexS_upload$rs1805008_T<-gsub('T/C', '1', HIrisPlexS_upload$rs1805008_T)
HIrisPlexS_upload$rs1805008_T<-gsub('TC', '1', HIrisPlexS_upload$rs1805008_T)
HIrisPlexS_upload$rs1805008_T<-gsub('C', 'NA', HIrisPlexS_upload$rs1805008_T)
HIrisPlexS_upload$rs1805008_T<-gsub('T', 'NA', HIrisPlexS_upload$rs1805008_T)
HIrisPlexS_upload$rs1805008_T<-gsub('G', 'NA', HIrisPlexS_upload$rs1805008_T)
HIrisPlexS_upload$rs1805008_T<-gsub('A', 'NA', HIrisPlexS_upload$rs1805008_T)
HIrisPlexS_upload$rs1805008_T<-gsub('./.', 'NA', HIrisPlexS_upload$rs1805008_T)
HIrisPlexS_upload$rs1805005_T<-gsub('G/G', '0', HIrisPlexS_upload$rs1805005_T)
HIrisPlexS_upload$rs1805005_T<-gsub('GG', '0', HIrisPlexS_upload$rs1805005_T)
HIrisPlexS_upload$rs1805005_T<-gsub('T/T', '2', HIrisPlexS_upload$rs1805005_T)
HIrisPlexS_upload$rs1805005_T<-gsub('TT', '2', HIrisPlexS_upload$rs1805005_T)
HIrisPlexS_upload$rs1805005_T<-gsub('G/T', '1', HIrisPlexS_upload$rs1805005_T)
HIrisPlexS_upload$rs1805005_T<-gsub('GT', '1', HIrisPlexS_upload$rs1805005_T)
HIrisPlexS_upload$rs1805005_T<-gsub('T/G', '1', HIrisPlexS_upload$rs1805005_T)
HIrisPlexS_upload$rs1805005_T<-gsub('TG', '1', HIrisPlexS_upload$rs1805005_T)
HIrisPlexS_upload$rs1805005_T<-gsub('G', 'NA', HIrisPlexS_upload$rs1805005_T)
HIrisPlexS_upload$rs1805005_T<-gsub('T', 'NA', HIrisPlexS_upload$rs1805005_T)
HIrisPlexS_upload$rs1805005_T<-gsub('C', 'NA', HIrisPlexS_upload$rs1805005_T)
HIrisPlexS_upload$rs1805005_T<-gsub('A', 'NA', HIrisPlexS_upload$rs1805005_T)
HIrisPlexS_upload$rs1805005_T<-gsub('./.', 'NA', HIrisPlexS_upload$rs1805005_T)
HIrisPlexS_upload$rs1805006_A<-gsub('C/C', '0', HIrisPlexS_upload$rs1805006_A)
HIrisPlexS_upload$rs1805006_A<-gsub('CC', '0', HIrisPlexS_upload$rs1805006_A)
HIrisPlexS_upload$rs1805006_A<-gsub('A/A', '2', HIrisPlexS_upload$rs1805006_A)
HIrisPlexS_upload$rs1805006_A<-gsub('AA', '2', HIrisPlexS_upload$rs1805006_A)
HIrisPlexS_upload$rs1805006_A<-gsub('C/A', '1', HIrisPlexS_upload$rs1805006_A)
HIrisPlexS_upload$rs1805006_A<-gsub('CA', '1', HIrisPlexS_upload$rs1805006_A)
HIrisPlexS_upload$rs1805006_A<-gsub('A/C', '1', HIrisPlexS_upload$rs1805006_A)
HIrisPlexS_upload$rs1805006_A<-gsub('AC', '1', HIrisPlexS_upload$rs1805006_A)
HIrisPlexS_upload$rs1805006_A<-gsub('C', 'NA', HIrisPlexS_upload$rs1805006_A)
HIrisPlexS_upload$rs1805006_A<-gsub('A', 'NA', HIrisPlexS_upload$rs1805006_A)
HIrisPlexS_upload$rs1805006_A<-gsub('G', 'NA', HIrisPlexS_upload$rs1805006_A)
HIrisPlexS_upload$rs1805006_A<-gsub('T', 'NA', HIrisPlexS_upload$rs1805006_A)
HIrisPlexS_upload$rs1805006_A<-gsub('./.', 'NA', HIrisPlexS_upload$rs1805006_A)
HIrisPlexS_upload$rs1805007_T<-gsub('C/C', '0', HIrisPlexS_upload$rs1805007_T)
HIrisPlexS_upload$rs1805007_T<-gsub('CC', '0', HIrisPlexS_upload$rs1805007_T)
HIrisPlexS_upload$rs1805007_T<-gsub('T/T', '2', HIrisPlexS_upload$rs1805007_T)
HIrisPlexS_upload$rs1805007_T<-gsub('TT', '2', HIrisPlexS_upload$rs1805007_T)
HIrisPlexS_upload$rs1805007_T<-gsub('C/T', '1', HIrisPlexS_upload$rs1805007_T)
HIrisPlexS_upload$rs1805007_T<-gsub('CT', '1', HIrisPlexS_upload$rs1805007_T)
HIrisPlexS_upload$rs1805007_T<-gsub('T/C', '1', HIrisPlexS_upload$rs1805007_T)
HIrisPlexS_upload$rs1805007_T<-gsub('TC', '1', HIrisPlexS_upload$rs1805007_T)
HIrisPlexS_upload$rs1805007_T<-gsub('C', 'NA', HIrisPlexS_upload$rs1805007_T)
HIrisPlexS_upload$rs1805007_T<-gsub('T', 'NA', HIrisPlexS_upload$rs1805007_T)
HIrisPlexS_upload$rs1805007_T<-gsub('G', 'NA', HIrisPlexS_upload$rs1805007_T)
HIrisPlexS_upload$rs1805007_T<-gsub('A', 'NA', HIrisPlexS_upload$rs1805007_T)
HIrisPlexS_upload$rs1805007_T<-gsub('./.', 'NA', HIrisPlexS_upload$rs1805007_T)
HIrisPlexS_upload$rs1805009_C<-gsub('G/G', '0', HIrisPlexS_upload$rs1805009_C)
HIrisPlexS_upload$rs1805009_C<-gsub('GG', '0', HIrisPlexS_upload$rs1805009_C)
HIrisPlexS_upload$rs1805009_C<-gsub('C/C', '2', HIrisPlexS_upload$rs1805009_C)
HIrisPlexS_upload$rs1805009_C<-gsub('CC', '2', HIrisPlexS_upload$rs1805009_C)
HIrisPlexS_upload$rs1805009_C<-gsub('G/C', '1', HIrisPlexS_upload$rs1805009_C)
HIrisPlexS_upload$rs1805009_C<-gsub('GC', '1', HIrisPlexS_upload$rs1805009_C)
HIrisPlexS_upload$rs1805009_C<-gsub('C/G', '1', HIrisPlexS_upload$rs1805009_C)
HIrisPlexS_upload$rs1805009_C<-gsub('CG', '1', HIrisPlexS_upload$rs1805009_C)
HIrisPlexS_upload$rs1805009_C<-gsub('G', 'NA', HIrisPlexS_upload$rs1805009_C)
HIrisPlexS_upload$rs1805009_C<-gsub('C', 'NA', HIrisPlexS_upload$rs1805009_C)
HIrisPlexS_upload$rs1805009_C<-gsub('A', 'NA', HIrisPlexS_upload$rs1805009_C)
HIrisPlexS_upload$rs1805009_C<-gsub('T', 'NA', HIrisPlexS_upload$rs1805009_C)
HIrisPlexS_upload$rs1805009_C<-gsub('./.', 'NA', HIrisPlexS_upload$rs1805009_C)
HIrisPlexS_upload$rs201326893_A<-gsub('C/C', '0', HIrisPlexS_upload$rs201326893_A)
HIrisPlexS_upload$rs201326893_A<-gsub('CC', '0', HIrisPlexS_upload$rs201326893_A)
HIrisPlexS_upload$rs201326893_A<-gsub('A/A', '2', HIrisPlexS_upload$rs201326893_A)
HIrisPlexS_upload$rs201326893_A<-gsub('AA', '2', HIrisPlexS_upload$rs201326893_A)
HIrisPlexS_upload$rs201326893_A<-gsub('C/A', '1', HIrisPlexS_upload$rs201326893_A)
HIrisPlexS_upload$rs201326893_A<-gsub('CA', '1', HIrisPlexS_upload$rs201326893_A)
HIrisPlexS_upload$rs201326893_A<-gsub('A/C', '1', HIrisPlexS_upload$rs201326893_A)
HIrisPlexS_upload$rs201326893_A<-gsub('AC', '1', HIrisPlexS_upload$rs201326893_A)
HIrisPlexS_upload$rs201326893_A<-gsub('C', 'NA', HIrisPlexS_upload$rs201326893_A)
HIrisPlexS_upload$rs201326893_A<-gsub('A', 'NA', HIrisPlexS_upload$rs201326893_A)
HIrisPlexS_upload$rs201326893_A<-gsub('G', 'NA', HIrisPlexS_upload$rs201326893_A)
HIrisPlexS_upload$rs201326893_A<-gsub('T', 'NA', HIrisPlexS_upload$rs201326893_A)
HIrisPlexS_upload$rs201326893_A<-gsub('./.', 'NA', HIrisPlexS_upload$rs201326893_A)
HIrisPlexS_upload$rs2228479_A<-gsub('G/G', '0', HIrisPlexS_upload$rs2228479_A)
HIrisPlexS_upload$rs2228479_A<-gsub('GG', '0', HIrisPlexS_upload$rs2228479_A)
HIrisPlexS_upload$rs2228479_A<-gsub('A/A', '2', HIrisPlexS_upload$rs2228479_A)
HIrisPlexS_upload$rs2228479_A<-gsub('AA', '2', HIrisPlexS_upload$rs2228479_A)
HIrisPlexS_upload$rs2228479_A<-gsub('G/A', '1', HIrisPlexS_upload$rs2228479_A)
HIrisPlexS_upload$rs2228479_A<-gsub('GA', '1', HIrisPlexS_upload$rs2228479_A)
HIrisPlexS_upload$rs2228479_A<-gsub('A/G', '1', HIrisPlexS_upload$rs2228479_A)
HIrisPlexS_upload$rs2228479_A<-gsub('AG', '1', HIrisPlexS_upload$rs2228479_A)
HIrisPlexS_upload$rs2228479_A<-gsub('G', 'NA', HIrisPlexS_upload$rs2228479_A)
HIrisPlexS_upload$rs2228479_A<-gsub('A', 'NA', HIrisPlexS_upload$rs2228479_A)
HIrisPlexS_upload$rs2228479_A<-gsub('C', 'NA', HIrisPlexS_upload$rs2228479_A)
HIrisPlexS_upload$rs2228479_A<-gsub('T', 'NA', HIrisPlexS_upload$rs2228479_A)
HIrisPlexS_upload$rs2228479_A<-gsub('./.', 'NA', HIrisPlexS_upload$rs2228479_A)
HIrisPlexS_upload$rs1110400_C<-gsub('T/T', '0', HIrisPlexS_upload$rs1110400_C)
HIrisPlexS_upload$rs1110400_C<-gsub('TT', '0', HIrisPlexS_upload$rs1110400_C)
HIrisPlexS_upload$rs1110400_C<-gsub('C/C', '2', HIrisPlexS_upload$rs1110400_C)
HIrisPlexS_upload$rs1110400_C<-gsub('CC', '2', HIrisPlexS_upload$rs1110400_C)
HIrisPlexS_upload$rs1110400_C<-gsub('T/C', '1', HIrisPlexS_upload$rs1110400_C)
HIrisPlexS_upload$rs1110400_C<-gsub('TC', '1', HIrisPlexS_upload$rs1110400_C)
HIrisPlexS_upload$rs1110400_C<-gsub('C/T', '1', HIrisPlexS_upload$rs1110400_C)
HIrisPlexS_upload$rs1110400_C<-gsub('CT', '1', HIrisPlexS_upload$rs1110400_C)
HIrisPlexS_upload$rs1110400_C<-gsub('T', 'NA', HIrisPlexS_upload$rs1110400_C)
HIrisPlexS_upload$rs1110400_C<-gsub('C', 'NA', HIrisPlexS_upload$rs1110400_C)
HIrisPlexS_upload$rs1110400_C<-gsub('A', 'NA', HIrisPlexS_upload$rs1110400_C)
HIrisPlexS_upload$rs1110400_C<-gsub('G', 'NA', HIrisPlexS_upload$rs1110400_C)
HIrisPlexS_upload$rs1110400_C<-gsub('./.', 'NA', HIrisPlexS_upload$rs1110400_C)
HIrisPlexS_upload$rs28777_C<-gsub('A/A', '0', HIrisPlexS_upload$rs28777_C)
HIrisPlexS_upload$rs28777_C<-gsub('AA', '0', HIrisPlexS_upload$rs28777_C)
HIrisPlexS_upload$rs28777_C<-gsub('C/C', '2', HIrisPlexS_upload$rs28777_C)
HIrisPlexS_upload$rs28777_C<-gsub('CC', '2', HIrisPlexS_upload$rs28777_C)
HIrisPlexS_upload$rs28777_C<-gsub('A/C', '1', HIrisPlexS_upload$rs28777_C)
HIrisPlexS_upload$rs28777_C<-gsub('AC', '1', HIrisPlexS_upload$rs28777_C)
HIrisPlexS_upload$rs28777_C<-gsub('C/A', '1', HIrisPlexS_upload$rs28777_C)
HIrisPlexS_upload$rs28777_C<-gsub('CA', '1', HIrisPlexS_upload$rs28777_C)
HIrisPlexS_upload$rs28777_C<-gsub('A', 'NA', HIrisPlexS_upload$rs28777_C)
HIrisPlexS_upload$rs28777_C<-gsub('C', 'NA', HIrisPlexS_upload$rs28777_C)
HIrisPlexS_upload$rs28777_C<-gsub('T', 'NA', HIrisPlexS_upload$rs28777_C)
HIrisPlexS_upload$rs28777_C<-gsub('G', 'NA', HIrisPlexS_upload$rs28777_C)
HIrisPlexS_upload$rs28777_C<-gsub('./.', 'NA', HIrisPlexS_upload$rs28777_C)
HIrisPlexS_upload$rs16891982_C<-gsub('G/G', '0', HIrisPlexS_upload$rs16891982_C)
HIrisPlexS_upload$rs16891982_C<-gsub('GG', '0', HIrisPlexS_upload$rs16891982_C)
HIrisPlexS_upload$rs16891982_C<-gsub('C/C', '2', HIrisPlexS_upload$rs16891982_C)
HIrisPlexS_upload$rs16891982_C<-gsub('CC', '2', HIrisPlexS_upload$rs16891982_C)
HIrisPlexS_upload$rs16891982_C<-gsub('G/C', '1', HIrisPlexS_upload$rs16891982_C)
HIrisPlexS_upload$rs16891982_C<-gsub('GC', '1', HIrisPlexS_upload$rs16891982_C)
HIrisPlexS_upload$rs16891982_C<-gsub('C/G', '1', HIrisPlexS_upload$rs16891982_C)
HIrisPlexS_upload$rs16891982_C<-gsub('CG', '1', HIrisPlexS_upload$rs16891982_C)
HIrisPlexS_upload$rs16891982_C<-gsub('G', 'NA', HIrisPlexS_upload$rs16891982_C)
HIrisPlexS_upload$rs16891982_C<-gsub('C', 'NA', HIrisPlexS_upload$rs16891982_C)
HIrisPlexS_upload$rs16891982_C<-gsub('A', 'NA', HIrisPlexS_upload$rs16891982_C)
HIrisPlexS_upload$rs16891982_C<-gsub('T', 'NA', HIrisPlexS_upload$rs16891982_C)
HIrisPlexS_upload$rs16891982_C<-gsub('./.', 'NA', HIrisPlexS_upload$rs16891982_C)
HIrisPlexS_upload$rs12821256_G<-gsub('T/T', '0', HIrisPlexS_upload$rs12821256_G)
HIrisPlexS_upload$rs12821256_G<-gsub('TT', '0', HIrisPlexS_upload$rs12821256_G)
HIrisPlexS_upload$rs12821256_G<-gsub('C/C', '2', HIrisPlexS_upload$rs12821256_G)
HIrisPlexS_upload$rs12821256_G<-gsub('CC', '2', HIrisPlexS_upload$rs12821256_G)
HIrisPlexS_upload$rs12821256_G<-gsub('T/C', '1', HIrisPlexS_upload$rs12821256_G)
HIrisPlexS_upload$rs12821256_G<-gsub('TC', '1', HIrisPlexS_upload$rs12821256_G)
HIrisPlexS_upload$rs12821256_G<-gsub('C/T', '1', HIrisPlexS_upload$rs12821256_G)
HIrisPlexS_upload$rs12821256_G<-gsub('CT', '1', HIrisPlexS_upload$rs12821256_G)
HIrisPlexS_upload$rs12821256_G<-gsub('T', 'NA', HIrisPlexS_upload$rs12821256_G)
HIrisPlexS_upload$rs12821256_G<-gsub('C', 'NA', HIrisPlexS_upload$rs12821256_G)
HIrisPlexS_upload$rs12821256_G<-gsub('A', 'NA', HIrisPlexS_upload$rs12821256_G)
HIrisPlexS_upload$rs12821256_G<-gsub('G', 'NA', HIrisPlexS_upload$rs12821256_G)
HIrisPlexS_upload$rs12821256_G<-gsub('./.', 'NA', HIrisPlexS_upload$rs12821256_G)
HIrisPlexS_upload$rs4959270_A<-gsub('C/C', '0', HIrisPlexS_upload$rs4959270_A)
HIrisPlexS_upload$rs4959270_A<-gsub('CC', '0', HIrisPlexS_upload$rs4959270_A)
HIrisPlexS_upload$rs4959270_A<-gsub('A/A', '2', HIrisPlexS_upload$rs4959270_A)
HIrisPlexS_upload$rs4959270_A<-gsub('AA', '2', HIrisPlexS_upload$rs4959270_A)
HIrisPlexS_upload$rs4959270_A<-gsub('C/A', '1', HIrisPlexS_upload$rs4959270_A)
HIrisPlexS_upload$rs4959270_A<-gsub('CA', '1', HIrisPlexS_upload$rs4959270_A)
HIrisPlexS_upload$rs4959270_A<-gsub('A/C', '1', HIrisPlexS_upload$rs4959270_A)
HIrisPlexS_upload$rs4959270_A<-gsub('AC', '1', HIrisPlexS_upload$rs4959270_A)
HIrisPlexS_upload$rs4959270_A<-gsub('C', 'NA', HIrisPlexS_upload$rs4959270_A)
HIrisPlexS_upload$rs4959270_A<-gsub('A', 'NA', HIrisPlexS_upload$rs4959270_A)
HIrisPlexS_upload$rs4959270_A<-gsub('G', 'NA', HIrisPlexS_upload$rs4959270_A)
HIrisPlexS_upload$rs4959270_A<-gsub('T', 'NA', HIrisPlexS_upload$rs4959270_A)
HIrisPlexS_upload$rs4959270_A<-gsub('./.', 'NA', HIrisPlexS_upload$rs4959270_A)
HIrisPlexS_upload$rs12203592_T<-gsub('C/C', '0', HIrisPlexS_upload$rs12203592_T)
HIrisPlexS_upload$rs12203592_T<-gsub('CC', '0', HIrisPlexS_upload$rs12203592_T)
HIrisPlexS_upload$rs12203592_T<-gsub('T/T', '2', HIrisPlexS_upload$rs12203592_T)
HIrisPlexS_upload$rs12203592_T<-gsub('TT', '2', HIrisPlexS_upload$rs12203592_T)
HIrisPlexS_upload$rs12203592_T<-gsub('C/T', '1', HIrisPlexS_upload$rs12203592_T)
HIrisPlexS_upload$rs12203592_T<-gsub('CT', '1', HIrisPlexS_upload$rs12203592_T)
HIrisPlexS_upload$rs12203592_T<-gsub('T/C', '1', HIrisPlexS_upload$rs12203592_T)
HIrisPlexS_upload$rs12203592_T<-gsub('TC', '1', HIrisPlexS_upload$rs12203592_T)
HIrisPlexS_upload$rs12203592_T<-gsub('C', 'NA', HIrisPlexS_upload$rs12203592_T)
HIrisPlexS_upload$rs12203592_T<-gsub('T', 'NA', HIrisPlexS_upload$rs12203592_T)
HIrisPlexS_upload$rs12203592_T<-gsub('G', 'NA', HIrisPlexS_upload$rs12203592_T)
HIrisPlexS_upload$rs12203592_T<-gsub('A', 'NA', HIrisPlexS_upload$rs12203592_T)
HIrisPlexS_upload$rs12203592_T<-gsub('./.', 'NA', HIrisPlexS_upload$rs12203592_T)
HIrisPlexS_upload$rs1042602_T<-gsub('C/C', '0', HIrisPlexS_upload$rs1042602_T)
HIrisPlexS_upload$rs1042602_T<-gsub('CC', '0', HIrisPlexS_upload$rs1042602_T)
HIrisPlexS_upload$rs1042602_T<-gsub('A/A', '2', HIrisPlexS_upload$rs1042602_T)
HIrisPlexS_upload$rs1042602_T<-gsub('AA', '2', HIrisPlexS_upload$rs1042602_T)
HIrisPlexS_upload$rs1042602_T<-gsub('C/A', '1', HIrisPlexS_upload$rs1042602_T)
HIrisPlexS_upload$rs1042602_T<-gsub('CA', '1', HIrisPlexS_upload$rs1042602_T)
HIrisPlexS_upload$rs1042602_T<-gsub('A/C', '1', HIrisPlexS_upload$rs1042602_T)
HIrisPlexS_upload$rs1042602_T<-gsub('AC', '1', HIrisPlexS_upload$rs1042602_T)
HIrisPlexS_upload$rs1042602_T<-gsub('C', 'NA', HIrisPlexS_upload$rs1042602_T)
HIrisPlexS_upload$rs1042602_T<-gsub('A', 'NA', HIrisPlexS_upload$rs1042602_T)
HIrisPlexS_upload$rs1042602_T<-gsub('G', 'NA', HIrisPlexS_upload$rs1042602_T)
HIrisPlexS_upload$rs1042602_T<-gsub('T', 'NA', HIrisPlexS_upload$rs1042602_T)
HIrisPlexS_upload$rs1042602_T<-gsub('./.', 'NA', HIrisPlexS_upload$rs1042602_T)
HIrisPlexS_upload$rs1800407_A<-gsub('C/C', '0', HIrisPlexS_upload$rs1800407_A)
HIrisPlexS_upload$rs1800407_A<-gsub('CC', '0', HIrisPlexS_upload$rs1800407_A)
HIrisPlexS_upload$rs1800407_A<-gsub('T/T', '2', HIrisPlexS_upload$rs1800407_A)
HIrisPlexS_upload$rs1800407_A<-gsub('TT', '2', HIrisPlexS_upload$rs1800407_A)
HIrisPlexS_upload$rs1800407_A<-gsub('C/T', '1', HIrisPlexS_upload$rs1800407_A)
HIrisPlexS_upload$rs1800407_A<-gsub('CT', '1', HIrisPlexS_upload$rs1800407_A)
HIrisPlexS_upload$rs1800407_A<-gsub('T/C', '1', HIrisPlexS_upload$rs1800407_A)
HIrisPlexS_upload$rs1800407_A<-gsub('TC', '1', HIrisPlexS_upload$rs1800407_A)
HIrisPlexS_upload$rs1800407_A<-gsub('C', 'NA', HIrisPlexS_upload$rs1800407_A)
HIrisPlexS_upload$rs1800407_A<-gsub('T', 'NA', HIrisPlexS_upload$rs1800407_A)
HIrisPlexS_upload$rs1800407_A<-gsub('G', 'NA', HIrisPlexS_upload$rs1800407_A)
HIrisPlexS_upload$rs1800407_A<-gsub('A', 'NA', HIrisPlexS_upload$rs1800407_A)
HIrisPlexS_upload$rs1800407_A<-gsub('./.', 'NA', HIrisPlexS_upload$rs1800407_A)
HIrisPlexS_upload$rs2402130_G<-gsub('A/A', '0', HIrisPlexS_upload$rs2402130_G)
HIrisPlexS_upload$rs2402130_G<-gsub('AA', '0', HIrisPlexS_upload$rs2402130_G)
HIrisPlexS_upload$rs2402130_G<-gsub('G/G', '2', HIrisPlexS_upload$rs2402130_G)
HIrisPlexS_upload$rs2402130_G<-gsub('GG', '2', HIrisPlexS_upload$rs2402130_G)
HIrisPlexS_upload$rs2402130_G<-gsub('G/A', '1', HIrisPlexS_upload$rs2402130_G)
HIrisPlexS_upload$rs2402130_G<-gsub('GA', '1', HIrisPlexS_upload$rs2402130_G)
HIrisPlexS_upload$rs2402130_G<-gsub('A/G', '1', HIrisPlexS_upload$rs2402130_G)
HIrisPlexS_upload$rs2402130_G<-gsub('AG', '1', HIrisPlexS_upload$rs2402130_G)
HIrisPlexS_upload$rs2402130_G<-gsub('A', 'NA', HIrisPlexS_upload$rs2402130_G)
HIrisPlexS_upload$rs2402130_G<-gsub('G', 'NA', HIrisPlexS_upload$rs2402130_G)
HIrisPlexS_upload$rs2402130_G<-gsub('T', 'NA', HIrisPlexS_upload$rs2402130_G)
HIrisPlexS_upload$rs2402130_G<-gsub('C', 'NA', HIrisPlexS_upload$rs2402130_G)
HIrisPlexS_upload$rs2402130_G<-gsub('./.', 'NA', HIrisPlexS_upload$rs2402130_G)
HIrisPlexS_upload$rs12913832_T<-gsub('G/G', '0', HIrisPlexS_upload$rs12913832_T)
HIrisPlexS_upload$rs12913832_T<-gsub('GG', '0', HIrisPlexS_upload$rs12913832_T)
HIrisPlexS_upload$rs12913832_T<-gsub('A/A', '2', HIrisPlexS_upload$rs12913832_T)
HIrisPlexS_upload$rs12913832_T<-gsub('AA', '2', HIrisPlexS_upload$rs12913832_T)
HIrisPlexS_upload$rs12913832_T<-gsub('G/A', '1', HIrisPlexS_upload$rs12913832_T)
HIrisPlexS_upload$rs12913832_T<-gsub('GA', '1', HIrisPlexS_upload$rs12913832_T)
HIrisPlexS_upload$rs12913832_T<-gsub('A/G', '1', HIrisPlexS_upload$rs12913832_T)
HIrisPlexS_upload$rs12913832_T<-gsub('AG', '1', HIrisPlexS_upload$rs12913832_T)
HIrisPlexS_upload$rs12913832_T<-gsub('G', 'NA', HIrisPlexS_upload$rs12913832_T)
HIrisPlexS_upload$rs12913832_T<-gsub('A', 'NA', HIrisPlexS_upload$rs12913832_T)
HIrisPlexS_upload$rs12913832_T<-gsub('C', 'NA', HIrisPlexS_upload$rs12913832_T)
HIrisPlexS_upload$rs12913832_T<-gsub('T', 'NA', HIrisPlexS_upload$rs12913832_T)
HIrisPlexS_upload$rs12913832_T<-gsub('./.', 'NA', HIrisPlexS_upload$rs12913832_T)
HIrisPlexS_upload$rs2378249_C<-gsub('A/A', '0', HIrisPlexS_upload$rs2378249_C)
HIrisPlexS_upload$rs2378249_C<-gsub('AA', '0', HIrisPlexS_upload$rs2378249_C)
HIrisPlexS_upload$rs2378249_C<-gsub('G/G', '2', HIrisPlexS_upload$rs2378249_C)
HIrisPlexS_upload$rs2378249_C<-gsub('GG', '2', HIrisPlexS_upload$rs2378249_C)
HIrisPlexS_upload$rs2378249_C<-gsub('A/G', '1', HIrisPlexS_upload$rs2378249_C)
HIrisPlexS_upload$rs2378249_C<-gsub('AG', '1', HIrisPlexS_upload$rs2378249_C)
HIrisPlexS_upload$rs2378249_C<-gsub('G/A', '1', HIrisPlexS_upload$rs2378249_C)
HIrisPlexS_upload$rs2378249_C<-gsub('GA', '1', HIrisPlexS_upload$rs2378249_C)
HIrisPlexS_upload$rs2378249_C<-gsub('A', 'NA', HIrisPlexS_upload$rs2378249_C)
HIrisPlexS_upload$rs2378249_C<-gsub('G', 'NA', HIrisPlexS_upload$rs2378249_C)
HIrisPlexS_upload$rs2378249_C<-gsub('T', 'NA', HIrisPlexS_upload$rs2378249_C)
HIrisPlexS_upload$rs2378249_C<-gsub('C', 'NA', HIrisPlexS_upload$rs2378249_C)
HIrisPlexS_upload$rs2378249_C<-gsub('./.', 'NA', HIrisPlexS_upload$rs2378249_C)
HIrisPlexS_upload$rs12896399_T<-gsub('G/G', '0', HIrisPlexS_upload$rs12896399_T)
HIrisPlexS_upload$rs12896399_T<-gsub('GG', '0', HIrisPlexS_upload$rs12896399_T)
HIrisPlexS_upload$rs12896399_T<-gsub('T/T', '2', HIrisPlexS_upload$rs12896399_T)
HIrisPlexS_upload$rs12896399_T<-gsub('TT', '2', HIrisPlexS_upload$rs12896399_T)
HIrisPlexS_upload$rs12896399_T<-gsub('G/T', '1', HIrisPlexS_upload$rs12896399_T)
HIrisPlexS_upload$rs12896399_T<-gsub('GT', '1', HIrisPlexS_upload$rs12896399_T)
HIrisPlexS_upload$rs12896399_T<-gsub('T/G', '1', HIrisPlexS_upload$rs12896399_T)
HIrisPlexS_upload$rs12896399_T<-gsub('TG', '1', HIrisPlexS_upload$rs12896399_T)
HIrisPlexS_upload$rs12896399_T<-gsub('G', 'NA', HIrisPlexS_upload$rs12896399_T)
HIrisPlexS_upload$rs12896399_T<-gsub('T', 'NA', HIrisPlexS_upload$rs12896399_T)
HIrisPlexS_upload$rs12896399_T<-gsub('C', 'NA', HIrisPlexS_upload$rs12896399_T)
HIrisPlexS_upload$rs12896399_T<-gsub('A', 'NA', HIrisPlexS_upload$rs12896399_T)
HIrisPlexS_upload$rs12896399_T<-gsub('./.', 'NA', HIrisPlexS_upload$rs12896399_T)
HIrisPlexS_upload$rs1393350_T<-gsub('G/G', '0', HIrisPlexS_upload$rs1393350_T)
HIrisPlexS_upload$rs1393350_T<-gsub('GG', '0', HIrisPlexS_upload$rs1393350_T)
HIrisPlexS_upload$rs1393350_T<-gsub('A/A', '2', HIrisPlexS_upload$rs1393350_T)
HIrisPlexS_upload$rs1393350_T<-gsub('AA', '2', HIrisPlexS_upload$rs1393350_T)
HIrisPlexS_upload$rs1393350_T<-gsub('G/A', '1', HIrisPlexS_upload$rs1393350_T)
HIrisPlexS_upload$rs1393350_T<-gsub('GA', '1', HIrisPlexS_upload$rs1393350_T)
HIrisPlexS_upload$rs1393350_T<-gsub('A/G', '1', HIrisPlexS_upload$rs1393350_T)
HIrisPlexS_upload$rs1393350_T<-gsub('AG', '1', HIrisPlexS_upload$rs1393350_T)
HIrisPlexS_upload$rs1393350_T<-gsub('G', 'NA', HIrisPlexS_upload$rs1393350_T)
HIrisPlexS_upload$rs1393350_T<-gsub('A', 'NA', HIrisPlexS_upload$rs1393350_T)
HIrisPlexS_upload$rs1393350_T<-gsub('C', 'NA', HIrisPlexS_upload$rs1393350_T)
HIrisPlexS_upload$rs1393350_T<-gsub('T', 'NA', HIrisPlexS_upload$rs1393350_T)
HIrisPlexS_upload$rs1393350_T<-gsub('./.', 'NA', HIrisPlexS_upload$rs1393350_T)
HIrisPlexS_upload$rs683_G<-gsub('A/A', '0', HIrisPlexS_upload$rs683_G)
HIrisPlexS_upload$rs683_G<-gsub('AA', '0', HIrisPlexS_upload$rs683_G)
HIrisPlexS_upload$rs683_G<-gsub('C/C', '2', HIrisPlexS_upload$rs683_G)
HIrisPlexS_upload$rs683_G<-gsub('CC', '2', HIrisPlexS_upload$rs683_G)
HIrisPlexS_upload$rs683_G<-gsub('C/A', '1', HIrisPlexS_upload$rs683_G)
HIrisPlexS_upload$rs683_G<-gsub('CA', '1', HIrisPlexS_upload$rs683_G)
HIrisPlexS_upload$rs683_G<-gsub('A/C', '1', HIrisPlexS_upload$rs683_G)
HIrisPlexS_upload$rs683_G<-gsub('AC', '1', HIrisPlexS_upload$rs683_G)
HIrisPlexS_upload$rs683_G<-gsub('A', 'NA', HIrisPlexS_upload$rs683_G)
HIrisPlexS_upload$rs683_G<-gsub('C', 'NA', HIrisPlexS_upload$rs683_G)
HIrisPlexS_upload$rs683_G<-gsub('T', 'NA', HIrisPlexS_upload$rs683_G)
HIrisPlexS_upload$rs683_G<-gsub('G', 'NA', HIrisPlexS_upload$rs683_G)
HIrisPlexS_upload$rs683_G<-gsub('./.', 'NA', HIrisPlexS_upload$rs683_G)
HIrisPlexS_upload$rs3114908_T<-gsub('C/C', '0', HIrisPlexS_upload$rs3114908_T)
HIrisPlexS_upload$rs3114908_T<-gsub('CC', '0', HIrisPlexS_upload$rs3114908_T)
HIrisPlexS_upload$rs3114908_T<-gsub('T/T', '2', HIrisPlexS_upload$rs3114908_T)
HIrisPlexS_upload$rs3114908_T<-gsub('TT', '2', HIrisPlexS_upload$rs3114908_T)
HIrisPlexS_upload$rs3114908_T<-gsub('C/T', '1', HIrisPlexS_upload$rs3114908_T)
HIrisPlexS_upload$rs3114908_T<-gsub('CT', '1', HIrisPlexS_upload$rs3114908_T)
HIrisPlexS_upload$rs3114908_T<-gsub('T/C', '1', HIrisPlexS_upload$rs3114908_T)
HIrisPlexS_upload$rs3114908_T<-gsub('TC', '1', HIrisPlexS_upload$rs3114908_T)
HIrisPlexS_upload$rs3114908_T<-gsub('C', 'NA', HIrisPlexS_upload$rs3114908_T)
HIrisPlexS_upload$rs3114908_T<-gsub('T', 'NA', HIrisPlexS_upload$rs3114908_T)
HIrisPlexS_upload$rs3114908_T<-gsub('G', 'NA', HIrisPlexS_upload$rs3114908_T)
HIrisPlexS_upload$rs3114908_T<-gsub('A', 'NA', HIrisPlexS_upload$rs3114908_T)
HIrisPlexS_upload$rs3114908_T<-gsub('./.', 'NA', HIrisPlexS_upload$rs3114908_T)
HIrisPlexS_upload$rs1800414_C<-gsub('T/T', '0', HIrisPlexS_upload$rs1800414_C)
HIrisPlexS_upload$rs1800414_C<-gsub('TT', '0', HIrisPlexS_upload$rs1800414_C)
HIrisPlexS_upload$rs1800414_C<-gsub('C/C', '2', HIrisPlexS_upload$rs1800414_C)
HIrisPlexS_upload$rs1800414_C<-gsub('CC', '2', HIrisPlexS_upload$rs1800414_C)
HIrisPlexS_upload$rs1800414_C<-gsub('T/C', '1', HIrisPlexS_upload$rs1800414_C)
HIrisPlexS_upload$rs1800414_C<-gsub('TC', '1', HIrisPlexS_upload$rs1800414_C)
HIrisPlexS_upload$rs1800414_C<-gsub('C/T', '1', HIrisPlexS_upload$rs1800414_C)
HIrisPlexS_upload$rs1800414_C<-gsub('CT', '1', HIrisPlexS_upload$rs1800414_C)
HIrisPlexS_upload$rs1800414_C<-gsub('T', 'NA', HIrisPlexS_upload$rs1800414_C)
HIrisPlexS_upload$rs1800414_C<-gsub('C', 'NA', HIrisPlexS_upload$rs1800414_C)
HIrisPlexS_upload$rs1800414_C<-gsub('A', 'NA', HIrisPlexS_upload$rs1800414_C)
HIrisPlexS_upload$rs1800414_C<-gsub('G', 'NA', HIrisPlexS_upload$rs1800414_C)
HIrisPlexS_upload$rs1800414_C<-gsub('./.', 'NA', HIrisPlexS_upload$rs1800414_C)
HIrisPlexS_upload$rs10756819_G<-gsub('A/A', '0', HIrisPlexS_upload$rs10756819_G)
HIrisPlexS_upload$rs10756819_G<-gsub('AA', '0', HIrisPlexS_upload$rs10756819_G)
HIrisPlexS_upload$rs10756819_G<-gsub('G/G', '2', HIrisPlexS_upload$rs10756819_G)
HIrisPlexS_upload$rs10756819_G<-gsub('GG', '2', HIrisPlexS_upload$rs10756819_G)
HIrisPlexS_upload$rs10756819_G<-gsub('A/G', '1', HIrisPlexS_upload$rs10756819_G)
HIrisPlexS_upload$rs10756819_G<-gsub('AG', '1', HIrisPlexS_upload$rs10756819_G)
HIrisPlexS_upload$rs10756819_G<-gsub('G/A', '1', HIrisPlexS_upload$rs10756819_G)
HIrisPlexS_upload$rs10756819_G<-gsub('GA', '1', HIrisPlexS_upload$rs10756819_G)
HIrisPlexS_upload$rs10756819_G<-gsub('A', 'NA', HIrisPlexS_upload$rs10756819_G)
HIrisPlexS_upload$rs10756819_G<-gsub('G', 'NA', HIrisPlexS_upload$rs10756819_G)
HIrisPlexS_upload$rs10756819_G<-gsub('T', 'NA', HIrisPlexS_upload$rs10756819_G)
HIrisPlexS_upload$rs10756819_G<-gsub('C', 'NA', HIrisPlexS_upload$rs10756819_G)
HIrisPlexS_upload$rs10756819_G<-gsub('./.', 'NA', HIrisPlexS_upload$rs10756819_G)
HIrisPlexS_upload$rs2238289_C<-gsub('A/A', '0', HIrisPlexS_upload$rs2238289_C)
HIrisPlexS_upload$rs2238289_C<-gsub('AA', '0', HIrisPlexS_upload$rs2238289_C)
HIrisPlexS_upload$rs2238289_C<-gsub('G/G', '2', HIrisPlexS_upload$rs2238289_C)
HIrisPlexS_upload$rs2238289_C<-gsub('GG', '2', HIrisPlexS_upload$rs2238289_C)
HIrisPlexS_upload$rs2238289_C<-gsub('A/G', '1', HIrisPlexS_upload$rs2238289_C)
HIrisPlexS_upload$rs2238289_C<-gsub('AG', '1', HIrisPlexS_upload$rs2238289_C)
HIrisPlexS_upload$rs2238289_C<-gsub('G/A', '1', HIrisPlexS_upload$rs2238289_C)
HIrisPlexS_upload$rs2238289_C<-gsub('GA', '1', HIrisPlexS_upload$rs2238289_C)
HIrisPlexS_upload$rs2238289_C<-gsub('A', 'NA', HIrisPlexS_upload$rs2238289_C)
HIrisPlexS_upload$rs2238289_C<-gsub('G', 'NA', HIrisPlexS_upload$rs2238289_C)
HIrisPlexS_upload$rs2238289_C<-gsub('T', 'NA', HIrisPlexS_upload$rs2238289_C)
HIrisPlexS_upload$rs2238289_C<-gsub('C', 'NA', HIrisPlexS_upload$rs2238289_C)
HIrisPlexS_upload$rs2238289_C<-gsub('./.', 'NA', HIrisPlexS_upload$rs2238289_C)
HIrisPlexS_upload$rs17128291_C<-gsub('A/A', '0', HIrisPlexS_upload$rs17128291_C)
HIrisPlexS_upload$rs17128291_C<-gsub('AA', '0', HIrisPlexS_upload$rs17128291_C)
HIrisPlexS_upload$rs17128291_C<-gsub('G/G', '2', HIrisPlexS_upload$rs17128291_C)
HIrisPlexS_upload$rs17128291_C<-gsub('GG', '2', HIrisPlexS_upload$rs17128291_C)
HIrisPlexS_upload$rs17128291_C<-gsub('A/G', '1', HIrisPlexS_upload$rs17128291_C)
HIrisPlexS_upload$rs17128291_C<-gsub('AG', '1', HIrisPlexS_upload$rs17128291_C)
HIrisPlexS_upload$rs17128291_C<-gsub('G/A', '1', HIrisPlexS_upload$rs17128291_C)
HIrisPlexS_upload$rs17128291_C<-gsub('GA', '1', HIrisPlexS_upload$rs17128291_C)
HIrisPlexS_upload$rs17128291_C<-gsub('A', 'NA', HIrisPlexS_upload$rs17128291_C)
HIrisPlexS_upload$rs17128291_C<-gsub('G', 'NA', HIrisPlexS_upload$rs17128291_C)
HIrisPlexS_upload$rs17128291_C<-gsub('T', 'NA', HIrisPlexS_upload$rs17128291_C)
HIrisPlexS_upload$rs17128291_C<-gsub('C', 'NA', HIrisPlexS_upload$rs17128291_C)
HIrisPlexS_upload$rs17128291_C<-gsub('./.', 'NA', HIrisPlexS_upload$rs17128291_C)
HIrisPlexS_upload$rs6497292_C<-gsub('A/A', '0', HIrisPlexS_upload$rs6497292_C)
HIrisPlexS_upload$rs6497292_C<-gsub('AA', '0', HIrisPlexS_upload$rs6497292_C)
HIrisPlexS_upload$rs6497292_C<-gsub('G/G', '2', HIrisPlexS_upload$rs6497292_C)
HIrisPlexS_upload$rs6497292_C<-gsub('GG', '2', HIrisPlexS_upload$rs6497292_C)
HIrisPlexS_upload$rs6497292_C<-gsub('A/G', '1', HIrisPlexS_upload$rs6497292_C)
HIrisPlexS_upload$rs6497292_C<-gsub('AG', '1', HIrisPlexS_upload$rs6497292_C)
HIrisPlexS_upload$rs6497292_C<-gsub('G/A', '1', HIrisPlexS_upload$rs6497292_C)
HIrisPlexS_upload$rs6497292_C<-gsub('GA', '1', HIrisPlexS_upload$rs6497292_C)
HIrisPlexS_upload$rs6497292_C<-gsub('A', 'NA', HIrisPlexS_upload$rs6497292_C)
HIrisPlexS_upload$rs6497292_C<-gsub('G', 'NA', HIrisPlexS_upload$rs6497292_C)
HIrisPlexS_upload$rs6497292_C<-gsub('T', 'NA', HIrisPlexS_upload$rs6497292_C)
HIrisPlexS_upload$rs6497292_C<-gsub('C', 'NA', HIrisPlexS_upload$rs6497292_C)
HIrisPlexS_upload$rs6497292_C<-gsub('./.', 'NA', HIrisPlexS_upload$rs6497292_C)
HIrisPlexS_upload$rs1129038_G<-gsub('T/T', '0', HIrisPlexS_upload$rs1129038_G)
HIrisPlexS_upload$rs1129038_G<-gsub('TT', '0', HIrisPlexS_upload$rs1129038_G)
HIrisPlexS_upload$rs1129038_G<-gsub('C/C', '2', HIrisPlexS_upload$rs1129038_G)
HIrisPlexS_upload$rs1129038_G<-gsub('CC', '2', HIrisPlexS_upload$rs1129038_G)
HIrisPlexS_upload$rs1129038_G<-gsub('T/C', '1', HIrisPlexS_upload$rs1129038_G)
HIrisPlexS_upload$rs1129038_G<-gsub('TC', '1', HIrisPlexS_upload$rs1129038_G)
HIrisPlexS_upload$rs1129038_G<-gsub('C/T', '1', HIrisPlexS_upload$rs1129038_G)
HIrisPlexS_upload$rs1129038_G<-gsub('CT', '1', HIrisPlexS_upload$rs1129038_G)
HIrisPlexS_upload$rs1129038_G<-gsub('T', 'NA', HIrisPlexS_upload$rs1129038_G)
HIrisPlexS_upload$rs1129038_G<-gsub('C', 'NA', HIrisPlexS_upload$rs1129038_G)
HIrisPlexS_upload$rs1129038_G<-gsub('A', 'NA', HIrisPlexS_upload$rs1129038_G)
HIrisPlexS_upload$rs1129038_G<-gsub('G', 'NA', HIrisPlexS_upload$rs1129038_G)
HIrisPlexS_upload$rs1129038_G<-gsub('./.', 'NA', HIrisPlexS_upload$rs1129038_G)
HIrisPlexS_upload$rs1667394_C<-gsub('T/T', '0', HIrisPlexS_upload$rs1667394_C)
HIrisPlexS_upload$rs1667394_C<-gsub('TT', '0', HIrisPlexS_upload$rs1667394_C)
HIrisPlexS_upload$rs1667394_C<-gsub('C/C', '2', HIrisPlexS_upload$rs1667394_C)
HIrisPlexS_upload$rs1667394_C<-gsub('CC', '2', HIrisPlexS_upload$rs1667394_C)
HIrisPlexS_upload$rs1667394_C<-gsub('T/C', '1', HIrisPlexS_upload$rs1667394_C)
HIrisPlexS_upload$rs1667394_C<-gsub('TC', '1', HIrisPlexS_upload$rs1667394_C)
HIrisPlexS_upload$rs1667394_C<-gsub('C/T', '1', HIrisPlexS_upload$rs1667394_C)
HIrisPlexS_upload$rs1667394_C<-gsub('CT', '1', HIrisPlexS_upload$rs1667394_C)
HIrisPlexS_upload$rs1667394_C<-gsub('T', 'NA', HIrisPlexS_upload$rs1667394_C)
HIrisPlexS_upload$rs1667394_C<-gsub('C', 'NA', HIrisPlexS_upload$rs1667394_C)
HIrisPlexS_upload$rs1667394_C<-gsub('A', 'NA', HIrisPlexS_upload$rs1667394_C)
HIrisPlexS_upload$rs1667394_C<-gsub('G', 'NA', HIrisPlexS_upload$rs1667394_C)
HIrisPlexS_upload$rs1667394_C<-gsub('./.', 'NA', HIrisPlexS_upload$rs1667394_C)
HIrisPlexS_upload$rs1126809_A<-gsub('G/G', '0', HIrisPlexS_upload$rs1126809_A)
HIrisPlexS_upload$rs1126809_A<-gsub('GG', '0', HIrisPlexS_upload$rs1126809_A)
HIrisPlexS_upload$rs1126809_A<-gsub('A/A', '2', HIrisPlexS_upload$rs1126809_A)
HIrisPlexS_upload$rs1126809_A<-gsub('AA', '2', HIrisPlexS_upload$rs1126809_A)
HIrisPlexS_upload$rs1126809_A<-gsub('G/A', '1', HIrisPlexS_upload$rs1126809_A)
HIrisPlexS_upload$rs1126809_A<-gsub('GA', '1', HIrisPlexS_upload$rs1126809_A)
HIrisPlexS_upload$rs1126809_A<-gsub('A/G', '1', HIrisPlexS_upload$rs1126809_A)
HIrisPlexS_upload$rs1126809_A<-gsub('AG', '1', HIrisPlexS_upload$rs1126809_A)
HIrisPlexS_upload$rs1126809_A<-gsub('G', 'NA', HIrisPlexS_upload$rs1126809_A)
HIrisPlexS_upload$rs1126809_A<-gsub('A', 'NA', HIrisPlexS_upload$rs1126809_A)
HIrisPlexS_upload$rs1126809_A<-gsub('C', 'NA', HIrisPlexS_upload$rs1126809_A)
HIrisPlexS_upload$rs1126809_A<-gsub('T', 'NA', HIrisPlexS_upload$rs1126809_A)
HIrisPlexS_upload$rs1126809_A<-gsub('./.', 'NA', HIrisPlexS_upload$rs1126809_A)
HIrisPlexS_upload$rs1470608_A<-gsub('G/G', '0', HIrisPlexS_upload$rs1470608_A)
HIrisPlexS_upload$rs1470608_A<-gsub('GG', '0', HIrisPlexS_upload$rs1470608_A)
HIrisPlexS_upload$rs1470608_A<-gsub('T/T', '2', HIrisPlexS_upload$rs1470608_A)
HIrisPlexS_upload$rs1470608_A<-gsub('TT', '2', HIrisPlexS_upload$rs1470608_A)
HIrisPlexS_upload$rs1470608_A<-gsub('G/T', '1', HIrisPlexS_upload$rs1470608_A)
HIrisPlexS_upload$rs1470608_A<-gsub('GT', '1', HIrisPlexS_upload$rs1470608_A)
HIrisPlexS_upload$rs1470608_A<-gsub('T/G', '1', HIrisPlexS_upload$rs1470608_A)
HIrisPlexS_upload$rs1470608_A<-gsub('TG', '1', HIrisPlexS_upload$rs1470608_A)
HIrisPlexS_upload$rs1470608_A<-gsub('G', 'NA', HIrisPlexS_upload$rs1470608_A)
HIrisPlexS_upload$rs1470608_A<-gsub('T', 'NA', HIrisPlexS_upload$rs1470608_A)
HIrisPlexS_upload$rs1470608_A<-gsub('C', 'NA', HIrisPlexS_upload$rs1470608_A)
HIrisPlexS_upload$rs1470608_A<-gsub('A', 'NA', HIrisPlexS_upload$rs1470608_A)
HIrisPlexS_upload$rs1470608_A<-gsub('./.', 'NA', HIrisPlexS_upload$rs1470608_A)
HIrisPlexS_upload$rs1426654_G<-gsub('A/A', '0', HIrisPlexS_upload$rs1426654_G)
HIrisPlexS_upload$rs1426654_G<-gsub('AA', '0', HIrisPlexS_upload$rs1426654_G)
HIrisPlexS_upload$rs1426654_G<-gsub('G/G', '2', HIrisPlexS_upload$rs1426654_G)
HIrisPlexS_upload$rs1426654_G<-gsub('GG', '2', HIrisPlexS_upload$rs1426654_G)
HIrisPlexS_upload$rs1426654_G<-gsub('G/A', '1', HIrisPlexS_upload$rs1426654_G)
HIrisPlexS_upload$rs1426654_G<-gsub('GA', '1', HIrisPlexS_upload$rs1426654_G)
HIrisPlexS_upload$rs1426654_G<-gsub('A/G', '1', HIrisPlexS_upload$rs1426654_G)
HIrisPlexS_upload$rs1426654_G<-gsub('AG', '1', HIrisPlexS_upload$rs1426654_G)
HIrisPlexS_upload$rs1426654_G<-gsub('A', 'NA', HIrisPlexS_upload$rs1426654_G)
HIrisPlexS_upload$rs1426654_G<-gsub('G', 'NA', HIrisPlexS_upload$rs1426654_G)
HIrisPlexS_upload$rs1426654_G<-gsub('T', 'NA', HIrisPlexS_upload$rs1426654_G)
HIrisPlexS_upload$rs1426654_G<-gsub('C', 'NA', HIrisPlexS_upload$rs1426654_G)
HIrisPlexS_upload$rs1426654_G<-gsub('./.', 'NA', HIrisPlexS_upload$rs1426654_G)
HIrisPlexS_upload$rs6119471_C<-gsub('C/C', '0', HIrisPlexS_upload$rs6119471_C)
HIrisPlexS_upload$rs6119471_C<-gsub('CC', '0', HIrisPlexS_upload$rs6119471_C)
HIrisPlexS_upload$rs6119471_C<-gsub('G/G', '2', HIrisPlexS_upload$rs6119471_C)
HIrisPlexS_upload$rs6119471_C<-gsub('GG', '2', HIrisPlexS_upload$rs6119471_C)
HIrisPlexS_upload$rs6119471_C<-gsub('C/G', '1', HIrisPlexS_upload$rs6119471_C)
HIrisPlexS_upload$rs6119471_C<-gsub('CG', '1', HIrisPlexS_upload$rs6119471_C)
HIrisPlexS_upload$rs6119471_C<-gsub('G/C', '1', HIrisPlexS_upload$rs6119471_C)
HIrisPlexS_upload$rs6119471_C<-gsub('GC', '1', HIrisPlexS_upload$rs6119471_C)
HIrisPlexS_upload$rs6119471_C<-gsub('C', 'NA', HIrisPlexS_upload$rs6119471_C)
HIrisPlexS_upload$rs6119471_C<-gsub('G', 'NA', HIrisPlexS_upload$rs6119471_C)
HIrisPlexS_upload$rs6119471_C<-gsub('A', 'NA', HIrisPlexS_upload$rs6119471_C)
HIrisPlexS_upload$rs6119471_C<-gsub('T', 'NA', HIrisPlexS_upload$rs6119471_C)
HIrisPlexS_upload$rs6119471_C<-gsub('./.', 'NA', HIrisPlexS_upload$rs6119471_C)
HIrisPlexS_upload$rs1545397_T<-gsub('A/A', '0', HIrisPlexS_upload$rs1545397_T)
HIrisPlexS_upload$rs1545397_T<-gsub('AA', '0', HIrisPlexS_upload$rs1545397_T)
HIrisPlexS_upload$rs1545397_T<-gsub('T/T', '2', HIrisPlexS_upload$rs1545397_T)
HIrisPlexS_upload$rs1545397_T<-gsub('TT', '2', HIrisPlexS_upload$rs1545397_T)
HIrisPlexS_upload$rs1545397_T<-gsub('A/T', '1', HIrisPlexS_upload$rs1545397_T)
HIrisPlexS_upload$rs1545397_T<-gsub('AT', '1', HIrisPlexS_upload$rs1545397_T)
HIrisPlexS_upload$rs1545397_T<-gsub('T/A', '1', HIrisPlexS_upload$rs1545397_T)
HIrisPlexS_upload$rs1545397_T<-gsub('TA', '1', HIrisPlexS_upload$rs1545397_T)
HIrisPlexS_upload$rs1545397_T<-gsub('A', 'NA', HIrisPlexS_upload$rs1545397_T)
HIrisPlexS_upload$rs1545397_T<-gsub('T', 'NA', HIrisPlexS_upload$rs1545397_T)
HIrisPlexS_upload$rs1545397_T<-gsub('C', 'NA', HIrisPlexS_upload$rs1545397_T)
HIrisPlexS_upload$rs1545397_T<-gsub('G', 'NA', HIrisPlexS_upload$rs1545397_T)
HIrisPlexS_upload$rs1545397_T<-gsub('./.', 'NA', HIrisPlexS_upload$rs1545397_T)
HIrisPlexS_upload$rs6059655_T<-gsub('G/G', '0', HIrisPlexS_upload$rs6059655_T)
HIrisPlexS_upload$rs6059655_T<-gsub('GG', '0', HIrisPlexS_upload$rs6059655_T)
HIrisPlexS_upload$rs6059655_T<-gsub('A/A', '2', HIrisPlexS_upload$rs6059655_T)
HIrisPlexS_upload$rs6059655_T<-gsub('AA', '2', HIrisPlexS_upload$rs6059655_T)
HIrisPlexS_upload$rs6059655_T<-gsub('G/A', '1', HIrisPlexS_upload$rs6059655_T)
HIrisPlexS_upload$rs6059655_T<-gsub('GA', '1', HIrisPlexS_upload$rs6059655_T)
HIrisPlexS_upload$rs6059655_T<-gsub('A/G', '1', HIrisPlexS_upload$rs6059655_T)
HIrisPlexS_upload$rs6059655_T<-gsub('AG', '1', HIrisPlexS_upload$rs6059655_T)
HIrisPlexS_upload$rs6059655_T<-gsub('G', 'NA', HIrisPlexS_upload$rs6059655_T)
HIrisPlexS_upload$rs6059655_T<-gsub('A', 'NA', HIrisPlexS_upload$rs6059655_T)
HIrisPlexS_upload$rs6059655_T<-gsub('C', 'NA', HIrisPlexS_upload$rs6059655_T)
HIrisPlexS_upload$rs6059655_T<-gsub('T', 'NA', HIrisPlexS_upload$rs6059655_T)
HIrisPlexS_upload$rs6059655_T<-gsub('./.', 'NA', HIrisPlexS_upload$rs6059655_T)
HIrisPlexS_upload$rs12441727_A<-gsub('G/G', '0', HIrisPlexS_upload$rs12441727_A)
HIrisPlexS_upload$rs12441727_A<-gsub('GG', '0', HIrisPlexS_upload$rs12441727_A)
HIrisPlexS_upload$rs12441727_A<-gsub('A/A', '2', HIrisPlexS_upload$rs12441727_A)
HIrisPlexS_upload$rs12441727_A<-gsub('AA', '2', HIrisPlexS_upload$rs12441727_A)
HIrisPlexS_upload$rs12441727_A<-gsub('G/A', '1', HIrisPlexS_upload$rs12441727_A)
HIrisPlexS_upload$rs12441727_A<-gsub('GA', '1', HIrisPlexS_upload$rs12441727_A)
HIrisPlexS_upload$rs12441727_A<-gsub('A/G', '1', HIrisPlexS_upload$rs12441727_A)
HIrisPlexS_upload$rs12441727_A<-gsub('AG', '1', HIrisPlexS_upload$rs12441727_A)
HIrisPlexS_upload$rs12441727_A<-gsub('G', 'NA', HIrisPlexS_upload$rs12441727_A)
HIrisPlexS_upload$rs12441727_A<-gsub('A', 'NA', HIrisPlexS_upload$rs12441727_A)
HIrisPlexS_upload$rs12441727_A<-gsub('C', 'NA', HIrisPlexS_upload$rs12441727_A)
HIrisPlexS_upload$rs12441727_A<-gsub('T', 'NA', HIrisPlexS_upload$rs12441727_A)
HIrisPlexS_upload$rs12441727_A<-gsub('./.', 'NA', HIrisPlexS_upload$rs12441727_A)
HIrisPlexS_upload$rs3212355_A<-gsub('C/C', '0', HIrisPlexS_upload$rs3212355_A)
HIrisPlexS_upload$rs3212355_A<-gsub('CC', '0', HIrisPlexS_upload$rs3212355_A)
HIrisPlexS_upload$rs3212355_A<-gsub('T/T', '2', HIrisPlexS_upload$rs3212355_A)
HIrisPlexS_upload$rs3212355_A<-gsub('TT', '2', HIrisPlexS_upload$rs3212355_A)
HIrisPlexS_upload$rs3212355_A<-gsub('C/T', '1', HIrisPlexS_upload$rs3212355_A)
HIrisPlexS_upload$rs3212355_A<-gsub('CT', '1', HIrisPlexS_upload$rs3212355_A)
HIrisPlexS_upload$rs3212355_A<-gsub('T/C', '1', HIrisPlexS_upload$rs3212355_A)
HIrisPlexS_upload$rs3212355_A<-gsub('TC', '1', HIrisPlexS_upload$rs3212355_A)
HIrisPlexS_upload$rs3212355_A<-gsub('C', 'NA', HIrisPlexS_upload$rs3212355_A)
HIrisPlexS_upload$rs3212355_A<-gsub('T', 'NA', HIrisPlexS_upload$rs3212355_A)
HIrisPlexS_upload$rs3212355_A<-gsub('G', 'NA', HIrisPlexS_upload$rs3212355_A)
HIrisPlexS_upload$rs3212355_A<-gsub('A', 'NA', HIrisPlexS_upload$rs3212355_A)
HIrisPlexS_upload$rs3212355_A<-gsub('./.', 'NA', HIrisPlexS_upload$rs3212355_A)
HIrisPlexS_upload$rs8051733_C<-gsub('A/A', '0', HIrisPlexS_upload$rs8051733_C)
HIrisPlexS_upload$rs8051733_C<-gsub('AA', '0', HIrisPlexS_upload$rs8051733_C)
HIrisPlexS_upload$rs8051733_C<-gsub('G/G', '2', HIrisPlexS_upload$rs8051733_C)
HIrisPlexS_upload$rs8051733_C<-gsub('GG', '2', HIrisPlexS_upload$rs8051733_C)
HIrisPlexS_upload$rs8051733_C<-gsub('G/A', '1', HIrisPlexS_upload$rs8051733_C)
HIrisPlexS_upload$rs8051733_C<-gsub('GA', '1', HIrisPlexS_upload$rs8051733_C)
HIrisPlexS_upload$rs8051733_C<-gsub('A/G', '1', HIrisPlexS_upload$rs8051733_C)
HIrisPlexS_upload$rs8051733_C<-gsub('AG', '1', HIrisPlexS_upload$rs8051733_C)
HIrisPlexS_upload$rs8051733_C<-gsub('A', 'NA', HIrisPlexS_upload$rs8051733_C)
HIrisPlexS_upload$rs8051733_C<-gsub('G', 'NA', HIrisPlexS_upload$rs8051733_C)
HIrisPlexS_upload$rs8051733_C<-gsub('T', 'NA', HIrisPlexS_upload$rs8051733_C)
HIrisPlexS_upload$rs8051733_C<-gsub('C', 'NA', HIrisPlexS_upload$rs8051733_C)
HIrisPlexS_upload$rs8051733_C<-gsub('./.', 'NA', HIrisPlexS_upload$rs8051733_C)

# Save transposed_row to CSV file
write.csv(HIrisPlexS_upload, file = paste0(output_dir, "HIrisPlexS_Estonian_upload.csv"), row.names = FALSE, quote = FALSE)

# Print confirmation message
cat("File 'HIrisPlexS_upload.csv' has been saved successfully in", output_dir, "\n")
cat("You can upload your file to https://hirisplex.erasmusmc.nl/ to get HIrisPlex-S predictions.\n")
