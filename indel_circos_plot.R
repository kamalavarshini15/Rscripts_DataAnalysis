library(vcfR)
library(BioCircos)
setwd("E:/CAG/frga_4_datasets/POONGARIRGC_28611-1")

# Load the VCF file
vcf_file <- read.vcfR("ERR606295_InDels.vcf")

# Extract the fixed data (e.g., CHROM, POS, REF, ALT)
vcf_data <- as.data.frame(vcf_file@fix)


# Add columns to identify the type and length of the indel
vcf_data$type <- ifelse(nchar(vcf_data$REF) < nchar(vcf_data$ALT), "Insertion", 
                        ifelse(nchar(vcf_data$REF) > nchar(vcf_data$ALT), "Deletion", NA))

vcf_data$length <- abs(nchar(vcf_data$REF) - nchar(vcf_data$ALT))

# Filter to retain only rows with indels
indels_data <- vcf_data[!is.na(vcf_data$type), ]

indels_data_df <- as.data.frame(indels_data)

write.csv(indels_data_df,"indels_data_ERR606295.csv",row.names = FALSE)

indel_chromosomes <- as.character(indels_data$CHROM)

indel_positions <- as.numeric(indels_data$POS)


indel_types <- indels_data$type


indel_lengths <- indels_data$length

#mapping chromosome 
indel_chr_updated <- indel_chromosomes
indel_chr_updated[indel_chromosomes == "FRIGIDA_Like"] <- "FRIGIDA_Like protein"
indel_chr_updated[indel_chromosomes == "Fibronectin_type-III_domain-containing_protein"] <- "Fibronectin_type-III_domain-containing_protein"
indel_chr_updated[indel_chromosomes == "Fibronectin"] <- "Fibronectin type III domain containing protein, expressed"
indel_chr_updated[indel_chromosomes == "Photoperiod"] <- "Photoperiod independent early flowering1"
indel_chr_updated[indel_chromosomes == "GPCR-type"] <- "GPCR-type G protein COLD1"
indel_chr_updated[indel_chromosomes == "CURLY"] <- "CURLY FLAG LEAF 1"
indel_chr_updated[indel_chromosomes == "Hd3a"] <- "Hd3a"
indel_chr_updated[indel_chromosomes == "RFT1"] <- "RFT1"
indel_chr_updated[indel_chromosomes == "OsMADS14"] <- "OsMADS14"
indel_chr_updated[indel_chromosomes == "OsMADS15"] <- "OsMADS15"
indel_chr_updated[indel_chromosomes == "Hd1"] <- "Hd1"
indel_chr_updated[indel_chromosomes == "Ehd1"] <- "Ehd1"
indel_chr_updated[indel_chromosomes == "Ghd7"] <- "Ghd7"
indel_chr_updated[indel_chromosomes == "Ghd8"] <- "Ghd8"
indel_chr_updated[indel_chromosomes == "Hd6"] <- "Hd6"
indel_chr_updated[indel_chromosomes == "Hd16"] <- "Hd16"
indel_chr_updated[indel_chromosomes == "OsMADS50"] <- "OsMADS50"
indel_chr_updated[indel_chromosomes == "OsMADS56"] <- "OsMADS56"
indel_chr_updated[indel_chromosomes == "OsMADS51"] <- "OsMADS51"
indel_chr_updated[indel_chromosomes == "Gigantea"] <- "Gigantea"
indel_chr_updated[indel_chromosomes == "OsPRR37"] <- "OsPRR37"
indel_chr_updated[indel_chromosomes == "Hd17"] <- "Hd17"
indel_chr_updated[indel_chromosomes == "OsPhyA"] <- "OsPhyA"
indel_chr_updated[indel_chromosomes == "Ehd3"] <- "Ehd3"


CCGenome = list("FRIGIDA_Like protein" = 3537,
                "Fibronectin_type-III_domain-containing_protein" = 4365,
                "Fibronectin type III domain containing protein, expressed" = 3352,
                "Photoperiod independent early flowering1" = 5886,
                "GPCR-type G protein COLD1" = 4785,
                "CURLY FLAG LEAF 1" = 2753,
                "Hd3a" = 2449,
                "RFT1" = 1652,
                "OsMADS14" = 10424,
                "OsMADS15" = 6482,
                "Hd1" = 2285,
                "Ehd1" = 1882,
                "Ghd7" = 2784,
                "Ghd8" = 1718,
                "Hd6" = 5710,
                "Hd16" = 7397,
                "OsMADS50" = 5384,
                "OsMADS56" = 10846,
                "OsMADS51" = 20256,
                "Gigantea" = 2286,
                "OsPRR37" = 12519,
                "Hd17" = 5044,
                "OsPhyA" = 7971,
                "Ehd3" = 4125
)


# Assign values based on type (Insertion = 1, Deletion = -1)
#indel_values <- ifelse(indel_types == "Insertion", 1, -1
# Assign colors based on the type
#indel_colors <- ifelse(indel_types == "Insertion", "green", "red")
points_chromosomes <- indel_chr_updated

points_coordinates <- indel_positions

points_values = 0:4

indel_track <- BioCircosSNPTrack("indelTrack", 
                                 #indel_chromosomes, 
                                 #indel_positions, 
                                 points_chromosomes,
                                 points_coordinates,
                                 points_values, 
                                 colors = c("green", "red"), 
                                 minRadius = 0.5, 
                                 maxRadius = 0.9)

tracklist = indel_track + BioCircosBackgroundTrack("myBackgroundTrack", 
                                                   minRadius = 0.5, maxRadius = 0.9,
                                                   borderColors = "#AAAAAA", borderSize = 0.6, fillColors = "#B3E6FF") 


BioCircos(tracklist, genome = CCGenome, genomeFillColor = "PuOr",
          chrPad = 0.05, displayGenomeBorder = FALSE, yChr =  FALSE, genomeLabelOrientation = 90,
          genomeTicksDisplay = FALSE,  genomeLabelTextSize = 10, genomeLabelDy = 0)