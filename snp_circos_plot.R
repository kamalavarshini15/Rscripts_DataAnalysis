setwd("E:/CAG/frga_4_datasets/POONGARIRGC_28611-1")

library(vcfR)

#Biocricos plot new 
#update.packages("VariantAnnotation")
#install.packages("RSQLite")

# Load required libraries
library(BioCircos)
library(Biostrings)
library(VariantAnnotation)




#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#BiocManager::install("VariantAnnotation")



vcf_file_snp <- read.vcfR("ERR606282_SNPs.vcf")


vcf_data_snp <- as.data.frame(vcf_file_snp@fix)

snps_data_df <- as.data.frame(vcf_data_snp)

write.csv(snps_data_df,"SNPS_data_ERR606282.csv",row.names = FALSE)



# Plotting 



# Process reference genome
fasta_file <- "E:/CAG/24gene_index/24gene_sequence.fasta"

# Read the FASTA file to get chromosome lengths
genome <- readDNAStringSet(fasta_file)
chr_lengths <- width(genome)
names(chr_lengths) <- names(genome)# Chromosome names
print(names(chr_lengths) <- names(genome))

# Prepare genome data for BioCircos
custom_genome <- list()
for (i in seq_along(chr_lengths)) {
  custom_genome[[names(chr_lengths)[i]]] <- chr_lengths[i]
}
custom_genome <- as.list(chr_lengths)  # Convert named vector to a named list


#Process SNP data
vcf_file <- "E:/CAG/frga_4_datasets/POONGARIRGC_28611-1/ERR606295_SNPs.vcf"

# Read the VCF file
vcf <- readVcf(vcf_file, genome = "Rice")
snp_data <- rowRanges(vcf)

snp_data_df <- as.data.frame(snp_data)




# Extract chromosome and position
snp_chr <- as.character(seqnames(snp_data))

snp_pos <- start(snp_data)




#mapping chromosome 
snp_chr_updated <- snp_chr
snp_chr_updated[snp_chr == "FRIGIDA_Like"] <- "FRIGIDA_Like protein"
snp_chr_updated[snp_chr == "Fibronectin_type-III_domain-containing_protein"] <- "Fibronectin_type-III_domain-containing_protein"
snp_chr_updated[snp_chr == "Fibronectin"] <- "Fibronectin type III domain containing protein, expressed"
snp_chr_updated[snp_chr == "Photoperiod"] <- "Photoperiod independent early flowering1"
snp_chr_updated[snp_chr == "GPCR-type"] <- "GPCR-type G protein COLD1"
snp_chr_updated[snp_chr == "CURLY"] <- "CURLY FLAG LEAF 1"
snp_chr_updated[snp_chr == "Hd3a"] <- "Hd3a"
snp_chr_updated[snp_chr == "RFT1"] <- "RFT1"
snp_chr_updated[snp_chr == "OsMADS14"] <- "OsMADS14"
snp_chr_updated[snp_chr == "OsMADS15"] <- "OsMADS15"
snp_chr_updated[snp_chr == "Hd1"] <- "Hd1"
snp_chr_updated[snp_chr == "Ehd1"] <- "Ehd1"
snp_chr_updated[snp_chr == "Ghd7"] <- "Ghd7"
snp_chr_updated[snp_chr == "Ghd8"] <- "Ghd8"
snp_chr_updated[snp_chr == "Hd6"] <- "Hd6"
snp_chr_updated[snp_chr == "Hd16"] <- "Hd16"
snp_chr_updated[snp_chr == "OsMADS50"] <- "OsMADS50"
snp_chr_updated[snp_chr == "OsMADS56"] <- "OsMADS56"
snp_chr_updated[snp_chr == "OsMADS51"] <- "OsMADS51"
snp_chr_updated[snp_chr == "Gigantea"] <- "Gigantea"
snp_chr_updated[snp_chr == "OsPRR37"] <- "OsPRR37"
snp_chr_updated[snp_chr == "Hd17"] <- "Hd17"
snp_chr_updated[snp_chr == "OsPhyA"] <- "OsPhyA"
snp_chr_updated[snp_chr == "Ehd3"] <- "Ehd3"




#print(length(names(chr_lengths) <- names(genome)))# Chromosome names

#table(snp_chr_updated) 
#chromosome_of_interest <- "Ehd3"

# Filter SNPs for the selected chromosome
#snp_positions_chr <- snp_pos[snp_chr_updated == chromosome_of_interest]

# Display the SNP coordinates
#print(snp_positions_chr)



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


#BioCircos(genome = CCGenome, genomeFillColor = c("tomato2", "darkblue"),
#genomeTicksScale = 4e+3)

# Chromosomes on which the points should be displayed

points_chromosomes <- snp_chr_updated

points_coordinates <- snp_pos

print(length(snp_pos))
print(length(snp_chr_updated))

points_values = 0:4

tracklist = BioCircosSNPTrack('SNPTrack', points_chromosomes, points_coordinates, 
                              points_values, colors = c("tomato2", "darkblue"), minRadius = 0.5, maxRadius = 0.9)



# Background are always placed below other tracks
tracklist = tracklist + BioCircosBackgroundTrack("myBackgroundTrack", 
                                                 minRadius = 0.5, maxRadius = 0.9,
                                                 borderColors = "#AAAAAA", borderSize = 0.6, fillColors = "#B3E6FF")  

BioCircos(tracklist, genome = CCGenome, genomeFillColor = "PuOr",
          chrPad = 0.05, displayGenomeBorder = FALSE, yChr =  FALSE, genomeLabelOrientation = 90,
          genomeTicksDisplay = FALSE,  genomeLabelTextSize = 10, genomeLabelDy = 0)
