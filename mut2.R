#######################
# Script: Mutational Signatures LIHC
# Author: Leonardo Daniel González López
# Date: 21/05/2026
#######################

setwd("/mnt/data/bioinfo-estadistica-2/lgonzalez/dani_clase/LIHC2/")
figdir <- "/mnt/data/bioinfo-estadistica-2/lgonzalez/dani_clase/LIHC2/graphs/"

# Load libraries
library(tidyverse) #version 2.0.0
cat("\nPackage tidyverse initiated.\n")
library(ggplot2)
.libPaths("/home/lgonzalez/R/x86_64-pc-linux-gnu-library/4.4")
.libPaths()
library(reticulate) #version 1.46.0
cat("\nPackage reticulate initiated.\n")
library(devtools) #version 2.5.0
cat("\nPackage devtools initiated.\n")
library(SigProfilerMatrixGeneratorR) #version 1.3.6
cat("\nPackage SigProfilerMatrixGeneratorR initiated.\n")
library(SigProfilerPlottingR) #version 1.4.3
cat("\nPackage SigProfilerPlottingR initiated.\n")
library(SigProfilerExtractorR) #version 1.2.6
cat("\nPackage SigProfilerExtractorR initiated.\n")
library(SigProfilerAssignmentR) #version 1.1.3
cat("\nPackage SigProfilerAssignmentR initiated.\n")

# --- Loafing files ---
# Read MAF file
maf_LIHC = read.delim('/mnt/data/bioinfo-estadistica-2/lgonzalez/datos_dani/cancer_genomics/data/LIHC/data_mutations_maf.txt')

cat("\n****************\nLoaded MAF.\n****************\n")

# --- 
# Selection of columns for SigProfiler
maf_SigProfiler = maf_LIHC %>%
  select(Hugo_Symbol, Entrez_Gene_Id, Center, NCBI_Build, Chromosome,
         Start_Position, End_Position, Strand, Variant_Classification,
         Variant_Type, Reference_Allele, Tumor_Seq_Allele1,
         Tumor_Seq_Allele2, dbSNP_RS, dbSNP_Val_Status, Tumor_Sample_Barcode)

# Only mutational signatures of SNPs
maf_SigProfiler = maf_SigProfiler %>%
  filter(Variant_Type == 'SNP')

# Save the filtered MAF file of LIHC 
#write.table(maf_SigProfiler, 'results/SPMG/data_mutations.maf', quote = F, row.names = F, sep = '\t')

cat("\n****************\nMAF of only SNPs saved.\n****************\n")

# --- Create conda environment ---
#This was done in bash before running this code because it was needed to proceed.
#A new conda environment was created in the home directory /home/lgonzalez/.conda/envs/mutational_signatures.
#This new conda environment contains all the SigProfiler packages in python so
#the script can be run.
#
#Code:
#conda config --add envs_dirs /mnt/data/bioinfo-estadistica-2/lgonzalez/conda_envs
#conda create -n mutational_signatures python=3.10 -y
#conda env list     #to check if the environment was created
#
#````
# conda environments:
#
#base                 * /cm/shared/apps/anaconda3/2025.06
#cutadapt191            /cm/shared/apps/anaconda3/2025.06/envs/cutadapt191
#deeptools-2.5.3        /cm/shared/apps/anaconda3/2025.06/envs/deeptools-2.5.3
#deeptools-3.5.6        /cm/shared/apps/anaconda3/2025.06/envs/deeptools-3.5.6
#kneaddata-0.12.4       /cm/shared/apps/anaconda3/2025.06/envs/kneaddata-0.12.4
#macs2-2.1.0            /cm/shared/apps/anaconda3/2025.06/envs/macs2-2.1.0
#multiqc-1.5            /cm/shared/apps/anaconda3/2025.06/envs/multiqc-1.5
#rsat                   /cm/shared/apps/anaconda3/2025.06/envs/rsat
#mutational_signatures   /mnt/data/bioinfo-estadistica-2/lgonzalez/conda_envs/mutational_signatures   <- Yes it was created
#````
#
#conda activate mutational_signatures
#pip install SigProfilerMatrixGenerator
#pip install SigProfilerPlotting
#pip install SigProfilerExtractor
#pip install SigProfilerAssignment
#SigProfilerAssignment      1.1.3
#SigProfilerExtractor       1.2.6
#SigProfilerMatrixGenerator 1.3.6
#sigProfilerPlotting        1.4.3
#
# #checar que se hayan instalado
#pip list | grep -i sigprofiler
#
#````
#SigProfilerAssignment      1.1.3
#igProfilerExtractor       1.2.6
#SigProfilerMatrixGenerator 1.3.6
#sigProfilerPlotting        1.4.3
#````

# --- SigPofiler MatrixGenerator ---

# Use conda environment
use_condaenv('mutational_signatures')

# Install the 37 version of the human genome reference
install('GRCh37', rsync=FALSE, bash=TRUE)

# Mutational profiles
mutational_profiles <- SigProfilerMatrixGeneratorR(project = "LIHC",
                                                   genome = "GRCh37",
                                                   matrix_path = "./results/SPMG",
                                                   plot = F,
                                                   exome = T)

cat("\nMutational profiles generated.\n")

# --- SigProfiler Plotting ---

# -- Average mutational profiles --
#plotSBS(matrix_path = 'results/SPMG/output/SBS/LIHC.SBS96.exome',
#        output_path = 'graphs/',
#        project = 'LIHC',
#        plot_type = '96',
#        percentage = FALSE)

cat("\n****************\nAverage mutational profiles plotted.\n****************\n")

# -- Average mutational profiles --
av_mut_profiles = mutational_profiles[['96']]

# Relative mutational matrix
relative_mut_matrix = apply(av_mut_profiles, 2, prop.table)

# Average mutational matrix
av_mut_matrix = rowMeans(relative_mut_matrix)
av_mut_matrix = data.frame(Average_LIHC = av_mut_matrix)

# Add row names as column and print
av_mut_matrix_print = cbind(rownames(av_mut_matrix),
                            av_mut_matrix)
colnames(av_mut_matrix_print)[1] = 'MutationType'
#write.table(av_mut_matrix_print, 'results/avg_LIHC.all',
#            quote = F, row.names = F, sep = '\t')

# Plot average mutational profiles
#plotSBS(matrix_path = 'results/avg_LIHC.all',
#        output_path = 'graphs/',
#        project = 'avg_LIHC',
#        plot_type = '96',
#        percentage = TRUE)

cat("\n****************\nAverage mutational matrix generated and plotted.\n****************\n")

# -- Average mutational profiles per subgroup --

# Load metadata
metadata = read.delim('/mnt/data/bioinfo-estadistica-2/lgonzalez/datos_dani/cancer_genomics/data/LIHC/lihc_tcga_pan_can_atlas_2018_clinical_data_mtsign.tsv')

# Samples with mutation info
metadata = metadata %>%
  filter(Sample.ID %in% maf_SigProfiler$Tumor_Sample_Barcode)

#save(metadata, file = "results/metadata.Rdata")

# Get samples from group
LIHC_subgroup1 = metadata %>%
  filter(Subtype == 'LIHC') %>% # <- There are no soubgrups found in this dataset, so we will only consider for LIHC as it is the only one found in the Subgroup column besides N.A. 
  pull(Sample.ID)

# Select group samples from main matrix and get average
group1_samples_av = rowMeans(relative_mut_matrix[,LIHC_subgroup1])
group1_samples_av = data.frame(group1_samples_av)

# Add row names as column and print
group1_samples_print = cbind(rownames(group1_samples_av),group1_samples_av)
colnames(group1_samples_print) = c('MutationType', 'LIHC')
#write.table(group1_samples_print,
#            'results/avg_LICH_sbgrp1.SBS96.all',
#            quote = F, row.names = F, sep = '\t')

# Plot average mutational profiles (note the percentage parameter now)
#plotSBS(matrix_path = 'results/avg_LICH_sbgrp1.SBS96.all',
#        output_path = 'graphs/',
#        project = 'avg_LIHC_sbgrp1.SBS96.all',
#        plot_type = '96',
#        percentage=TRUE)

cat("\n****************\nPer goup mutational profiles finished.\n****************\n")

# --- Extracting mutational signatures (SigProfiler Extractor) ---

#sigprofilerextractor(input_type = 'matrix',
#                     output = 'results/SPE/',
#                     input_data = 'results/SPMG/output/SBS/LIHC.SBS96.exome',
#                     nmf_replicates = 100,   # <- for stable results
#                     minimum_signatures = 1,
#                     maximum_signatures = 5, # <- Because we have > 200 samples
#                     exome = T)

cat("\n****************\nExtracting of mutational signatures completed.\n****************\n")

# --- Assessing reference mutational signatures ---

# Assignment analysis
#cosmic_fit(samples = 'results/SPMG/output/SBS/LIHC.SBS96.exome',
#           output = 'results/SPA',
#           input_type='matrix',
#           exome = T)

cat("\n****************\nAssessing reference of mutational signatures completed.\n****************\n")

# --- Downstream analysis of signature assignment results ---

# -- Assessing the accuracy of the signature activities --

stats <- read.delim('results/SPA/Assignment_Solution/Solution_Stats/Assignment_Solution_Samples_Stats.txt')

accuary_sig_act_plot <- ggplot(stats) +
  aes(x=Cosine.Similarity) +
  labs(x='')+
  geom_histogram(aes(y = after_stat(density))) +
  geom_density(col = 4, lwd = 1.5) +
  geom_vline(aes(xintercept = 0.8),
             col = 2, lwd = 1.5) +
  labs(x = 'Cosine Similarity') +
  theme_bw()

#ggsave(filename = paste0(figdir, "Accuary_signature_activities_plot.png"), plot = accuary_sig_act_plot)

cat("\n****************\nAccuary of the signature activities plot completed.\n****************\n")

# -- Visualizing signature activities from SigProfilerExtractor  --

# Activities matrix
act_matrix = read.delim('results/SPE/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities.txt')

# Average activities per signature
avg_act_per_signature = colMeans(act_matrix[,-1])
cat("\n****************\nActivities per signature\n****************\n")
avg_act_per_signature
cat("\n****************\n")

# Save plot because the cluster can not renderize images
#avg_act_barplot <- ggplot(data = data.frame(Signature = names(avg_act_per_signature),
#                                            Avg = avg_act_per_signature)) +
#  aes(x = reorder(Signature, Avg), y = Avg) +
#  geom_bar(stat = 'identity', fill = 'steelblue') +
#  coord_flip() +
#  labs(x = 'Signature', y = 'Average mutations') +
#  theme_bw()

#ggsave(filename = paste0(figdir, "Average_activities_barplot.png"), plot = avg_act_barplot)

cat("\n****************\nAverage activities per signature barplot completed.\n****************\n")

# Reformat dataframe to use ggplot
acts_reformat = act_matrix %>%
  pivot_longer(cols = !Samples,
               names_to = 'Signature',
               values_to = 'Mutations')

# Generate percent stacked barplot
#stacked_barplot <- ggplot(acts_reformat) +
#  aes(x = Samples, y = Mutations, fill = Signature) +
#  geom_bar(position = 'fill', stat = 'identity') +
#  theme_bw() +
#  theme(axis.title.x=element_blank(),
#        axis.text.x=element_blank(),
#        axis.ticks.x=element_blank())

#ggsave(filename = paste0(figdir, "stacked_barplot.png"), plot = stacked_barplot)

cat("\n****************\nAccuary of the signature activities percent stacked plot completed.\n****************\n")

# Calculate number of mutations per sample
num_mutations_per_sample = rowSums(act_matrix[,-1])

# Selecting the activities of only the top 10 mutated cases
top_10_mutated_samples = act_matrix[order(num_mutations_per_sample,
                                          decreasing = T)[1:10],]

# Reformatting and plotting
top_10_mut_samples_plot <- top_10_mutated_samples %>%
  pivot_longer(cols = !Samples,
               names_to = 'Signature',
               values_to = 'Mutations') %>%
  ggplot() +
  aes(x = reorder(Samples, Mutations), y = Mutations, fill = Signature) +
  geom_bar(position = 'fill', stat = 'identity') +
  theme_bw() +
  labs(x = 'Samples')  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

#ggsave(filename = paste0(figdir, "top_10_mut_samples.png"), plot = top_10_mut_samples_plot)

cat("\n****************\nAccuary of the signature activities percent stacked plot of teh top 10 mutated samples completed.\n****************\n")

# -- Associating signature activities with specific metadata --

# Merge activities and metadata tables (The samples column needs to be renamed in one of them)
acts_and_metadata = acts_reformat %>%
  rename(Sample.ID = Samples) %>%
  left_join(metadata)

# Calculate average activities per subtype
acts_per_subgroup = acts_and_metadata %>%
  group_by(Subtype, Signature) %>%
  summarise(Avg_mutations = mean(Mutations))

cat("\n****************\nAverage activities per subtypes\n****************\n")
head(acts_per_subgroup)
cat("\n****************\n")

# Selecting only LIHC subtypes
##### THERE ARE NO SUBTYPES IN THIS DATASET
#acts_per_subgroup = acts_per_subgroup %>%
#  filter(grepl('LIHC', Subtype))

# Plotting stacked barplot per subtype
#stacked_per_subtype_barplot <- ggplot(acts_per_subgroup) +
#  aes(x = reorder(Subtype, Avg_mutations), y = Avg_mutations, fill = Signature) +
#  geom_bar(position = 'fill', stat = 'identity') +
#  theme_bw() +
#  labs(x = 'Liver Hepatocarcinoma Subtype',
#       y = 'Average number of mutations')

#ggsave(filename = paste0(figdir, "stacked_per_subtype_plot.png"), plot = stacked_per_subtype_barplot)

# --- Associating signature activities with specific metadata ---

# -- CTNNB1 --
# From the initial dataset we only pick variants in CTNNB1
maf_LIHC_CTNNB1 = maf_LIHC %>%
  filter(Hugo_Symbol == 'CTNNB1')

# We record into Samples_with_CTNNB1_var the Tumor_Samples_Barcode of samples with variants in the most mutated region of CTNNB1
Samples_with_CTNNB1_var = maf_LIHC_CTNNB1 %>% 
  filter(between(Protein_position, 32, 45)) %>% #Position range obtained from COSMIC Gene View
  distinct(Tumor_Sample_Barcode)

# Now, we need to know what samples of the cohort had Signatures Activities related to CTNNB1 
SBS96_Activities = read.delim('results/SPE/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities.txt')

# Selectes the 4 most common signatures
SBS96_Activities_CTNNB1 = SBS96_Activities %>% 
  select(Samples, SBS5, SBS40a, SBS22a)

# Here we are print a list of samples that have a "variant in the most mutated region of CTNNB1" and it's activities in the signatures selected
SBS96_Activities_CTNNB1 %>%
  filter(Samples %in% Samples_with_CTNNB1_var$Tumor_Sample_Barcode)

# Filter samples with CTNNB1 variants and add total mutations column
SBS96_Activities_CTNNB1_filtered = SBS96_Activities_CTNNB1 %>%
  filter(Samples %in% Samples_with_CTNNB1_var$Tumor_Sample_Barcode) %>%
  mutate(Total = SBS5 + SBS40a + SBS22a)

# now we are plotting the mutational burden (quantity of variants per sample) of the samples that contain a varianat in the most mutated region of CTNNB1 gene
maf_LIHC_var_CTNNB1_sign = maf_LIHC %>% 
  filter(Tumor_Sample_Barcode %in% Samples_with_CTNNB1_var$Tumor_Sample_Barcode)

CTNNB1_sign <- ggplot(data = maf_LIHC_var_CTNNB1_sign) +
  geom_bar(mapping = aes(x = Tumor_Sample_Barcode)) +
  theme(axis.text.x = element_text(angle = 90))

CTNNB1_total_plot <- SBS96_Activities_CTNNB1_filtered %>%
  pivot_longer(cols = c(SBS5, SBS40a, SBS22a),
               names_to = 'Signature',
               values_to = 'Mutations') %>%
  ggplot() +
  aes(x = reorder(Samples, -Total), y = Mutations, fill = Signature) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = c(
    'SBS5'  = 'cornflowerblue', 
    'SBS40a' = 'firebrick1', 
    'SBS22a' = 'plum2'  
  )) +
  labs(x = 'Sample', y = 'Number of mutations', title = 'CTNNB1 - Signature Activities') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

#ggsave(filename = paste0(figdir, "CTNNB1_sign.png"), CTNNB1_sign)
#write.table(x = maf_LIHC_var_CTNNB1_sign, file = "CTNNB1_signature_activity.tsv", row.names = F, sep = '\t')
ggsave(filename = paste0(figdir, "CTNNB1_signature_activity.png"), plot = CTNNB1_total_plot)


# -- NFE2L2 --
# From the initial dataset we only pick variants in NFE2L2
maf_LIHC_NFE2L2 = maf_LIHC %>%
  filter(Hugo_Symbol == 'NFE2L2')

# We record into Samples_with_NFE2L2_var the Tumor_Samples_Barcode of samples with variants in the most mutated region of NFE2L2
Samples_with_NFE2L2_var = maf_LIHC_NFE2L2 %>% 
  filter(between(Protein_position, 77, 82)) %>% #Position range obtained from COSMIC Gene View
  distinct(Tumor_Sample_Barcode)

# Now, we need to know what samples of the cohort had Signatures Activities related to NFE2L2 
SBS96_Activities = read.delim('results/SPE/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities.txt')

# Selectes the 4 most common signatures
SBS96_Activities_NFE2L2 = SBS96_Activities %>% 
  select(Samples, SBS5, SBS40a, SBS22a)

# Here we are print a list of samples that have a "variant in the most mutated region of NFE2L2" and it's activities in the signatures selected
SBS96_Activities_NFE2L2 %>%
  filter(Samples %in% Samples_with_NFE2L2_var$Tumor_Sample_Barcode)

# Filter samples with NFE2L2 variants and add total mutations column
SBS96_Activities_NFE2L2_filtered = SBS96_Activities_NFE2L2 %>%
  filter(Samples %in% Samples_with_NFE2L2_var$Tumor_Sample_Barcode) %>%
  mutate(Total = SBS5 + SBS40a + SBS22a)

# now we are plotting the mutational burden (quantity of variants per sample) of the samples that contain a varianat in the most mutated region of NFE2L2 gene
maf_LIHC_var_NFE2L2_sign = maf_LIHC %>% 
  filter(Tumor_Sample_Barcode %in% Samples_with_NFE2L2_var$Tumor_Sample_Barcode)

NFE2L2_sign <- ggplot(data = maf_LIHC_var_NFE2L2_sign) +
  geom_bar(mapping = aes(x = Tumor_Sample_Barcode)) +
  theme(axis.text.x = element_text(angle = 90))

NFE2L2_total_plot <- SBS96_Activities_NFE2L2_filtered %>%
  pivot_longer(cols = c(SBS5, SBS40a, SBS22a),
               names_to = 'Signature',
               values_to = 'Mutations') %>%
  ggplot() +
  aes(x = reorder(Samples, -Total), y = Mutations, fill = Signature) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = c(
    'SBS5'  = 'cornflowerblue', 
    'SBS40a' = 'firebrick1', 
    'SBS22a' = 'plum2'  
  )) +
  labs(x = 'Sample', y = 'Number of mutations', title = 'NFE2L2 - Signature Activities') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

#ggsave(filename = paste0(figdir, "NFE2L2_sign.png"), NFE2L2_sign)
#write.table(x = maf_LIHC_var_NFE2L2_sign, file = "NFE2L2_signature_activity.tsv", row.names = F, sep = '\t')
ggsave(filename = paste0(figdir, "NFE2L2_signature_activity.png"), plot = NFE2L2_total_plot)

# -- CRIP3 --
# From the initial dataset we only pick variants in CRIP3
maf_LIHC_CRIP3 = maf_LIHC %>%
  filter(Hugo_Symbol == 'CRIP3')

# We record into Samples_with_CRIP3_var the Tumor_Samples_Barcode of samples with variants in the most mutated region of CRIP3
Samples_with_CRIP3_var = maf_LIHC_CRIP3 %>% 
  filter(between(Protein_position, 120, 195)) %>% #Position range obtained from COSMIC Gene View and is the zinc-finger2 domain
  distinct(Tumor_Sample_Barcode)

# Now, we need to know what samples of the cohort had Signatures Activities related to CRIP3 
SBS96_Activities = read.delim('results/SPE/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities.txt')

# Selectes the 4 most common signatures
SBS96_Activities_CRIP3 = SBS96_Activities %>% 
  select(Samples, SBS5, SBS40a, SBS22a)

# Here we are print a list of samples that have a "variant in the most mutated region of CRIP3" and it's activities in the signatures selected
SBS96_Activities_CRIP3 %>%
  filter(Samples %in% Samples_with_CRIP3_var$Tumor_Sample_Barcode)

# Filter samples with CRIP3 variants and add total mutations column
SBS96_Activities_CRIP3_filtered = SBS96_Activities_CRIP3 %>%
  filter(Samples %in% Samples_with_CRIP3_var$Tumor_Sample_Barcode) %>%
  mutate(Total = SBS5 + SBS40a + SBS22a)

# now we are plotting the mutational burden (quantity of variants per sample) of the samples that contain a varianat in the most mutated region of CRIP3 gene
maf_LIHC_var_CRIP3_sign = maf_LIHC %>% 
  filter(Tumor_Sample_Barcode %in% Samples_with_CRIP3_var$Tumor_Sample_Barcode)

CRIP3_sign <- ggplot(data = maf_LIHC_var_CRIP3_sign) +
  geom_bar(mapping = aes(x = Tumor_Sample_Barcode)) +
  theme(axis.text.x = element_text(angle = 90))

CRIP3_total_plot <- SBS96_Activities_CRIP3_filtered %>%
  pivot_longer(cols = c(SBS5, SBS40a, SBS22a),
               names_to = 'Signature',
               values_to = 'Mutations') %>%
  ggplot() +
  aes(x = reorder(Samples, -Total), y = Mutations, fill = Signature) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = c(
    'SBS5'  = 'cornflowerblue', 
    'SBS40a' = 'firebrick1', 
    'SBS22a' = 'plum2'  
  )) +
  labs(x = 'Sample', y = 'Number of mutations', title = 'CRIP3 - Signature Activities') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

#ggsave(filename = paste0(figdir, "CRIP3_sign.png"), CRIP3_sign)
#write.table(x = maf_LIHC_var_CRIP3_sign, file = "CRIP3_signature_activity.tsv", row.names = F, sep = '\t')
ggsave(filename = paste0(figdir, "CRIP3_signature_activity.png"), plot = CRIP3_total_plot)

# -- AXIN1 --
# From the initial dataset we only pick variants in AXIN1
maf_LIHC_AXIN1 = maf_LIHC %>%
  filter(Hugo_Symbol == 'AXIN1')

# We record into Samples_with_AXIN1_var the Tumor_Samples_Barcode of samples with variants in the most mutated region of AXIN1
Samples_with_AXIN1_var = maf_LIHC_AXIN1 %>% 
  filter(between(Protein_position, 88, 211)) %>% #Position range obtained from COSMIC Gene View and is the zinc-finger2 domain
  distinct(Tumor_Sample_Barcode)

# Now, we need to know what samples of the cohort had Signatures Activities related to AXIN1 
SBS96_Activities = read.delim('results/SPE/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities.txt')

# Selectes the 4 most common signatures
SBS96_Activities_AXIN1 = SBS96_Activities %>% 
  select(Samples, SBS5, SBS40a, SBS22a)

# Here we are print a list of samples that have a "variant in the most mutated region of AXIN1" and it's activities in the signatures selected
SBS96_Activities_AXIN1 %>%
  filter(Samples %in% Samples_with_AXIN1_var$Tumor_Sample_Barcode)

# Filter samples with ACIN1 variants and add total mutations column
SBS96_Activities_AXIN1_filtered = SBS96_Activities_AXIN1 %>%
  filter(Samples %in% Samples_with_AXIN1_var$Tumor_Sample_Barcode) %>%
  mutate(Total = SBS5 + SBS40a + SBS22a)

# now we are plotting the mutational burden (quantity of variants per sample) of the samples that contain a varianat in the most mutated region of AXIN1 gene
maf_LIHC_var_AXIN1_sign = maf_LIHC %>% 
  filter(Tumor_Sample_Barcode %in% Samples_with_AXIN1_var$Tumor_Sample_Barcode)

AXIN1_sign <- ggplot(data = maf_LIHC_var_AXIN1_sign) +
  geom_bar(mapping = aes(x = Tumor_Sample_Barcode)) +
  theme(axis.text.x = element_text(angle = 90))

AXIN1_total_plot <- SBS96_Activities_AXIN1_filtered %>%
  pivot_longer(cols = c(SBS5, SBS40a, SBS22a),
               names_to = 'Signature',
               values_to = 'Mutations') %>%
  ggplot() +
  aes(x = reorder(Samples, -Total), y = Mutations, fill = Signature) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = c(
    'SBS5'  = 'cornflowerblue', 
    'SBS40a' = 'firebrick1', 
    'SBS22a' = 'plum2'  
  )) +
  labs(x = 'Sample', y = 'Number of mutations', title = 'AXIN1 - Signature Activities') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

#ggsave(filename = paste0(figdir, "AXIN1_sign.png"), AXIN1_sign)
#write.table(x = maf_LIHC_var_AXIN1_sign, file = "AXIN1_signature_activity.tsv", row.names = F, sep = '\t')
ggsave(filename = paste0(figdir, "AXIN1_signature_activity.png"), plot = AXIN1_total_plot)

# -- TP53 --
# From the initial dataset we only pick variants in TP53
maf_LIHC_TP53 = maf_LIHC %>%
  filter(Hugo_Symbol == 'TP53')

# We record into Samples_with_TP53_var the Tumor_Samples_Barcode of samples with variants in the most mutated region of TP53
Samples_with_TP53_var = maf_LIHC_TP53 %>% 
  distinct(Tumor_Sample_Barcode)

# Now, we need to know what samples of the cohort had Signatures Activities related to TP53 
SBS96_Activities = read.delim('results/SPE/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities.txt')

# Selectes the 4 most common signatures
SBS96_Activities_TP53 = SBS96_Activities %>% 
  select(Samples, SBS5, SBS40a, SBS22a)

# Here we are print a list of samples that have a "variant in the most mutated region of TP53" and it's activities in the signatures selected
SBS96_Activities_TP53 %>%
  filter(Samples %in% Samples_with_TP53_var$Tumor_Sample_Barcode)

# Filter samples with TP53 variants and add total mutations column
SBS96_Activities_TP53_filtered = SBS96_Activities_TP53 %>%
  filter(Samples %in% Samples_with_TP53_var$Tumor_Sample_Barcode) %>%
  mutate(Total = SBS5 + SBS40a + SBS22a)

# now we are plotting the mutational burden (quantity of variants per sample) of the samples that contain a varianat in the most mutated region of TP53 gene
maf_LIHC_var_TP53_sign = maf_LIHC %>% 
  filter(Tumor_Sample_Barcode %in% Samples_with_TP53_var$Tumor_Sample_Barcode)

TP53_sign <- ggplot(data = maf_LIHC_var_TP53_sign) +
  geom_bar(mapping = aes(x = Tumor_Sample_Barcode)) +
  theme(axis.text.x = element_text(angle = 90))

TP53_total_plot <- SBS96_Activities_TP53_filtered %>%
  pivot_longer(cols = c(SBS5, SBS40a, SBS22a),
               names_to = 'Signature',
               values_to = 'Mutations') %>%
  ggplot() +
  aes(x = reorder(Samples, -Total), y = Mutations, fill = Signature) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = c(
    'SBS5'  = 'cornflowerblue', 
    'SBS40a' = 'firebrick1', 
    'SBS22a' = 'plum2'  
  )) +
  labs(x = 'Sample', y = 'Number of mutations', title = 'TP53 - Signature Activities') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

#ggsave(filename = paste0(figdir, "TP53_sign.png"), TP53_sign)
#write.table(x = maf_LIHC_var_TP53_sign, file = "TP53_signature_activity.tsv", row.names = F, sep = '\t')
ggsave(filename = paste0(figdir, "TP53_signature_activity.png"), plot = TP53_total_plot)


