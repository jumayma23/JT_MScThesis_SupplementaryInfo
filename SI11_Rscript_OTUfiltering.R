# Load necessary package
library(tidyverse)

# Read the OTU table
otu <- read.table("otutab.0.95.tsv", header=TRUE, sep="\t") %>% as_tibble()

# Read the BLAST results
blast <- read.table("all_samples_sorted.cluster_0.95.non-chim.fasta.tufDB_v2.0.80.tab", header=FALSE, sep="\t") %>% as_tibble()
colnames(blast) <- c("OTU.ID", "taxonomy", "identity", "alignment_length", "mismatches", "gap_opens", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

# Read the sequence table
seq <- read.table("cleaned_sequences.tab", header=FALSE, sep="\t") %>% as_tibble()
colnames(seq) <- c("OTU.ID", "seq")

# Create a unique, sorted list of OTU IDs from the OTU table
otu_ids <- unique(otu$OTU.ID)

# Create a mapping from old OTU IDs to new OTU IDs
new_ids <- sprintf("OTU_%05d", seq_along(otu_ids))
otu_id_mapping <- tibble(OTU.ID = otu_ids, New.OTU.ID = new_ids)

# Apply the mapping to the OTU table
otu <- otu %>% left_join(otu_id_mapping, by = "OTU.ID") %>% 
  mutate(OTU.ID = New.OTU.ID) %>% 
  select(-New.OTU.ID)

# Apply the mapping to the BLAST results
blast <- blast %>% left_join(otu_id_mapping, by = "OTU.ID") %>% 
  mutate(OTU.ID = New.OTU.ID) %>% 
  select(-New.OTU.ID)

# Apply the mapping to the sequence table
seq <- seq %>% left_join(otu_id_mapping, by = "OTU.ID") %>% 
  mutate(OTU.ID = New.OTU.ID) %>% 
  select(-New.OTU.ID)

# Join the tables
full_table <- inner_join(blast, seq, by="OTU.ID") %>% inner_join(otu, by="OTU.ID")

# Proceed with the rest of your analysis
full_table <- full_table %>% select(OTU.ID, everything())



# Check length (i.e., how many rows) of the joined table
cat("Number of rows in the full table:", nrow(full_table), "\n")

# 5. Split the taxonomy into separate columns
full_table <- full_table %>%
  mutate(original_taxonomy = taxonomy) %>%
  separate(taxonomy, into = c("ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";|\\|", extra = "merge")

# 6. Remove rows with NA values in any of the taxonomy columns
full_table <- full_table %>% drop_na()

# Check length (i.e., how many rows) after removing NAs
cat("Number of rows after removing NAs:", nrow(full_table), "\n")

# 7. Filter columns with 0 or 1 sequences in the OTU counts
cols_to_remove <- full_table %>% select(starts_with("REV")) %>% colSums() %>% .[. <= 1] %>% names()
full_table <- full_table %>% select(-all_of(cols_to_remove))

# 8. Remove Bacteria from the Domain column
full_table <- full_table %>% filter(!grepl("Bacteria", Domain))

# 9. Check the taxonomy at various levels
cat("Taxonomy summary at the Domain level:\n")
print(table(full_table$Domain))
cat("\nTaxonomy summary at the Phylum level:\n")
print(table(full_table$Phylum))
cat("\nTaxonomy summary at the Class level:\n")
print(table(full_table$Class))
cat("\nTaxonomy summary at the Order level:\n")
print(table(full_table$Order))
cat("\nTaxonomy summary at the Species level:\n")
print(table(full_table$Species))

# 10. Write the final table to a file
#write_tsv(full_table, "full_table.tsv")

# Optionally, you can also save the sequences in a separate FASTA file
#write_lines(paste0(">", full_table$OTU.ID, "\n", full_table$seq), "filtered_sequences.fasta")





# Filter the rows for Phylum Chlorophyta
#chlorophyta_table <- full_table %>% filter(Phylum == "Chlorophyta")

# Check the number of rows to verify
#cat("Number of rows for Chlorophyta:", nrow(chlorophyta_table), "\n")

# Save the Chlorophyta table to a new TSV file
#write_tsv(chlorophyta_table, "chlorophyta_table.tsv")

# Optionally, save the sequences for Chlorophyta in a FASTA file
#write_lines(paste0(">", chlorophyta_table$OTU.ID, "\n", chlorophyta_table$seq), "chlorophyta_sequences.fasta")



# Read the OTU table
otu_table <- read.table("otutab.0.95.tsv", header=TRUE, sep="\t") %>% as_tibble()
otu_table <- otu_table %>% 
  mutate(total_size = rowSums(select(., starts_with("REV"))))

# Move the "total_size" column to the second position for better readability
otu_table <- otu_table %>% select(OTU.ID, total_size, everything())

# Remove low-abundance OTUs (i.e., remove OTUs with total_size less than 10)
# Adjust the threshold as needed
otu_table <- otu_table %>% filter(total_size >= 5)

# Read the BLAST results
blast_results <- read.table("all_samples_sorted.cluster_0.95.non-chim.fasta.tufDB_v2.0.80.tab", header=FALSE, sep="\t") %>% as_tibble()
colnames(blast_results) <- c("OTU.ID", "taxonomy", "identity", "alignment_length", "mismatches", "gap_opens", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

# Read the sequence table
seq_table <- read.table("cleaned_sequences.tab", header=FALSE, sep="\t") %>% as_tibble()
colnames(seq_table) <- c("OTU.ID", "sequence")

# Create a unique, sorted list of OTU IDs from the OTU table
unique_OTU.IDs <- unique(otu_table$OTU.ID)

# Create a mapping from old OTU IDs to new OTU IDs
new_OTU.IDs <- sprintf("OTU_%05d", seq_along(unique_OTU.IDs))
otu_mapping <- tibble(OTU.ID = unique_OTU.IDs, New_OTU.ID = new_OTU.IDs)

# Apply the mapping to the OTU table
otu_table <- otu_table %>% left_join(otu_mapping, by = "OTU.ID") %>% 
  mutate(OTU.ID = New_OTU.ID) %>% 
  select(-New_OTU.ID)

# Apply the mapping to the BLAST results
blast_results <- blast_results %>% left_join(otu_mapping, by = "OTU.ID") %>% 
  mutate(OTU.ID = New_OTU.ID) %>% 
  select(-New_OTU.ID)

# Apply the mapping to the sequence table
seq_table <- seq_table %>% left_join(otu_mapping, by = "OTU.ID") %>% 
  mutate(OTU.ID = New_OTU.ID) %>% 
  select(-New_OTU.ID)




#BLAST RESULTS TABLES FOR CHLOROPHYTA AND RHODOPHYTA (excluding kappa & ulva c.)
# Split the taxonomy into separate columns in the BLAST results
blast_results <- blast_results %>%
  mutate(original_taxonomy = taxonomy) %>%
  separate(taxonomy, into = c("ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";|\\|", extra = "merge")

# Create a table for Chlorophyta excluding Ulva_compressa
chlorophyta_blast_results <- blast_results %>%
  filter(Phylum == "Chlorophyta") %>%
  filter(Species != "Ulva_compressa")

# Remove rows with NA values in any of the taxonomy columns
chlorophyta_blast_results <- chlorophyta_blast_results %>% drop_na()

# Save the Chlorophyta BLAST results to a TSV file
write_tsv(chlorophyta_blast_results, "chlorophyta_blast_results_excluding_ulva_compressa.tsv")

# Optionally, print the number of rows to verify
cat("Number of rows for Chlorophyta (excluding Ulva_compressa):", nrow(chlorophyta_blast_results), "\n")


# Create a table for Rhodophyta excluding Kappaphycus_alvarezii and Chondrus_crispus
rhodophyta_blast_results <- blast_results %>%
  filter(Phylum == "Rhodophyta") %>%
  filter(Species != "Kappaphycus_alvarezii" & Species != "Chondrus_crispus")

# Remove rows with NA values in any of the taxonomy columns
rhodophyta_blast_results <- rhodophyta_blast_results %>% drop_na()


# Save the Rhodophyta BLAST results to a TSV file
write_tsv(rhodophyta_blast_results, "rhodophyta_blast_results_excluding_kappaphycus_chondrus.tsv")

# Optionally, print the number of rows to verify
cat("Number of rows for Rhodophyta (excluding Kappaphycus_alvarezii and Chondrus_crispus):", nrow(rhodophyta_blast_results), "\n")



#CHLORO AND RHODO SELECTED SPECIES< NO ULVA COMPRESSA/K.ALVAREZII
# Filter the OTU table for Chlorophyta OTUs
chlorophyta_otu_table <- otu_table %>%
  filter(OTU.ID %in% chlorophyta_blast_results$OTU.ID)

# Filter the OTU table for Rhodophyta OTUs
rhodophyta_otu_table <- otu_table %>%
  filter(OTU.ID %in% rhodophyta_blast_results$OTU.ID)

# Combine the Chlorophyta and Rhodophyta OTU tables
combined_otu_table <- bind_rows(chlorophyta_otu_table, rhodophyta_otu_table)

# Merge the combined OTU table with the corresponding BLAST results
# This will allow you to keep the taxonomy and sequence information as well
chlorophyta_blast_merged <- chlorophyta_blast_results %>%
  left_join(chlorophyta_otu_table, by = "OTU.ID")

rhodophyta_blast_merged <- rhodophyta_blast_results %>%
  left_join(rhodophyta_otu_table, by = "OTU.ID")

# Select the columns of interest
chlorophyta_small_table <- chlorophyta_blast_merged %>%
  select(OTU.ID, total_size, starts_with("REV"), Domain, Phylum, Class, Order, Family, Genus, Species)

rhodophyta_small_table <- rhodophyta_blast_merged %>%
  select(OTU.ID, total_size, starts_with("REV"), Domain, Phylum, Class, Order, Family, Genus, Species)

# Combine both small tables into one
combined_small_table <- bind_rows(chlorophyta_small_table, rhodophyta_small_table)

combined_table <-bind_rows(chlorophyta_blast_merged, rhodophyta_blast_merged)

# Save the final combined small table to a TSV file
write_tsv(combined_small_table, "combined_chlorophyta_rhodophyta_small_table.tsv")

# Save the final combined small table to a TSV file
write_tsv(combined_table, "combined_chlorophyta_rhodophytatable.tsv")



#WITH HOST SPECIES
# Create a table for Rhodophyta including Kappaphycus_alvarezii and Chondrus_crispus
rhodophyta_blast_results2 <- blast_results %>%
  filter(Phylum == "Rhodophyta") 

# Remove rows with NA values in any of the taxonomy columns
rhodophyta_blast_results2 <- rhodophyta_blast_results2 %>% drop_na()

# Save the Rhodophyta BLAST results to a TSV file
write_tsv(rhodophyta_blast_results2, "rhodophyta_blast_results.tsv")


# Filter the OTU table for Rhodophyta OTUs
rhodophyta_otu_table <- otu_table %>%
  filter(OTU.ID %in% rhodophyta_blast_results2$OTU.ID)

# Merge the Rhodophyta OTU table with the corresponding BLAST results
rhodophyta_blast_merged <- rhodophyta_blast_results2 %>%
  left_join(rhodophyta_otu_table, by = "OTU.ID")

# Select the columns of interest
rhodophyta_small_table <- rhodophyta_blast_merged %>%
  select(OTU.ID, total_size, starts_with("REV"), Domain, Phylum, Class, Order, Family, Genus, Species)

# Save the final combined small table to a TSV file
write_tsv(rhodophyta_small_table, "rhodophyta_combined_small_table.tsv")


# Filter the OTU table for Rhodophyta OTUs
rhodophyta_otu_table <- otu_table %>%
  filter(OTU.ID %in% rhodophyta_blast_results2$OTU.ID)

# Merge the Rhodophyta OTU table with the corresponding BLAST results
rhodophyta_blast_merged <- rhodophyta_blast_results2 %>%
  left_join(rhodophyta_otu_table, by = "OTU.ID")

# Save the final combined table with all info to a TSV file
write_tsv(rhodophyta_blast_merged, "rhodophyta_combined_all_info.tsv")




#WITH ULVA COMPRESSA

# Create a table for Chlorophyta including Ulva_compressa
chlorophyta_blast_results2 <- blast_results %>%
  filter(Phylum == "Chlorophyta") 

# Remove rows with NA values in any of the taxonomy columns
chlorophyta_blast_results2 <- chlorophyta_blast_results2 %>% drop_na()

# Save the Chlorophyta BLAST results to a TSV file
#write_tsv(chlorophyta_blast_results2, "chlorophyta_blast_results.tsv")

# Filter the OTU table for Chlorophyta OTUs
chlorophyta_otu_table <- otu_table %>%
  filter(OTU.ID %in% chlorophyta_blast_results2$OTU.ID)

# Merge the Chlorophyta OTU table with the corresponding BLAST results
chlorophyta_blast_merged <- chlorophyta_blast_results2 %>%
  left_join(chlorophyta_otu_table, by = "OTU.ID")

# Save the final combined table with all info to a TSV file
write_tsv(chlorophyta_blast_merged, "chlorophyta_combined_all_info.tsv")






# Join the tables
merged_table <- inner_join(blast_results, seq_table, by="OTU.ID") %>% inner_join(otu_table, by="OTU.ID")

# Proceed with the rest of your analysis
merged_table <- merged_table %>% select(OTU.ID, everything())

# Check the number of rows in the joined table
cat("Number of rows in the merged table:", nrow(merged_table), "\n")

# Split the taxonomy into separate columns
merged_table <- merged_table %>%
  mutate(original_taxonomy = taxonomy) %>%
  separate(taxonomy, into = c("ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";|\\|", extra = "merge")

# Remove rows with NA values in any of the taxonomy columns
merged_table <- merged_table %>% drop_na()

# Check the number of rows after removing NAs
cat("Number of rows after removing NAs:", nrow(merged_table), "\n")

# Filter out columns with 5 or fewer sequences in the OTU counts
columns_to_filter <- merged_table %>% select(starts_with("REV")) %>% colSums() %>% .[. <= 5] %>% names()
filtered_table <- merged_table %>% select(-all_of(columns_to_filter))

# Remove Bacteria from the Domain column
filtered_table <- filtered_table %>% filter(!grepl("Bacteria", Domain))

# Check the taxonomy at various levels
cat("Taxonomy summary at the Domain level:\n")
print(table(filtered_table$Domain))
cat("\nTaxonomy summary at the Phylum level:\n")
print(table(filtered_table$Phylum))
cat("\nTaxonomy summary at the Class level:\n")
print(table(filtered_table$Class))
cat("\nTaxonomy summary at the Order level:\n")
print(table(filtered_table$Order))
cat("\nTaxonomy summary at the Species level:\n")
print(table(filtered_table$Species))

# Write the final table to a file
write_tsv(filtered_table, "filtered_table_lt10.tsv")

# Optionally, save the sequences in a separate FASTA file
write_lines(paste0(">", filtered_table$OTU.ID, "\n", filtered_table$sequence), "filtered_sequenceslt10.fasta")

# Filter the rows for Phylum Chlorophyta
chlorophyta_filtered_table <- filtered_table %>% filter(Phylum == "Chlorophyta")

# Check the number of rows to verify
cat("Number of rows for Chlorophyta:", nrow(chlorophyta_filtered_table), "\n")

# Save the Chlorophyta table to a new TSV file
write_tsv(chlorophyta_filtered_table, "chlorophyta_filtered_tablelt10.tsv")

# Optionally, save the sequences for Chlorophyta in a FASTA file
write_lines(paste0(">", chlorophyta_filtered_table$OTU.ID, "\n", chlorophyta_filtered_table$sequence), "chlorophyta_filtered_sequenceslt10.fasta")




# Filter the rows for Phylum Rhodophyta
rhodophyta_filtered_table <- filtered_table %>% filter(Phylum == "Rhodophyta" & Species != "Kappaphycus_striatus" & Species != "Kappaphycus_alvarezii")

# Check the number of rows to verify
cat("Number of rows for Rhodophyta:", nrow(rhodophyta_filtered_table), "\n")

# Save the Rhodophyta table to a new TSV file
write_tsv(rhodophyta_filtered_table, "rhodophyta_filtered_tablelt10.tsv")

# Optionally, save the sequences for Rhodophyta in a FASTA file
write_lines(paste0(">", rhodophyta_filtered_table$OTU.ID, "\n", rhodophyta_filtered_table$sequence), "rhodophyta_filtered_sequenceslt10.fasta")


















# Filter the BLAST results for Bacteria with identity >= 90.0
bacteria_blast_results <- blast_results %>%
  filter(Domain == "Bacteria") %>%
  filter(identity >= 90.0)

# Remove rows with NA values in any of the taxonomy columns
bacteria_blast_results <- bacteria_blast_results %>% drop_na()

# Save the Bacteria BLAST results to a TSV file
write_tsv(bacteria_blast_results, "bacteria_blast_results_filtered.tsv")

# Optionally, print the number of rows to verify
cat("Number of rows for Bacteria with identity >= 90.0:", nrow(bacteria_blast_results), "\n")

# Filter the OTU table for Bacteria OTUs
bacteria_otu_table <- otu_table %>%
  filter(OTU.ID %in% bacteria_blast_results$OTU.ID)

# Merge the Bacteria OTU table with the corresponding BLAST results
bacteria_blast_merged <- bacteria_blast_results %>%
  left_join(bacteria_otu_table, by = "OTU.ID")

# Select the columns of interest
bacteria_small_table <- bacteria_blast_merged %>%
  select(OTU.ID, total_size, starts_with("REV"), Domain, Phylum, Class, Order, Family, Genus, Species)

# Save the final Bacteria small table to a TSV file
write_tsv(bacteria_small_table, "bacteria_small_table_filtered.tsv")

# Save the full merged Bacteria table to a TSV file
write_tsv(bacteria_blast_merged, "bacteria_full_table_filtered.tsv")


