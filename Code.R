#install :  remotes::install_github("fkeck/bioseq")


library(bioseq)
library(tidyverse)



#First steps

x <- dna(Seq_1 = "ACCTAG", Seq_2 = "GGTATATACC", Seq_3 = "AGTC")
is_dna(x)

x




x[c("Seq_3", "Seq_1")]
x[2]
x[c(FALSE, FALSE, TRUE)]




y <- dna("?AcGF")
y




#Operations on sequences

x_dna <- dna("ATGTCACCACAAACAGAGACT")
x_dna

x_rna <- seq_transcribe(x_dna)
x_rna

x_aa <- seq_translate(x_rna)
x_aa




dna_from_rna <- seq_rev_transcribe(x_rna)
dna_from_rna

dna_from_aa <- seq_rev_translate(x_aa)
dna_from_aa

x_dna_comp <- seq_complement(x_dna)
x_dna_comp_rev <- seq_reverse(x_dna_comp)


dna(x_dna, x_dna_comp, x_dna_comp_rev)

#String operations

x <- dna("CTGAAAACTG", "ATGAAAACTG", "CTGCTG")

#Detection and selection

x[seq_detect_pattern(x, "AAAA")]

x[seq_detect_pattern(x, "A{4}")]

# This works
x[seq_detect_pattern(x, dna("AAAA"))]


# This fails because x is a DNA vector and pattern is an amino acid vector
#x[seq_detect_pattern(x, aa("AAAA"))]


# This works because W can be A or T.
x[seq_detect_pattern(x, dna("WAWA"))]


seq_disambiguate_IUPAC(dna("WAWA"))



#Remove and replace


seq_remove_pattern(x, "A{4}")


seq_replace_pattern(x,
                    pattern = dna("AAAA"),
                    replacement = dna("----"))



x <- seq_remove_pattern(x, "A{4}")
seq_replace_position(x, 4, 6,
                     replacement = dna("CCC"))


x
seq_replace_position(x, 1:3, 6,
                     replacement = dna("-", "--", "---"))





#Cleaning and exploring NCBI data with the bioseq package


library(bioseq)
library(tidyverse)



data(fragilaria, package = "bioseq")
fra_data <- read_fasta(fragilaria)

fra_data



seq_nchar(fra_data) %>% range()


fra_data <- tibble(label = names(fra_data), sequence = fra_data)


fra_data <- fra_data %>% 
  mutate(genbank_id = str_extract(label, "([^\\s]+)"),
         taxa = str_extract(label, "(?<= ).*")) %>% 
  select(genbank_id, taxa, sequence)

fra_data <- fra_data %>% 
  mutate(n_base = seq_nchar(sequence))

fra_data


#Cropping sequences

FWD <- dna("AGGTGAAGTAAAAGGTTCWTACTTAAA",
           "AGGTGAAGTTAAAGGTTCWTAYTTAAA",
           "AGGTGAAACTAAAGGTTCWTACTTAAA")

REV <- dna("CAGTWGTWGGTAAATTAGAAGG",
           "CTGTTGTWGGTAAATTAGAAGG")



seq_disambiguate_IUPAC(FWD)


fra_data <- fra_data %>% 
  mutate(barcode = seq_crop_pattern(sequence,
                                    pattern_in = list(FWD),
                                    pattern_out = list(REV)))

fra_data


fra_data <- fra_data %>% 
  filter(seq_nchar(barcode) == 312)

fra_data


#Consensus sequences and phylogeny

fra_consensus <- fra_data %>% 
  group_by(taxa) %>% 
  summarise(consensus_barcode = seq_consensus(barcode))

fra_consensus


fra_consensus %>% 
  as_DNAbin(consensus_barcode, taxa) %>% 
  ape::dist.dna() %>% 
  ape::bionj() %>% 
  plot()


#Clustering sequences


duplicated(fra_consensus$consensus_barcode)

fra_consensus <- 
  fra_consensus %>% 
  mutate(cluster = seq_cluster(consensus_barcode,
                               threshold = 0.001))
fra_consensus


fra_consensus <-
  fra_consensus %>% 
  group_by(cluster) %>% 
  summarise(taxa_group = paste(taxa, collapse = "/"),
            consensus_barcode = seq_consensus(consensus_barcode))

fra_consensus %>% 
  as_DNAbin(consensus_barcode, taxa_group) %>% 
  ape::dist.dna() %>% 
  ape::bionj() %>% 
  plot()


#Exporting data


fra_consensus %>% 
  select(taxa_group, consensus_barcode) %>% 
  deframe() %>% 
  write_fasta("my_sequences.fasta")
  
  
  
  #https://fkeck.github.io/bioseq/
