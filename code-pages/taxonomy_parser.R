#so you want to push the available known taxonomic levels as specific as possible?

library(tidyverse)

#example: I have a taxonomic assignment string that looks like this
#Bacteria;Cyanobacteria;Cyanobacteria;SubsectionIII;Unassigned;Pseudanabaena
#or like this:
#Bacteria;Proteobacteria;Deltaproteobacteria;Myxococcales
#where anythign more specific than the order level is absent because it couldn't classify farther

#dplyr::coalesce function: 
#take two inputs X and Y
#if x is NA, then coalesce y into x
#otherwise, keep x


#first, make sure your taxonomy strings are separated into columns
#then replace any "Unassigned" or troublesome assignments

taxonomy_df <- data.frame(Taxonomy = c("Bacteria;Proteobacteria;Deltaproteobacteria;Myxococcales",
                                       "Bacteria;Cyanobacteria;Cyanobacteria;SubsectionIII;Unassigned;Pseudanabaena",
                                       "Bacteria;Proteobacteria;Epsilonproteobacteria",
                                       # "Bacteria;Proteobacteria;Epsilonproteobacteria;Campylobacterales;Helicobacteraceae;Sulfuricurvum",
                                       "Bacteria;Proteobacteria;Gammaproteobacteria;Thiotrichales;Thiotrichaceae;Beggiatoa")) %>%
  tidyr::separate_wider_delim(Taxonomy, delim = ";", 
                              names = c("domain", "phylum", "class", "order", "family", "genus"),
                              too_few = "align_start", cols_remove = FALSE) %>%
  dplyr::relocate(Taxonomy) %>%
  dplyr::mutate(across(c(domain:genus), ~gsub("Unassigned", NA, .x))) #replace any "Unassigned" with NA

new_taxonomy_df <- taxonomy_df %>%
  dplyr::mutate(class = coalesce(class, phylum)) %>%
  dplyr::mutate(order = coalesce(order, class)) %>%
  dplyr::mutate(family = coalesce(family, order)) %>%
  dplyr::mutate(genus = coalesce(genus, family)) %>%
  droplevels
