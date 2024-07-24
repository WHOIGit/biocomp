# 00z_heatmap_example.R

#this script plots the expression of MAG bins


# load packages -----------------------------------------------------------

library(tidyverse)
library(data.table)
library(viridis)
library(BiocManager)
library(BiocParallel)
# library(DESeq2)
library(cli)
library(furrr)
library(progressr)
library(formattable)
library(pals)
library(pheatmap)
library(gplots)
library(dendextend)
# library(ggpmisc)
library(ggplotify)

# library(batchtools)
# library(future)
# library(future.apply)

# Plan for resource allocation --------------------------------------------
if(grepl("arch64", Sys.getenv("R_PLATFORM"))){
  print("Using parallelly...")
  nthreads <- parallelly::availableCores() - 1
  future::plan(multisession, workers = nthreads)
  options(future.globals.maxSize = 10e9)
} else {
  print("Using data.table")
  nthreads <- data.table::getDTthreads()
  future::plan(sequential)
  options(future.globals.maxSize = 10e9)
}

# cluster <- multidplyr::new_cluster(n = nthreads)
if(nthreads > 4){
  bpparam_multi <- BiocParallel::MulticoreParam(timeout = 100,
                                                workers = nthreads - 2,
                                                stop.on.error = TRUE,
                                                RNGseed = 48105,
                                                progressbar = TRUE)
} else {
  bpparam_multi <- BiocParallel::MulticoreParam(timeout = 100,
                                                workers = nthreads,
                                                stop.on.error = TRUE,
                                                RNGseed = 48105,
                                                progressbar = TRUE)
}

register(bpparam_multi)
cluster <- multidplyr::new_cluster(n = nthreads)
multidplyr::cluster_library(cluster, "dplyr")


if(grepl("arch64", Sys.getenv("R_PLATFORM")) | grepl("shiny-fiesta", getwd())){
  projectpath <- getwd()
} else {
  if(grepl("user", getwd())){
    projectpath <- "/user/sharon.grim/projects/sabrina/"
  } else {
    projectpath <- "/proj/omics/huber/sgrim/projects/sabrina/" 
  }
}

# Set up custom functions -------------------------------------------------
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
no_progress <- function() {} #your placeholder function for progress reporting
# Import sample metadata --------------------------------------------------

mtg_example_metadata <- tribble(
    ~sample_source, ~sample_type, ~Vent, ~year,
    "FS903", "vent", "marker113", "2013",
    "FS904", "vent", "marker33", "2013",
    "FS906", "vent", "marker113", "2014",
    "LVWS2", "vent", "anemone", "2013",
    "FS907", "vent", "anemone", "2014",
    "FS908", "vent", "marker33", "2014",
    "AnemonePlume", "bg", "plume", "2015",
    "CTDBack", "bg", "CTD", "2015",
    "FS914", "vent", "anemone", "2015",
    "FS915", "vent", "marker113", "2015",
    "FS917", "vent", "marker33", "2015"
  ) %>%
  droplevels %>%
  dplyr::mutate(Vent = factor(Vent, levels = c("marker113", "anemone", "marker33", "plume", "CTD")),
                year = factor(year),
                sample_type = factor(sample_type, levels = c("vent", "bg"))) %>%
  dplyr::arrange(sample_type, Vent, sample_source, year) %>%
  # tibble::column_to_rownames(var = "sample_source")
  droplevels


#make a fake dataframe:

#if contig has "A", it is most abundant in Anemone vent
#if contig has "B" it is most abundant in Marker33
#if contig has "C" it is most abundant in marker113
#if contig has "D" it ias most abundant in CTD
#if contig has "E" it is most abundant in plume

# mtg_example_tbl <- matrix(data = rnorm(1100, mean = 50, sd = 45),
#                           nrow = 1000, 
#                           ncol = length(mtg_example_metadata[["sample_source"]]),
#                           dimnames = list(c(paste0("Coassembly_contig_", LETTERS[1:5], "_", seq(1:1000))),
#                                                  c(mtg_example_metadata[["sample_source"]]))) %>%
#   abs(.) %>%
#   round(., digits = 0)

mtg_example_tbl <- data.frame(contig = grep("_A_", paste0("Coassembly_contig_", LETTERS[1:5], "_", seq(1:1000)), value = TRUE),
                              LVWS2 = c(sample(rnorm(1100, mean = 500, sd = 45), 200)) %>%
                                abs(.) %>%
                                round(., digits = 0),
                              FS907 = c(sample(rnorm(1100, mean = 500, sd = 45), 200)) %>%
                                abs(.) %>%
                                round(., digits = 0),
                              FS914 = c(sample(rnorm(1100, mean = 500, sd = 45), 200)) %>%
                                abs(.) %>%
                                round(., digits = 0)) %>%
  bind_rows(., (data.frame(contig = grep("_C_", paste0("Coassembly_contig_", LETTERS[1:5], "_", seq(1:1000)), value = TRUE),
                           FS903 = c(sample(rnorm(1100, mean = 500, sd = 45), 200)) %>%
                             abs(.) %>%
                             round(., digits = 0),
                           FS906 = c(sample(rnorm(1100, mean = 500, sd = 45), 200)) %>%
                             abs(.) %>%
                             round(., digits = 0),
                           FS915 = c(sample(rnorm(1100, mean = 500, sd = 45), 200)) %>%
                             abs(.) %>%
                             round(., digits = 0))),
            (data.frame(contig = grep("_B_", paste0("Coassembly_contig_", LETTERS[1:5], "_", seq(1:1000)), value = TRUE),
                        FS904 = c(sample(rnorm(1100, mean = 500, sd = 45), 200)) %>%
                          abs(.) %>%
                          round(., digits = 0),
                        FS908 = c(sample(rnorm(1100, mean = 500, sd = 45), 200)) %>%
                          abs(.) %>%
                          round(., digits = 0),
                        FS917 = c(sample(rnorm(1100, mean = 500, sd = 45), 200)) %>%
                          abs(.) %>%
                          round(., digits = 0))),
            (data.frame(contig = grep("_D_", paste0("Coassembly_contig_", LETTERS[1:5], "_", seq(1:1000)), value = TRUE),
                        CTDBack = c(sample(rnorm(1100, mean = 500, sd = 45), 200)) %>%
                          abs(.) %>%
                          round(., digits = 0))),
            (data.frame(contig = grep("_E_", paste0("Coassembly_contig_", LETTERS[1:5], "_", seq(1:1000)), value = TRUE),
                        AnemonePlume = c(sample(rnorm(1100, mean = 500, sd = 45), 200)) %>%
                          abs(.) %>%
                          round(., digits = 0)))) %>%
  tidyr::pivot_longer(., cols = !c("contig"),
                      names_to = "sample_source",
                      values_to = "counts") %>%
  dplyr::rowwise(.) %>%
  dplyr::mutate(counts = dplyr::case_when(is.na(counts) ~ round(abs(rnorm(1, mean = 10, sd = 10)), digits = 0),
                                          .default = counts)) %>%
  droplevels %>%
  dplyr::arrange(desc(counts)) %>%
  tidyr::pivot_wider(., id_cols = "contig",
                     names_from = "sample_source",
                     values_from = "counts") %>%
  tibble::column_to_rownames(var = "contig") %>%
  dplyr::select(sort(mtg_example_metadata[["sample_source"]])) #make it disordered :)

# dend_mtg_example_contig <- data.frame(contig = paste0("Coassembly_contig_", LETTERS[1:5], "_", seq(1:1000))) %>%
#   dplyr::mutate(group = dplyr::case_when(grepl("_A_", contig) ~ "anemone",
#                                          grepl("_B_", contig) ~ "marker33",
#                                          grepl("_C_", contig) ~ "marker113",
#                                          grepl("_D_", contig) ~ "CTD",
#                                          grepl("_E_", contig) ~ "plume")) %>%
#   dplyr::mutate(across(everything(), factor)) %>%
#   tibble::column_to_rownames(var = "contig") %>% 
#   droplevels %>%
#   dplyr::mutate(across(everything(), as.numeric)) %>%
#   # droplevels
#   dist(t(.), method = "euclidean") %>%
#   hclust(method = "ward.D2") %>%
#   as.dendrogram

dend_mtg_example_contig <- mtg_example_tbl %>%
  # t(.) %>%
  dist(., method = "euclidean") %>%
  hclust(method = "ward.D2") %>%
  as.dendrogram

dend_mtg_example_abund <- mtg_example_tbl %>%
  t(.) %>%
  dist(., method = "euclidean") %>%
  hclust(method = "ward.D2") %>%
  as.dendrogram

#maybe a factor like Year is more important for visualizing data:
dend_mtg_example_sample <- mtg_example_metadata %>%
  tibble::column_to_rownames(var = "sample_source") %>%
  dplyr::mutate(across(everything(), factor)) %>%
  dplyr::mutate(across(everything(), as.numeric)) %>%
  dplyr::select(year, sample_type, Vent) %>%
  dplyr::mutate(year = year*1000,
                sample_type = sample_type*100,
                Vent = Vent) %>%
  dist(., method = "euclidean") %>%
  hclust(method = "ward.D2") %>%
  as.dendrogram

mtg_example_sample_colorbar <- mtg_example_metadata %>%
  dplyr::arrange(sample_type, Vent, year) %>%
  tibble::column_to_rownames(var = "sample_source")

mtg_example_contig_colorbar <- data.frame(contig = paste0("Coassembly_contig_", LETTERS[1:5], "_", seq(1:1000))) %>%
  dplyr::mutate(Vent = dplyr::case_when(grepl("_A_", contig) ~ "anemone",
                                         grepl("_B_", contig) ~ "marker33",
                                         grepl("_C_", contig) ~ "marker113",
                                         grepl("_D_", contig) ~ "CTD",
                                         grepl("_E_", contig) ~ "plume")) %>%
  dplyr::mutate(across(everything(), factor)) %>%
  tibble::column_to_rownames(var = "contig")

example_color_names <- c(
  list(Vent = c("navy", "pink", "yellow", "lavender", "cornflowerblue") %>%
         setNames(., c("marker113", "anemone", "marker33", "plume", "CTD")),
       sample_type = c("forestgreen", "purple") %>%
         setNames(c("bg", "vent")),
       year = c("brown", "khaki", "azure") %>%
         setNames(c("2013", "2014", "2015"))))

#without any ordering:

pheatmap(mtg_example_tbl,
         border_color = "black",
         scale = "none",
         na_col = "black",
         cluster_cols = F,
         cluster_rows = F,
         annotation_row = mtg_example_contig_colorbar,
         annotation_col = mtg_example_sample_colorbar,
         annotation_colors = example_color_names,
         angle_col = 90,
         main = paste0("Per-contig coverage across mtg samples"),
         fontsize = 8,
         legend = TRUE,
         drop_levels = TRUE,
         silent = TRUE,
         show_rownames = F,
         show_colnames = T,
) %>%
  as.ggplot(.)


#to order columns, 3 options:
#1. rearrange the matrix so that the columns are in the order you want plotted, and 'cluster_cols = F'
mtg_example_tbl_v1 <- mtg_example_tbl[, rownames(mtg_example_sample_colorbar)]
pheatmap(mtg_example_tbl_v1,
                  border_color = "black",
                  scale = "none",
                  na_col = "black",
                  cluster_cols = F,
                  cluster_rows = F,
                  annotation_row = mtg_example_contig_colorbar,
                  annotation_col = mtg_example_sample_colorbar,
                  annotation_colors = example_color_names,
                  angle_col = 90,
                  main = paste0("Per-contig coverage across mtg samples"),
                  fontsize = 8,
                  legend = TRUE,
                  drop_levels = TRUE,
                  silent = TRUE,
                  show_rownames = F,
                  show_colnames = T,
) %>%
  as.ggplot(.)


#2. order the columns by the dendrogram you made, cluster_cols = F :
mtg_example_tbl_v2 <- mtg_example_tbl[, labels(dend_mtg_example_sample)]

pheatmap(mtg_example_tbl_v2,
         border_color = "black",
         scale = "none",
         na_col = "black",
         cluster_cols = F,
         cluster_rows = F,
         annotation_row = mtg_example_contig_colorbar,
         annotation_col = mtg_example_sample_colorbar,
         annotation_colors = example_color_names,
         # cutree_rows = length(namevar_idx),
         angle_col = 90,
         main = paste0("Per-contig coverage across mtg samples"),
         fontsize = 8,
         # breaks = temp_legend,
         # legend_breaks = temp_legendBreaks,
         # legend_labels = temp_legendLabels,
         # legend = FALSE,
         legend = TRUE,
         # annotation_legend = FALSE,
         drop_levels = TRUE,
         silent = TRUE,
         show_rownames = F,
         show_colnames = T,
) %>%
  as.ggplot(.)


#order the rows now:
mtg_example_tbl_v3 <- mtg_example_tbl[labels(dend_mtg_example_contig), ]

pheatmap(mtg_example_tbl_v3,
         border_color = "black",
         scale = "none",
         na_col = "black",
         cluster_rows = F,
         cluster_cols = F,
         annotation_row = mtg_example_contig_colorbar,
         annotation_col = mtg_example_sample_colorbar,
         annotation_colors = example_color_names,
         angle_col = 90,
         main = paste0("Per-contig coverage across mtg samples"),
         fontsize = 8,
         legend = TRUE,
         drop_levels = TRUE,
         silent = TRUE,
         show_rownames = F,
         show_colnames = T,
) %>%
  as.ggplot(.)

#or use 'cluster_rows = T'
pheatmap(mtg_example_tbl_v3,
         border_color = "black",
         scale = "none",
         na_col = "black",
         cluster_rows = T,
         cluster_cols = F,
         annotation_row = mtg_example_contig_colorbar,
         annotation_col = mtg_example_sample_colorbar,
         annotation_colors = example_color_names,
         angle_col = 90,
         main = paste0("Per-contig coverage across mtg samples"),
         fontsize = 8,
         legend = TRUE,
         drop_levels = TRUE,
         silent = TRUE,
         show_rownames = F,
         show_colnames = T,
) %>%
  as.ggplot(.)



#putting it together, if you have cluster dendrograms for your data:
#use the Year-clustered matrix 
pheatmap(mtg_example_tbl_v2,
         border_color = "black",
         scale = "none",
         na_col = "black",
         cluster_rows = as.hclust(dend_mtg_example_contig),
         cluster_cols = F,
         annotation_row = mtg_example_contig_colorbar,
         annotation_col = mtg_example_sample_colorbar,
         annotation_colors = example_color_names,
         angle_col = 90,
         main = paste0("Per-contig coverage across mtg samples"),
         fontsize = 8,
         legend = TRUE,
         drop_levels = TRUE,
         silent = TRUE,
         show_rownames = F,
         show_colnames = T,
) %>%
  as.ggplot(.)

#use the cluster dendrogram from the sample abundances that are related to sample source
pheatmap(mtg_example_tbl,
         border_color = "black",
         scale = "none",
         na_col = "black",
         cluster_rows = as.hclust(dend_mtg_example_contig),
         cluster_cols = as.hclust(dend_mtg_example_abund),
         annotation_row = mtg_example_contig_colorbar,
         annotation_col = mtg_example_sample_colorbar,
         annotation_colors = example_color_names,
         # cutree_rows = length(namevar_idx),
         angle_col = 90,
         main = paste0("Per-contig coverage across mtg samples"),
         fontsize = 8,
         # breaks = temp_legend,
         # legend_breaks = temp_legendBreaks,
         # legend_labels = temp_legendLabels,
         # legend = FALSE,
         legend = TRUE,
         # annotation_legend = FALSE,
         drop_levels = TRUE,
         silent = TRUE,
         show_rownames = F,
         show_colnames = T,
) %>%
  as.ggplot(.)


#fancy labeling for legend:
# #if you have a log scale:
# temp_breaks <- c(quantile(mtg_example_tbl, probs = seq(0, 1, 0.05), names = TRUE, na.rm = TRUE)) %>% #granular
#   signif(., digits = 1) %>%
#   unique(.) 

#if not:
temp_breaks <- c(seq(from = 0, to = ceiling(max(mtg_example_tbl, na.rm = TRUE)), length = 20)) %>%
  unique(.)

hm_legendBreaks <- c(quantile(temp_breaks, probs = seq(0, 1, 0.2), names = TRUE))
hm_legend <- temp_breaks
# hm_legendLabels <- c(100*(as.numeric(hm_legendBreaks))) %>% #if you rescaled your data but want to present the legend unscaled:
hm_legendLabels <- c((as.numeric(hm_legendBreaks))) %>%
  scales::number(., big.mark = ",", decimal.mark = ".")

pheatmap(mtg_example_tbl,
         color = colorRampPalette(pals::ocean.thermal(n = 3))(length(temp_breaks)),
         border_color = "black",
         scale = "none",
         na_col = "black",
         cluster_rows = as.hclust(dend_mtg_example_contig),
         cluster_cols = as.hclust(dend_mtg_example_abund),
         annotation_row = mtg_example_contig_colorbar,
         annotation_col = mtg_example_sample_colorbar,
         annotation_colors = example_color_names,
         cutree_cols = 4, #fancy!
         cutree_rows = 5, #whoa
         angle_col = 90,
         main = paste0("Per-contig coverage across mtg samples"),
         fontsize = 8,
         breaks = hm_legend,
         legend_breaks = hm_legendBreaks,
         legend_labels = hm_legendLabels,
         legend = TRUE,
         drop_levels = TRUE,
         silent = TRUE,
         show_rownames = F,
         show_colnames = T,
) %>%
  as.ggplot(.)
