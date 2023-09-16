#This is an example of Kruskal-Wallis testing between two (or more) groups of observations
#This snippet sets up a matrix for the final results of KW testing
# then prepares the vectors of observations for testing

#In here, df is a dataframe of observations of abundances of taxa
#ASV  SampleName1 SampleName2 SampleName3
#ASV1 5 20  8
#ASV2 18  300 4975
#etc.


start_vector <- sample(seq(1:1000), size = 1000, replace = TRUE) %>% sort()

df <- data.frame(Type_A = sample(start_vector, size = 10000, replace = TRUE, prob = c(rep(1, times = 100), rep(10, times = 500), rep(100, times = 400))),
                 Type_B = sample(start_vector, size = 10000, replace = TRUE, prob = c(rep(1, times = 100), rep(10, times = 500), rep(100, times = 400))),
                 Type_C = sample(start_vector, size = 10000, replace = TRUE, prob = c(rep(1, times = 300), rep(10, times = 500), rep(100, times = 200))),
                 Type_D = sample(start_vector, size = 10000, replace = TRUE, prob = c(rep(1, times = 300), rep(10, times = 500), rep(100, times = 200))),
                 Type_E = sample(start_vector, size = 10000, replace = TRUE, prob = c(rep(1, times = 300), rep(10, times = 500), rep(100, times = 200))),
                 Type_F = sample(start_vector, size = 10000, replace = TRUE, prob = c(rep(1, times = 100), rep(10, times = 500), rep(100, times = 400))),
                 Type_G = sample(start_vector, size = 10000, replace = TRUE, prob = c(rep(1, times = 100), rep(10, times = 500), rep(100, times = 400))),
                 ASV = paste0("ASV", seq(1:10000))) %>%
  # t() %>%
  tibble::column_to_rownames(var = "ASV") %>%
  as.data.frame()

#set up character vectors with the names of your samples in the groups of interest
pigments_high_chl <- c("Type_A", "Type_B", "Type_F", "Type_G")
pigments_low_chl <- c("Type_C", "Type_D", "Type_E")
pigments_high_cyano <- c("Type_A", "Type_B")
pigments_high_chlorobi <- c("Type_C", "Type_D", "Type_E")
pigments_high_diatom <- c("Type_E", "Type_F", "Type_G")

#make the dataframe where you'll put in the p-values from KW-testing of the abundances of ASVs between groups
#nrows = number of taxa/types of observations
#ncols = number of comparisons you want to conduct

kw.test <- matrix(nrow = nrow(df), ncol = 2)
colnames(kw.test) <- c("Between chlorophyll levels", "Between photosynthetic taxa")
rownames(kw.test) <- rownames(df)

for(i in 1:nrow(kw.test)){
  a <- na.omit(as.numeric(df[i, pigments_high_chl]))
  b <- na.omit(as.numeric(df[i, pigments_low_chl])) #if you have more groups to include in this comparison, make other vectors like this one (d, e, etc.)
  ifelse(length(a[!is.na(a)]) > 0 & length(b[!is.na(b)]) > 0, #test to make sure the length of the vector of observations in the groups is nonzero
         kw.test[i,1] <- kruskal.test(list(a[!is.na(a)], b[!is.na(b)]))$p.value,
         kw.test[i,1] <- NA)
  if(!is.na(kw.test[i,1]))
  {if(kw.test[i,1] > 0.05)
  {(kw.test[i,1] <- NA)}
  }
}

for(i in 1:nrow(kw.test)){
  a <- na.omit(as.numeric(df[i, pigments_high_cyano]))
  b <- na.omit(as.numeric(df[i, pigments_high_diatom]))
  c <- na.omit(as.numeric(df[i, pigments_high_chlorobi]))
  ifelse(length(a[!is.na(a)]) > 0 & length(b[!is.na(b)]) > 0 & length(c[!is.na(c)]) > 0,
         kw.test[i,2] <- kruskal.test(list(a[!is.na(a)], b[!is.na(b)], c[!is.na(c)]))$p.value,
         kw.test[i,2] <- NA)
  if(!is.na(kw.test[i,2]))
  {if(kw.test[i,2] > 0.05)
  {(kw.test[i,2] <- NA)}
  }
}
rm(a)
rm(b)
rm(c)
