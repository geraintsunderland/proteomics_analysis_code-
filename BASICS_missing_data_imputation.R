library("DEP")
library(tidyverse)
library("tidyr")
library("purrr")
library("ggplot2")
library("SummarizedExperiment")


replicates = 28
replicates1 = 19
replicates2 = 9




Prot1 <-read.csv('BASICS_data_sample_revisions.csv', row.names = NULL)
rownames(Prot1) <- paste(Prot1[,2])

#Filter out all proteins not present in any samples
Prot1 <-  Prot1[rowSums(Prot1[,3:ncol(Prot1)])>0,]


head(Prot1)

Prot1_metadata <- read_csv('Revisions_Phenotype.csv')

Prot1_unique_names <- DEP::make_unique(Prot1, "Protein_ID", "Gene_Name", delim = ";")
Prot1_unique_names <- distinct(Prot1_unique_names, Protein_ID, Gene_Name, .keep_all = TRUE)

Clean <- grep("Clean", colnames(Prot1_unique_names))
Infected <- grep("Infected", colnames(Prot1_unique_names))



#Experimental design



experimental_design <- data_frame(
  label = colnames(Prot1)[!colnames(Prot1) %in% c("Protein_ID", "Gene_Name", "Number")],
  condition = c(rep("Infected", replicates2),
                rep("Clean", replicates1)),
  replicate = c((1:9), (1:19)))

Prot1se <- DEP::make_se(Prot1_unique_names, c(Infected, Clean), experimental_design)

#Need to decide regarding filtering specific proteins, all missing proteins or just those with too many missing values.
#Plot a barplot of protein quantification overlap between samples

plot_frequency(Prot1se)
#Shows that vast majority of proteins are represented in all 28 samples (19 Uninfected, 9 Infected)
plot_numbers(Prot1se)
#Barplot of number of identified proteins per sample
plot_coverage(Prot1se)
#Barplot of protein identification overlap between samples


#Filter options -
#1.Not filter any at all
no_filter<-Prot1se
#2. Filter proteins quantified in all replicates of at least one condition
condition_filter <- filter_proteins(Prot1se, "condition", thr = 0)
#3. Filter proteins that have no missing values only
complete_cases<- filter_proteins(Prot1se, "complete")
#4. Filter for proteins quantified on at least a specified proportion of samples e.g. 80%
frac_filtered <- filter_proteins(Prot1se, "fraction", min = 0.8)

single_sample <- filter_proteins(Prot1se, "fraction", min = 0.036)

#Produces a summarized experiment object. To extract data use assays operatorm e.g. for frac_filtered (proteins presebt in 80% samples)
Prot_filtered <- data.frame(assays(frac_filtered))
write.csv(Prot_filtered, 'Prot_filtered.csv')



#Scale and variance stabilise (normalise) data
# Scale and variance stabilize
no_filter <- normalize_vsn(Prot1se)
condition_filter <- normalize_vsn(condition_filter)
complete_cases <- normalize_vsn(complete_cases)
frac_filtered <- normalize_vsn(frac_filtered)

#Plot mean vs sd scatter
meanSdPlot(no_filter)
meanSdPlot(condition_filter)
meanSdPlot(complete_cases)
meanSdPlot(frac_filtered)

#Next need to consider whether to impute the missing values
#Need to do exploratory analysis
#heatmap:
plot_missval(no_filter)

plot_missval(frac_filtered)

#Check whether missing values are biases towards lower intensity proteins can plot density plot with/without missing values.
plot_detect(no_filter)

plot_detect(frac_filtered)
#look for differences in the 2 distributions.

#No imputation
no_imputation <- no_filter

#Use MinProb method - random draws from gaussian distribution centered around minimal value (for MNAR)
#leave fun blank and error will show all available imputation methods
MinProb_imputation <- impute(no_filter, fun = " ")

MinProb_imputation <- impute(frac_filtered, fun = "MinProb", q=0.01)
# Impute missing data using random draws from a
# manually defined left-shifted Gaussian distribution (for MNAR)
manual_imputation <- impute(frac_filtered, fun = "man", shift = 1.8, scale = 0.3)

# Impute missing data using the k-nearest neighbour approach (for MAR)
knn_imputation <- impute(frac_filtered, fun = "knn", rowmax = 0.9)

#Plot densiry plots to show effect of different methods
plot_imputation(no_imputation, frac_filtered, MinProb_imputation, manual_imputation, knn_imputation)

#extract filtered data as df

BASICS_rev_filtered <- data.frame(assays(frac_filtered))
#replace NAs
BASICS_rev_filtered <- BASICS_rev_filtered %>% mutate(across(everything(), ~replace_na(.x, 0)))%>% dplyr::select(-'group', -'group_name')
rownames(BASICS_rev_filtered) = paste(frac_filtered@NAMES)

#Mixed imputation methods
# Extract protein names with missing values in all replicates of at least one condition
proteins_MNAR <- get_df_long(no_filter) %>%
  group_by(name, condition) %>%
  summarize(NAs = all(is.na(intensity))) %>%
  filter(NAs) %>%
  pull(name) %>%
  unique()

# Get a logical vector
MNAR <- names(frac_filtered) %in% proteins_MNAR

# Perform a mixed imputation
mixed_imputation <- impute(
  frac_filtered,
  fun = "mixed",
  randna = !MNAR, # we have to define MAR which is the opposite of MNAR
  mar = "knn", # imputation function for MAR
  mnar = "zero") # imputation function for MNAR

plot_imputation(frac_filtered, mixed_imputation)

#Test for differential expression - to compare performance of imputation methods
# Function that wraps around test_diff, add_rejections and get_results functions
DE_analysis <- function(Prot1se) {
  Prot1se %>%
    test_diff(., type = "control", control = "Clean") %>%
    add_rejections(., alpha = 0.1, lfc = 0) %>%
    get_results()
}

# DE analysis on no, knn, MinProb and mixed imputation
no_imputation_results <- DE_analysis(no_imputation)
knn_imputation_results <- DE_analysis(knn_imputation)
MinProb_imputation_results <- DE_analysis(MinProb_imputation)
mixed_imputation_results <- DE_analysis(mixed_imputation)
# Function to extract number of DE proteins
DE_prots <- function(results) {
  data_frame(Dataset = gsub("_results", "", results),
             significant_proteins = get(results) %>%
               filter(significant) %>%
               nrow())
}

objects <- c("no_imputation_results",
             "knn_imputation_results",
             "MinProb_imputation_results",
             "mixed_imputation_results")

map_df(objects, DE_prots)

#ROC curves used to visualise results of different methods for comparison

# Function to obtain ROC data
get_ROC_df <- function(results) {
  get(results) %>%
    select(name, Infected_vs_Clean_p.val, significant) %>%
    mutate(
      DE = grepl("DE", name),
      BG = grepl("bg", name)) %>%
    arrange(Infected_vs_Clean_p.val) %>%
    mutate(
      TPR = cumsum(as.numeric(DE)) / 300,
      FPR = cumsum(as.numeric(BG)) / 3000,
      method = results)
}

# Get ROC data for no, knn, MinProb and mixed imputation
ROC_df <- map_df(objects, get_ROC_df)

# Plot ROC curves
ggplot(ROC_df, aes(FPR, TPR, col = method)) +
  geom_line() +
  theme_DEP1() +
  ggtitle("ROC-curve")

#Then useful to compare the number of true positives and false positives that the different methods are generating.
#Which proteins are not identified as differentially expressed proteins in the datasets with no or knn imputation? And which proteins are specifically for mixed imputation? We look at both true and false positive hits as well as the missing values.

# Function to obtain summary data
get_rejected_proteins <- function(results) {
  get(results) %>%
    filter(significant) %>%
    left_join(., select(sim, name, MAR, MNAR), by = "name") %>%
    mutate(
      DE = grepl("DE", name),
      BG = grepl("bg", name),
      method = results)
}

# Get summary data for no, knn, MinProb and mixed imputation
objects <- c("no_imputation_results",
             "knn_imputation_results",
             "MinProb_imputation_results",
             "mixed_imputation_results")

summary_df <- map_df(objects, get_rejected_proteins)

# Plot number of DE proteins by technique (True and False)
no_imputation_results %>%
  group_by(method) %>%
  summarize(TP = sum(DE), FP = sum(BG)) %>%
  gather(category, number, -method) %>%
  mutate(method = gsub("_results", "", method)) %>%
  ggplot(aes(method, number, fill = category)) +
  geom_col(position = position_dodge()) +
  theme_DEP2() +
  labs(title = "True and False Hits",
       x = "",
       y = "Number of DE proteins",
       fill = "False or True")
