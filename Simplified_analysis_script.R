BASICS_filtered <- read.csv('BASICS_hithru_filt_comb_no_outliers.csv')
library(tidyverse)
library(ggrepel)

dimnames(BASICS_filtered)[[1]] <- BASICS_filtered[,2]
BASICS_filtered <- BASICS_filtered[,-c(1,2)]
BASICS_filtered <- BASICS_filtered +1

# Use the t-test function to calculate the difference between poor and good outcome.

colnames(BASICS_filtered)

t.test(BASICS_filtered[1,1:8], BASICS_filtered[1,9:26])

y <- t.test(BASICS_filtered[2,1:8], BASICS_filtered[2, 9:26])

dim(BASICS_filtered)

ttestBASICS_filtered <- function(df, grp1, grp2) {
  x <- df[grp1]
  y <- df[grp2]
  x <- as.numeric(x)
  y <- as.numeric(y)
  results <- t.test(x, y)
  results$p.value
}
rawpvalue <- apply(BASICS_filtered, 1, ttestBASICS_filtered, grp1 <- c(1:8), grp2 <- c(9:26))

# Plot a histogram of the p-values.

hist(rawpvalue)


# Log2 the data, calculate the mean for each gene per group. Then calculate the fold change between the groups.

##transform data into log2 base.
BASICS_filtered <- log2(BASICS_filtered)

#calculate the mean of each gene per control group
control <- apply(BASICS_filtered[,9:26], 1, mean)

#calculate the mean of each gene per test group
test <- apply(BASICS_filtered[, 1:8], 1, mean)

#confirming that we have a vector of numbers
class(control)
## [1] "numeric"
#confirming we have a vector of numbers
class(test)
## [1] "numeric"
#because our data is already log2 transformed, we can take the difference between the means.  And this is our log2 Fold Change or log2 ratio <-<- log2(control / test)
foldchange <- test - control
# Plot a histogram of the fold change values.

hist(foldchange, xlab = "log2 Fold Change (Control vs Test)")


# Transform the p-value (-1*log(p-value)) and create a volcano plot using ggplot2.

results <- cbind(foldchange, rawpvalue)
results <- as.data.frame(results)
results$probename = rownames(results)

library(ggplot2)
volcano <- ggplot(data = results, aes(x = foldchange, y = -1*log10(rawpvalue)))
volcano + geom_point()


# add a column of NAs
results$diffexpressed = "NO"
# if logFC > 1 and pvalue < 0.05, set as "UP"
results$diffexpressed[results$foldchange > 1 & results$rawpvalue < 0.5] = "UP"
# if logFC < 1 and pvalue < 0.05, set as "DOWN"
results$diffexpressed[results$foldchange < -1 & results$rawpvalue < 0.5] = "DOWN"
results$label <- NA
results$label[results$diffexpressed != "NO"] <- results$probename[results$diffexpressed != "NO"]

write.csv(results, 'BASICS_mod_filtered_univariate_results.csv')

top <- 20
top_genes = bind_rows(
  results %>%
    filter(diffexpressed == 'UP') %>%
    arrange(desc(foldchange), desc(abs(rawpvalue))) %>%
    head(top),
  results %>%
    filter(diffexpressed == 'DOWN') %>%
    arrange(foldchange, desc(abs(rawpvalue))) %>%
    head(top))

ggplot(data=results, aes(x=foldchange, y=-log10(rawpvalue), col=diffexpressed, label=label))+
  geom_point() +
  scale_color_manual(values = c("#4575B4", "grey60", "#D73027"), name = "Expression Level",
                     labels = c("Down-regulated\n(p\u22640.05)\n", "Unchanged\n(p>0.05)\n", "Up-regulated\n(p\u22640.05)"))+
  geom_label_repel(data = top_genes,
                   mapping = aes(x=foldchange, y= -log10(rawpvalue), label = probename),
                   size = 2, show.legend = "none", max.overlaps = 100)+
  geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.5), col="red")+
  xlab("log2 Fold Change")

results %>% count(diffexpressed)
results %>% count()


