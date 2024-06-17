library(dplyr)
library(ggplot2)
library(tidyr)
meta <- read.table(file = "/data/sample_meta_data.txt",
                   sep = "\t", header = TRUE)

dir <- "/data/dapars_output/"

files <- list.files(path = dir, pattern = "\\.txt$", full.names = TRUE)

dfs <- list()

# Loop through each file, read it using read.table, and store it in the list
for (file in files) {
  df <- read.table(file, header = TRUE)  # Adjust header argument as per your file structure
  dfs <- c(dfs, list(df))
}

# Combine all data frames into a single data frame
combined_df <- do.call(rbind, dfs)
rm(dir, file, files, dfs, df)
df <- combined_df
new_id <- meta$new_id
colnames(combined_df)[5:71] <- new_id

rf <- df %>% rowwise(Gene, fit_value, Predicted_Proximal_APA, Loci)
rf[is.na(rf)] <- 0

tmp <- rf %>% mutate(mean = mean(c_across(Synthego_G1_1_30_S2_PDUI:Synthego_Wildtype_Peroxide_Treated_3hour_4_S12_PDUI)),
                     sd = sd(c_across(Synthego_G1_1_30_S2_PDUI:Synthego_Wildtype_Peroxide_Treated_3hour_4_S12_PDUI)))
colnames(tmp)[5:71] <- new_id
colnames(rf)[5:71] <- new_id

pdui_long <- pivot_longer(rf, cols = -c("Gene", "fit_value",
                                        "Predicted_Proximal_APA", "Loci"),
                          names_to = "sample", values_to = "PDUI")

groups <- data.frame(groups = unique(pdui_long$sample))
group_stats <- pdui_long %>%
  group_by(Gene, fit_value, Predicted_Proximal_APA, Loci, sample) %>%
  summarise(mean = mean(PDUI),
            sd = sd(PDUI),
            median = median(PDUI))

rm(df)

library(dplyr)

# Function to perform Wilcoxon test for each gene and pair of samples
# Apply the function, perform significance testing for each comparison, you will have to change the
# here is for IL6 4 hours vs WT serum starved
{perform_wilcox_tests <- function(df) {
  df %>%
    group_by(Gene, Loci) %>%
    summarize(
      p_value = wilcox.test(PDUI[sample == "WT_IL6_4"], PDUI[sample == "WT_Serum_Starved"])$p.value
    )
}
  ntc_vs_four <- perform_wilcox_tests(pdui_long)
  ntc_vs_four$comparison <- "NTC_vs_IL6_4_hours"
}

# Function to perform Wilcoxon test for each gene and pair of samples
# here is for 24 hours vs WT serum starved
{perform_wilcox_tests <- function(df) {
  df %>%
    group_by(Gene, Loci) %>%
    summarize(
      p_value = wilcox.test(PDUI[sample == "WT_IL6_24"], PDUI[sample == "WT_Serum_Starved"])$p.value
    )
}
  ntc_vs_24 <- perform_wilcox_tests(pdui_long)
  ntc_vs_24$comparison <- "NTC_vs_IL6_24_hours"
}

# Function to perform Wilcoxon test for each gene and pair of samples
# here is comparison for 4 hours vs 24 hours
{perform_wilcox_tests <- function(df) {
  df %>%
    group_by(Gene, Loci) %>%
    summarize(
      p_value = wilcox.test(PDUI[sample == "WT_IL6_4"], PDUI[sample == "WT_IL6_24"])$p.value
    )
}
  four_vs_24 <- perform_wilcox_tests(pdui_long)
  four_vs_24$comparison <- "IL6_4_vs_IL6_24"
}
# four_vs_24$sample <- "WT_IL6_4"
# now we will merge significant results for all comparisons
ntc_vs_24$sample <- "WT_IL6_24"
ntc_vs_four$sample <- "WT_IL6_4"
df <- rbind(ntc_vs_24, ntc_vs_four)

tmp <- pdui_long %>%
  filter(sample %in% c("WT_IL6_4","WT_IL6_24"))
tmp <- tmp %>%
  group_by(Gene, fit_value, Predicted_Proximal_APA, Loci, sample) %>%
  summarise(mean = mean(PDUI),
            sd = sd(PDUI),
            median = median(PDUI))

tmp2 <- merge(tmp, df, by = c("Gene", "Loci", "sample"))
tmp2 <- tmp2 %>%
  filter(p_value < .05)
wt <- group_stats %>%
  filter(sample == "WT_Serum_Starved")

wt <- wt[,c(1:4, 6)]
test <- merge(tmp2, wt, by = c("Gene","fit_value","Predicted_Proximal_APA","Loci"))
test$dPDUI <- test$mean.x - test$mean.y

two_v_four <- merge(tmp, four_vs_24, by = c("Gene", "Loci", "sample"))
two_v_four <- two_v_four %>%
  filter(p_value < .05)
t4 <- group_stats %>%
  filter(sample == "WT_IL6_24")
t4 <- t4[,c(1:4, 6)]
two_v_four <- merge(two_v_four, t4, by = c("Gene", "fit_value", "fit_value",
                                           "Predicted_Proximal_APA", "Loci"))
two_v_four$dPDUI <- two_v_four$mean.x - two_v_four$mean.y

df <- rbind(two_v_four, test)

# now we write out the significant results 
write.table(df, file = "data/dPDUI_Significant_Results.txt", sep = "\t", 
            row.names = F, quote = F)

write.table(group_stats, file = "data/group_stats_reduced.txt", sep = "\t", row.names = F, quote = F)




