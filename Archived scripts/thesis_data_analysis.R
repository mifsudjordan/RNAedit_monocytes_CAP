library(tidyverse)
library(psych)
library(readr)
library(ggplot2)
library(tidyr)
library(DESeq2)
library(apeglm)
library(biomaRt)
library(gplots)
library(gtools)
library(ggrepel)
library(Hmisc)
library(MASS) # load it only when needed
ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", 
                      dataset = "hsapiens_gene_ensembl", 
                      host = "www.ensembl.org")


# Get data
thesis_data_norm <- read_csv("thesis_data_norm.csv")
thesis_data_strict <- read_csv("thesis_data_strict.csv")



recovery_com_var <-read_csv("recovery_com_var.csv")
hosp_com_var <-read_csv("hosp_com_var.csv")
controls_com_var <-read_csv("controls_com_var.csv")
gene_counts <- read_csv("gene_counts.csv")

modified_intersect_atoi_all <- read_delim("modified_intersect_atoi_all.csv", 
                                          delim = "\t", escape_double = FALSE, 
                                          trim_ws = TRUE)

modified_intersect_ctou_all <- read_delim("modified_intersect_ctou_all.csv", 
                                          delim = "\t", escape_double = FALSE, 
                                          trim_ws = TRUE)

atoi<-as.data.frame(modified_intersect_atoi_all)
ctou<-as.data.frame(modified_intersect_ctou_all)
all_variants<-rbind(atoi, ctou)
rm(atoi)
rm(ctou)

write.table(all_variants, file="all_variants.csv", sep = "\t", 
            row.names = FALSE, quote = FALSE, col.names = TRUE)

# Cleaning before final filtration strategy
rm(modified_intersect_atoi_all)
rm(modified_intersect_ctou_all)
rm(all_variants)




# Create a table with the counts of C to U common variants in all groups

hosp_com_var <- as.data.frame(hosp_com_var)
recovery_com_var <- as.data.frame(recovery_com_var)
controls_com_var <- as.data.frame(controls_com_var)
acute_var_count <- hosp_com_var$chr %>% 
  table()
recovery_var_count <-recovery_com_var$chr %>% 
  table()
controls_var_count <- controls_com_var$chr %>% 
  table()

tables_list <- list(
  acute = acute_var_count,
  controls = controls_var_count,
  recovery = recovery_var_count
)
str(tables_list)
hosp_com_var
# Combine tables into a data frame
ctou_var_counts <- do.call(cbind, tables_list)
ctou_var_counts <-as.data.frame(ctou_var_counts)


# Select important variables

tdf <- thesis_datav4 %>% 
  select(-BioSample, -Bases, -Bytes, -Experiment, -GEO_Accession, -create_date)

tdf <- as.data.frame(tdf)
     
#Find the mean variant counts
mean_var_counts <- tdf %>% 
  select(Time_point, CtoU, AtoI_TOTALmatches) %>% 
  group_by(Time_point) %>% 
  summarise(Mean_CtoU = mean(CtoU), Mean_AtoImatched = mean(AtoI_TOTALmatches))

mean_var_counts

variance_var_counts <- tdf %>% 
  select(Time_point, CtoU, AtoI_TOTALmatches) %>% 
  group_by(Time_point) %>% 
  summarise(Var_CtoU = var(CtoU), Var_AtoImatched = var(AtoI_TOTALmatches))

variance_var_counts

# Reshape data to long format
mean_var_counts_long <- pivot_longer(mean_var_counts, cols = c(Mean_CtoU, Mean_AtoImatched),
                                     names_to = "Variable", values_to = "Mean_Value")

# Plot mean C to U and A to I in the same plot
ggplot(mean_var_counts_long, aes(x = Time_point, y = Mean_Value, fill = Variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = round(Mean_Value, 2)),
            position = position_dodge(width = 0.9),
            vjust = -0.1) +
  labs(title = "Mean C to U & A to I variant counts",
       x = "Time Point",
       y = "Mean count of variants") +
  facet_wrap(~ Variable, scales = "free_y", ncol = 1)

#Histograms to visualise distributions
hist(tdf$CtoU)
hist(tdf$AtoI_TOTALmatches)
shapiro.test(tdf$CtoU) # results indicate data is not normally distributed.
shapiro.test(tdf$AtoI_TOTALmatches) # results indicate data is not normally distributed.

CtoU_controls <- tdf %>%
  select(age, Time_point, CtoU) %>% 
  filter(Time_point == "control")

CtoU_acute <- tdf %>%
  select(age, Time_point, CtoU) %>%
  filter(Time_point == "acute")

CtoU_rec <- tdf %>% 
  select(age, Time_point, CtoU) %>%
  filter(Time_point == "recovery")

AtoI_acute <- tdf %>% 
  select(Time_point, AtoI_TOTALmatches) %>%
  filter(Time_point == "acute")

AtoI_recovery <- tdf %>% 
  select(Time_point, AtoI_TOTALmatches) %>%
  filter(Time_point == "recovery")

AtoI_control <- tdf %>% 
  select(Time_point, AtoI_TOTALmatches) %>%
  filter(Time_point == "control")

hist(CtoU_controls$CtoU)
hist(CtoU_acute$CtoU)
hist(CtoU_rec$CtoU)
hist(AtoI_acute$AtoI_TOTALmatches)
hist(AtoI_recovery$AtoI_TOTALmatches)
hist(AtoI_control$AtoI_TOTALmatches)

# Q-Q plot for CtoU
qqnorm(tdf$CtoU)
qqline(tdf$CtoU)

# Q-Q plot for AtoI
qqnorm(tdf$AtoI_TOTALmatches)
qqline(tdf$AtoI_TOTALmatches)


#Negative binomial test results
nb_test_result_CtoU <- glm.nb(CtoU ~ Time_point, data = tdf)
nb_test_result_AtoI <- glm.nb(AtoI_TOTALmatches ~ Time_point, data = tdf)
summary(nb_test_result_CtoU)
summary(nb_test_result_AtoI)


#box plots showing C to U and A to I variants by Time_point
ggplot(tdf, aes(x = Time_point, y = CtoU, color = Time_point)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75), alpha = 0.7) +
  geom_point(position = position_jitter(width = 0.1), size = 1) +
  labs(title = "Boxplot and Scatter Plot of C to U by Time_point",
       x = "Time_point",
       y = "count of C to U variants") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(tdf, aes(x = Time_point, y = AtoI_TOTALmatches, color = Time_point)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75), alpha = 0.7) +
  geom_point(position = position_jitter(width = 0.1), size = 1) +
  labs(title = "Boxplot and Scatter Plot of A to I matched to ADAR db by Time_point",
       x = "Time_point",
       y = "count of A to I variants") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(tdf, aes(x = Time_point, y = AtoI, color = Time_point)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75), alpha = 0.7) +
  geom_point(position = position_jitter(width = 0.1), size = 1) +
  labs(title = "Boxplot and Scatter Plot of A to I by Time_point",
       x = "Time_point",
       y = "count of A to I variants") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# C to U and A to I count differences in controls, paired acute and paired rec
paired_data_controls <- tdf %>% 
  select(Participant_ID, Paired_single, Time_point, CtoU, AtoI_TOTALmatches) %>% 
  filter(Paired_single =="paired") %>% 
  filter(!Participant_ID =="1090") %>% #outlier
  filter(!Participant_ID == "3008") #outlier
paired_data_controls

paired_data_controls$Participant_ID <- as.factor(paired_data_controls$Participant_ID)
paired_data_controls

# Paired bar plots of Acute vs Recovery 
ggplot(paired_data_controls, aes(x = Participant_ID, y = CtoU, fill = Time_point)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  labs(title = "Paired Bar Plot of C to U Counts for Acute and Recovery",
       x = "Participant ID",
       y = "C to U Counts",
       fill = "Time Point") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(paired_data_controls, aes(x = Participant_ID, y = AtoI_TOTALmatches, fill = Time_point)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  labs(title = "Paired Bar Plot of A to I Counts for Acute and Recovery",
       x = "Participant ID",
       y = "A to I Counts",
       fill = "Time Point") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# negative binomial test for paired data acute vs recovery
nb_test_result_CtoU_paired <- glm.nb(CtoU ~ Time_point, data = paired_data_controls)
summary(nb_test_result_CtoU_paired)

nb_test_result_AtoI_paired <- glm.nb(AtoI_TOTALmatches ~ Time_point, data = paired_data_controls)
summary(nb_test_result_AtoI_paired)

paired_data_acute <- paired_data_controls %>% 
  filter(Time_point == "acute")

paired_data_rec <- paired_data_controls %>% 
  filter(Time_point == "recovery")

delta_data <- paired_data_controls %>%
  group_by(Participant_ID) %>%
  summarise(delta_CtoU = diff(CtoU))
delta_data

delta_dataAtoI <- paired_data_controls %>%
  group_by(Participant_ID) %>%
  summarise(delta_AtoI = diff(AtoI_TOTALmatches))

delta_dataAtoI
mean(delta_dataAtoI$delta_AtoI)
plot(delta_dataAtoI$delta_AtoI)

delta_dataAtoI <- delta_dataAtoI %>% 
  mutate(percentage_difference = (delta_dataAtoI$delta_AtoI / paired_data_acute$AtoI_TOTALmatches) * 100)

plot(delta_dataAtoI$percentage_difference)
hist(delta_dataAtoI$delta_AtoI)
hist(delta_dataAtoI$percentage_difference)
mean(delta_dataAtoI$percentage_difference)

delta_data <- delta_data %>% 
  mutate(percentage_difference = (delta_data$delta_CtoU / paired_data_acute$CtoU) * 100)

delta_data

mean(delta_data$delta_CtoU)
mean(delta_data$percentage_difference)

hist(delta_data$delta_CtoU)
hist(delta_dataAtoI$delta_AtoI)
boxplot(delta_data$percentage_difference)

controls_data <- tdf %>% 
  select(Participant_ID, Paired_single, Time_point, CtoU, AtoI_TOTALmatches) %>% 
  filter(Time_point =="control")

paired_data_controls <- rbind(paired_data_controls, controls_data)
paired_data_controls

paired_data_controls$Time_point <- factor(
  paired_data_controls$Time_point,
  levels = c("acute", "recovery", "control")
)

ggplot(paired_data_controls, aes(x = Time_point, y = CtoU)) +
  geom_point() +
  labs(title = "Scatter Plot of C to U by group in paired data",
       x = "Group",
       y = "count of C to U variants")

summary(paired_data_controls$CtoU)
var(paired_data_controls$CtoU)

# Scatter plot of C to U changes vs Age

# Define the age breaks
age_breaks <- seq(20, 100, by = 10)

# Create age groups using cut
age_groups <- cut(tdf$age, breaks = c(-Inf, age_breaks, Inf), 
                  labels = FALSE, include.lowest = TRUE)

age_groups
hist(tdf$age)

# Count the number of individuals in each age group
table(age_groups)
table(tdf$age)

ggplot(tdf, aes(x = age, y = , color = age)) +
  geom_point() +
  labs(title = "Scatter Plot of C to U by Age",
       x = "Age",
       y = "count of C to U variants") +
  scale_color_gradient(low = "lightblue", high = "darkblue")

#boxplot everyone C to U vs age groups

tdf$age_group <- cut(tdf$age, breaks = seq(20, 100, by = 10),
                            labels = FALSE)
age_groups_names <- c("31-40","41-50","51-60","61-70","71-80","81-90", "91-100")

ggplot(tdf, aes(x = as.factor(age_group), y = CtoU, fill = as.factor(age_group))) +
  geom_boxplot() +
  labs(title = "Boxplot of C to U by Age Group 
       in all samples",
       x = "Age Group",
       y = "Count of C to U variants") +
  scale_x_discrete(labels = age_groups_names) +  # Set custom labels for x-axis
  scale_fill_discrete(labels = age_groups_names,name = "Age ranges")


ggplot(tdf, aes(x = as.factor(age_group), y = AtoI_TOTALmatches, fill = as.factor(age_group))) +
  geom_boxplot() +
  labs(title = "Boxplot of A to I by Age Group 
       in all samples",
       x = "Age Group",
       y = "Count of A to I variants") +
  scale_x_discrete(labels = age_groups_names) +  # Set custom labels for x-axis
  scale_fill_discrete(labels = age_groups_names,name = "Age ranges")


ggplot(AtoI_acute, aes(x = as.factor(age_group), y = AtoI_TOTALmatches, fill = as.factor(age_group))) +
  geom_boxplot() +
  labs(title = "Boxplot of A to I by Age Group 
       in Acute sample",
       x = "Age Group",
       y = "Count of A to I variants") +
  scale_x_discrete(labels = age_groups_names) +  # Set custom labels for x-axis
  scale_fill_discrete(labels = age_groups_names,name = "Age ranges")

summary(tdf$AtoI_TOTALmatches)
summary(tdf$CtoU)

ggplot(CtoU_acute, aes(x = age, y = CtoU, color = age)) +
  geom_point() +
  labs(title = "Scatter Plot of C to U in Acute by Age",
       x = "Age",
       y = "count of C to U variants") +
  scale_color_gradient(low = "lightblue", high = "darkblue")

# grouping individuals by age group

CtoU_acute$age_group <- cut(CtoU_acute$age, breaks = seq(20, 100, by = 10),
                            labels = FALSE)

CtoU_controls$age_group <- cut(CtoU_controls$age, breaks = seq(20, 100, by = 10),
                            labels = FALSE)

CtoU_rec$age_group <- cut(CtoU_rec$age, breaks = seq(20, 100, by = 10),
                            labels = FALSE)

# Plotting scatter and box plots

ggplot(CtoU_acute, aes(x = as.factor(age_group), y = CtoU, fill = as.factor(age_group))) +
  geom_boxplot() +
  labs(title = "Boxplot of C to U by Age Group 
       in Acute samples",
       x = "Age Group",
       y = "Count of C to U variants")

ggplot(CtoU_controls, aes(x = age, y = CtoU, color = age)) +
  geom_point() +
  labs(title = "Scatter Plot of CtoU at Controls by Age",
       x = "Age",
       y = "count of C to U variants") +
  scale_color_gradient(low = "lightblue", high = "darkblue")

ggplot(CtoU_controls, aes(x = as.factor(age_group), y = CtoU, fill = as.factor(age_group))) +
  geom_boxplot() +
  labs(title = "Boxplot of C to U by Age Group 
       in Control samples",
       x = "Age Group",
       y = "Count of C to U variants")

ggplot(CtoU_rec, aes(x = age, y = CtoU, color = age)) +
  geom_point() +
  labs(title = "Scatter Plot of CtoU at Recovery by Age",
       x = "Age",
       y = "count of C to U variants") +
  scale_color_gradient(low = "lightblue", high = "darkblue")

ggplot(CtoU_rec, aes(x = as.factor(age_group), y = CtoU, fill = as.factor(age_group))) +
  geom_boxplot() +
  labs(title = "Boxplot of C to U by Age Group 
       in Recovery samples",
       x = "Age Group",
       y = "Count of C to U variants")

# trying to combine them into a single plot:

# Combine the three data frames into a single data frame with a new 'Group' column
CtoU_acute$Group <- "Acute"
CtoU_controls$Group <- "Control"
CtoU_rec$Group <- "Recovery"

combined_data_ctou <- rbind(CtoU_acute, CtoU_controls, CtoU_rec)

# Plot the combined data using facet_wrap
ggplot(combined_data_ctou, aes(x = as.factor(age_group), y = CtoU, fill = as.factor(age_group))) +
  geom_boxplot() +
  facet_wrap(~ Group, scales = "free_y") +
  labs(title = "Boxplot of C to U by Age Group",
       x = "Age Group",
       y = "Count of C to U variants") + 
  scale_x_discrete(labels = age_groups_names) +  # Set custom labels for x-axis
  scale_fill_discrete(labels = age_groups_names,name = "Age ranges")

ggplot(CtoU_acute, aes(x = age, y = CtoU, color = age)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear trend line
  labs(title = "Scatter Plot of C to U by Age in Acute Group",
       x = "Age",
       y = "Count of C to U variants") +
  scale_color_gradient(low = "lightblue", high = "darkblue")

lm_ctoU_acute <- lm(CtoU ~ age, data = CtoU_acute)
summary(lm_ctoU_acute)

CtoU_acute$age_group <- as.factor(CtoU_acute$age_group)

# Kruskal-Wallis test
kruskal.test(CtoU ~ age_group, data = CtoU_acute)


# Sex

summary(as.factor(tdf$sex))

ggplot(tdf, aes(x = as.factor(sex), y = CtoU)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = as.factor(Time_point)), position = position_jitter(width = 0.1), size = 2) +
  labs(title = "Boxplot of C to U by sex",
       x = "Sex",
       y = "Count of C to U variants")

ggplot(tdf, aes(x = as.factor(sex), y = AtoI_TOTALmatches)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = as.factor(Time_point)), position = position_jitter(width = 0.1), size = 2) +
  labs(title = "Boxplot of A to I by sex",
       x = "Sex",
       y = "Count of A to I variants")

#Negative binomial test results
nb_test_result_CtoU_sex <- glm.nb(CtoU ~ sex, data = tdf)
nb_test_result_AtoI_sex <- glm.nb(AtoI_TOTALmatches ~ sex, data = tdf)
summary(nb_test_result_CtoU_sex)
summary(nb_test_result_AtoI_sex)

# Mann-whitney U test
wilcox.test(CtoU ~ sex, data = tdf)
wilcox.test(AtoI_TOTALmatches ~ sex, data = tdf)

# COPD
summary(as.factor(tdf$copd))
table(tdf$Time_point, tdf$copd)

CtoU_acute_copd <- tdf %>%
  filter(Time_point == "acute")

ggplot(CtoU_acute_copd, aes(x = as.factor(copd), y = CtoU)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1), size = 2) +
  labs(title = "Boxplot of C to U by COPD in acute samples",
       x = "COPD",
       y = "Count of C to U variants")

ggplot(CtoU_acute_copd, aes(x = as.factor(copd), y = AtoI_TOTALmatches)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1), size = 2) +
  labs(title = "Boxplot of A to I by COPD in acute samples",
       x = "COPD",
       y = "Count of A to I variants")

# Chromosomal analysis of variant counts of chromosome 16
ggplot(tdf, aes(x = as.factor(Time_point), y = AtoI_Ch16matches)) +
  geom_boxplot() +
  labs(title = "Boxplot of A to I counts in Chromosome 16",
       x = "Group",
       y = "Count of A to I variants")

# Negative binomial tests and kruskal tests
nb_test_result_chr16 <- glm.nb(AtoI_Ch16matches ~ Time_point, data = tdf)
summary(nb_test_result_chr16)

library(broom)

# Columns to consider
ctou_chr_col <- grep("^CtoU_", names(thesis_datav4), value = TRUE)

# Initialize an empty data frame to store the results
ctou_kw_results <- data.frame()

# Loop over each chromosome column
for (col in ctou_chr_col) {
  # Create the formula for the test
  formula <- as.formula(paste(col, "~ Time_point"))
  
  # Perform the test
  test_result <- kruskal.test(formula, data = thesis_datav4)
  
  # Extract the p-value
  p_value <- tidy(test_result)$p.value
  
  # Add the result to the data frame
  ctou_kw_results <- rbind(ctou_kw_results, data.frame(Column = col, P_Value = p_value))
}

# Adjust the p-values using the Bonferroni method
ctou_kw_results$Adjusted_P_Value <- p.adjust(ctou_kw_results$P_Value, method = "BH")

# Print the results
print(ctou_kw_results)

chr_atoi_long <- thesis_datav4 %>%
  select(Time_point, starts_with("AtoI_Ch")) %>% 
  pivot_longer(cols = starts_with("AtoI_Ch"), names_to = "Chromosome", values_to = "Count")

chr_ctou_long <- thesis_datav4 %>%
  select(Time_point, starts_with("CtoU_")) %>% 
  pivot_longer(cols = starts_with("CtoU_"), names_to = "Chromosome", values_to = "Count")

#post-hoc analysis
library(FSA)
dunnTest(CtoU_chr18 ~ Time_point, data=thesis_datav4, method="bh")
dunnTest(CtoU_chrX ~ Time_point, data=thesis_datav4, method="bh")
dunnTest(CtoU_chr19 ~ Time_point, data=thesis_datav4, method="bh")
dunnTest(CtoU_chr14 ~ Time_point, data=thesis_datav4, method="bh")
dunnTest(CtoU_chr9 ~ Time_point, data=thesis_datav4, method="bh")
dunnTest(CtoU_chr4 ~ Time_point, data=thesis_datav4, method="bh")
dunnTest(CtoU_chr8 ~ Time_point, data=thesis_datav4, method="bh")



# Finding the mean count of A to I or C to U variants on each chromosome
# Ran the following block for each group
mean_ctou_recovery <- thesis_datav4 %>% 
  filter(Time_point == "recovery") %>% 
  select(starts_with("CtoU_"))

mean <-colMeans(mean_ctou_recovery)

mean_ctou_acute <- as.data.frame(mean)
mean_ctou_controls <-as.data.frame(mean)
mean_ctou_recovery <-as.data.frame(mean)

mean_ctou_acute$source <-"mean_ctou_acute"
mean_ctou_controls$source <-"mean_ctou_controls"
mean_ctou_recovery$source <-"mean_ctou_recovery"

mean_ctou <- bind_rows(mean_ctou_acute, mean_ctou_recovery, mean_ctou_controls)


mean_ctou$chromosome <- rownames(mean_ctou)
str(mean_ctou)

mean_ctou <- data.frame(
  chromosome = rep(c(1:22, "X", "Y", "MT"), times = 75),
  group = rep(c("acute", "recovery", "control"), each = 25),
  mean = mean_ctou$mean)

mean_ctou <- mean_ctou[1:75,]

mean_ctou

# Plot the bar graph
ggplot(mean_ctou,
       aes(x = factor(chromosome, levels = mixedsort(unique(chromosome))),
           y = mean, fill = factor(group))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Mean count of C to U Variants by Chromosome and Group",
       x = "Chromosome",
       y = "Mean count", 
       fill = "Group") +
  theme(axis.text.x = element_text(hjust = 1))
mean_ctou

# Create a box plot with facets
ggplot(chr_ctou_long, aes(x = as.factor(Time_point), y = Count, fill = as.factor(Time_point))) +
  geom_boxplot() +
  facet_wrap(~Chromosome, scales = "free_y", ncol = 5) +
  labs(title = "Boxplot of C to U counts by Chromosome and Time_point",
       x = "Time_point",
       y = "Count of C to U variants") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(chr_atoi_long, aes(x = as.factor(Time_point), y = Count, fill = as.factor(Time_point))) +
  geom_boxplot() +
  facet_wrap(~Chromosome, scales = "free_y", ncol = 5) +
  labs(title = "Boxplot of A to I counts by Chromosome and Time_point",
       x = "Time_point",
       y = "Count of A to I variants") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Negative binomial tests and kruskal tests
nb_test_result_chrMT <- glm.nb(AtoI_ChMTmatches ~ Time_point, data = tdf)
summary(nb_test_result_chrMT)

kruskal.test(AtoI_ChMTmatches ~ Time_point, data = tdf)

plot(tdf$AtoI_ChMTmatches)

# Common C to U variant count on each chromosome in the groups:
# The following block plots the diveristy of C to U variants on
# each chromosome per group
ctou_var_counts <- ctou_var_counts %>%
  rownames_to_column(var = "Chromosome")

# Reshape the data for ggplot
ctou_var_counts_long <- ctou_var_counts %>%
  pivot_longer(cols = -Chromosome, names_to = "Group", values_to = "Count") %>%
  arrange(factor(Chromosome, levels = mixedsort(unique(ctou_var_counts$Chromosome)))) %>% 
  arrange(factor(Group, levels = c("acute", "recovery", "controls")))

# Plot the bar graph
ggplot(ctou_var_counts_long,
       aes(x = factor(Chromosome, levels = mixedsort(unique(ctou_var_counts$Chromosome))),
           y = Count, fill = factor(Group, levels = c("acute", "recovery", "controls")))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Diversity of C to U Variants on Chromosomes",
       x = "Chromosome",
       y = "Count", 
       fill = "Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ctou_var_counts <- ctou_var_counts %>% 
  mutate(difacutecontrol = ((acute - controls) / acute) * 100) %>% 
  mutate(difacuterec = ((acute - recovery) / acute) * 100) %>% 
  mutate(difacuterecabs = acute - recovery) %>% 
  mutate(difacutecontrolabs = acute - controls)
ctou_var_counts

max(ctou_var_counts$difacutecontrol)
max(ctou_var_counts$difacuterec)

table(acute_var_unique$chr)
table(acute_varNotInControls$chr)
table(acute_varNotInRec$chr)


#------------------------#
#Common Variant analysis
hosp_com_var <- read_csv("hosp_com_var.csv")
controls_com_var <- read_csv("controls_com_var.csv")
recovery_com_var <- read_csv("recovery_com_var.csv")

# Renaming a variable and cleaning data from tailing "," from it's records
hosp_com_var <- hosp_com_var %>% 
  rename(found_in = "Found in")

# Remove the last comma from the "found_in" column
hosp_com_var$found_in <- substring(hosp_com_var$found_in, 1, 
                                     nchar(hosp_com_var$found_in) - 1)

controls_com_var <- controls_com_var %>% 
  rename(found_in = "Found in")

controls_com_var$found_in <- substring(controls_com_var$found_in, 1, 
                                   nchar(controls_com_var$found_in) - 1)

recovery_com_var <- recovery_com_var %>% 
  rename(found_in = "Found in")

recovery_com_var$found_in <- substring(recovery_com_var$found_in, 1, 
                                   nchar(recovery_com_var$found_in) - 1)

# Finding variants found only in hospital but not in other groups and show 
# their % occurance.

acute_varNotInControls <- anti_join(hosp_com_var, controls_com_var, by = "unique_var_id")
controls_var_NotInAcute <- anti_join(controls_com_var, hosp_com_var, by = "unique_var_id")
acute_varNotInRec <- anti_join(hosp_com_var, recovery_com_var, by = "unique_var_id")
acute_var_unique <- anti_join(acute_varNotInControls, recovery_com_var, by =  "unique_var_id")

output_file1 <- "acute_varNotInControls.csv"
output_file2 <- "controls_var_NotInAcute.csv"
output_file3 <- "acute_varNotInRec.csv"
output_file4 <- "acute_var_unique.csv"

# Export tables as tab-delimited text files.
write.table(acute_varNotInControls, file = output_file1, sep = "\t", 
            row.names = FALSE, quote = FALSE, col.names = TRUE)

write.table(controls_var_NotInAcute, file = output_file2, sep = "\t", 
            row.names = FALSE, quote = FALSE, col.names = TRUE)

write.table(acute_varNotInRec, file = output_file3, sep = "\t", 
            row.names = FALSE, quote = FALSE, col.names = TRUE)

write.table(acute_var_unique, file = output_file4, sep = "\t", 
            row.names = FALSE, quote = FALSE, col.names = TRUE)

# Split the comma-delimited "found_in" values into separate rows
hosp_com_var <- hosp_com_var %>%
  separate_rows(found_in, sep = ",")

controls_com_var <- controls_com_var %>%
  separate_rows(found_in, sep = ",")

recovery_com_var <- recovery_com_var %>%
  separate_rows(found_in, sep = ",")

# Merging all variant tables:
all_variants <- bind_rows(hosp_com_var, controls_com_var, recovery_com_var)

all_variants <- all_variants %>% 
  select(-Time_point.y)


# Filter paired, hospital admission and not recovery
paired_variants <- all_variants %>% 
  select(unique_var_id, ref, alt, matches, perc_hosp, Time_point.x, 
         Paired_single, Participant_ID) %>% 
  filter(Time_point.x == "Hospital admission" | Time_point.x == "After 1 month recovery", 
         Paired_single == "paired")

output_file <- "paired_variants.txt"

# Export the dataframe as tab-delimited text
write.table(paired_variants, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)


# Perform a left join to get the corresponding data from the main df
hosp_com_var <- hosp_com_var %>%
  left_join(thesis_datav2, by = c("found_in" = "Run"))

controls_com_var <- controls_com_var %>%
  left_join(thesis_datav2, by = c("found_in" = "Run"))

recovery_com_var <- recovery_com_var %>%
  left_join(thesis_datav2, by = c("found_in" = "Run"))




# Creating a df with unique variant ids and their respective % occurance in groups

controls_per <- controls_com_var %>% 
  select(unique_var_id, perc_controls)

hosp_per <- hosp_com_var %>% 
  select(unique_var_id, perc_hosp)

recovery_per <- recovery_com_var %>% 
  select(unique_var_id, perc_rec)

combined_df <- full_join(hosp_per, controls_per, by = "unique_var_id") %>%
  full_join(., recovery_per, by = "unique_var_id")

#Change all NA values in df to "0"
combined_df <- combined_df %>%
  mutate_all(~ifelse(is.na(.), 0, .))

combined_df <- combined_df %>%
  mutate(diff_hosp_controls = abs(perc_hosp - perc_controls),
         diff_hosp_rec = abs(perc_hosp - perc_rec),
         diff_controls_rec = abs(perc_controls - perc_rec))

unique_var_id_with_biggest_diff <- combined_df %>%
  arrange(desc(pmax(diff_hosp_controls, diff_hosp_rec, diff_controls_rec))) %>%
  select(unique_var_id, diff_hosp_controls, diff_hosp_rec, diff_controls_rec)



#---------------Gene expression analysis--------------------#
#Filter-out 0 counts

gene_counts <- gene_counts[which(rowSums(gene_counts[, -1]) > 0), ]


# getting gene names
BiocManager::install("org.Hs.eg.db")
BiocManager::install("apeglm")
library("org.Hs.eg.db")

# convert gene_counts to regular dataframe
gene_counts <- as.data.frame(gene_counts)
str(gene_counts)


# check if the genes are all different
any(duplicated(gene_counts$Gene))

# Exclude the first column (gene IDs)
gene_counts_without_gene_id <- gene_counts[, -1]

str(gene_counts_without_gene_id)

# setting the row names
row.names(gene_counts_without_gene_id) <- gene_counts$Gene

#maybe read: https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#multi-factor-designs

# Creating a list of factors indicating the group of each sample
condition <- factor(c(thesis_datav4$Time_point))
coldata <- data.frame(row.names = colnames(gene_counts_without_gene_id), condition)
coldata

# Construct DESEQDataSet Object
dds2 <- DESeqDataSetFromMatrix(countData = gene_counts_without_gene_id,
                              colData = coldata,
                              design = ~condition)

model.matrix(~0+~condition, data=colData(dds2))

#no. of genes in dds2 object
nrow(dds2)

#check number of genes with expression of less than 10
genes_less_10_reads <- dds2[rowSums(counts(dds2)) < 10, ]
nrow(genes_less_10_reads)

# estimate size factors
dds2_sizefactors <- estimateSizeFactors(dds2)
dds2_sizefactors$sizeFactor

plot(dds2_sizefactors$sizeFactor)

# Running the DEseq function
DEdds2 <- DESeq(dds2)

# Variance stabilizing transformation 
# to stabilize the variance across the mean. 
# In RNA-seq data analysis, the variance of counts tends to increase 
# with the mean.
vsdata2 <- vst(dds2, blind = FALSE)

#PCA plot
plotPCA(vsdata2, intgroup = "condition")

# Dispersion Estimates plot of DEdds2
plotDispEsts(DEdds2)

# results
res2 <- results(DEdds2, contrast = c("condition", "acute", "control"))

write.csv(res2, file="res_acutevscontrol.csv", row.names = TRUE)

res2_rec <- results(DEdds2, contrast = c("condition", "acute", "recovery"))
write.csv(res2_rec, file="res_acutevsrecovery.csv", row.names = TRUE)

res2_recvcon <- results(DEdds2, contrast = c("condition", "recovery", "control"))
write.csv(res2_recvcon, file="res_recvcon.csv", row.names = TRUE)

# significant results
sig_acutevscontrol <-na.omit(res2)
sig_acutevscontrol <- sig_acutevscontrol[sig_acutevscontrol$padj < 0.05,]
sig_acutevsrecovery <-na.omit(res2_rec)
sig_acutevsrecovery <- sig_acutevsrecovery[sig_acutevsrecovery$padj < 0.05,]

# Convert Ensemblb IDs to gene name and place gene_name as first column
# acutevscontrol
sig_acutevscontrol.df <- as.data.frame(sig_acutevscontrol)

sig_acutevscontrol.df <- data.frame(gene_name = mapIds(org.Hs.eg.db, 
                                                       keys = rownames(sig_acutevscontrol.df),
                                                       keytype = "ENSEMBL",
                                                       column = "SYMBOL"),
                                    sig_acutevscontrol.df)

sig_acutevscontrol.df
# acutevsrecovery
sig_acutevsrecovery.df <- as.data.frame(sig_acutevsrecovery)

sig_acutevsrecovery.df <- data.frame(gene_name = mapIds(org.Hs.eg.db, 
                                                       keys = rownames(sig_acutevsrecovery.df),
                                                       keytype = "ENSEMBL",
                                                       column = "SYMBOL"),
                                    sig_acutevsrecovery.df)

# gene names output to file
writeLines(unlist(gene_name), con = "gene_names.txt")

# target gene list
target_genes = c("APOBEC4", "APOBEC3A", "APOBEC3B", "APOBEC3B-AS1", 
                 "APOBEC3C", "APOBEC3D", "APOBEC3F", "APOBEC3G", 
                 "APOBEC3H", "APOBEC2", "ADAR", "ADARB2", "ADARB1")

# Subsetting results of RNA editing genes
subset_res <- sig_acutevscontrol.df[sig_acutevscontrol.df$gene_name %in% target_genes, ]
subset_res_recovery <- sig_acutevsrecovery.df[sig_acutevsrecovery.df$gene_name %in% target_genes, ]

# Bar plots
ggplot(subset_res, aes(x = gene_name, y = log2FoldChange)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Significant Log2fold change of APOBEC and ADAR Genes in Acute vs Control",
       x = "Gene Name",
       y = "log2FoldChange") +
  theme_minimal()

ggplot(subset_res_recovery, aes(x = gene_name, y = log2FoldChange)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Significant Log2fold change of APOBEC and ADAR Genes in Acute vs Recovery",
       x = "Gene Name",
       y = "log2FoldChange") +
  theme_minimal()

# Extracting gene expression level of ENSG00000100298 (APOBEC3H)
# in the three groups 

normalized_counts <- counts(DEdds2, normalized = TRUE)
apobec3h <- normalized_counts["ENSG00000100298", ]


# Transpose data for plotting
apobec3h <- t(t(apobec3h))
apobec3h <- as.data.frame(apobec3h)
apobec3h
colnames(apobec3h)[colnames(apobec3h) == "V1"] <- "Normalised_counts" # rename

# Adding cytokines to apobec3h data frame, starting il-10

apobec3h$Il_10 <- thesis_datav4$LPS_IL_10
apobec3h$ctou <- thesis_datav4$CtoU
apobec3h$atoi <- thesis_datav4$AtoI_TOTALmatches
apobec3h$condition <- thesis_datav4$Time_point
apobec3h


# Plotting
# adding 1 to avoid  taking the log of zero

# Box plots and scatter
ggplot(apobec3h, aes(x = condition, y = log2(Normalised_counts + 1))) +
  geom_boxplot() +  
  geom_jitter(width = 0.15) +          # Add jitter to better show individual points
  labs(title = "Gene Expression of APOBEC3H across Conditions",
       x = "Condition",
       y = "Log2(Normalized Counts + 1)")

apobec3h <-na.omit(apobec3h)

# Box plots showing secretion of Il-10 across groups
ggplot(apobec3h, aes(x = condition, y = Il_10)) +
  geom_boxplot() +
  labs(title = "Il-10 secretion across Groups",
       x = "Group",
       y = "Il-10 section")

ggplot(apobec3h, aes(x = Il_10, fill = condition)) +
  geom_histogram(binwidth = 250, color = "black", alpha = 0.7) +
  geom_bar(position = "stack", binwidth = 250, alpha = 0.5) +
  scale_x_continuous(breaks = seq(0, max(apobec3h$Il_10), by = 500)) +
  labs(title = "Histogram with Stacked Bars of Conditions",
       x = "LPS_IL_10",
       y = "Frequency")


apobec3h$LPS_IL10_groups <- cut(apobec3h$Il_10, 
                                breaks = c(-100, 300, 600, 1000, 2000,  Inf),
                                labels = c("up to 300", "300-600", "600-1000", "1000-2000","over 2000"),
                                include.lowest = TRUE)

ggplot(apobec3h, aes(x = condition, y = log2(Normalised_counts + 1))) +
  geom_boxplot() +  
  geom_jitter(aes(color = as.factor(LPS_IL10_groups)), width = 0.15) +  # Color points by LPS_IL_10
  labs(title = "Gene Expression of APOBEC3H across Conditions",
       x = "Condition",
       y = "Log2(Normalized Counts + 1)")

ggplot(apobec3h, aes(x = as.factor(LPS_IL10_groups), y = log2(Normalised_counts + 1))) +
  geom_boxplot() +  
  labs(title = "Gene Expression of APOBEC3H vs Il-10 secretion",
       x = "delta Il-10 secretion",
       y = "Log2(Normalized Counts + 1)")

# C to U counts vs section of il-10
ggplot(apobec3h, aes(x = as.factor(LPS_IL10_groups), y = ctou)) +
  geom_boxplot() +  
  labs(title = "C to U count vs Il-10 secretion",
       x = "delta Il-10 secretion",
       y = "C to U count")

matching_indices <- match(row.names(apobec3h), tdf$Run)

apobec3h$AtoI_Totalmatch <- tdf$AtoI_TOTALmatches[matching_indices]

# box plots of A to I counts vs Il-10 secretion
ggplot(apobec3h, aes(x = as.factor(LPS_IL10_groups), y = AtoI_Totalmatch)) +
  geom_boxplot() +  
  labs(title = "A to I count vs Il-10 secretion",
       x = "delta Il-10 secretion",
       y = "A to I count")

# box plot of apobec3h expression across groups with some labelled outliers
ggplot(apobec3h, aes(x = condition, y = log2(Normalised_counts + 1))) +
  geom_boxplot() +  
  geom_jitter(width = 0.15) +
  geom_text_repel(data = subset(apobec3h, log2(Normalised_counts + 1) > 5),
                  aes(label = id),  # label points with id
                  box.padding = 0.2, point.padding = 0.2) +  # adjust padding as needed
  labs(title = "Gene Expression of APOBEC3H across Conditions",
       x = "Condition",
       y = "Log2(Normalized Counts + 1)")

# Il-1b
apobec3h <- normalized_counts["ENSG00000100298", ]
# Preparing dataframe for plotting
apobec3h <- t(t(apobec3h))
apobec3h <- as.data.frame(apobec3h)
colnames(apobec3h)[colnames(apobec3h) == "V1"] <- "Normalised_counts"
apobec3h$Il_1b <- thesis_datav4$LPS_IL_1b
apobec3h$ctou <- thesis_datav4$CtoU
apobec3h$atoi <- thesis_datav4$AtoI_TOTALmatches
apobec3h$condition <- thesis_datav4$Time_point

# Box plots showing secretion of Il-1b across groups
plot(apobec3h$Il_1b)
apobec3h <- apobec3h[rownames(apobec3h) != "SRR12926276", , drop = FALSE] # remove outlier
ggplot(apobec3h, aes(x = condition, y = log2(Il_1b+1))) +
  geom_boxplot() +
  labs(title = "Il-1b secretion across Groups",
       x = "Group",
       y = "log 2 of Il-1b section")

apobec3h$Il_1b
summary(apobec3h$Il_1b)
ggplot(apobec3h, aes(x = Il_1b, fill = condition)) +
  geom_histogram(binwidth = 1000, color = "black", alpha = 0.7) +
  scale_x_continuous(breaks = seq(-80512, 130790, by = 40000)) +  # Adjust the range based on your data
  labs(title = "Histogram with Stacked Bars of Conditions",
       x = "Il_1b",
       y = "Frequency") +
  theme_minimal()

apobec3h <-na.omit(apobec3h)

apobec3h$LPS_IL1b_groups <- cut(apobec3h$Il_1b, 
                                breaks = c(-Inf, 5000, 10000, 15000, 20000,  Inf),
                                labels = c("up to 5000", "5000-10000",
                                           "10000-15000", "15000-20000", "over 20000"),
                                include.lowest = TRUE)


ggplot(apobec3h, aes(x = condition, y = log2(Normalised_counts + 1))) +
  geom_boxplot() +  
  geom_jitter(aes(color = as.factor(LPS_IL1b_groups)), width = 0.15) +  # Color points by LPS_IL_10
  labs(title = "Gene Expression of APOBEC3H across Conditions",
       x = "Condition",
       y = "Log2(Normalized Counts + 1)")

ggplot(apobec3h, aes(x = as.factor(LPS_IL1b_groups), y = log2(Normalised_counts + 1))) +
  geom_boxplot() +  
  labs(title = "Gene Expression of APOBEC3H across Il-1b groups",
       x = "delta Il-1b secretion",
       y = "Log2(Normalized Counts + 1)")

# C to U counts vs section of il-1b
ggplot(apobec3h, aes(x = as.factor(LPS_IL1b_groups), y = ctou)) +
  geom_boxplot() +  
  labs(title = "C to U count vs Il-1b secretion",
       x = "delta Il-1b secretion",
       y = "C to U count")


# box plots of A to I counts vs Il-1b secretion
ggplot(apobec3h, aes(x = as.factor(LPS_IL1b_groups), y = ati)) +
  geom_boxplot() +  
  labs(title = "A to I count vs Il-1b secretion",
       x = "delta Il-10 secretion",
       y = "A to I count")

# TNF-alpha
apobec3h <- normalized_counts["ENSG00000100298", ]
# Preparing dataframe for plotting
apobec3h <- t(t(apobec3h))
apobec3h <- as.data.frame(apobec3h)
colnames(apobec3h)[colnames(apobec3h) == "V1"] <- "Normalised_counts"
apobec3h$tnf <- thesis_datav4$LPS_TNF_alpha
apobec3h$ctou <- thesis_datav4$CtoU
apobec3h$atoi <- thesis_datav4$AtoI_TOTALmatches
apobec3h$condition <- thesis_datav4$Time_point

# Box plots showing secretion of tnf-alpha across groups
plot(apobec3h$tnf)

ggplot(apobec3h, aes(x = condition, y = log2(tnf))) +
  geom_boxplot() +
  labs(title = "TNF-alpha secretion across Groups",
       x = "Group",
       y = "TNF secretion")



ggplot(apobec3h, aes(x = tnf, fill = condition)) +
  geom_histogram(binwidth = 1000, color = "black", alpha = 0.7) +
 # scale_x_continuous(breaks = seq(-80512, 130790, by = 40000)) +  # Adjust the range based on your data
  labs(title = "Histogram with Stacked Bars of Conditions",
       x = "TNF-alpha",
       y = "Frequency") +
  theme_minimal()

apobec3h <-na.omit(apobec3h)

apobec3h$tnf_groups <- cut(apobec3h$tnf, 
                                breaks = c(-Inf, 3000, 6000, 9000, 12000, 18000,  Inf),
                                labels = c("up to 3000", "3000-6000",
                                           "6000-9000", "9000-12000", "12000-18000", "over 18000"),
                                include.lowest = TRUE)


ggplot(apobec3h, aes(x = condition, y = log2(Normalised_counts + 1))) +
  geom_boxplot() +  
  geom_jitter(aes(color = as.factor(tnf_groups)), width = 0.15) +  # Color points by LPS_IL_10
  labs(title = "Gene Expression of APOBEC3H across Conditions",
       x = "Condition",
       y = "Log2(Normalized Counts + 1)")

ggplot(apobec3h, aes(x = as.factor(tnf_groups), y = log2(Normalised_counts + 1))) +
  geom_boxplot() +  
  labs(title = "Gene Expression of APOBEC3H across TNF-alpha groups",
       x = "delta TNF-alpha secretion",
       y = "Log2(Normalized Counts + 1)")

# C to U counts vs section of il-1b
ggplot(apobec3h, aes(x = as.factor(tnf_groups), y = ctou)) +
  geom_boxplot() +  
  labs(title = "C to U count vs Il-1b secretion",
       x = "delta Il-1b secretion",
       y = "C to U count")


# box plots of A to I counts vs Il-1b secretion
ggplot(apobec3h, aes(x = as.factor(tnf_groups), y = atoi)) +
  geom_boxplot() +  
  labs(title = "A to I count vs Il-1b secretion",
       x = "delta Il-10 secretion",
       y = "A to I count")

# Il-6

apobec3h <- normalized_counts["ENSG00000100298", ]
# Preparing dataframe for plotting
apobec3h <- t(t(apobec3h))
apobec3h <- as.data.frame(apobec3h)
colnames(apobec3h)[colnames(apobec3h) == "V1"] <- "Normalised_counts"
apobec3h$il6 <- thesis_datav4$LPS_IL_6
apobec3h$ctou <- thesis_datav4$CtoU
apobec3h$atoi <- thesis_datav4$AtoI_TOTALmatches
apobec3h$condition <- thesis_datav4$Time_point

# Box plots showing secretion of tnf-alpha across groups
plot(apobec3h$il6)

ggplot(apobec3h, aes(x = condition, y = il6)) +
  geom_boxplot() +
  labs(title = "IL-6 secretion across Groups",
       x = "Group",
       y = "IL-6 secretion")
apobec3h$il6


ggplot(apobec3h, aes(x = il6, fill = condition)) +
  geom_histogram(binwidth = 1000, color = "black", alpha = 0.7) +
  # scale_x_continuous(breaks = seq(-80512, 130790, by = 40000)) +  # Adjust the range based on your data
  labs(title = "Histogram with Stacked Bars of Conditions",
       x = "Il-6",
       y = "Frequency") +
  theme_minimal()

apobec3h <-na.omit(apobec3h)

apobec3h$il6_groups <- cut(apobec3h$il6, 
                           breaks = c(-Inf, 10000, 17500, 25000, 32000, 40000,  Inf),
                           labels = c("up to 10000", "10000-17500",
                                      "17500-25000", "25000-32000", "32000-40000", "over 40000"),
                           include.lowest = TRUE)


ggplot(apobec3h, aes(x = condition, y = log2(Normalised_counts + 1))) +
  geom_boxplot() +  
  geom_jitter(aes(color = as.factor(il6_groups)), width = 0.15) +  # Color points by LPS_IL_10
  labs(title = "Gene Expression of APOBEC3H across Conditions",
       x = "Condition",
       y = "Log2(Normalized Counts + 1)")

ggplot(apobec3h, aes(x = as.factor(il6_groups), y = log2(Normalised_counts + 1))) +
  geom_boxplot() +  
  labs(title = "Gene Expression of APOBEC3H across Il6 groups",
       x = "delta Il-6 secretion",
       y = "Log2(Normalized Counts + 1)")



# C to U counts vs section of il-1b
ggplot(apobec3h, aes(x = as.factor(tnf_groups), y = ctou)) +
  geom_boxplot() +  
  labs(title = "C to U count vs tnf secretion",
       x = "delta Il-1b secretion",
       y = "C to U count")


# box plots of A to I counts vs Il-1b secretion
ggplot(apobec3h, aes(x = as.factor(tnf_groups), y = atoi)) +
  geom_boxplot() +  
  labs(title = "A to I count vs tnf secretion",
       x = "delta Il-10 secretion",
       y = "A to I count")

# Filter data for the acute condition
acute_exp_data <- subset(apobec3h, condition == "acute")

# Create a scatter plot of all conditions together (change to acute_exp_data for 
# just the acute data)
ggplot(apobec3h, aes(x = ctou, y = log2(Normalised_counts + 1), color = age)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = "Relationship between log2Norm_counts of apobec3h and ctou variants)",
       x = "ctou",
       y = "log2Normalized Counts")

# Compute Spearman correlation
rcorr(apobec3h$Normalised_counts, apobec3h$ctou, type = "spearman")

acute_ctou_data <- subset(thesis_datav4, Time_point == "acute")
plot(acute_ctou_data$CtoU, acute_ctou_data$age)
rcorr(acute_ctou_data$CtoU, acute_ctou_data$age, type = "spearman")


#getting normalised counts of all significant rna editing genes in acute vs control
rna_edit_genes_apobec <- normalized_counts[c("ENSG00000128383", "ENSG00000244509", "ENSG00000100298"), ]

# Transpose data for plotting
rna_edit_genes_apobec <- t(rna_edit_genes_apobec)

rna_edit_genes_apobec <- as.data.frame(rna_edit_genes_apobec)
rna_edit_genes_apobec$condition <-coldata$condition
rna_edit_genes_apobec$id <- thesis_datav4$Participant_ID
rna_edit_genes_apobec$age <- thesis_datav4$age
rna_edit_genes_apobec$ctou <- thesis_datav4$CtoU
str(rna_edit_genes_apobec)

selected_columns_ctou <- c("ENSG00000128383", "ENSG00000244509", "ENSG00000100298", "ctou")
selected_columns_age <- c("ENSG00000128383", "ENSG00000244509", "ENSG00000100298", "age")

# Spearman correlation matrix between counts of rna_editing genes and c to u changes
# and age.
cor_matrix_ctou <- rcorr(as.matrix(rna_edit_genes_apobec[, selected_columns_ctou]), type = "spearman")
cor_matrix_age <- rcorr(as.matrix(rna_edit_genes_apobec[, selected_columns_age]), type = "spearman")

# Extract the correlation matrix
cor_matrix_ctou$r
cor_matrix_age$r

# Filter data for the acute condition
acute_exp_data_apobec <- subset(rna_edit_genes_apobec, condition == "acute")

# Spearman correlation matrix between counts of rna editing genes and 
# c to u and age only in acute group
acute_cor_matrix_ctou <- rcorr(as.matrix(acute_exp_data_apobec[, selected_columns_ctou]), type = "spearman")
acute_cor_matrix_age <- rcorr(as.matrix(acute_exp_data_apobec[, selected_columns_age]), type = "spearman")

# Extract the correlation matrix
acute_cor_matrix_ctou$r
acute_cor_matrix_age$r

# Comparing the gex of ENSG00000100298 (apobec3h) in runs in which a variant in 
# RELL1 (RELT-like protein 1) were detected.

# Define the sample names for acute condition
rell1_var_samples <- c("SRR12926256", "SRR12926262", "SRR12926269", "SRR12926284",
                   "SRR12926292", "SRR12926339", "SRR12926347", "SRR12926351")

rell1_var_samples_paired <- c("SRR12926256", "SRR12926262", "SRR12926269",
                       "SRR12926292", "SRR12926339", "SRR12926347")


rell1_var_samples_recovery <- c("SRR12926257", "SRR12926263", "SRR12926270", "SRR12926293",
                                "SRR12926340", "SRR12926348")

# Extract the normalised counts for ENSG00000100298 in acute condition
rell1_counts <- normalized_counts["ENSG00000100298", rell1_var_samples]
rell1_counts_paired <- normalized_counts["ENSG00000100298", rell1_var_samples_paired]
rell1_counts_rec <-normalized_counts["ENSG00000100298", rell1_var_samples_recovery]


# Extract the normalised counts for ENSG00000100298 in all other samples
no_rell1_counts <- normalized_counts["ENSG00000100298", !colnames(normalized_counts) %in% rell1_var_samples]



boxplot(log2(rell1_counts), log2(no_rell1_counts))
boxplot(log2(rell1_counts), log2(rell1_counts_rec))
boxplot(log2(rell1_counts_paired), log2(rell1_counts_rec))


# Create a data frame for each set of values
df_rell1 <- data.frame(
  Condition = rep("Acute", length(rell1_counts)),
  Log2Counts = log2(rell1_counts)
)

df_rell1_paired <- data.frame(
  Condition = rep("Acute", length(rell1_counts_paired)),
  Log2Counts = rell1_counts_paired
)

df_no_rell1 <- data.frame(
  Condition = rep("Other", length(no_rell1_counts)),
  Log2Counts = log2(no_rell1_counts)
)

df_rell1_counts_rec <- data.frame(
  Condition = rep("recovery", length(rell1_counts_rec)),
  Log2Counts = rell1_counts_rec
)

# Combine the data frames
combined_df <- rbind(df_rell1, df_no_rell1)
combined_df <- rbind(df_rell1, df_rell1_counts_rec)
combined_df <- rbind(df_rell1_paired, df_rell1_counts_rec)

combined_df

ggplot(combined_df, aes(x = Condition, y = Log2Counts, fill = Condition)) +
  geom_boxplot() +
  labs(title = "apobec3h gene expression in individuals with RELL1 variant",
       x = "Condition",
       y = "Normalised counts of Apobec3h")

combined_df$Run <- rownames(combined_df)

# Merge the data frames based on the "Run" column
merged_df <- merge(combined_df, thesis_datav4, by = "Run", all.x = TRUE)

ggplot(merged_df, aes(x = Condition, y = CtoU)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.75)) +
  geom_text(aes(label = Participant_ID), position = position_dodge(width = 0.75), vjust = -0.5) +  # Add labels
  geom_jitter(position = position_dodge(width = 0.75), alpha = 0.5) +
  labs(title = "Comparison of condition vs count of C to U variants",
       x = "Condition",
       y = "C to U variants")

ctou_in_paired_samples <- merged_df %>%
  select(Condition, CtoU) %>% 
  group_by(Condition) %>% 
  summarise(mean_ctou = mean(CtoU))

# displaying mean counts
ctou_in_paired_samples

ctou_in_paired_samples_age <- merged_df %>%
  select(Condition, CtoU, age) %>% 
  group_by(Condition) %>% 
  summarise(mean_age = mean(age))
ctou_in_paired_samples_age

# mean ages of acute, recovery and controls
ages<-thesis_datav4 %>% 
  select(Time_point, age) %>% 
  group_by(Time_point) %>% 
  summarise(mean_age = mean(age))
ages

# The 2:44175023 G>T variant in PPM1B 
PPM1B_var_samples <- c(
  "SRR12926207",
  "SRR12926212",
  "SRR12926224",
  "SRR12926264",
  "SRR12926266",
  "SRR12926271",
  "SRR12926278",
  "SRR12926349",
  "SRR12926352",
  "SRR12926354",
  "SRR12926366"
)

PPM1B_var_samples_paired<- c(
  "SRR12926212",
  "SRR12926264",
  "SRR12926271",
  "SRR12926278",
  "SRR12926349",
  "SRR12926352",
  "SRR12926366"
)

PPM1B_var_samples_paired_rec<- c(
  "SRR12926213",
  "SRR12926265",
  "SRR12926272",
  "SRR12926279",
  "SRR12926350",
  "SRR12926353",
  "SRR12926367"
)


# Extract the normalised counts for apobec3h in acute and rec condition
PPM1B_counts <- normalized_counts["ENSG00000100298", PPM1B_var_samples]
PPM1B_counts_paired <- normalized_counts["ENSG00000100298", PPM1B_var_samples_paired]
PPM1B_counts_paired_rec <- normalized_counts["ENSG00000100298", PPM1B_var_samples_paired_rec]

# Extract the normalised counts for apobec3h in all other samples
no_PPM1B_counts <- normalized_counts["ENSG00000100298", !colnames(normalized_counts) %in% PPM1B_var_samples]

boxplot(PPM1B_counts, no_PPM1B_counts)
boxplot(PPM1B_counts_paired, PPM1B_counts_paired_rec)

PPM1B_counts_combined <- data.frame(
  Sample = rep(c(names(PPM1B_counts_paired), names(PPM1B_counts_paired_rec)), each = 1),
  Counts = c(PPM1B_counts_paired, PPM1B_counts_paired_rec),
  Condition = rep(c("Acute", "Recovery"), each = length(PPM1B_counts_paired))
)

PPM1B_counts_combined

# Plotting of apobec3h exp and C to U variants counts in PPM1B variant individuals
ggplot(PPM1B_counts_combined, aes(x = Condition, y = Counts, fill = Condition)) +
  geom_boxplot() +
  labs(title = "apobec3h gene expression in individuals with PPM1B variant",
       x = "Condition",
       y = "Normalised counts of Apobec3h ")

PPM1B_counts_combined <- merge(PPM1B_counts_combined, thesis_datav4, by.x = "Sample", by.y = "Run", all.x = TRUE)

ggplot(PPM1B_counts_combined, aes(x = Condition, y = CtoU, fill = Condition)) +
  geom_boxplot() +
  labs(title = "C to U variants in individuals with PPM1B variant",
       x = "Condition",
       y = "Count of C to U variants")

# Finding mean C to U variant count
ctou_in_paired_samples <- PPM1B_counts_combined %>%
  select(Condition, CtoU) %>% 
  group_by(Condition) %>% 
  summarise(mean_ctou = mean(CtoU))

# displaying mean counts
ctou_in_paired_samples

# Removing samples without PPM1B and RELL1 variants

samples_to_remove <- c(
  "SRR12926256", "SRR12926262", "SRR12926269", "SRR12926284",
  "SRR12926292", "SRR12926339", "SRR12926347", "SRR12926351",
  "SRR12926212", "SRR12926264", "SRR12926271", "SRR12926278",
  "SRR12926349", "SRR12926352", "SRR12926366"
)

samples <- c(
  "SRR12926202", "SRR12926204", "SRR12926207", "SRR12926208",
  "SRR12926209", "SRR12926212", "SRR12926214", "SRR12926216",
  "SRR12926218", "SRR12926220", "SRR12926222", "SRR12926224",
  "SRR12926225", "SRR12926227", "SRR12926229", "SRR12926230",
  "SRR12926233", "SRR12926235", "SRR12926237", "SRR12926239",
  "SRR12926240", "SRR12926242", "SRR12926244", "SRR12926246",
  "SRR12926249", "SRR12926251", "SRR12926253", "SRR12926255",
  "SRR12926256", "SRR12926258", "SRR12926260", "SRR12926262",
  "SRR12926264", "SRR12926266", "SRR12926267", "SRR12926269",
  "SRR12926271", "SRR12926273", "SRR12926275", "SRR12926276",
  "SRR12926278", "SRR12926280", "SRR12926281", "SRR12926282",
  "SRR12926283", "SRR12926284", "SRR12926285", "SRR12926287",
  "SRR12926289", "SRR12926291", "SRR12926292", "SRR12926335",
  "SRR12926337", "SRR12926339", "SRR12926341", "SRR12926343",
  "SRR12926345", "SRR12926347", "SRR12926349", "SRR12926351",
  "SRR12926352", "SRR12926354", "SRR12926355", "SRR12926357",
  "SRR12926359", "SRR12926361", "SRR12926363", "SRR12926365",
  "SRR12926366"
)

remaining_samples <- setdiff(samples, samples_to_remove)
# Extract the normalised counts for apobec3h in acute and rec condition
apobec3h_counts_non_variants <- normalized_counts["ENSG00000100298", remaining_samples]
# Convert to data frame
apobec3h_counts_non_variants <- data.frame(
  Sample = names(apobec3h_counts_non_variants),
  Counts = as.numeric(apobec3h_counts_non_variants),
  row.names = NULL
)

# Add rest of data
apobec3h_counts_non_variants <- merge(apobec3h_counts_non_variants, thesis_datav4, by.x = "Sample", by.y = "Run", all.x = TRUE) %>%
  filter(!grepl("single", Paired_single))

id_to_extract <- apobec3h_counts_non_variants %>% 
  select(Participant_ID)

samples_to_extract <- thesis_datav4 %>% 
  select(Run, Participant_ID, Time_point) %>%
  filter(grepl("recovery", Time_point)) %>% 
  filter(Participant_ID %in% id_to_extract$Participant_ID)

samples_to_extract <- samples_to_extract %>% 
  select(Run)

samples_to_extract <- c(samples_to_extract$Run)
samples_to_extract

# Extract the normalised counts for apobec3h in all other samples
no_variants_counts_rec <- normalized_counts["ENSG00000100298", colnames(normalized_counts) %in% samples_to_extract]
str(no_variants_counts_rec)

#convert to dataframe
no_variants_counts_rec <- data.frame(
  Sample = names(no_variants_counts_rec),
  Counts = as.numeric(no_variants_counts_rec),
  row.names = NULL
)

no_variants_counts_rec
apobec3h_counts_non_variants

apobec3h_non_variants <-rbind(no_variants_counts_rec, apobec3h_counts_non_variants)

# Add rest of data
apobec3h_non_variants <- merge(apobec3h_non_variants, thesis_datav4, by.x = "Sample", by.y = "Run", all.x = TRUE) %>%
  filter(!grepl("single", Paired_single))

ggplot(apobec3h_non_variants, aes(x = Time_point, y = Counts, fill = Time_point)) +
  geom_boxplot() +
  labs(title = "apobec3h gene expression in individuals without PPM1B oR RELL1 variant",
       x = "Condition",
       y = "Normalised counts of Apobec3h ")

ggplot(apobec3h_non_variants, aes(x = Time_point, y = CtoU, fill = Time_point)) +
  geom_boxplot() +
  labs(title = "C to U variants in individuals without PPM1B oR RELL1 variant",
       x = "Condition",
       y = "Count of C to U variants")

#------------------------------------------------------------#
# Combining gene expression levels with the rest of the data

thesis_datav4_df <- as.data.frame(thesis_datav4)
normalized_counts_DEdds2 <- t(counts(DEdds2, normalized = TRUE))
str(normalized_counts_DEdds2)

# Extract the normalised counts data from the matrix
matrix_data <- as.matrix(normalized_counts_DEdds2)

# Create a data frame
matrix_dataframe <- as.data.frame(matrix_data)

# Set row names from the first list in dimnames
rownames(matrix_dataframe) <- dimnames(matrix_data)[[1]]

# Set column names from the second list in dimnames
colnames(matrix_dataframe) <- dimnames(matrix_data)[[2]]

# Convert row names of matrix_dataframe to a column named "Run"
matrix_dataframe$Run <- rownames(matrix_dataframe)

# Merged everything!!!
merged_df <- merge(thesis_datav4_df, matrix_dataframe, by = "Run", all.x = TRUE)


# plotting box plots of cytokine release data
ggplot(merged_df, aes(x = Time_point, y = log2(LPS_IFN_gamma), fill = Time_point)) +
  geom_boxplot() +
  labs(title = "LPS_IFN_gamma",
       x = "Condition",
       y = "log2 of delta LPS_IFN_gamma")

ggplot(merged_df, aes(x = Time_point, y = LPS_IL_10, fill = Time_point)) +
  geom_boxplot() +
  labs(title = "LPS_IL_10",
       x = "Condition",
       y = "delta LPS_IL_10")

ggplot(merged_df, aes(x = Time_point, y = log2(LPS_IL_27), fill = Time_point)) +
  geom_boxplot() +
  labs(title = "log 2 of LPS_IL_27",
       x = "Condition",
       y = "log 2 delta LPS_IL_27")

ggplot(merged_df, aes(x = Time_point, y = log2(LPS_IL_1b), fill = Time_point)) +
  geom_boxplot() +
  labs(title = "log 2 of LPS_IL_1b",
       x = "Condition",
       y = "log 2 of delta LPS_IL_1b")

ggplot(merged_df, aes(x = Time_point, y = LPS_IL_6, fill = Time_point)) +
  geom_boxplot() +
  labs(title = "LPS_IL_6",
       x = "Condition",
       y = "delta LPS_IL_6")

ggplot(merged_df, aes(x = Time_point, y = log2(LPS_IL_12p70), fill = Time_point)) +
  geom_boxplot() +
  labs(title = "log2 LPS_IL_12p70",
       x = "Condition",
       y = "log2 delta LPS_IL_12p70")

ggplot(merged_df, aes(x = Time_point, y = (log2(LPS_TNF_alpha)+1), fill = Time_point)) +
  geom_boxplot() +
  labs(title = "log2 LPS_TNF_alpha",
       x = "Condition",
       y = "log2 delta LPS_TNF_alpha")

#plots of acute important cytokine-release data
cytokines_acute <- merged_df %>%
  filter(Time_point == "acute") %>%
  select(Participant_ID, LPS_TNF_alpha, LPS_IL_1b, LPS_IL_6, LPS_IL_10)

cytokines_acute
plot(cytokines_acute$LPS_TNF_alpha)
plot(merged_df$LPS_TNF_alpha, col = merged_df$Time_point)
str(merged_df)
plot(cytokines_acute$LPS_IL_1b)
plot(cytokines_acute$LPS_IL_6)
plot(cytokines_acute$LPS_IL_10)
summary(cytokines_acute[,-1])


ggplot(merged_df, aes(x = Time_point, y = log2(IFN_gamma), fill = Time_point)) +
  geom_boxplot() +
  labs(title = "log2 Kp IFN_gamma",
       x = "Condition",
       y = "log2 delta Kp IFN_gamma")

ggplot(merged_df, aes(x = Time_point, y = log2(IL_10), fill = Time_point)) +
  geom_boxplot() +
  labs(title = "log2 Kp IL_10",
       x = "Condition",
       y = "log2 delta Kp IL_10")

ggplot(merged_df, aes(x = Time_point, y = log2(IL_27), fill = Time_point)) +
  geom_boxplot() +
  labs(title = "log2 Kp IL_27",
       x = "Condition",
       y = "log2 delta Kp IL_27")

ggplot(merged_df, aes(x = Time_point, y = log2(IL_1b), fill = Time_point)) +
  geom_boxplot() +
  labs(title = "log2 Kp IL_1b",
       x = "Condition",
       y = "log2 delta Kp IL_1b")

ggplot(merged_df, aes(x = Time_point, y = IL_6, fill = Time_point)) +
  geom_boxplot() +
  labs(title = "Kp IL_6",
       x = "Condition",
       y = "delta Kp IL_6")

ggplot(merged_df, aes(x = Time_point, y = log2(IL_12p70), fill = Time_point)) +
  geom_boxplot() +
  labs(title = "log2 Kp IL_12p70",
       x = "Condition",
       y = "log2 delta Kp IL_12p70")

ggplot(merged_df, aes(x = Time_point, y = log2(TNF_alpha), fill = Time_point)) +
  geom_boxplot() +
  labs(title = "log2 Kp TNF_alpha",
       x = "Condition",
       y = "log2 delta Kp TNF_alpha")

#Plots of the distribution of cytokines in the different groups:
vars_of_int <- merged_df %>% 
  select(LPS_TNF_alpha, ENSG00000100298 )
# Checking for correlation between TNF alpha release and CtoU variant count
# and Apobec3h normalised gene counts.
TNFvsApobec3h<-merged_df %>% 
  select(Time_point, Participant_ID, LPS_TNF_alpha, ENSG00000100298, CtoU) %>% 
  filter(Time_point == "acute") %>% 
  filter(LPS_TNF_alpha < 40000) # remove outliers

scatter.smooth(TNFvsApobec3h$CtoU, TNFvsApobec3h$LPS_TNF_alpha) # no relationship
scatter.smooth(TNFvsApobec3h$CtoU, TNFvsApobec3h$LPS_TNF_alpha)
rcorr(TNFvsApobec3h$CtoU, TNFvsApobec3h$LPS_TNF_alpha, type = "spearman")

merged_df %>% 
  select(Time_point, ENSG00000100298, LPS_TNF_alpha) %>% 
  filter(Time_point == "acute", LPS_TNF_alpha < 40000, ENSG00000100298 < 50) %>% 
  select(LPS_TNF_alpha, ENSG00000100298) %>% 
  plot()

merged_df %>% 
  select(Time_point, CtoU, LPS_TNF_alpha) %>% 
  filter(Time_point == "acute", LPS_TNF_alpha < 40000) %>% 
  select(LPS_TNF_alpha, CtoU) %>% 
  plot()

plot(merged_df$LPS_TNF_alpha)
plot(y=log2(merged_df$LPS_TNF_alpha), x=log2(merged_df$ENSG00000100298))

merged_df %>% 
  select(Time_point, ENSG00000100298, LPS_TNF_alpha) %>% 
  filter(Time_point == "acute", LPS_TNF_alpha < 40000, ENSG00000100298 < 50) %>% 
  select(LPS_TNF_alpha, ENSG00000100298) %>% 
  plot()

merged_df %>% 
  select(Time_point, CtoU, LPS_TNF_alpha) %>% 
  filter(Time_point == "acute", LPS_TNF_alpha < 40000) %>% 
  select(LPS_TNF_alpha, CtoU) %>% 
  plot()

# Is there a relationship between any of the APOBEC genes of interest 
# and the cytokine release data?
# ALL SAMPLES
covariates <- merged_df %>% 
  select(CtoU, AtoI_TOTALmatches, 
         ENSG00000128383, ENSG00000179750, ENSG00000244509,
         ENSG00000100298, LPS_TNF_alpha, LPS_IL_1b, 
         LPS_IL_6, LPS_IL_10)
covariates <-na.omit(covariates)
cor_matrix <- rcorr(as.matrix(covariates), type = "spearman")
cor_matrix

#Acute stage only:
covariates <- merged_df %>% 
  select(Time_point, CtoU, AtoI_TOTALmatches, 
         ENSG00000128383, ENSG00000179750, ENSG00000244509,
         ENSG00000100298, LPS_TNF_alpha, LPS_IL_1b, 
         LPS_IL_6, LPS_IL_10) %>% 
  filter(Time_point == "acute") %>% 
  select(CtoU, AtoI_TOTALmatches, 
         ENSG00000128383, ENSG00000179750, ENSG00000244509,
         ENSG00000100298, LPS_TNF_alpha, LPS_IL_1b, 
         LPS_IL_6, LPS_IL_10)

covariates <-na.omit(covariates)
cor_matrix <- rcorr(as.matrix(covariates), type = "spearman")
cor_matrix

plot(covariates$CtoU, covariates$AtoI_TOTALmatches)

# Is there a relationship between any of the ADAR genes of interest 
# and the cytokine release data?
# ALL SAMPLES

covariates <- merged_df %>% 
  select(CtoU, AtoI_TOTALmatches, 
         ENSG00000160710,
         ENSG00000197381,
         LPS_TNF_alpha, LPS_IL_1b, 
         LPS_IL_6, LPS_IL_10)
covariates <-na.omit(covariates)
cor_matrix <- rcorr(as.matrix(covariates), type = "spearman")
cor_matrix

#Acute stage only:
covariates <- merged_df %>% 
  select(Time_point, CtoU, AtoI_TOTALmatches, 
         ENSG00000160710,
         ENSG00000197381, LPS_TNF_alpha, LPS_IL_1b, 
         LPS_IL_6, LPS_IL_10) %>% 
  filter(Time_point == "acute") %>% 
  select(CtoU, AtoI_TOTALmatches, 
         ENSG00000160710,
         ENSG00000197381, LPS_TNF_alpha, LPS_IL_1b, 
         LPS_IL_6, LPS_IL_10)
covariates <-na.omit(covariates)
cor_matrix <- rcorr(as.matrix(covariates), type = "spearman")
cor_matrix

# Box plots of cytokine release 
# working with LPS data
# Acute Individuals with and without variants of interest
# rell1

with_rell1_var <- merged_df %>%
  select(Time_point, LPS_TNF_alpha, LPS_IL_1b, 
         LPS_IL_6, LPS_IL_10, Run) %>%
  filter(Time_point == "acute")
with_rell1_var <- with_rell1_var[with_rell1_var$Run %in% rell1_var_samples, ]
without_rell1_var <- merged_df %>%
  select(Time_point, LPS_TNF_alpha, LPS_IL_1b, 
         LPS_IL_6, LPS_IL_10, Run) %>%
  filter(Time_point == "acute")
without_rell1_var <- without_rell1_var[!(without_rell1_var$Run %in% rell1_var_samples), ]
without_rell1_var <- na.omit(without_rell1_var)
without_rell1_var

# Getting cytokine release data of the other groups as well

cytokines_controls <- merged_df %>%
  select(Time_point, LPS_TNF_alpha, LPS_IL_1b, 
         LPS_IL_6, LPS_IL_10) %>%
  filter(Time_point == "control")
cytokines_controls <- na.omit(cytokines_controls)

cytokines_recovery <- merged_df %>%
  select(Time_point, LPS_TNF_alpha, LPS_IL_1b, 
         LPS_IL_6, LPS_IL_10) %>%
  filter(Time_point == "recovery")
cytokines_recovery <- na.omit(cytokines_recovery)

boxplot(with_rell1_var$LPS_IL_10, without_rell1_var$LPS_IL_10, 
        cytokines_recovery$LPS_IL_10, cytokines_controls$LPS_IL_10,  
        main = "Comparison of IL10 release
        with RELL1 variant and the rest", 
        xlab = "Group", ylab = "delta LPS_IL_10", 
        names = c("Var", "!Var","Rec","Cont")) # IL-10

t_test_result <- t.test(with_rell1_var$LPS_IL_10, without_rell1_var$LPS_IL_10)
print(t_test_result)

wilcox_test_result <- wilcox.test(with_rell1_var$LPS_IL_10, without_rell1_var$LPS_IL_10)
print(wilcox_test_result)

boxplot(log2(with_rell1_var$LPS_IL_1b), log2(without_rell1_var$LPS_IL_1b), 
        log2(cytokines_recovery$LPS_IL_1b), log2(cytokines_controls$LPS_IL_1b),  
        main = "Comparison of IL-1b release
        with RELL1 variant and the rest", 
        xlab = "Group", ylab = "log2(delta LPS_IL_1b)", 
        names = c("Var", "!Var","Rec","Cont")) # IL-1b

boxplot(with_rell1_var$LPS_IL_6, without_rell1_var$LPS_IL_6, 
        cytokines_recovery$LPS_IL_6, cytokines_controls$LPS_IL_6,  
        main = "Comparison of IL6 release
        with RELL1 variant and the rest", 
        xlab = "Group", ylab = "delta LPS_IL_6", 
        names = c("Var", "!Var","Rec","Cont")) # IL-6

t_test_result <- t.test(with_rell1_var$LPS_IL_6, without_rell1_var$LPS_IL_6)
print(t_test_result)

wilcox_test_result <- wilcox.test(with_rell1_var$LPS_IL_6, without_rell1_var$LPS_IL_6)
print(wilcox_test_result)

boxplot(log2(with_rell1_var$LPS_TNF_alpha), log2(without_rell1_var$LPS_TNF_alpha), 
        log2(cytokines_recovery$LPS_TNF_alpha), log2(cytokines_controls$LPS_TNF_alpha),  
        main = "Comparison of TNF-a release
        with RELL1 variant and the rest", 
        xlab = "Group", ylab = "log 2 delta TNF-a", 
        names = c("Var", "!Var","Rec","Cont")) # TNF-a

t_test_result <- t.test(with_rell1_var$LPS_TNF_alpha, without_rell1_var$LPS_TNF_alpha)
print(t_test_result)

wilcox_test_result <- wilcox.test(with_rell1_var$LPS_TNF_alpha, without_rell1_var$LPS_TNF_alpha)
print(wilcox_test_result)

# PPM1B (PPM1B_var_samples)
with_PPM1B_var <- merged_df %>%
  select(Time_point, LPS_TNF_alpha, LPS_IL_1b, 
         LPS_IL_6, LPS_IL_10, Run) %>%
  filter(Time_point == "acute")
with_PPM1B_var <- with_PPM1B_var[with_PPM1B_var$Run %in% PPM1B_var_samples, ]
with_PPM1B_var <- na.omit(with_PPM1B_var)

without_PPM1B_var <- merged_df %>%
  select(Time_point, LPS_TNF_alpha, LPS_IL_1b, 
         LPS_IL_6, LPS_IL_10, Run) %>%
  filter(Time_point == "acute")
without_PPM1B_var <- without_PPM1B_var[!(without_PPM1B_var$Run %in% PPM1B_var_samples), ]
without_PPM1B_var <- na.omit(without_PPM1B_var)
without_PPM1B_var

#removing a extremely high Il-10 outlier from the acute group
without_PPM1B_var <- without_PPM1B_var[without_PPM1B_var$Run != "SRR12926276", ]

boxplot(with_PPM1B_var$LPS_IL_10, without_PPM1B_var$LPS_IL_10, 
        cytokines_recovery$LPS_IL_10, cytokines_controls$LPS_IL_10,  
        main = "Comparison of IL10 release
        with PPM1B variant and the rest", 
        xlab = "Group", ylab = "delta LPS_IL_10", 
        names = c("Var", "!Var","Rec","Cont")) # IL-10

boxplot(with_PPM1B_var$LPS_IL_1b, without_PPM1B_var$LPS_IL_1b, 
        cytokines_recovery$LPS_IL_1b, cytokines_controls$LPS_IL_1b,  
        main = "Comparison of IL1b release
        with PPM1B variant and the rest", 
        xlab = "Group", ylab = "delta LPS_IL_1b", 
        names = c("Var", "!Var","Rec","Cont")) # IL-1b

boxplot(with_PPM1B_var$LPS_IL_6, without_PPM1B_var$LPS_IL_6, 
        cytokines_recovery$LPS_IL_6, cytokines_controls$LPS_IL_6,  
        main = "Comparison of IL6 release
        with PPM1B variant and the rest", 
        xlab = "Group", ylab = "delta LPS_IL_6", 
        names = c("Var", "!Var","Rec","Cont")) # IL-6

t_test_result <- t.test(with_PPM1B_var$LPS_IL_6, without_PPM1B_var$LPS_IL_6)
print(t_test_result)

wilcox_test_result <- wilcox.test(with_PPM1B_var$LPS_IL_6, without_PPM1B_var$LPS_IL_6)
print(wilcox_test_result)

boxplot(log2(with_PPM1B_var$LPS_TNF_alpha), log2(without_PPM1B_var$LPS_TNF_alpha), 
        log2(cytokines_recovery$LPS_TNF_alpha), log2(cytokines_controls$LPS_TNF_alpha),  
        main = "Comparison of TNF-a release
        with PPM1B variant and the rest", 
        xlab = "Group", ylab = "log2 delta TNF-a", 
        names = c("Var", "!Var","Rec","Cont")) # TNF-a

t_test_result <- t.test(with_PPM1B_var$LPS_TNF_alpha, without_PPM1B_var$LPS_TNF_alpha)
print(t_test_result)

# Age

plot(merged_df$age, merged_df$AtoIandCtoU, main = "No of C to U variants vs Age", 
     xlab = "Age", ylab = "CtoU")

# No. of C to U vs age

plot(merged_df$age, merged_df$CtoU, main = "No of C to U variants vs Age", 
     xlab = "Age", ylab = "CtoU")

# Add a linear trend line
abline(lm(merged_df$CtoU ~ merged_df$age), col = "red")
rcorr(merged_df$CtoU, merged_df$age, type = "spearman")


merged_df %>%
  filter(Time_point == "acute") %>%
  ggplot(aes(x = age, y = CtoU)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, col = "blue") +
  labs(title = "Scatter Plot with Trend Line",
       x = "Age",
       y = "CtoU")

#------------------ Responders and non-responders --------------#
cytokine_column_names <- colnames(merged_df[,10:23])


cytokine_responses <- merged_df %>% 
  select(Participant_ID, Time_point, Paired_single, all_of(cytokine_column_names)) %>% 
  filter(Paired_single == "paired")

cytokine_responders <- cytokine_responses %>%
  group_by(Participant_ID) %>%
  mutate(across(all_of(cytokine_column_names), 
                ~ ifelse(.[Time_point == "acute"] < .[Time_point == "recovery"], "lower", "higher")))

cytokine_responders <- cytokine_responders %>% 
  filter(Time_point == "acute")

cytokine_responders <- as.data.frame(cytokine_responders)

#counting the number of "lower" and "higher" cytokine releases:
cytokine_responders$lower_count <- apply(cytokine_responders == "lower", 1, function(x) sum(x, na.rm = TRUE))
cytokine_responders$higher_count <- apply(cytokine_responders == "higher", 1, function(x) sum(x, na.rm = TRUE))

#renaming cytokine response type columns
names(cytokine_responders) <- gsub("^LPS", "LPS_resp", names(cytokine_responders))
names(cytokine_responders) <- gsub("^IL", "IL_resp", names(cytokine_responders))
cytokine_responders$TNF_resp_alpha <- cytokine_responders$TNF_alpha
cytokine_responders <- cytokine_responders %>% 
  select(-TNF_alpha)

# Getting runs SRR numbers
acute_paired_runs <- thesis_datav4_df %>%
  select(Run, Time_point, Paired_single) %>% 
  filter(Time_point == "acute") %>%
  filter(Paired_single == "paired") %>% 
  select(Run)

cytokine_responders$Run <-acute_paired_runs$Run

#Move Run as the first column
cytokine_responders <- cytokine_responders %>% relocate(Run)

